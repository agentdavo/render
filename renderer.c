// renderer.c
#include "renderer.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>
#include <math.h>

// External definitions (from other modules)
extern Light lights[MAX_LIGHTS];
extern int num_lights;

// Define framebuffer and zbuffer as global variables, partitioned per thread
ThreadBuffers thread_buffers[MAX_THREADS];

// Initialize framebuffer and zbuffer partitions
void initialize_buffers(ThreadBuffers buffers[MAX_THREADS]) {
    for(int i = 0; i < MAX_THREADS; i++) {
        buffers[i].framebuffer = aligned_alloc(16, sizeof(Color) * WIDTH * HEIGHT / MAX_THREADS * SSAA_FACTOR * SSAA_FACTOR);
        if(!buffers[i].framebuffer) {
            fprintf(stderr, "Error: Failed to allocate framebuffer for thread %d\n", i);
            exit(EXIT_FAILURE);
        }

        buffers[i].zbuffer = aligned_alloc(16, sizeof(float) * WIDTH * HEIGHT / MAX_THREADS * SSAA_FACTOR * SSAA_FACTOR);
        if(!buffers[i].zbuffer) {
            fprintf(stderr, "Error: Failed to allocate zbuffer for thread %d\n", i);
            exit(EXIT_FAILURE);
        }

        // Initialize zbuffer with maximum float values
        for(int j = 0; j < (WIDTH * HEIGHT / MAX_THREADS) * SSAA_FACTOR * SSAA_FACTOR; j++) {
            buffers[i].zbuffer[j] = FLT_MAX;
        }

        // Initialize framebuffer with transparent black
        for(int j = 0; j < (WIDTH * HEIGHT / MAX_THREADS) * SSAA_FACTOR * SSAA_FACTOR; j++) {
            buffers[i].framebuffer[j].r = 0;
            buffers[i].framebuffer[j].g = 0;
            buffers[i].framebuffer[j].b = 0;
            buffers[i].framebuffer[j].a = 0;
        }
    }
}

// Cleanup buffers
void cleanup_buffers(ThreadBuffers buffers[MAX_THREADS]) {
    for(int i = 0; i < MAX_THREADS; i++) {
        free(buffers[i].framebuffer);
        free(buffers[i].zbuffer);
    }
}

// Reference Counting Functions

RefCounted* retain_resource(RefCounted* resource) {
    if(resource) {
        atomic_fetch_add(&(resource->ref_count), 1);
    }
    return resource;
}

void release_resource(RefCounted* resource) {
    if(resource) {
        if(atomic_fetch_sub(&(resource->ref_count), 1) == 1) {
            switch(resource->type) {
                case RESOURCE_BSP_NODE:
                    free_bsp_tree_recursively((BSPNode*)resource);
                    break;
                case RESOURCE_TEXTURE:
                    // Implement texture-specific cleanup if necessary
                    free(resource);
                    break;
                default:
                    fprintf(stderr, "Error: Unknown resource type in release_resource\n");
                    break;
            }
        }
    }
}

// Mathematical utility functions

Vec3 add_vec3(const Vec3* a, const Vec3* b) {
    Vec3 result = {a->x + b->x, a->y + b->y, a->z + b->z};
    return result;
}

Vec3 subtract_vec3(const Vec3* a, const Vec3* b) {
    Vec3 result = {a->x - b->x, a->y - b->y, a->z - b->z};
    return result;
}

Vec3 multiply_vec3(const Vec3* a, float scalar) {
    Vec3 result = {a->x * scalar, a->y * scalar, a->z * scalar};
    return result;
}

float dot_product(const Vec3* a, const Vec3* b) {
    return (a->x * b->x + a->y * b->y + a->z * b->z);
}

void normalize(Vec3* v) {
    float length = sqrtf(v->x * v->x + v->y * v->y + v->z * v->z);
    if(length > 1e-6f) {
        v->x /= length;
        v->y /= length;
        v->z /= length;
    }
}

Mat4 create_look_at_matrix(const Vec3* eye, const Vec3* target, const Vec3* up) {
    Vec3 forward = subtract_vec3(target, eye);
    normalize(&forward);

    Vec3 side;
    // Compute cross product forward x up
    side.x = forward.y * up->z - forward.z * up->y;
    side.y = forward.z * up->x - forward.x * up->z;
    side.z = forward.x * up->y - forward.y * up->x;
    normalize(&side);

    Vec3 true_up;
    // Compute cross product side x forward
    true_up.x = side.y * forward.z - side.z * forward.y;
    true_up.y = side.z * forward.x - side.x * forward.z;
    true_up.z = side.x * forward.y - side.y * forward.x;

    Mat4 view = {0};
    view.m[0][0] = side.x;
    view.m[0][1] = side.y;
    view.m[0][2] = side.z;
    view.m[0][3] = -dot_product(&side, eye);

    view.m[1][0] = true_up.x;
    view.m[1][1] = true_up.y;
    view.m[1][2] = true_up.z;
    view.m[1][3] = -dot_product(&true_up, eye);

    view.m[2][0] = -forward.x;
    view.m[2][1] = -forward.y;
    view.m[2][2] = -forward.z;
    view.m[2][3] = dot_product(&forward, eye);

    view.m[3][3] = 1.0f;

    return view;
}

Mat4 create_perspective_matrix(float fov, float aspect, float near_plane, float far_plane) {
    float f = 1.0f / tanf((fov / 2.0f) * (M_PI / 180.0f));
    Mat4 proj = {0};
    proj.m[0][0] = f / aspect;
    proj.m[1][1] = f;
    proj.m[2][2] = (far_plane + near_plane) / (near_plane - far_plane);
    proj.m[2][3] = (2.0f * far_plane * near_plane) / (near_plane - far_plane);
    proj.m[3][2] = -1.0f;
    return proj;
}

Mat4 multiply_matrices(const Mat4* a, const Mat4* b) {
    Mat4 result = {0};
    for(int row = 0; row < 4; row++) {
        for(int col = 0; col < 4; col++) {
            for(int k = 0; k < 4; k++) {
                result.m[row][col] += a->m[row][k] * b->m[k][col];
            }
        }
    }
    return result;
}

Vec3 transform_vector(const Mat4* mat, const Vec3* vec) {
    Vec3 result;
    result.x = mat->m[0][0] * vec->x + mat->m[0][1] * vec->y + mat->m[0][2] * vec->z + mat->m[0][3];
    result.y = mat->m[1][0] * vec->x + mat->m[1][1] * vec->y + mat->m[1][2] * vec->z + mat->m[1][3];
    result.z = mat->m[2][0] * vec->x + mat->m[2][1] * vec->y + mat->m[2][2] * vec->z + mat->m[2][3];
    float w = mat->m[3][0] * vec->x + mat->m[3][1] * vec->y + mat->m[3][2] * vec->z + mat->m[3][3];
    if(w != 0.0f) {
        result.x /= w;
        result.y /= w;
        result.z /= w;
    }
    return result;
}

// Function to load ASCII STL models with improved error handling
int load_ascii_stl(const char* filename, Triangle* triangles, int* num_triangles) {
    FILE* file = fopen(filename, "r");
    if (!file) {
        fprintf(stderr, "Error: Unable to open STL file %s\n", filename);
        return -1;
    }

    char line[256];
    int triangle_count = 0;
    Triangle current_tri;
    int vertex_index = 0;

    // Default material initialization
    current_tri.material.diffuse = (Color){200, 200, 200, 255};
    current_tri.material.specular = (Color){255, 255, 255, 255};
    current_tri.material.emissive = (Color){0, 0, 0, 255};
    current_tri.material.shininess = 32.0f;

    while (fgets(line, sizeof(line), file)) {
        if (strncmp(line, "facet normal", 12) == 0) {
            float nx, ny, nz;
            if (sscanf(line, "facet normal %f %f %f", &nx, &ny, &nz) != 3) {
                fprintf(stderr, "Error: Malformed normal in STL file %s\n", filename);
                fclose(file);
                return -1;
            }
            current_tri.normal.x = FLOAT_TO_FIXED16_16(nx);
            current_tri.normal.y = FLOAT_TO_FIXED16_16(ny);
            current_tri.normal.z = FLOAT_TO_FIXED16_16(nz);
            vertex_index = 0; // Reset for new facet
        }
        else if (strncmp(line, "vertex", 6) == 0) {
            float x, y, z;
            if (sscanf(line, "vertex %f %f %f", &x, &y, &z) != 3) {
                fprintf(stderr, "Error: Malformed vertex in STL file %s\n", filename);
                fclose(file);
                return -1;
            }
            switch(vertex_index) {
                case 0:
                    current_tri.v0.x = FLOAT_TO_FIXED16_16(x);
                    current_tri.v0.y = FLOAT_TO_FIXED16_16(y);
                    current_tri.v0.z = FLOAT_TO_FIXED16_16(z);
                    current_tri.u0 = (x + 1.0f) / 2.0f;
                    current_tri.v0_tex = (y + 1.0f) / 2.0f;
                    break;
                case 1:
                    current_tri.v1.x = FLOAT_TO_FIXED16_16(x);
                    current_tri.v1.y = FLOAT_TO_FIXED16_16(y);
                    current_tri.v1.z = FLOAT_TO_FIXED16_16(z);
                    current_tri.u1 = (x + 1.0f) / 2.0f;
                    current_tri.v1_tex = (y + 1.0f) / 2.0f;
                    break;
                case 2:
                    current_tri.v2.x = FLOAT_TO_FIXED16_16(x);
                    current_tri.v2.y = FLOAT_TO_FIXED16_16(y);
                    current_tri.v2.z = FLOAT_TO_FIXED16_16(z);
                    current_tri.u2 = (x + 1.0f) / 2.0f;
                    current_tri.v2_tex = (y + 1.0f) / 2.0f;
                    break;
                default:
                    // More than 3 vertices per facet is invalid
                    fprintf(stderr, "Error: More than 3 vertices in facet in STL file %s\n", filename);
                    fclose(file);
                    return -1;
            }
            vertex_index++;
        }
        else if (strncmp(line, "endfacet", 8) == 0) {
            if(vertex_index != 3) {
                fprintf(stderr, "Error: Incomplete facet in STL file %s\n", filename);
                fclose(file);
                return -1;
            }

            // Normalize the normal vector
            Vec3 normal_float = {
                FIXED16_16_TO_FLOAT(current_tri.normal.x),
                FIXED16_16_TO_FLOAT(current_tri.normal.y),
                FIXED16_16_TO_FLOAT(current_tri.normal.z)
            };
            normalize(&normal_float);
            current_tri.normal.x = FLOAT_TO_FIXED16_16(normal_float.x);
            current_tri.normal.y = FLOAT_TO_FIXED16_16(normal_float.y);
            current_tri.normal.z = FLOAT_TO_FIXED16_16(normal_float.z);

            // Add to triangles array
            if (triangle_count < MAX_TRIANGLES) {
                triangles[triangle_count++] = current_tri;
            } else {
                fprintf(stderr, "Warning: Maximum number of triangles (%d) reached.\n", MAX_TRIANGLES);
                break;
            }
        }
    }

    fclose(file);
    *num_triangles = triangle_count;
    return 0;
}

// Comparator function for qsort_r
int compare_triangles(void *context, const void *a, const void *b) {
    int axis = *(int *)context;
    const Triangle *tri_a = (const Triangle *)a;
    const Triangle *tri_b = (const Triangle *)b;

    float centroid_a, centroid_b;
    switch (axis) {
        case 0:
            centroid_a = (FIXED16_16_TO_FLOAT(tri_a->v0.x) + FIXED16_16_TO_FLOAT(tri_a->v1.x) + FIXED16_16_TO_FLOAT(tri_a->v2.x)) / 3.0f;
            centroid_b = (FIXED16_16_TO_FLOAT(tri_b->v0.x) + FIXED16_16_TO_FLOAT(tri_b->v1.x) + FIXED16_16_TO_FLOAT(tri_b->v2.x)) / 3.0f;
            break;
        case 1:
            centroid_a = (FIXED16_16_TO_FLOAT(tri_a->v0.y) + FIXED16_16_TO_FLOAT(tri_a->v1.y) + FIXED16_16_TO_FLOAT(tri_a->v2.y)) / 3.0f;
            centroid_b = (FIXED16_16_TO_FLOAT(tri_b->v0.y) + FIXED16_16_TO_FLOAT(tri_b->v1.y) + FIXED16_16_TO_FLOAT(tri_b->v2.y)) / 3.0f;
            break;
        case 2:
            centroid_a = (FIXED16_16_TO_FLOAT(tri_a->v0.z) + FIXED16_16_TO_FLOAT(tri_a->v1.z) + FIXED16_16_TO_FLOAT(tri_a->v2.z)) / 3.0f;
            centroid_b = (FIXED16_16_TO_FLOAT(tri_b->v0.z) + FIXED16_16_TO_FLOAT(tri_b->v1.z) + FIXED16_16_TO_FLOAT(tri_b->v2.z)) / 3.0f;
            break;
        default:
            return 0;
    }

    if (centroid_a < centroid_b) return -1;
    if (centroid_a > centroid_b) return 1;
    return 0;
}

// Function to build BSP tree recursively with memory safety improvements
BSPNode* build_bsp_tree(Triangle* triangles, int num_triangles, int depth) {
    if(num_triangles == 0) return NULL;

    int axis = depth % 3;

    // Sort triangles based on centroid along the selected axis using qsort_r
    qsort_r(triangles, num_triangles, sizeof(Triangle), compare_triangles, &axis);

    int median = num_triangles / 2;
    Triangle pivot = triangles[median];

    // Calculate centroid of the pivot triangle
    float centroid_x = (FIXED16_16_TO_FLOAT(pivot.v0.x) + FIXED16_16_TO_FLOAT(pivot.v1.x) + FIXED16_16_TO_FLOAT(pivot.v2.x)) / 3.0f;
    float centroid_y = (FIXED16_16_TO_FLOAT(pivot.v0.y) + FIXED16_16_TO_FLOAT(pivot.v1.y) + FIXED16_16_TO_FLOAT(pivot.v2.y)) / 3.0f;
    float centroid_z = (FIXED16_16_TO_FLOAT(pivot.v0.z) + FIXED16_16_TO_FLOAT(pivot.v1.z) + FIXED16_16_TO_FLOAT(pivot.v2.z)) / 3.0f;

    // Compute plane distance 'd' using the plane equation: ax + by + cz + d = 0
    Vec3 normal_float = {
        FIXED16_16_TO_FLOAT(pivot.normal.x),
        FIXED16_16_TO_FLOAT(pivot.normal.y),
        FIXED16_16_TO_FLOAT(pivot.normal.z)
    };
    normalize(&normal_float);
    float d = -(normal_float.x * centroid_x + normal_float.y * centroid_y + normal_float.z * centroid_z);

    // Allocate memory for BSPNode using posix_memalign
    BSPNode* node;
    if(posix_memalign((void**)&node, 16, sizeof(BSPNode)) != 0) {
        fprintf(stderr, "Error: Failed to allocate memory for BSPNode\n");
        exit(EXIT_FAILURE);
    }

    node->normal = (Vec3){ FLOAT_TO_FIXED16_16(normal_float.x),
                            FLOAT_TO_FIXED16_16(normal_float.y),
                            FLOAT_TO_FIXED16_16(normal_float.z) };
    node->d = d;
    node->front = node->back = NULL;
    node->triangles = NULL;
    node->num_triangles = 0;
    node->base.ref_count = 1;
    node->base.type = RESOURCE_BSP_NODE;

    // Partition triangles into front and back sets
    int front_count = 0, back_count = 0;
    for(int i = 0; i < num_triangles; i++) {
        float centroid;
        switch(axis) {
            case 0:
                centroid = (FIXED16_16_TO_FLOAT(triangles[i].v0.x) + FIXED16_16_TO_FLOAT(triangles[i].v1.x) + FIXED16_16_TO_FLOAT(triangles[i].v2.x)) / 3.0f;
                break;
            case 1:
                centroid = (FIXED16_16_TO_FLOAT(triangles[i].v0.y) + FIXED16_16_TO_FLOAT(triangles[i].v1.y) + FIXED16_16_TO_FLOAT(triangles[i].v2.y)) / 3.0f;
                break;
            case 2:
                centroid = (FIXED16_16_TO_FLOAT(triangles[i].v0.z) + FIXED16_16_TO_FLOAT(triangles[i].v1.z) + FIXED16_16_TO_FLOAT(triangles[i].v2.z)) / 3.0f;
                break;
            default:
                centroid = 0.0f;
        }

        float distance = normal_float.x * centroid + normal_float.y * centroid + normal_float.z * centroid + d;

        if(distance >= 0.0f) {
            front_count++;
        } else {
            back_count++;
        }
    }

    // Allocate memory for front and back triangles
    Triangle* front_tris = NULL;
    Triangle* back_tris = NULL;
    if(front_count > 0) {
        front_tris = malloc(sizeof(Triangle) * front_count);
        if(!front_tris) {
            fprintf(stderr, "Error: Failed to allocate memory for front triangles\n");
            free(node);
            exit(EXIT_FAILURE);
        }
    }
    if(back_count > 0) {
        back_tris = malloc(sizeof(Triangle) * back_count);
        if(!back_tris) {
            fprintf(stderr, "Error: Failed to allocate memory for back triangles\n");
            free(front_tris);
            free(node);
            exit(EXIT_FAILURE);
        }
    }

    // Partition the triangles
    front_count = back_count = 0;
    for(int i = 0; i < num_triangles; i++) {
        float centroid;
        switch(axis) {
            case 0:
                centroid = (FIXED16_16_TO_FLOAT(triangles[i].v0.x) + FIXED16_16_TO_FLOAT(triangles[i].v1.x) + FIXED16_16_TO_FLOAT(triangles[i].v2.x)) / 3.0f;
                break;
            case 1:
                centroid = (FIXED16_16_TO_FLOAT(triangles[i].v0.y) + FIXED16_16_TO_FLOAT(triangles[i].v1.y) + FIXED16_16_TO_FLOAT(triangles[i].v2.y)) / 3.0f;
                break;
            case 2:
                centroid = (FIXED16_16_TO_FLOAT(triangles[i].v0.z) + FIXED16_16_TO_FLOAT(triangles[i].v1.z) + FIXED16_16_TO_FLOAT(triangles[i].v2.z)) / 3.0f;
                break;
            default:
                centroid = 0.0f;
        }

        float distance = normal_float.x * centroid + normal_float.y * centroid + normal_float.z * centroid + d;

        if(distance >= 0.0f) {
            if(front_tris)
                front_tris[front_count++] = triangles[i];
        } else {
            if(back_tris)
                back_tris[back_count++] = triangles[i];
        }
    }

    // Recursively build front and back subtrees
    node->front = build_bsp_tree(front_tris, front_count, depth + 1);
    node->back = build_bsp_tree(back_tris, back_count, depth + 1);

    // Free temporary triangle arrays
    free(front_tris);
    free(back_tris);

    return node;
}

// Traverse BSP tree and collect visible triangles
void traverse_bsp_tree(const BSPNode* node, const Camera* camera, Triangle** visible_tris, int* count) {
    if(node == NULL) return;

    Vec3 view_dir = subtract_vec3(&camera->target, &camera->position);
    normalize(&view_dir);

    float distance = FIXED16_16_TO_FLOAT(node->normal.x) * camera->position.x +
                     FIXED16_16_TO_FLOAT(node->normal.y) * camera->position.y +
                     FIXED16_16_TO_FLOAT(node->normal.z) * camera->position.z +
                     node->d;

    if(distance >= 0.0f) {
        traverse_bsp_tree(node->front, camera, visible_tris, count);
        for(int i = 0; i < node->num_triangles; i++) {
            visible_tris[(*count)++] = &node->triangles[i];
        }
        traverse_bsp_tree(node->back, camera, visible_tris, count);
    } else {
        traverse_bsp_tree(node->back, camera, visible_tris, count);
        for(int i = 0; i < node->num_triangles; i++) {
            visible_tris[(*count)++] = &node->triangles[i];
        }
        traverse_bsp_tree(node->front, camera, visible_tris, count);
    }
}

// Free BSP tree recursively with safety checks
void free_bsp_tree_recursively(BSPNode* node) {
    if(node == NULL) return;
    free_bsp_tree_recursively(node->front);
    free_bsp_tree_recursively(node->back);
    if(node->triangles) free(node->triangles);
    free(node);
}

// Shadow testing function
int is_in_shadow(const Vec3* frag_pos, const Vec3* normal, const Light* light, const ShadowMap* shadow_map) {
    // Transform fragment position into light's clip space
    Vec3 frag_pos_float = {
        FIXED16_16_TO_FLOAT(frag_pos->x),
        FIXED16_16_TO_FLOAT(frag_pos->y),
        FIXED16_16_TO_FLOAT(frag_pos->z)
    };

    Vec3 frag_light_space = transform_vector(&shadow_map->view_projection, &frag_pos_float);

    // Perform perspective divide and ensure we avoid division by zero
    float w = frag_light_space.z;
    if (fabsf(w) < 1e-5f) {
        return 0; // Avoid divide-by-zero; assume not in shadow
    }

    float x = frag_light_space.x / w;
    float y = frag_light_space.y / w;
    float z = frag_light_space.z / w;

    // Convert to [0, WIDTH] and [0, HEIGHT] for shadow map sampling
    float u = (x + 1.0f) * 0.5f * WIDTH;
    float v = (y + 1.0f) * 0.5f * HEIGHT;

    int tex_x = clamp((int)floorf(u), 0, WIDTH - 1);
    int tex_y = clamp((int)floorf(v), 0, HEIGHT - 1);

    // Retrieve depth from shadow map
    float shadow_depth = shadow_map->depth[tex_y * WIDTH + tex_x];

    // Bias to prevent shadow acne
    float bias = fmaxf(0.005f * (1.0f - dot_product(normal, &light->direction)), 0.001f);

    // If the fragment is behind the shadow map depth, it is in shadow
    return (z - bias > shadow_depth) ? 1 : 0;
}

// Shadow mapping thread function
void* generate_shadow_map_thread(void* arg) {
    ShadowMapTask* task = (ShadowMapTask*)arg;
    int thread_id = task->thread_id;
    int total_threads = task->total_threads;

    for(int i = thread_id; i < task->num_lights; i += total_threads) {
        const Light* light = &task->lights[i];
        ShadowMap* shadow_map = &task->shadow_maps[i];

        if(light->type == DIRECTIONAL_LIGHT || light->type == SPOTLIGHT) {
            Camera light_camera;
            if(light->type == DIRECTIONAL_LIGHT) {
                // Position the camera far in the opposite direction of the light
                light_camera.position = (Vec3){
                    -light->direction.x * 100.0f,
                    -light->direction.y * 100.0f,
                    -light->direction.z * 100.0f
                };
                light_camera.target = (Vec3){0.0f, 0.0f, 0.0f};
                light_camera.up = (Vec3){0.0f, 1.0f, 0.0f};
            } else { // SPOTLIGHT
                light_camera.position = light->position;
                Vec3 target_float = {
                    light->position.x + light->direction.x,
                    light->position.y + light->direction.y,
                    light->position.z + light->direction.z
                };
                light_camera.target = (Vec3){target_float.x, target_float.y, target_float.z};
                light_camera.up = (Vec3){0.0f, 1.0f, 0.0f};
            }

            Mat4 light_view = create_look_at_matrix(&light_camera.position, &light_camera.target, &light_camera.up);
            Mat4 light_projection = create_perspective_matrix(90.0f, 1.0f, 1.0f, 100.0f);
            Mat4 light_view_projection = multiply_matrices(&light_projection, &light_view);

            // Initialize shadow map depth buffer with FLT_MAX
            for(int j = 0; j < WIDTH * HEIGHT; j++) {
                shadow_map->depth[j] = FLT_MAX;
            }
            shadow_map->view_projection = light_view_projection;

            // Traverse BSP tree to get all triangles visible to the light
            Triangle* all_triangles[MAX_TRIANGLES];
            int triangle_count = 0;
            traverse_bsp_tree(task->bsp_root, &light_camera, all_triangles, &triangle_count);

            // Rasterize each triangle into the shadow map
            for(int t = 0; t < triangle_count; t++) {
                rasterize_triangle(all_triangles[t], &light_view, &light_projection, NULL, 0, NULL, 0, shadow_map, 1);
            }

            printf("Thread %d: Generated shadow map for light %d\n", thread_id, i);
        }
    }

    return NULL;
}

// Generate shadow maps in parallel
void generate_shadow_maps_parallel(const BSPNode* bsp_root, const Light* lights, int num_lights, ShadowMap* shadow_maps, int num_threads) {
    pthread_t threads[num_threads];
    ShadowMapTask tasks[num_threads];

    for(int i = 0; i < num_threads; i++) {
        tasks[i].shadow_maps = shadow_maps;
        tasks[i].lights = lights;
        tasks[i].num_lights = num_lights;
        tasks[i].bsp_root = bsp_root;
        tasks[i].thread_id = i;
        tasks[i].total_threads = num_threads;

        if(pthread_create(&threads[i], NULL, generate_shadow_map_thread, &tasks[i]) != 0) {
            fprintf(stderr, "Error: Failed to create shadow map thread %d\n", i);
            exit(EXIT_FAILURE);
        }
    }

    for(int i = 0; i < num_threads; i++) {
        pthread_join(threads[i], NULL);
    }
}

// Function to rasterize a triangle (with shadow mapping)
void rasterize_triangle(const Triangle* tri, const Mat4* view, const Mat4* projection, const Light* lights, int num_lights, Texture* textures, int num_textures, ShadowMap* shadow_maps, int num_shadow_maps) {
    // Transform vertices to clip space
    Vec3 v0_clip = transform_vector(projection, &tri->v0);
    Vec3 v1_clip = transform_vector(projection, &tri->v1);
    Vec3 v2_clip = transform_vector(projection, &tri->v2);

    // Clipping against frustum (simplified)
    Triangle clipped_tris[10];
    int clipped_count = clip_triangle_all_planes(tri, clipped_tris);
    if(clipped_count == 0) return;

    for(int i = 0; i < clipped_count; i++) {
        Triangle current_tri = clipped_tris[i];

        // Convert to screen space
        float x0 = (current_tri.v0.x + 1.0f) * 0.5f * WIDTH * SSAA_FACTOR;
        float y0 = (1.0f - current_tri.v0.y) * 0.5f * HEIGHT * SSAA_FACTOR;
        float z0 = current_tri.v0.z;

        float x1 = (current_tri.v1.x + 1.0f) * 0.5f * WIDTH * SSAA_FACTOR;
        float y1 = (1.0f - current_tri.v1.y) * 0.5f * HEIGHT * SSAA_FACTOR;
        float z1 = current_tri.v1.z;

        float x2 = (current_tri.v2.x + 1.0f) * 0.5f * WIDTH * SSAA_FACTOR;
        float y2 = (1.0f - current_tri.v2.y) * 0.5f * HEIGHT * SSAA_FACTOR;
        float z2 = current_tri.v2.z;

        // Bounding box
        int min_x = clamp((int)floorf(fminf(fminf(x0, x1), x2)), 0, WIDTH * SSAA_FACTOR - 1);
        int max_x = clamp((int)ceilf(fmaxf(fmaxf(x0, x1), x2)), 0, WIDTH * SSAA_FACTOR - 1);
        int min_y = clamp((int)floorf(fminf(fminf(y0, y1), y2)), 0, HEIGHT * SSAA_FACTOR - 1);
        int max_y = clamp((int)ceilf(fmaxf(fmaxf(y0, y1), y2)), 0, HEIGHT * SSAA_FACTOR - 1);

        // Edge functions
        float edge0_a = y1 - y0;
        float edge0_b = x0 - x1;
        float edge0_c = x1 * y0 - x0 * y1;

        float edge1_a = y2 - y1;
        float edge1_b = x1 - x2;
        float edge1_c = x2 * y1 - x1 * y2;

        float edge2_a = y0 - y2;
        float edge2_b = x2 - x0;
        float edge2_c = x0 * y2 - x2 * y0;

        // Area of the triangle
        float area = edge0_a * x2 + edge0_b * y2 + edge0_c;
        if (fabsf(area) < 1e-5f) continue; // Degenerate triangle

        // Iterate over the bounding box
        for(int y = min_y; y <= max_y; y++) {
            for(int x = min_x; x <= max_x; x++) {
                // Compute edge functions
                float f0 = edge0_a * x + edge0_b * y + edge0_c;
                float f1 = edge1_a * x + edge1_b * y + edge1_c;
                float f2 = edge2_a * x + edge2_b * y + edge2_c;

                // Check if inside triangle
                if((f0 >= 0 && f1 >= 0 && f2 >= 0) || (f0 <= 0 && f1 <= 0 && f2 <= 0)) {
                    // Barycentric coordinates
                    float w0 = f1 / area;
                    float w1 = f2 / area;
                    float w2 = 1.0f - w0 - w1;

                    // Interpolate depth
                    float z = w0 * z2 + w1 * z0 + w2 * z1;

                    // Determine which thread's buffer to write to based on y coordinate
                    int thread_id = y / (HEIGHT / MAX_THREADS);
                    if(thread_id >= MAX_THREADS) thread_id = MAX_THREADS - 1;

                    int idx = (y % (HEIGHT / MAX_THREADS)) * WIDTH * SSAA_FACTOR + x;

                    // Depth test
                    if(z < thread_buffers[thread_id].zbuffer[idx]) {
                        thread_buffers[thread_id].zbuffer[idx] = z;

                        // Interpolate texture coordinates
                        float u = w0 * current_tri.u2 + w1 * current_tri.u0 + w2 * current_tri.u1;
                        float v = w0 * current_tri.v2_tex + w1 * current_tri.v0_tex + w2 * current_tri.v1_tex;

                        // Select texture (for simplicity, using first texture)
                        Color tex_color = {255, 255, 255, 255}; // Default white
                        if(textures && num_textures > 0) {
                            // Nearest-neighbor sampling
                            int tex_x = clamp((int)(u * (TEXTURE_SIZE - 1)), 0, TEXTURE_SIZE - 1);
                            int tex_y = clamp((int)(v * (TEXTURE_SIZE - 1)), 0, TEXTURE_SIZE - 1);
                            tex_color.r = textures[0].data[tex_y][tex_x][0];
                            tex_color.g = textures[0].data[tex_y][tex_x][1];
                            tex_color.b = textures[0].data[tex_y][tex_x][2];
                            tex_color.a = 255;
                        }

                        // Compute fragment position (approximation)
                        Vec3 frag_pos = {
                            w0 * FIXED16_16_TO_FLOAT(current_tri.v2.x) + w1 * FIXED16_16_TO_FLOAT(current_tri.v0.x) + w2 * FIXED16_16_TO_FLOAT(current_tri.v1.x),
                            w0 * FIXED16_16_TO_FLOAT(current_tri.v2.y) + w1 * FIXED16_16_TO_FLOAT(current_tri.v0.y) + w2 * FIXED16_16_TO_FLOAT(current_tri.v1.y),
                            z
                        };

                        // Compute view direction (camera assumed at origin)
                        Vec3 view_dir = {0.0f - frag_pos.x,
                                         0.0f - frag_pos.y,
                                         0.0f - frag_pos.z};
                        normalize(&view_dir);

                        // Apply shadow mapping
                        int shadow = 0;
                        for(int sm = 0; sm < num_shadow_maps; sm++) {
                            shadow += is_in_shadow(&frag_pos, &current_tri.normal, &lights[sm], &shadow_maps[sm]);
                        }

                        // Compute lighting
                        Color lighting = compute_lighting(&current_tri.normal, &view_dir, lights, num_lights, &current_tri.material, &frag_pos);
                        if (shadow) {
                            lighting.r = clamp((int)(lighting.r * 0.5f), 0, 255);
                            lighting.g = clamp((int)(lighting.g * 0.5f), 0, 255);
                            lighting.b = clamp((int)(lighting.b * 0.5f), 0, 255);
                        }

                        // Fog calculation based on distance
                        float distance = sqrtf(frag_pos.x * frag_pos.x +
                                               frag_pos.y * frag_pos.y +
                                               frag_pos.z * frag_pos.z);

                        Color fog_color = (Color){100, 100, 150, 255}; // Fog color
                        float fog_density = 0.05f; // Adjust for effect
                        Color final_color = apply_fog(&lighting, distance, fog_density, &fog_color);

                        // Combine texture color with lighting and fog
                        final_color.r = clamp((int)(tex_color.r * (final_color.r / 255.0f)), 0, 255);
                        final_color.g = clamp((int)(tex_color.g * (final_color.g / 255.0f)), 0, 255);
                        final_color.b = clamp((int)(tex_color.b * (final_color.b / 255.0f)), 0, 255);
                        final_color.a = 255;

                        // Write to framebuffer
                        thread_buffers[thread_id].framebuffer[idx] = final_color;

                        // Increment fill rate counters
                        increment_pixels_filled();
                        increment_textures_filled();
                    }
                }
            }
        }
    }
}

// Lighting computation function (Example implementation using Phong shading)
Color compute_lighting(const Vec3* normal, const Vec3* view_dir, const Light* lights, int num_lights, const Material* material, const Vec3* frag_pos) {
    Color result = {0, 0, 0, 255};
    for(int i = 0; i < num_lights; i++) {
        const Light* light = &lights[i];
        Vec3 light_dir;
        if(light->type == DIRECTIONAL_LIGHT) {
            light_dir = multiply_vec3(&light->direction, -1.0f);
        } else { // POINT_LIGHT or SPOTLIGHT
            Vec3 light_pos = {light->position.x, light->position.y, light->position.z};
            Vec3 frag_pos_float = {FIXED16_16_TO_FLOAT(frag_pos->x),
                                   FIXED16_16_TO_FLOAT(frag_pos->y),
                                   FIXED16_16_TO_FLOAT(frag_pos->z)};
            Vec3 to_light = subtract_vec3(&light_pos, &frag_pos_float);
            normalize(&to_light);
            light_dir = to_light;
        }

        // Diffuse shading
        float diff = fmaxf(dot_product(normal, &light_dir), 0.0f);
        Color diffuse = (Color){
            (uint8_t)(material->diffuse.r * light->color.r / 255.0f * diff * light->intensity),
            (uint8_t)(material->diffuse.g * light->color.g / 255.0f * diff * light->intensity),
            (uint8_t)(material->diffuse.b * light->color.b / 255.0f * diff * light->intensity),
            255
        };

        // Specular shading
        Vec3 temp = multiply_vec3(normal, 2.0f * dot_product(normal, &light_dir));
        Vec3 reflect_dir = subtract_vec3(&temp, &light_dir);
        float spec = powf(fmaxf(dot_product(&reflect_dir, view_dir), 0.0f), material->shininess);
        Color specular = (Color){
            (uint8_t)(material->specular.r * light->color.r / 255.0f * spec * light->intensity),
            (uint8_t)(material->specular.g * light->color.g / 255.0f * spec * light->intensity),
            (uint8_t)(material->specular.b * light->color.b / 255.0f * spec * light->intensity),
            255
        };

        // Ambient shading (optional)
        Color ambient = (Color){
            (uint8_t)(material->diffuse.r * 0.1f),
            (uint8_t)(material->diffuse.g * 0.1f),
            (uint8_t)(material->diffuse.b * 0.1f),
            255
        };

        // Accumulate colors
        result.r = clamp(result.r + diffuse.r + specular.r + ambient.r, 0, 255);
        result.g = clamp(result.g + diffuse.g + specular.g + ambient.g, 0, 255);
        result.b = clamp(result.b + diffuse.b + specular.b + ambient.b, 0, 255);
    }
    return result;
}

// Function to load a texture from BMP file (Stub Implementation)
int load_texture_from_bmp(const char* filename, Texture* tex) {
    // Implement BMP loading logic here
    // For now, return 0 to indicate success
    // In a real implementation, you would parse the BMP file and fill tex->data
    return 0;
}

// Initialize checkerboard texture
void init_checkerboard_texture(Texture* tex) {
    for(int y = 0; y < TEXTURE_SIZE; y++) {
        for(int x = 0; x < TEXTURE_SIZE; x++) {
            if(((x / 16) % 2) == ((y / 16) % 2)) {
                tex->data[y][x][0] = 255; // Red
                tex->data[y][x][1] = 255; // Green
                tex->data[y][x][2] = 255; // Blue
            } else {
                tex->data[y][x][0] = 0;
                tex->data[y][x][1] = 0;
                tex->data[y][x][2] = 0;
            }
        }
    }
}

// Fill rate tracking functions (Example Implementation)
static atomic_int pixels_filled = 0;
static atomic_int textures_filled = 0;

void increment_pixels_filled() {
    atomic_fetch_add(&pixels_filled, 1);
}

void increment_textures_filled() {
    atomic_fetch_add(&textures_filled, 1);
}

// Clipping function (stub implementation)
int clip_triangle_all_planes(const Triangle* tri, Triangle* clipped_tris) {
    // Implement actual clipping against frustum planes
    // For simplicity, assume no clipping is necessary
    clipped_tris[0] = *tri;
    return 1;
}

// Fog application function
Color apply_fog(const Color* original, float distance, float fog_density, const Color* fog_color) {
    // Simple exponential fog
    float fog_factor = expf(-distance * fog_density);
    fog_factor = clamp(fog_factor, 0.0f, 1.0f);

    // Linearly interpolate between original color and fog color
    Color final_color;
    final_color.r = (uint8_t)(original->r * fog_factor + fog_color->r * (1.0f - fog_factor));
    final_color.g = (uint8_t)(original->g * fog_factor + fog_color->g * (1.0f - fog_factor));
    final_color.b = (uint8_t)(original->b * fog_factor + fog_color->b * (1.0f - fog_factor));
    final_color.a = original->a;  // Preserve alpha

    return final_color;
}

// Function to initialize font
void init_font(CharBitmap* font) {
    // Simple 8x8 font for ASCII characters 32-126
    // Initialize all characters to blank
    for(int i = 0; i < 95; i++) { // ASCII 32 to 126
        font[i].character = (char)(i + 32);
        memset(font[i].bitmap, 0, 8);
    }

    // Define 'A'
    font['A' - 32].bitmap[0] = 0b00011000;
    font['A' - 32].bitmap[1] = 0b00100100;
    font['A' - 32].bitmap[2] = 0b01000010;
    font['A' - 32].bitmap[3] = 0b01000010;
    font['A' - 32].bitmap[4] = 0b01111110;
    font['A' - 32].bitmap[5] = 0b01000010;
    font['A' - 32].bitmap[6] = 0b01000010;
    font['A' - 32].bitmap[7] = 0b01000010;

    // Define 'B'
    font['B' - 32].bitmap[0] = 0b01111100;
    font['B' - 32].bitmap[1] = 0b01000010;
    font['B' - 32].bitmap[2] = 0b01000010;
    font['B' - 32].bitmap[3] = 0b01111100;
    font['B' - 32].bitmap[4] = 0b01000010;
    font['B' - 32].bitmap[5] = 0b01000010;
    font['B' - 32].bitmap[6] = 0b01000010;
    font['B' - 32].bitmap[7] = 0b01111100;

    // Define 'C'
    font['C' - 32].bitmap[0] = 0b00111100;
    font['C' - 32].bitmap[1] = 0b01000010;
    font['C' - 32].bitmap[2] = 0b10000000;
    font['C' - 32].bitmap[3] = 0b10000000;
    font['C' - 32].bitmap[4] = 0b10000000;
    font['C' - 32].bitmap[5] = 0b10000000;
    font['C' - 32].bitmap[6] = 0b01000010;
    font['C' - 32].bitmap[7] = 0b00111100;

    // Define '0'
    font['0' - 32].bitmap[0] = 0b00111100;
    font['0' - 32].bitmap[1] = 0b01000010;
    font['0' - 32].bitmap[2] = 0b01000110;
    font['0' - 32].bitmap[3] = 0b01001010;
    font['0' - 32].bitmap[4] = 0b01010010;
    font['0' - 32].bitmap[5] = 0b01100010;
    font['0' - 32].bitmap[6] = 0b01000010;
    font['0' - 32].bitmap[7] = 0b00111100;
}

// Render text onto framebuffer
void render_text(const char* text, int x, int y, const Color* color, const CharBitmap* font) {
    int start_x = x;
    for(int i = 0; text[i] != '\0'; i++) {
        char c = text[i];
        if(c < 32 || c > 126) continue; // Unsupported character
        const CharBitmap* cb = &font[c - 32];
        for(int row = 0; row < 8; row++) {
            for(int col = 0; col < 8; col++) {
                if(cb->bitmap[row] & (1 << (7 - col))) {
                    int px = x + col;
                    int py = y + row;
                    if(px >= 0 && px < WIDTH * SSAA_FACTOR && py >= 0 && py < HEIGHT * SSAA_FACTOR) {
                        // Determine which thread's buffer to write to based on y coordinate
                        int thread_id = py / (HEIGHT / MAX_THREADS);
                        if(thread_id >= MAX_THREADS) thread_id = MAX_THREADS - 1;

                        int idx = (py % (HEIGHT / MAX_THREADS)) * WIDTH * SSAA_FACTOR + px;
                        thread_buffers[thread_id].framebuffer[idx].r = color->r;
                        thread_buffers[thread_id].framebuffer[idx].g = color->g;
                        thread_buffers[thread_id].framebuffer[idx].b = color->b;
                        thread_buffers[thread_id].framebuffer[idx].a = color->a;
                    }
                }
            }
        }
        x += 8; // Move to next character position
    }
}

// Combine thread buffers into a single framebuffer
void combine_framebuffers(ThreadBuffers buffers[MAX_THREADS], Color* final_framebuffer) {
    for(int t = 0; t < MAX_THREADS; t++) {
        // Assuming each thread's buffer corresponds to a horizontal partition
        // Adjust as needed based on actual partitioning strategy
        int partition_height = HEIGHT / MAX_THREADS;
        for(int y = 0; y < partition_height * SSAA_FACTOR; y++) {
            for(int x = 0; x < WIDTH * SSAA_FACTOR; x++) {
                int final_idx = (t * partition_height * SSAA_FACTOR + y) * WIDTH * SSAA_FACTOR + x;
                int thread_idx = y * WIDTH * SSAA_FACTOR + x;
                final_framebuffer[final_idx] = buffers[t].framebuffer[thread_idx];
            }
        }
    }
}

// Main rendering function
void render_scene(const BSPNode* bsp_root, const Camera* camera, const Light* lights, int num_lights, Texture* textures, int num_textures, ShadowMap* shadow_maps, int num_shadow_maps, Color* final_framebuffer) {
    // Initialize thread buffers
    initialize_buffers(thread_buffers);

    // Generate shadow maps in parallel
    generate_shadow_maps_parallel(bsp_root, lights, num_lights, shadow_maps, MAX_THREADS);

    // Traverse BSP tree to get visible triangles
    Triangle* visible_tris[MAX_TRIANGLES];
    int visible_count = 0;
    traverse_bsp_tree(bsp_root, camera, visible_tris, &visible_count);

    // Rasterize each visible triangle
    for(int i = 0; i < visible_count; i++) {
        rasterize_triangle(visible_tris[i], NULL, NULL, lights, num_lights, textures, num_textures, shadow_maps, num_shadow_maps);
    }

    // Render text (example usage)
    CharBitmap font[95];
    init_font(font);
    Color text_color = {255, 255, 255, 255}; // White color
    render_text("ABC 012", 10, 10, &text_color, font);

    // Combine thread buffers into the final framebuffer
    combine_framebuffers(thread_buffers, final_framebuffer);

    // Cleanup thread buffers
    cleanup_buffers(thread_buffers);
}

// Example usage (main function)
#ifdef RENDERER_MAIN
int main(int argc, char* argv[]) {
    if(argc < 2) {
        printf("Usage: %s <model.stl>\n", argv[0]);
        return -1;
    }

    // Load model
    Triangle triangles[MAX_TRIANGLES];
    int num_triangles = 0;
    if(load_ascii_stl(argv[1], triangles, &num_triangles) != 0) {
        return -1;
    }

    printf("Loaded %d triangles from %s\n", num_triangles, argv[1]);

    // Build BSP tree
    BSPNode* bsp_root = build_bsp_tree(triangles, num_triangles, 0);
    printf("BSP tree built.\n");

    // Initialize camera
    Camera camera = {
        .position = (Vec3){0.0f, 0.0f, -5.0f},
        .target = (Vec3){0.0f, 0.0f, 0.0f},
        .up = (Vec3){0.0f, 1.0f, 0.0f}
    };

    // Initialize lights
    Light scene_lights[MAX_LIGHTS];
    num_lights = 2;

    // Directional light
    scene_lights[0].type = DIRECTIONAL_LIGHT;
    scene_lights[0].direction = (Vec3){-1.0f, -1.0f, -1.0f};
    scene_lights[0].intensity = 1.0f;
    scene_lights[0].color = (Color){255, 255, 255, 255};

    // Point light
    scene_lights[1].type = POINT_LIGHT;
    scene_lights[1].position = (Vec3){2.0f, 2.0f, -3.0f};
    scene_lights[1].direction = (Vec3){0.0f, 0.0f, 0.0f}; // Not used for point light
    scene_lights[1].intensity = 0.8f;
    scene_lights[1].color = (Color){255, 200, 200, 255};

    // Initialize textures
    Texture textures[MAX_TEXTURES];
    // For simplicity, initialize checkerboard texture
    init_checkerboard_texture(&textures[0]);

    // Initialize shadow maps
    ShadowMap shadow_maps[MAX_LIGHTS];

    // Allocate final framebuffer
    Color* final_framebuffer = malloc(sizeof(Color) * WIDTH * HEIGHT * SSAA_FACTOR * SSAA_FACTOR);
    if(!final_framebuffer) {
        fprintf(stderr, "Error: Failed to allocate final framebuffer\n");
        return -1;
    }

    // Render the scene
    render_scene(bsp_root, &camera, scene_lights, num_lights, textures, 1, shadow_maps, num_lights, final_framebuffer);

    // TODO: Save or display the framebuffer as an image

    // Cleanup
    free(final_framebuffer);
    free_bsp_tree_recursively(bsp_root);

    printf("Rendering complete.\n");
    return 0;
}
#endif // RENDERER_MAIN
