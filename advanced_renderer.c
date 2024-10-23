// advanced_renderer.c
#define _POSIX_C_SOURCE 200112L

#define GL_SILENCE_DEPRECATION  // Silences OpenGL deprecation warnings on macOS

#include <GLFW/glfw3.h>
#include <pthread.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>

#include "renderer.h"

// Framebuffer and Z-buffer are defined in renderer.c
extern Color framebuffer[];
extern float zbuffer[];

// Texture structures
Texture textures[MAX_LIGHTS]; // Assuming one texture per material for simplicity

// Shadow maps
ShadowMap shadow_maps[MAX_SHADOW_MAPS];
int shadow_map_count = 0;

// Timer structure for profiling
typedef struct {
    double start_time;
    double end_time;
} Timer;

// Function to start the timer
void start_timer(Timer* timer) {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    timer->start_time = tv.tv_sec + tv.tv_usec / 1000000.0;
}

// Function to end the timer
void end_timer(Timer* timer) {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    timer->end_time = tv.tv_sec + tv.tv_usec / 1000000.0;
}

// Function to get elapsed time in seconds
double get_elapsed_time(const Timer* timer) {
    return timer->end_time - timer->start_time;
}

// Function to clear framebuffer and zbuffer
void clear_buffers() {
    memset(framebuffer, 0, sizeof(Color) * WIDTH * SSAA_FACTOR * HEIGHT * SSAA_FACTOR);
    memset(zbuffer, 0, sizeof(float) * WIDTH * SSAA_FACTOR * HEIGHT * SSAA_FACTOR);
}

// Downsample framebuffer for anti-aliasing
void downsample_framebuffer(Color* dest, const Color* src) {
    for(int y = 0; y < HEIGHT; y++) {
        for(int x = 0; x < WIDTH; x++) {
            int r = 0, g = 0, b = 0;
            for(int dy = 0; dy < SSAA_FACTOR; dy++) {
                for(int dx = 0; dx < SSAA_FACTOR; dx++) {
                    int src_idx = (y * SSAA_FACTOR + dy) * WIDTH * SSAA_FACTOR + (x * SSAA_FACTOR + dx);
                    r += src[src_idx].r;
                    g += src[src_idx].g;
                    b += src[src_idx].b;
                }
            }
            int samples = SSAA_FACTOR * SSAA_FACTOR;
            dest[y * WIDTH + x].r = r / samples;
            dest[y * WIDTH + x].g = g / samples;
            dest[y * WIDTH + x].b = b / samples;
            dest[y * WIDTH + x].a = 255;
        }
    }
}

// Function to create OpenGL texture from framebuffer
GLuint create_texture(const Color* framebuffer_downsampled) {
    GLuint texture;
    glGenTextures(1, &texture);
    if(texture == 0) {
        fprintf(stderr, "Error: Failed to generate OpenGL texture\n");
        exit(EXIT_FAILURE);
    }
    glBindTexture(GL_TEXTURE_2D, texture);

    // Set texture parameters
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

    // Upload framebuffer data
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, WIDTH, HEIGHT, 0, GL_RGBA, GL_UNSIGNED_BYTE, framebuffer_downsampled);

    return texture;
}

// Function to render fullscreen quad with texture
void render_fullscreen_quad(GLuint texture) {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D, texture);

    glBegin(GL_QUADS);
        glTexCoord2f(0.0f, 0.0f); glVertex2f(-1.0f, -1.0f);
        glTexCoord2f(1.0f, 0.0f); glVertex2f( 1.0f, -1.0f);
        glTexCoord2f(1.0f, 1.0f); glVertex2f( 1.0f,  1.0f);
        glTexCoord2f(0.0f, 1.0f); glVertex2f(-1.0f,  1.0f);
    glEnd();

    glDisable(GL_TEXTURE_2D);
}

// Structure to hold triangles per thread
typedef struct {
    int thread_id;
    int start_idx;
    int end_idx;
    Triangle* triangles;
    Light* lights;
    int num_lights;
    Texture* textures;
    int num_textures;
    ShadowMap* shadow_maps;
    int num_shadow_maps;
} ThreadTask;

// Function to rasterize triangles in a thread
void* thread_rasterize(void* arg) {
    ThreadTask* task = (ThreadTask*)arg;
    for(int i = task->start_idx; i < task->end_idx; i++) {
        rasterize_triangle(&task->triangles[i], NULL, NULL, task->lights, task->num_lights, task->textures, task->num_textures, task->shadow_maps, task->num_shadow_maps);
    }
    pthread_exit(NULL);
}

// Function to collect visible triangles with frustum culling
int collect_visible_triangles(const BSPNode* bsp, const Camera* camera, Triangle** visible_tris) {
    int count = 0;
    traverse_bsp_tree(bsp, camera, visible_tris, &count);
    return count;
}

// Function to handle user input for camera movement
void handle_input(GLFWwindow* window, Camera* camera, float delta_time) {
    const float speed = 2.0f * delta_time;
    Vec3 forward = {
        FIXED16_16_TO_FLOAT(camera->target.x) - FIXED16_16_TO_FLOAT(camera->position.x),
        FIXED16_16_TO_FLOAT(camera->target.y) - FIXED16_16_TO_FLOAT(camera->position.y),
        FIXED16_16_TO_FLOAT(camera->target.z) - FIXED16_16_TO_FLOAT(camera->position.z)
    };
    normalize(&forward);

    // Compute right vector
    Vec3 up = camera->up;
    Vec3 right = cross_product(&forward, &up);
    normalize(&right);

    // Move forward
    if(glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS) {
        camera->position.x += FLOAT_TO_FIXED16_16(forward.x * speed);
        camera->position.y += FLOAT_TO_FIXED16_16(forward.y * speed);
        camera->position.z += FLOAT_TO_FIXED16_16(forward.z * speed);
        camera->target.x += FLOAT_TO_FIXED16_16(forward.x * speed);
        camera->target.y += FLOAT_TO_FIXED16_16(forward.y * speed);
        camera->target.z += FLOAT_TO_FIXED16_16(forward.z * speed);
    }
    // Move backward
    if(glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS) {
        camera->position.x -= FLOAT_TO_FIXED16_16(forward.x * speed);
        camera->position.y -= FLOAT_TO_FIXED16_16(forward.y * speed);
        camera->position.z -= FLOAT_TO_FIXED16_16(forward.z * speed);
        camera->target.x -= FLOAT_TO_FIXED16_16(forward.x * speed);
        camera->target.y -= FLOAT_TO_FIXED16_16(forward.y * speed);
        camera->target.z -= FLOAT_TO_FIXED16_16(forward.z * speed);
    }
    // Move right
    if(glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS) {
        camera->position.x += FLOAT_TO_FIXED16_16(right.x * speed);
        camera->position.y += FLOAT_TO_FIXED16_16(right.y * speed);
        camera->position.z += FLOAT_TO_FIXED16_16(right.z * speed);
        camera->target.x += FLOAT_TO_FIXED16_16(right.x * speed);
        camera->target.y += FLOAT_TO_FIXED16_16(right.y * speed);
        camera->target.z += FLOAT_TO_FIXED16_16(right.z * speed);
    }
    // Move left
    if(glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS) {
        camera->position.x -= FLOAT_TO_FIXED16_16(right.x * speed);
        camera->position.y -= FLOAT_TO_FIXED16_16(right.y * speed);
        camera->position.z -= FLOAT_TO_FIXED16_16(right.z * speed);
        camera->target.x -= FLOAT_TO_FIXED16_16(right.x * speed);
        camera->target.y -= FLOAT_TO_FIXED16_16(right.y * speed);
        camera->target.z -= FLOAT_TO_FIXED16_16(right.z * speed);
    }
    // Move up
    if(glfwGetKey(window, GLFW_KEY_Q) == GLFW_PRESS) {
        camera->position.y += FLOAT_TO_FIXED16_16(speed);
        camera->target.y += FLOAT_TO_FIXED16_16(speed);
    }
    // Move down
    if(glfwGetKey(window, GLFW_KEY_E) == GLFW_PRESS) {
        camera->position.y -= FLOAT_TO_FIXED16_16(speed);
        camera->target.y -= FLOAT_TO_FIXED16_16(speed);
    }
}

int main(void) {
    // Initialize GLFW
    if (!glfwInit()) {
        fprintf(stderr, "Error: Failed to initialize GLFW\n");
        return EXIT_FAILURE;
    }

    // Create GLFW window
    GLFWwindow* window = glfwCreateWindow(WIDTH, HEIGHT, "Advanced Multithreaded Renderer", NULL, NULL);
    if (!window) {
        fprintf(stderr, "Error: Failed to create GLFW window\n");
        glfwTerminate();
        return EXIT_FAILURE;
    }

    // Make the OpenGL context current
    glfwMakeContextCurrent(window);

    // Initialize OpenGL settings
    glViewport(0, 0, WIDTH, HEIGHT);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(-1, 1, -1, 1, -1, 1);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glDisable(GL_DEPTH_TEST);
    glDisable(GL_LIGHTING);
    glDisable(GL_CULL_FACE);

    // Initialize textures
    init_checkerboard_texture(&textures[0]);
    // Load additional textures if needed
    /*
    if(load_texture_from_bmp("path_to_texture.bmp", &textures[1]) != 0) {
        fprintf(stderr, "Error: Failed to load BMP texture\n");
        glfwDestroyWindow(window);
        glfwTerminate();
        return EXIT_FAILURE;
    }
    */

    // Initialize triangles
    Triangle triangles_list[MAX_TRIANGLES];
    int num_triangles = 0;
    // Load an ASCII STL model
    if(load_ascii_stl("model.stl", triangles_list, &num_triangles) != 0) {
        fprintf(stderr, "Error: Failed to load STL model\n");
        glfwDestroyWindow(window);
        glfwTerminate();
        return EXIT_FAILURE;
    }
    printf("Loaded %d triangles from model.stl\n", num_triangles);

    // Build BSP tree
    Timer bsp_timer;
    start_timer(&bsp_timer);
    BSPNode* bsp_root = build_bsp_tree(triangles_list, num_triangles, 0);
    end_timer(&bsp_timer);
    printf("BSP Tree Construction Time: %.6f seconds\n", get_elapsed_time(&bsp_timer));

    // Initialize lights
    Light lights[MAX_LIGHTS];
    int num_lights = 0;
    // Add a directional light
    Vec3 dir_light_dir = { FLOAT_TO_FIXED16_16(-1.0f), FLOAT_TO_FIXED16_16(-1.0f), FLOAT_TO_FIXED16_16(-1.0f) };
    normalize(&dir_light_dir);
    lights[num_lights].type = DIRECTIONAL_LIGHT;
    lights[num_lights].direction = dir_light_dir;
    lights[num_lights].color = (Color){255, 255, 255, 255};
    lights[num_lights].intensity = 1.0f;
    lights[num_lights].specular_strength = 0.5f;
    lights[num_lights].cutoff_angle = 0.0f; // Not used for directional lights
    num_lights++;

    // Add a point light
    Vec3 point_light_pos = { FLOAT_TO_FIXED16_16(2.0f), FLOAT_TO_FIXED16_16(2.0f), FLOAT_TO_FIXED16_16(2.0f) };
    lights[num_lights].type = POINT_LIGHT;
    lights[num_lights].position = point_light_pos;
    lights[num_lights].color = (Color){255, 200, 200, 255};
    lights[num_lights].intensity = 0.8f;
    lights[num_lights].specular_strength = 0.3f;
    lights[num_lights].cutoff_angle = 0.0f; // Not used for point lights
    num_lights++;

    // Add a spotlight
    Vec3 spotlight_pos = { FLOAT_TO_FIXED16_16(-2.0f), FLOAT_TO_FIXED16_16(2.0f), FLOAT_TO_FIXED16_16(2.0f) };
    Vec3 spotlight_dir = { FLOAT_TO_FIXED16_16(1.0f), FLOAT_TO_FIXED16_16(-1.0f), FLOAT_TO_FIXED16_16(-1.0f) };
    normalize(&spotlight_dir);
    lights[num_lights].type = SPOTLIGHT;
    lights[num_lights].position = spotlight_pos;
    lights[num_lights].direction = spotlight_dir;
    lights[num_lights].color = (Color){200, 255, 200, 255};
    lights[num_lights].intensity = 0.9f;
    lights[num_lights].specular_strength = 0.4f;
    lights[num_lights].cutoff_angle = 30.0f;
    num_lights++;

    // Initialize camera with Y+ into the screen, Z+ vertical, X+ right
    Camera camera = {
        .position = { FLOAT_TO_FIXED16_16(0.0f), FLOAT_TO_FIXED16_16(0.0f), FLOAT_TO_FIXED16_16(5.0f) },
        .target = { FLOAT_TO_FIXED16_16(0.0f), FLOAT_TO_FIXED16_16(0.0f), FLOAT_TO_FIXED16_16(0.0f) },
        .up = { FLOAT_TO_FIXED16_16(0.0f), FLOAT_TO_FIXED16_16(1.0f), FLOAT_TO_FIXED16_16(0.0f) }
    };

    // Create projection matrix
    Mat4 projection = create_perspective_matrix(90.0f, (float)WIDTH / HEIGHT, 0.1f, 100.0f);

    // Generate shadow maps in parallel
    Timer shadow_timer;
    start_timer(&shadow_timer);
    generate_shadow_maps_parallel(bsp_root, lights, num_lights, shadow_maps, MAX_TEXTURE_THREADS);
    end_timer(&shadow_timer);
    printf("Shadow Map Generation Time: %.6f seconds\n", get_elapsed_time(&shadow_timer));
    shadow_map_count = num_lights < MAX_SHADOW_MAPS ? num_lights : MAX_SHADOW_MAPS;

    // Initialize fill rate tracker
    init_fillrate_tracker();

    // Initialize texture renderer
    if(init_texture_renderer() != 0) {
        fprintf(stderr, "Error: Failed to initialize texture renderer\n");
        release_resource((RefCounted*)bsp_root);
        glfwDestroyWindow(window);
        glfwTerminate();
        return EXIT_FAILURE;
    }

    // Initialize font
    CharBitmap font[95];
    init_font(font);

    // Variables for timing
    double last_time = glfwGetTime();
    float delta_time = 0.0f;
    int frame_count = 0;
    Timer frame_timer;
    Timer perf_timer;
    double total_pixels = 0.0f;
    start_timer(&perf_timer);

    // Main loop
    while (!glfwWindowShouldClose(window)) {
        // Calculate delta time
        double current_time = glfwGetTime();
        delta_time = (float)(current_time - last_time);
        last_time = current_time;
        frame_count++;

        // Handle user input
        handle_input(window, &camera, delta_time);

        // Create view matrix
        Mat4 view = create_look_at_matrix(&camera.position, &camera.target, &camera.up);

        // Clear buffers
        clear_buffers();

        // Start total frame timer
        start_timer(&frame_timer);

        // Collect visible triangles with frustum culling
        Triangle* visible_triangles[MAX_TRIANGLES];
        int visible_count = collect_visible_triangles(bsp_root, &camera, visible_triangles);
        printf("Visible Triangles after Frustum Culling: %d\n", visible_count);

        // Apply backface culling and collect front-facing triangles
        Triangle* front_facing_triangles[MAX_TRIANGLES];
        int front_facing_count = 0;
        for(int i = 0; i < visible_count; i++) {
            // Compute view direction
            Vec3 frag_dir = {
                FIXED16_16_TO_FLOAT(camera.target.x) - FIXED16_16_TO_FLOAT(camera.position.x),
                FIXED16_16_TO_FLOAT(camera.target.y) - FIXED16_16_TO_FLOAT(camera.position.y),
                FIXED16_16_TO_FLOAT(camera.target.z) - FIXED16_16_TO_FLOAT(camera.position.z)
            };
            normalize(&frag_dir);

            float dp = dot_product(&visible_triangles[i]->normal, &frag_dir);
            if(dp > 0.0f) { // Front-facing
                front_facing_triangles[front_facing_count++] = visible_triangles[i];
            }
        }
        printf("Front-Facing Triangles after Backface Culling: %d\n", front_facing_count);

        // Start threading timer
        Timer threading_timer;
        start_timer(&threading_timer);

        // Setup threading tasks for rasterization
        pthread_t threads[MAX_THREADS];
        ThreadTask tasks[MAX_THREADS];
        int triangles_per_thread = front_facing_count / MAX_THREADS;
        int remaining_triangles = front_facing_count % MAX_THREADS;

        for(int i = 0; i < MAX_THREADS; i++) {
            tasks[i].thread_id = i;
            tasks[i].start_idx = i * triangles_per_thread;
            tasks[i].end_idx = (i +1) * triangles_per_thread;
            if(i == MAX_THREADS -1) {
                tasks[i].end_idx += remaining_triangles; // Last thread handles remaining triangles
            }
            tasks[i].triangles = front_facing_triangles;
            tasks[i].lights = lights;
            tasks[i].num_lights = num_lights;
            tasks[i].textures = textures;
            tasks[i].num_textures = 1; // Assuming one texture for simplicity
            tasks[i].shadow_maps = shadow_maps;
            tasks[i].num_shadow_maps = shadow_map_count;

            if(pthread_create(&threads[i], NULL, thread_rasterize, &tasks[i]) != 0) {
                fprintf(stderr, "Error: Failed to create thread %d\n", i);
                // Cleanup and exit
                shutdown_texture_renderer();
                release_resource((RefCounted*)bsp_root);
                glfwDestroyWindow(window);
                glfwTerminate();
                return EXIT_FAILURE;
            }
        }

        // Join threads
        for(int i = 0; i < MAX_THREADS; i++) {
            pthread_join(threads[i], NULL);
        }

        end_timer(&threading_timer);
        printf("Multithreaded Rasterization Time: %.6f seconds\n", get_elapsed_time(&threading_timer));

        // Downsample framebuffer for anti-aliasing
        Timer downsample_timer;
        start_timer(&downsample_timer);
        Color downsampled_framebuffer[WIDTH * HEIGHT] __attribute__((aligned(16)));
        downsample_framebuffer(downsampled_framebuffer, framebuffer);
        end_timer(&downsample_timer);
        printf("Downsampling Time: %.6f seconds\n", get_elapsed_time(&downsample_timer));

        // Create texture from downsampled framebuffer
        GLuint texture = create_texture(downsampled_framebuffer);

        // Render the texture to a fullscreen quad
        Timer render_timer;
        start_timer(&render_timer);
        render_fullscreen_quad(texture);
        end_timer(&render_timer);
        printf("Rendering to Screen Time: %.6f seconds\n", get_elapsed_time(&render_timer));

        // Render text overlay
        // Calculate FPS
        static int frame_counter = 0;
        static double fps_timer = 0.0;
        static double last_fps_time = 0.0;
        frame_counter++;
        double current_fps_time = glfwGetTime();
        if(current_fps_time - last_fps_time >= 1.0) {
            fps_timer = frame_counter / (current_fps_time - last_fps_time);
            frame_counter = 0;
            last_fps_time = current_fps_time;
        }
        char fps_text[32];
        snprintf(fps_text, sizeof(fps_text), "FPS: %.2f", fps_timer);
        Color text_color = {255, 255, 255, 255};
        render_text(fps_text, 10, 10, &text_color, font);

        // Enqueue texture rendering task (example: invert the checkerboard texture every second)
        if((int)current_fps_time % 2 == 0) { // Every 2 seconds
            enqueue_texture_task(example_texture_task, &textures[0]);
        }

        // Swap buffers and poll events
        glfwSwapBuffers(window);
        glfwPollEvents();

        // Delete the texture
        glDeleteTextures(1, &texture);

        // End total frame timer and log
        end_timer(&frame_timer);
        printf("Total Frame Time: %.6f seconds\n\n", get_elapsed_time(&frame_timer));

        // Advanced profiling: Calculate megapixels per second and fill rates
        total_pixels += (double)(WIDTH * HEIGHT);
        if(get_elapsed_time(&perf_timer) >= 1.0f) { // Every second
            double elapsed = get_elapsed_time(&perf_timer);
            double mpixels = (total_pixels / 1e6) / elapsed;

            // Get and reset fill rates
            FillRate fr = get_and_reset_fillrate();

            printf("Performance: %.2f Megapixels/sec | Pixels Filled/sec: %lu | Textures Filled/sec: %lu\n",
                   mpixels, fr.pixels_filled, fr.textures_filled);

            total_pixels = 0.0f;
            start_timer(&perf_timer);
        }
    }

    // Cleanup and exit
    shutdown_texture_renderer();
    release_resource((RefCounted*)bsp_root);
    glfwDestroyWindow(window);
    glfwTerminate();
    return EXIT_SUCCESS;
}
