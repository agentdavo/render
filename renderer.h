// renderer.h
#ifndef RENDERER_H
#define RENDERER_H

#include <stdint.h>
#include <pthread.h>
#include <stdatomic.h>
#include <float.h>

#include "config.h"
#include "maths.h"
#include "texture.h"
#include "fillrate.h"

// Fixed-point conversion macros
#define FLOAT_TO_FIXED16_16(x) ((int32_t)((x) * 65536.0f))
#define FIXED16_16_TO_FLOAT(x) ((float)(x) / 65536.0f)
#define FLOAT_TO_FIXED32_32(x) ((int64_t)((x) * 4294967296.0))
#define FIXED32_32_TO_FLOAT(x) ((float)(x) / 4294967296.0f)

// Clamp macro
#define clamp(val, min, max) ((val) < (min) ? (min) : ((val) > (max) ? (max) : (val)))

// Resource Types for Reference Counting
typedef enum {
    RESOURCE_BSP_NODE,
    RESOURCE_TEXTURE,
    // Add other resource types as needed
} ResourceType;

// Basic Vector Structure
typedef struct {
    float x;
    float y;
    float z;
} Vec3;

// Basic Matrix Structure (4x4)
typedef struct {
    float m[4][4];
} Mat4;

// Camera Structure
typedef struct {
    Vec3 position;
    Vec3 target;
    Vec3 up;
} Camera;

// Color Structure
typedef struct {
    uint8_t r;
    uint8_t g;
    uint8_t b;
    uint8_t a;
} Color;

// Material Structure
typedef struct {
    Color diffuse;
    Color specular;
    Color emissive;
    float shininess;
} Material;

// Light Types
typedef enum {
    DIRECTIONAL_LIGHT,
    POINT_LIGHT,
    SPOTLIGHT
} LightType;

// Light Structure
typedef struct {
    LightType type;
    Vec3 position;      // For point lights and spotlights
    Vec3 direction;     // For directional lights and spotlights
    float intensity;
    Color color;
    float specular_strength;  // Add this field for specular reflection strength
    float cutoff_angle;       // Add this field for spotlights
} Light;

// Reference-counted Resource Base Structure
typedef struct {
    atomic_int ref_count;
    ResourceType type;
} RefCounted;

// Texture structure with reference counting
typedef struct {
    RefCounted base;
    uint8_t data[TEXTURE_SIZE][TEXTURE_SIZE][3];
} Texture;

// ShadowMap structure
typedef struct {
    float depth[WIDTH * HEIGHT]; // Depth buffer
    Mat4 view_projection;        // Combined view and projection matrix for the light
} ShadowMap;

// Triangle structure
typedef struct {
    Vec3 v0;
    Vec3 v1;
    Vec3 v2;
    Vec3 normal;
    float u0, v0_tex;
    float u1, v1_tex;
    float u2, v2_tex;
    Material material;
} Triangle;

// BSPNode structure with reference counting
typedef struct BSPNode BSPNode;
struct BSPNode {
    Vec3 normal;
    float d;
    BSPNode* front;
    BSPNode* back;
    Triangle* triangles;
    int num_triangles;
    RefCounted base; // For reference counting
};

// Font character bitmap structure
typedef struct {
    char character;
    uint8_t bitmap[8];
} CharBitmap;

// Shadow Mapping Task Structure
typedef struct {
    ShadowMap* shadow_maps;
    const Light* lights;
    int num_lights;
    const BSPNode* bsp_root;
    int thread_id;
    int total_threads;
} ShadowMapTask;

// Function declarations

// Reference Counting Functions
RefCounted* retain_resource(RefCounted* resource);
void release_resource(RefCounted* resource);

// Model loading
int load_ascii_stl(const char* filename, Triangle* triangles, int* num_triangles);

// BSP tree functions
BSPNode* build_bsp_tree(Triangle* triangles, int num_triangles, int depth);
void traverse_bsp_tree(const BSPNode* node, const Camera* camera, Triangle** visible_tris, int* count);
void free_bsp_tree_recursively(BSPNode* node);

// Texture functions
int load_texture_from_bmp(const char* filename, Texture* tex);
void init_checkerboard_texture(Texture* tex);

// Lighting functions
Color compute_lighting(const Vec3* normal, const Vec3* view_dir, const Light* lights, int num_lights, const Material* material, const Vec3* frag_pos);

// Clipping functions
int clip_triangle_all_planes(const Triangle* tri, Triangle* clipped_tris);

// Shadow mapping functions
void* generate_shadow_map_thread(void* arg);
void generate_shadow_maps_parallel(const BSPNode* bsp_root, const Light* lights, int num_lights, ShadowMap* shadow_maps, int num_threads);

// Shadow testing
int is_in_shadow(const Vec3* frag_pos, const Vec3* normal, const Light* light, const ShadowMap* shadow_map);

// Rasterization functions
void rasterize_triangle(const Triangle* tri, const Mat4* view, const Mat4* projection, const Light* lights, int num_lights, Texture* textures, int num_textures, ShadowMap* shadow_maps, int num_shadow_maps);

// Fog functions
Color apply_fog(const Color* color, float distance, float fog_density, const Color* fog_color);

// Text rendering functions
void init_font(CharBitmap* font);
void render_text(const char* text, int x, int y, const Color* color, const CharBitmap* font);

// Fill rate tracking functions (Assumed to be defined elsewhere)
void increment_pixels_filled();
void increment_textures_filled();

// Mathematical utility functions
Vec3 add_vec3(const Vec3* a, const Vec3* b);
Vec3 subtract_vec3(const Vec3* a, const Vec3* b);
Vec3 multiply_vec3(const Vec3* a, float scalar);
float dot_product(const Vec3* a, const Vec3* b);
void normalize(Vec3* v);
Mat4 create_look_at_matrix(const Vec3* eye, const Vec3* target, const Vec3* up);
Mat4 create_perspective_matrix(float fov, float aspect, float near_plane, float far_plane);
Mat4 multiply_matrices(const Mat4* a, const Mat4* b);
Vec3 transform_vector(const Mat4* mat, const Vec3* vec);

// Framebuffer and zbuffer structures
typedef struct {
    Color* framebuffer;
    float* zbuffer;
} ThreadBuffers;

// Initialize and cleanup functions for buffers
void initialize_buffers(ThreadBuffers buffers[MAX_THREADS]);
void cleanup_buffers(ThreadBuffers buffers[MAX_THREADS]);

// Function to combine thread buffers into a single framebuffer
void combine_framebuffers(ThreadBuffers buffers[MAX_THREADS], Color* final_framebuffer);

#endif // RENDERER_H
