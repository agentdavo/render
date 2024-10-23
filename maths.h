#ifndef MATHS_H
#define MATHS_H

#include <stdint.h>
#include <math.h>

#include "libfixmath/fix16.h"

// Fixed-point types
typedef int32_t fixed16_16;
typedef int64_t fixed32_32;

// 3D vector and 4D vector structures (aligned to 16 bytes for SIMD)
typedef struct Vec3 {
    fixed16_16 x, y, z;
} Vec3 __attribute__((aligned(16)));

typedef struct Vec4 {
    fixed16_16 x, y, z, w;
} Vec4 __attribute__((aligned(16)));

// 4x4 matrix for transformation and projection operations
typedef struct Mat4 {
    fixed16_16 m[4][4];
} Mat4 __attribute__((aligned(16)));

// Fixed-point conversion macros
#define INT_TO_FIXED16_16(x) ((fixed16_16)((x) << 16))
#define FLOAT_TO_FIXED16_16(x) ((fixed16_16)((x) * 65536.0f))
#define FIXED16_16_TO_FLOAT(x) ((float)(x) / 65536.0f)

#define INT_TO_FIXED32_32(x) ((fixed32_32)((x) << 32))
#define FLOAT_TO_FIXED32_32(x) ((fixed32_32)((x) * 4294967296.0f))
#define FIXED32_32_TO_FLOAT(x) ((float)(x) / 4294967296.0f)

// Vec3 operations
void normalize(Vec3* v);
fixed16_16 dot_product_vec3(const Vec3* a, const Vec3* b);
Vec3 cross_product(const Vec3* a, const Vec3* b);
Vec3 scale_vector(const Vec3* v, fixed16_16 scale);
Vec3 add_vectors(const Vec3* a, const Vec3* b);
Vec3 subtract_vectors(const Vec3* a, const Vec3* b);
fixed16_16 length_vec3(const Vec3* v);
fixed16_16 squared_length_vec3(const Vec3* v);

// Vec4 operations (useful for colors, lighting, and transformation)
Vec4 add_vec4(const Vec4* a, const Vec4* b);
Vec4 subtract_vec4(const Vec4* a, const Vec4* b);
Vec4 scale_vec4(const Vec4* v, fixed16_16 scale);
Vec4 normalize_vec4(Vec4* v);
fixed16_16 dot_product_vec4(const Vec4* a, const Vec4* b);
Vec4 lerp_color(const Vec4* c1, const Vec4* c2, fixed16_16 t);
Vec4 premultiply_alpha(const Vec4* color);
Vec4 blend_colors(const Vec4* c1, const Vec4* c2);
Vec4 transform_vector4(const Mat4* mat, const Vec4* v);

// Matrix operations
Mat4 multiply_matrices(const Mat4* a, const Mat4* b);
Vec3 transform_vector(const Mat4* mat, const Vec3* v);
Mat4 create_perspective_matrix(fixed16_16 fov, fixed16_16 aspect, fixed16_16 near, fixed16_16 far);
Mat4 create_look_at_matrix(const Vec3* position, const Vec3* target, const Vec3* up);
Mat4 create_identity_matrix();
Mat4 create_translation_matrix(const Vec3* translation);
Mat4 create_rotation_matrix_x(fixed16_16 angle);
Mat4 create_rotation_matrix_y(fixed16_16 angle);
Mat4 create_rotation_matrix_z(fixed16_16 angle);
Mat4 transpose_matrix(const Mat4* mat);

#endif // MATHS_H
