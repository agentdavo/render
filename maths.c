#include "maths.h"

#include <math.h>

// Vec3 operations
void normalize(Vec3* v) {
    fixed32_32 length_sq = fix16_mul(v->x, v->x) + fix16_mul(v->y, v->y) + fix16_mul(v->z, v->z);
    if (length_sq == 0) return;

    fixed16_16 inv_length = fix16_div(fix16_one, fix16_sqrt(length_sq));
    v->x = fix16_mul(v->x, inv_length);
    v->y = fix16_mul(v->y, inv_length);
    v->z = fix16_mul(v->z, inv_length);
}

fixed16_16 dot_product_vec3(const Vec3* a, const Vec3* b) {
    return fix16_add(fix16_add(fix16_mul(a->x, b->x), fix16_mul(a->y, b->y)), fix16_mul(a->z, b->z));
}

Vec3 cross_product(const Vec3* a, const Vec3* b) {
    Vec3 result;
    result.x = fix16_sub(fix16_mul(a->y, b->z), fix16_mul(a->z, b->y));
    result.y = fix16_sub(fix16_mul(a->z, b->x), fix16_mul(a->x, b->z));
    result.z = fix16_sub(fix16_mul(a->x, b->y), fix16_mul(a->y, b->x));
    return result;
}

Vec3 scale_vector(const Vec3* v, fixed16_16 scale) {
    Vec3 result;
    result.x = fix16_mul(v->x, scale);
    result.y = fix16_mul(v->y, scale);
    result.z = fix16_mul(v->z, scale);
    return result;
}

Vec3 add_vectors(const Vec3* a, const Vec3* b) {
    Vec3 result;
    result.x = fix16_add(a->x, b->x);
    result.y = fix16_add(a->y, b->y);
    result.z = fix16_add(a->z, b->z);
    return result;
}

Vec3 subtract_vectors(const Vec3* a, const Vec3* b) {
    Vec3 result;
    result.x = fix16_sub(a->x, b->x);
    result.y = fix16_sub(a->y, b->y);
    result.z = fix16_sub(a->z, b->z);
    return result;
}

fixed16_16 length_vec3(const Vec3* v) {
    return fix16_sqrt(fix16_mul(v->x, v->x) + fix16_mul(v->y, v->y) + fix16_mul(v->z, v->z));
}

fixed16_16 squared_length_vec3(const Vec3* v) {
    return fix16_mul(v->x, v->x) + fix16_mul(v->y, v->y) + fix16_mul(v->z, v->z);
}

// Vec4 operations
Vec4 add_vec4(const Vec4* a, const Vec4* b) {
    Vec4 result;
    result.x = fix16_add(a->x, b->x);
    result.y = fix16_add(a->y, b->y);
    result.z = fix16_add(a->z, b->z);
    result.w = fix16_add(a->w, b->w);
    return result;
}

Vec4 subtract_vec4(const Vec4* a, const Vec4* b) {
    Vec4 result;
    result.x = fix16_sub(a->x, b->x);
    result.y = fix16_sub(a->y, b->y);
    result.z = fix16_sub(a->z, b->z);
    result.w = fix16_sub(a->w, b->w);
    return result;
}

Vec4 scale_vec4(const Vec4* v, fixed16_16 scale) {
    Vec4 result;
    result.x = fix16_mul(v->x, scale);
    result.y = fix16_mul(v->y, scale);
    result.z = fix16_mul(v->z, scale);
    result.w = fix16_mul(v->w, scale);
    return result;
}

Vec4 normalize_vec4(Vec4* v) {
    fixed16_16 length = fix16_sqrt(fix16_add(fix16_add(fix16_mul(v->x, v->x), fix16_mul(v->y, v->y)), fix16_mul(v->z, v->z)));
    Vec4 result = {fix16_div(v->x, length), fix16_div(v->y, length), fix16_div(v->z, length), v->w};
    return result;
}

fixed16_16 dot_product_vec4(const Vec4* a, const Vec4* b) {
    return fix16_add(fix16_add(fix16_add(fix16_mul(a->x, b->x), fix16_mul(a->y, b->y)), fix16_mul(a->z, b->z)), fix16_mul(a->w, b->w));
}

Vec4 lerp_color(const Vec4* c1, const Vec4* c2, fixed16_16 t) {
    Vec4 result;
    result.x = fix16_add(fix16_mul(fix16_sub(c2->x, c1->x), t), c1->x);
    result.y = fix16_add(fix16_mul(fix16_sub(c2->y, c1->y), t), c1->y);
    result.z = fix16_add(fix16_mul(fix16_sub(c2->z, c1->z), t), c1->z);
    result.w = fix16_add(fix16_mul(fix16_sub(c2->w, c1->w), t), c1->w);
    return result;
}

Vec4 premultiply_alpha(const Vec4* color) {
    Vec4 result;
    result.x = fix16_mul(color->x, color->w);
    result.y = fix16_mul(color->y, color->w);
    result.z = fix16_mul(color->z, color->w);
    result.w = color->w;
    return result;
}

Vec4 blend_colors(const Vec4* c1, const Vec4* c2) {
    fixed16_16 one_minus_alpha1 = fix16_sub(fix16_one, c1->w);
    Vec4 result;
    result.x = fix16_add(fix16_mul(c1->x, c1->w), fix16_mul(c2->x, one_minus_alpha1));
    result.y = fix16_add(fix16_mul(c1->y, c1->w), fix16_mul(c2->y, one_minus_alpha1));
    result.z = fix16_add(fix16_mul(c1->z, c1->w), fix16_mul(c2->z, one_minus_alpha1));
    result.w = fix16_add(fix16_mul(c1->w, c1->w), fix16_mul(c2->w, one_minus_alpha1));
    return result;
}

Vec4 transform_vector4(const Mat4* mat, const Vec4* v) {
    Vec4 result;
    result.x = fix16_add(fix16_add(fix16_add(fix16_mul(mat->m[0][0], v->x), fix16_mul(mat->m[0][1], v->y)), fix16_mul(mat->m[0][2], v->z)), fix16_mul(mat->m[0][3], v->w));
    result.y = fix16_add(fix16_add(fix16_add(fix16_mul(mat->m[1][0], v->x), fix16_mul(mat->m[1][1], v->y)), fix16_mul(mat->m[1][2], v->z)), fix16_mul(mat->m[1][3], v->w));
    result.z = fix16_add(fix16_add(fix16_add(fix16_mul(mat->m[2][0], v->x), fix16_mul(mat->m[2][1], v->y)), fix16_mul(mat->m[2][2], v->z)), fix16_mul(mat->m[2][3], v->w));
    result.w = fix16_add(fix16_add(fix16_add(fix16_mul(mat->m[3][0], v->x), fix16_mul(mat->m[3][1], v->y)), fix16_mul(mat->m[3][2], v->z)), fix16_mul(mat->m[3][3], v->w));
    return result;
}

// Matrix operations
Mat4 multiply_matrices(const Mat4* a, const Mat4* b) {
    Mat4 result;
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            result.m[i][j] = fix16_add(fix16_add(fix16_add(fix16_mul(a->m[i][0], b->m[0][j]), fix16_mul(a->m[i][1], b->m[1][j])), fix16_mul(a->m[i][2], b->m[2][j])), fix16_mul(a->m[i][3], b->m[3][j]));
        }
    }
    return result;
}

Vec3 transform_vector(const Mat4* mat, const Vec3* v) {
    Vec3 result;
    fixed16_16 w = fix16_add(fix16_add(fix16_add(fix16_mul(mat->m[3][0], v->x), fix16_mul(mat->m[3][1], v->y)), fix16_mul(mat->m[3][2], v->z)), mat->m[3][3]);
    if (w != fix16_one) {
        result.x = fix16_div(fix16_add(fix16_add(fix16_mul(mat->m[0][0], v->x), fix16_mul(mat->m[0][1], v->y)), fix16_mul(mat->m[0][2], v->z)), w);
        result.y = fix16_div(fix16_add(fix16_add(fix16_mul(mat->m[1][0], v->x), fix16_mul(mat->m[1][1], v->y)), fix16_mul(mat->m[1][2], v->z)), w);
        result.z = fix16_div(fix16_add(fix16_add(fix16_mul(mat->m[2][0], v->x), fix16_mul(mat->m[2][1], v->y)), fix16_mul(mat->m[2][2], v->z)), w);
    } else {
        result.x = fix16_add(fix16_add(fix16_mul(mat->m[0][0], v->x), fix16_mul(mat->m[0][1], v->y)), fix16_mul(mat->m[0][2], v->z));
        result.y = fix16_add(fix16_add(fix16_mul(mat->m[1][0], v->x), fix16_mul(mat->m[1][1], v->y)), fix16_mul(mat->m[1][2], v->z));
        result.z = fix16_add(fix16_add(fix16_mul(mat->m[2][0], v->x), fix16_mul(mat->m[2][1], v->y)), fix16_mul(mat->m[2][2], v->z));
    }
    return result;
}

Mat4 create_identity_matrix() {
    Mat4 mat = {0};
    mat.m[0][0] = fix16_one;
    mat.m[1][1] = fix16_one;
    mat.m[2][2] = fix16_one;
    mat.m[3][3] = fix16_one;
    return mat;
}

Mat4 create_translation_matrix(const Vec3* translation) {
    Mat4 mat = create_identity_matrix();
    mat.m[0][3] = translation->x;
    mat.m[1][3] = translation->y;
    mat.m[2][3] = translation->z;
    return mat;
}

Mat4 create_rotation_matrix_x(fixed16_16 angle) {
    Mat4 mat = create_identity_matrix();
    fixed16_16 c = fix16_cos(angle);
    fixed16_16 s = fix16_sin(angle);
    mat.m[1][1] = c;
    mat.m[1][2] = fix16_neg(s);
    mat.m[2][1] = s;
    mat.m[2][2] = c;
    return mat;
}

Mat4 create_rotation_matrix_y(fixed16_16 angle) {
    Mat4 mat = create_identity_matrix();
    fixed16_16 c = fix16_cos(angle);
    fixed16_16 s = fix16_sin(angle);
    mat.m[0][0] = c;
    mat.m[0][2] = s;
    mat.m[2][0] = fix16_neg(s);
    mat.m[2][2] = c;
    return mat;
}

Mat4 create_rotation_matrix_z(fixed16_16 angle) {
    Mat4 mat = create_identity_matrix();
    fixed16_16 c = fix16_cos(angle);
    fixed16_16 s = fix16_sin(angle);
    mat.m[0][0] = c;
    mat.m[0][1] = fix16_neg(s);
    mat.m[1][0] = s;
    mat.m[1][1] = c;
    return mat;
}

Mat4 transpose_matrix(const Mat4* mat) {
    Mat4 result;
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            result.m[i][j] = mat->m[j][i];
        }
    }
    return result;
}
