// config.h
#ifndef CONFIG_H
#define CONFIG_H

// Screen and rendering configurations
#define WIDTH 800
#define HEIGHT 600
#define SSAA_FACTOR 2 // Supersampling factor for anti-aliasing

// Texture size
#define TEXTURE_SIZE 256

// Maximum counts
#define MAX_TRIANGLES 100000
#define MAX_LIGHTS 10
#define MAX_SHADOW_MAPS 10

// Maximum threads per process
#define MAX_THREADS 1
#define MAX_TEXTURE_THREADS 1

#endif // CONFIG_H
