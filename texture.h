// texture_renderer.h
#ifndef TEXTURE_H
#define TEXTURE_H

#include "config.h"

#include <pthread.h>
#include <stdbool.h>

// Structure for a texture rendering task
typedef struct TextureTask {
    void (*task_function)(void*); // Function pointer for the task
    void* task_arg;                // Argument to the task function
    struct TextureTask* next;      // Pointer to the next task in the queue
} TextureTask;

// Initialize the texture renderer thread pool
int init_texture_renderer();

// Enqueue a texture rendering task
int enqueue_texture_task(void (*task_function)(void*), void* task_arg);

// Shutdown the texture renderer thread pool
void shutdown_texture_renderer();

#endif // TEXTURE_RENDERER_H
