// texture_renderer.c
#include <stdlib.h>
#include <stdio.h>
#include <pthread.h>

#include "texture.h"

// Task queue
static TextureTask* task_queue_head = NULL;
static TextureTask* task_queue_tail = NULL;

// Mutex and condition variable for task queue synchronization
static pthread_mutex_t queue_mutex = PTHREAD_MUTEX_INITIALIZER;
static pthread_cond_t queue_cond = PTHREAD_COND_INITIALIZER;

// Flag to indicate if the renderer is shutting down
static bool shutdown_flag = false;

// Texture rendering threads
static pthread_t texture_threads[MAX_TEXTURE_THREADS];

// Worker function for texture rendering threads
static void* texture_worker(void* arg) {
    while (1) {
        pthread_mutex_lock(&queue_mutex);

        // Wait for tasks to be available or shutdown signal
        while (task_queue_head == NULL && !shutdown_flag) {
            pthread_cond_wait(&queue_cond, &queue_mutex);
        }

        // If shutting down and no tasks, exit
        if (shutdown_flag && task_queue_head == NULL) {
            pthread_mutex_unlock(&queue_mutex);
            break;
        }

        // Dequeue the next task
        TextureTask* task = task_queue_head;
        if (task_queue_head != NULL) {
            task_queue_head = task_queue_head->next;
            if (task_queue_head == NULL) {
                task_queue_tail = NULL;
            }
        }
        pthread_mutex_unlock(&queue_mutex);

        // Execute the task
        if (task != NULL) {
            task->task_function(task->task_arg);
            free(task);
        }
    }

    return NULL;
}

// Initialize the texture renderer thread pool
int init_texture_renderer() {
    shutdown_flag = false;
    task_queue_head = task_queue_tail = NULL;

    for (int i = 0; i < MAX_TEXTURE_THREADS; i++) {
        if (pthread_create(&texture_threads[i], NULL, texture_worker, NULL) != 0) {
            fprintf(stderr, "Error: Failed to create texture rendering thread %d\n", i);
            return -1;
        }
    }

    return 0;
}

// Enqueue a texture rendering task
int enqueue_texture_task(void (*task_function)(void*), void* task_arg) {
    TextureTask* new_task = (TextureTask*)malloc(sizeof(TextureTask));
    if (!new_task) {
        fprintf(stderr, "Error: Failed to allocate memory for texture task\n");
        return -1;
    }
    new_task->task_function = task_function;
    new_task->task_arg = task_arg;
    new_task->next = NULL;

    pthread_mutex_lock(&queue_mutex);
    if (task_queue_tail == NULL) {
        task_queue_head = task_queue_tail = new_task;
    } else {
        task_queue_tail->next = new_task;
        task_queue_tail = new_task;
    }
    pthread_cond_signal(&queue_cond);
    pthread_mutex_unlock(&queue_mutex);

    return 0;
}

// Shutdown the texture renderer thread pool
void shutdown_texture_renderer() {
    pthread_mutex_lock(&queue_mutex);
    shutdown_flag = true;
    pthread_cond_broadcast(&queue_cond);
    pthread_mutex_unlock(&queue_mutex);

    // Join all texture rendering threads
    for (int i = 0; i < MAX_TEXTURE_THREADS; i++) {
        pthread_join(texture_threads[i], NULL);
    }

    // Cleanup remaining tasks
    while (task_queue_head != NULL) {
        TextureTask* temp = task_queue_head;
        task_queue_head = task_queue_head->next;
        free(temp);
    }

    task_queue_tail = NULL;
}
