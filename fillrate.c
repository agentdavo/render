// fillrate.c
#include "fillrate.h"
#include <stdint.h>

// Using GCC built-in atomics for thread-safe counters
static volatile unsigned long pixels_filled = 0;
static volatile unsigned long textures_filled = 0;

// Initialize the fill rate tracker
void init_fillrate_tracker() {
    __sync_lock_test_and_set(&pixels_filled, 0);
    __sync_lock_test_and_set(&textures_filled, 0);
}

// Increment the pixel fill counter
void increment_pixels_filled() {
    __sync_fetch_and_add(&pixels_filled, 1);
}

// Increment the texture fill counter
void increment_textures_filled() {
    __sync_fetch_and_add(&textures_filled, 1);
}

// Get the current fill rate (pixels/sec, textures/sec) and reset the counters
FillRate get_and_reset_fillrate() {
    FillRate fr;
    fr.pixels_filled = __sync_lock_test_and_set(&pixels_filled, 0);
    fr.textures_filled = __sync_lock_test_and_set(&textures_filled, 0);
    return fr;
}
