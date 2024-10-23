// fillrate.h
#ifndef FILLRATE_H
#define FILLRATE_H

#include <stdint.h>

#include "config.h"

// Structure to hold fill rate information
typedef struct {
    unsigned long pixels_filled;
    unsigned long textures_filled;
} FillRate;

// Initialize the fill rate tracker
void init_fillrate_tracker();

// Increment the pixel fill counter
void increment_pixels_filled();

// Increment the texture fill counter
void increment_textures_filled();

// Get the current fill rate (pixels/sec, textures/sec) and reset the counters
FillRate get_and_reset_fillrate();

#endif // FILLRATE_H
