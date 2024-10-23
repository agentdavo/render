# Minimum CMake version requirement
cmake_minimum_required(VERSION 3.10)

# Project declaration
project(render C)

# Set the C standard to GNU99
set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS ON)

# =====================================================================
#                             Configuration
# =====================================================================

# =====================================================================
#                             Dependencies
# =====================================================================

find_package(OpenGL REQUIRED)
find_package(PkgConfig REQUIRED)
pkg_search_module(GLFW REQUIRED glfw3)
find_package(Threads REQUIRED)

include(FetchContent)
FetchContent_Declare(
    libfixmath
    GIT_REPOSITORY https://github.com/PetteriAimonen/libfixmath.git
    GIT_TAG d308e466e1a09118d03f677c52e5fbf402f6fdd0
)
FetchContent_MakeAvailable(libfixmath)

target_compile_definitions(libfixmath PRIVATE
    # FIXMATH_FAST_SIN
    # FIXMATH_NO_64BIT
    # FIXMATH_NO_CACHE
    # FIXMATH_NO_HARD_DIVISION
    # FIXMATH_NO_OVERFLOW
    # FIXMATH_NO_ROUNDING
    # FIXMATH_OPTIMIZE_8BIT
)

# =====================================================================
#                             Include Directories
# =====================================================================

include_directories(
    ${GLFW_INCLUDE_DIRS}
    ${PROJECT_SOURCE_DIR}  # To include headers in the root directory
)

# =====================================================================
#                             Source Files
# =====================================================================

# List all source files
set(SOURCES
    advanced_renderer.c
    renderer.c
    maths.c
    fillrate.c
    texture.c
)

# =====================================================================
#                             Executable
# =====================================================================

# Define the executable target
add_executable(${PROJECT_NAME} ${SOURCES})

# =====================================================================
#                             Compiler Definitions
# =====================================================================

# =====================================================================
#                             Linking Libraries
# =====================================================================

# Link against OpenGL, GLFW, Threads, and math library
target_link_libraries(${PROJECT_NAME}
    ${OPENGL_LIBRARIES}    # Links against OpenGL
    ${GLFW_LIBRARIES}      # Links against GLFW
    Threads::Threads       # Links against the Threads library (pthreads on Unix)
    m                      # Links against the math library (libm)
    libfixmath             # Links against the fixed math library (libfixmath)
)

# =====================================================================
#                             Platform-Specific Settings
# =====================================================================

# Linux-specific: Link with 'dl' for dynamic loading if necessary
if(UNIX AND NOT APPLE)
    target_link_libraries(${PROJECT_NAME} dl)
endif()

# Windows-specific: Link against 'opengl32' and GLFW libraries
if(WIN32)
    target_link_libraries(${PROJECT_NAME} opengl32)
    # If using precompiled GLFW binaries on Windows, specify the path or ensure they're in the library search path
endif()

# macOS-specific: Link against Cocoa, IOKit, and CoreVideo frameworks for GLFW
if(APPLE)
    find_library(COCOA_LIBRARY Cocoa REQUIRED)
    find_library(IOKIT_LIBRARY IOKit REQUIRED)
    find_library(COREVIDEO_LIBRARY CoreVideo REQUIRED)
    target_link_libraries(${PROJECT_NAME} ${COCOA_LIBRARY} ${IOKIT_LIBRARY} ${COREVIDEO_LIBRARY})
endif()

# =====================================================================
#                             RPATH Configuration
# =====================================================================

# Set RPATH for Unix-like systems to find shared libraries at runtime
if(UNIX AND NOT APPLE)
    set_target_properties(${PROJECT_NAME} PROPERTIES BUILD_WITH_INSTALL_RPATH TRUE)
    set_target_properties(${PROJECT_NAME} PROPERTIES INSTALL_RPATH "$ORIGIN")
endif()

# =====================================================================
#                             Installation
# =====================================================================
