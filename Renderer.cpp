//
// Created by goksu on 2/25/20.
//

#include <fstream>
#include "Scene.hpp"
#include "Renderer.hpp"

// ******* The following lines have to be added for stb_image_write to
// compile properly
#define STBI_MSC_SECURE_CRT
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define M_PI 3.14159265358979323846  /* pi */

inline float deg2rad(const float& deg) { return deg * M_PI / 180.0; }

const float EPSILON = 0.00001;



inline Vector3f castRayIterativeFunc(int i, int j, const Scene &scene, 
    float imageAspectRatio, float scale, Vector3f eye_pos)
{
    int num_of_samples = 16;
    Vector3f col = Vector3f(0.0f);
    for (int s = 0; s < num_of_samples; ++s) {
        // generate primary ray direction
        float x;
        float y;

        // generate primary ray direction
        float x_ndc = (i + get_random_float()) / scene.width;
        float y_ndc = (j + get_random_float()) / scene.height;

        //float x_screen = 2 * x_ndc - 1;
        //float y_screen = 1- 2 * y_ndc;

        float x_camera = (2 * x_ndc - 1) * imageAspectRatio * scale;
        float y_camera = (1 - 2 * y_ndc) * scale;
        // TODO: Find the x and y positions of the current pixel to get the
        // direction
        //  vector that passes through it.
        // Also, don't forget to multiply both of them with the variable
        // *scale*, and x (horizontal) variable with the *imageAspectRatio*

        // Don't forget to normalize this direction!
        Vector3f dir = Vector3f(x_camera, y_camera, -1);
        normalize(dir);


        col += scene.castRayIterative(Ray(eye_pos, dir), 0);
    }

    col = col / num_of_samples;
    return col;
}


inline Vector3f castRayRecusriveFunc(int i, int j, const Scene& scene,
    float imageAspectRatio, float scale, Vector3f eye_pos)
{
    // generate primary ray direction
    float x;
    float y;

    // generate primary ray direction
    float x_ndc = (i + get_random_float()) / scene.width;
    float y_ndc = (j + get_random_float()) / scene.height;

    //float x_screen = 2 * x_ndc - 1;
    //float y_screen = 1- 2 * y_ndc;

    float x_camera = (2 * x_ndc - 1) * imageAspectRatio * scale;
    float y_camera = (1 - 2 * y_ndc) * scale;
    // TODO: Find the x and y positions of the current pixel to get the
    // direction
    //  vector that passes through it.
    // Also, don't forget to multiply both of them with the variable
    // *scale*, and x (horizontal) variable with the *imageAspectRatio*

    // Don't forget to normalize this direction!
    Vector3f dir = Vector3f(x_camera, y_camera, -1);
    normalize(dir);


    return scene.castRay(Ray(eye_pos, dir), 0);
}

// The main render function. This where we iterate over all pixels in the image,
// generate primary rays and cast these rays into the scene. The content of the
// framebuffer is saved to a file.
void Renderer::Render(const Scene& scene)
{
    //framebuffer for output
    #define CHANNEL_NUM 3

    //framebuffer for the stb image
    uint8_t* stb_framebuffer = new uint8_t[scene.width * scene.height * CHANNEL_NUM];

    std::vector<Vector3f> framebuffer(scene.width * scene.height);

    float scale = tan(deg2rad(scene.fov * 0.5));
    float imageAspectRatio = scene.width / (float)scene.height;
    //Vector3f eye_pos(50, 30, 150);//-1, 5, 10
    Vector3f eye_pos(-1, 5, 10);
    int m = 0;
    int index = 0;
    for (uint32_t j = 0; j < scene.height; ++j) {
        for (uint32_t i = 0; i < scene.width; ++i) {
            
            /*int num_of_samples = 2;
            Vector3f col = Vector3f(0.0f);
            for (int s = 0; s < num_of_samples; ++s) {
                // generate primary ray direction
                float x;
                float y;

                // generate primary ray direction
                float x_ndc = (i + get_random_float()) / scene.width;
                float y_ndc = (j + get_random_float()) / scene.height;

                //float x_screen = 2 * x_ndc - 1;
                //float y_screen = 1- 2 * y_ndc;

                float x_camera = (2 * x_ndc - 1) * imageAspectRatio * scale;
                float y_camera = (1 - 2 * y_ndc) * scale;
                // TODO: Find the x and y positions of the current pixel to get the
                // direction
                //  vector that passes through it.
                // Also, don't forget to multiply both of them with the variable
                // *scale*, and x (horizontal) variable with the *imageAspectRatio*

                // Don't forget to normalize this direction!
                Vector3f dir = Vector3f(x_camera, y_camera, -1);
                normalize(dir);
                

                col += scene.castRayIterative(Ray(eye_pos, dir), 0);
            }
            col = col / num_of_samples;*/
            auto col = castRayIterativeFunc(i, j, scene, imageAspectRatio, scale, eye_pos);
            framebuffer[m++] = col;

            int ired = int(255.99 * col.x);
            int igreen = int(255.99 * col.y);
            int iblue = int(255.99 * col.z);

            stb_framebuffer[index++] = ired;
            stb_framebuffer[index++] = igreen;
            stb_framebuffer[index++] = iblue;
        }
        UpdateProgress(j / (float)scene.height);
    }
    UpdateProgress(1.f);

    // save framebuffer to file
    std::ofstream ofs;
    ofs.open("./out.ppm");
    ofs << "P6\n" << scene.width << " " << scene.height << "\n255\n";
    for (uint32_t i = 0; i < scene.height * scene.width; ++i) {
        auto r = (uint8_t)(255 * clamp(0, 1, framebuffer[i].x));
        auto g = (uint8_t)(255 * clamp(0, 1, framebuffer[i].y));
        auto b = (uint8_t)(255 * clamp(0, 1, framebuffer[i].z));
        ofs << r << g << b;
    }
    ofs.close();

    /* save out the image */
    stbi_write_jpg("output\\bunny_dragon_buddha.jpg", scene.width, scene.height, 3, stb_framebuffer, 100);
    /* */

    /* cleanup */
    delete[] stb_framebuffer;
}
