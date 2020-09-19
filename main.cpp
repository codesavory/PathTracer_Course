#include "Renderer.hpp"
#include "Scene.hpp"
#include "Triangle.hpp"
#include "Vector.hpp"
#include "global.hpp"
#include "Sphere.hpp"
#include <chrono>
#include "Material.hpp"

// In the main function of the program, we create the scene (create objects and
// lights) as well as set the options for the render (image width and height,
// maximum recursion depth, field-of-view, etc.). We then call the render
// function().
int main(int argc, char** argv)
{
    Scene scene(1500, 1500);
	Material* red_mat = new Material(MaterialType::DIFFUSE_AND_GLOSSY,
		Vector3f(0.913, 0.003, 0.078), Vector3f(0, 0, 0));
	red_mat->Kd = 0.6;
	red_mat->Ks = 0.6;//default 0
	red_mat->specularExponent = 2;

	Material* green_mat = new Material(MaterialType::DIFFUSE_AND_GLOSSY,
		Vector3f(0.301, 0.560, 0.329), Vector3f(0, 0, 0));
	green_mat->Kd = 0.6;
	green_mat->Ks = 0.6;//default 0
	green_mat->specularExponent = 2;//default 0

	Material* blue_mat = new Material(MaterialType::DIFFUSE_AND_GLOSSY,
		Vector3f(0.180, 0.305, 0.878), Vector3f(0, 0, 0));
	blue_mat->Kd = 0.6;
	blue_mat->Ks = 0.0;//default 0
	blue_mat->specularExponent = 0;//default 0
	//blue_mat->ior = 0.5;

	Material* grey_floor = new Material(MaterialType::DIFFUSE_AND_GLOSSY,
		Vector3f(0.819, 0.819, 0.819), Vector3f(0, 0, 0));
	grey_floor->Kd = 0.6;
	grey_floor->Ks = 0.8;//default 0
	grey_floor->specularExponent = 10;//default 0

	Material* ip_mat = new Material(MaterialType::DIFFUSE_AND_GLOSSY,
		Vector3f(0.913, 0.003, 0.078), Vector3f(0, 0, 0));
	ip_mat->Kd = 0.6;
	ip_mat->Ks = 0.6;//default 0
	ip_mat->specularExponent = 2;//default 0

	Material* white_mat = new Material(MaterialType::DIFFUSE_AND_GLOSSY,
		Vector3f(0.0), Vector3f(0, 0, 0));
	grey_floor->Kd = 0.6;
	grey_floor->Ks = 0.6;//default 0
	grey_floor->specularExponent = 2;//default 0

	//Cornell Box Scene
	/*objl::Loader loader;
	loader.LoadFile("models\\CornellBox-Original\\CB_with_light.obj");
	assert(loader.LoadedMeshes.size() >= 1);
	MeshTriangle floor(loader.LoadedMeshes[0], grey_floor);
	MeshTriangle ceiling(loader.LoadedMeshes[1], grey_floor);
	MeshTriangle backWall(loader.LoadedMeshes[2], grey_floor);
	MeshTriangle rightWall(loader.LoadedMeshes[3], grey_floor);
	MeshTriangle leftWall(loader.LoadedMeshes[4], grey_floor);
	MeshTriangle shortBox(loader.LoadedMeshes[5], green_mat);
	MeshTriangle tallBox(loader.LoadedMeshes[6], red_mat);
	MeshTriangle sphere_light(loader.LoadedMeshes[8], red_mat);
	scene.Add(&floor);
	scene.Add(&ceiling);
	scene.Add(&backWall);
	scene.Add(&rightWall);
	scene.Add(&leftWall);
	scene.Add(&shortBox);
	scene.Add(&tallBox);
	//scene.Add(&sphere_light);*/

	//load buddha-dragon-bunny scene
	objl::Loader loader;
	loader.LoadFile("models\\complex_scene.obj");
	assert(loader.LoadedMeshes.size() >= 1);
	MeshTriangle bunny(loader.LoadedMeshes[0], blue_mat);
	MeshTriangle buddha(loader.LoadedMeshes[1], red_mat);
	MeshTriangle dragon(loader.LoadedMeshes[2], green_mat);
	MeshTriangle plane(loader.LoadedMeshes[3], grey_floor);
	scene.Add(&bunny);
	scene.Add(&dragon);
	scene.Add(&buddha);
	scene.Add(&plane);

	/*objl::Loader loader;
	loader.LoadFile("D:/CG_Source/Rendering_Pjts/Path Tracer/Path_Tracer_Extension/The_Path_Tracer/models/bunny_floor.obj");
	assert(loader.LoadedMeshes.size() >= 1);
	MeshTriangle bunny(loader.LoadedMeshes[0], blue_mat);
	MeshTriangle plane(loader.LoadedMeshes[1], grey_floor);
	scene.Add(&bunny);
	scene.Add(&plane);*/

	Vector3f lightDir(1.0f, 1.5f, 1.0f);
	lightDir = -1 * normalize(lightDir);
	scene.Add(std::make_unique<Light>(lightDir, 2.0f));
	//scene.Add(std::make_unique<Light>(Vector3f(-1, 5, 10), 1));
    //scene.Add(std::make_unique<Light>(Vector3f(-20, 70, 20), 1));
    //scene.Add(std::make_unique<Light>(Vector3f(20, 70, 20), 1));
    //scene.Add(std::make_unique<AreaLight>(Vector3f(3, 5, -5), 0.6));
    scene.buildBVH();

    Renderer r;

    auto start = std::chrono::system_clock::now();
    r.Render(scene);
    auto stop = std::chrono::system_clock::now();

    std::cout << "Render complete: \n";
    std::cout << "Time taken: " << std::chrono::duration_cast<std::chrono::hours>(stop - start).count() << " hours\n";
    std::cout << "          : " << std::chrono::duration_cast<std::chrono::minutes>(stop - start).count() << " minutes\n";
    std::cout << "          : " << std::chrono::duration_cast<std::chrono::seconds>(stop - start).count() << " seconds\n";

    return 0;
}