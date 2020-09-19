#include "Scene.hpp"

void Scene::buildBVH() {
	printf(" - Generating BVH...\n\n");
	this->bvh = new BVHAccel(objects, 1, BVHAccel::SplitMethod::NAIVE);
}

Intersection Scene::intersect(const Ray &ray) const
{
	return this->bvh->Intersect(ray);
}

bool Scene::trace(
	const Ray &ray,
	const std::vector<Object*> &objects,
	float &tNear, uint32_t &index, Object **hitObject)
{
	*hitObject = nullptr;
	for (uint32_t k = 0; k < objects.size(); ++k) {
		float tNearK = kInfinity;
		uint32_t indexK;
		Vector2f uvK;
		if (objects[k]->intersect(ray, tNearK, indexK) && tNearK < tNear) {
			*hitObject = objects[k];
			tNear = tNearK;
			index = indexK;
		}
	}


	return (*hitObject != nullptr);
}

// Implementation of the Whitted-syle light transport algorithm (E [S*] (D|G) L)
//
// This function is the function that compute the color at the intersection point
// of a ray defined by a position and a direction. Note that thus function is recursive (it calls itself).
//
// If the material of the intersected object is either reflective or reflective and refractive,
// then we compute the reflection/refracton direction and cast two new rays into the scene
// by calling the castRay() function recursively. When the surface is transparent, we mix
// the reflection and refraction color using the result of the fresnel equations (it computes
// the amount of reflection and refractin depending on the surface normal, incident view direction
// and surface refractive index).
//
// If the surface is duffuse/glossy we use the Phong illumation model to compute the color
// at the intersection point.

std::default_random_engine generator;
std::uniform_real_distribution<float> distribution(0, 1);
//uniform_distribution = std::uniform_real_distribution<float>(0.0f, 1.0f);

void createCoordinateSystem(const Vector3f &N, Vector3f &Nt, Vector3f &Nb)
{
	if (std::fabs(N.x) > std::fabs(N.y))
		Nt = Vector3f(N.z, 0, -N.x) / sqrtf(N.x * N.x + N.z * N.z);
	else
		Nt = Vector3f(0, -N.z, N.y) / sqrtf(N.y * N.y + N.z * N.z);
	Nb = crossProduct(N, Nt);
}

//inversion method using polar coordinates theta,phi -> cartesian coordinates x,y,z
Vector3f uniformSampleHemisphere(const float &r1, const float &r2)
{
	// cos(theta) = u1 = y
	// cos^2(theta) + sin^2(theta) = 1 -> sin(theta) = srtf(1 - cos^2(theta))
	float sinTheta = sqrtf(1 - r1 * r1);
	float phi = 2 * M_PI * r2;
	float x = sinTheta * cosf(phi);
	float z = sinTheta * sinf(phi);
	return Vector3f(x, r1, z);
}

//cosine-weighted lambertian
Vector3f importanceSampleHemisphere(const float &r1, const float &r2)
{
	float phi = 2 * M_PI * r2;
	float SinTheta = sqrtf(1 -	r1);
	float x = SinTheta * cosf(phi);
	float z = SinTheta * sinf(phi);

	return Vector3f(x, r1, z);
}

Vector3f generateCosineHemisphere(const float& r1, const float& r2)
{
	float u = r1;
	float v = r2;
	float r = sqrtf(u);
	float phi = v * 2 * M_PI;
	float x = r * cosf(phi);
	float y = r * sinf(phi);
	float z = sqrtf(std::max(0.0f, 1.0f - u));
	return Vector3f(x, y, z); //prob = z*1/PI
}

Vector3f Scene::castRay(const Ray &ray, int depth) const
{
	if (depth > this->maxDepth) {
		return Vector3f(0.0, 0.0, 0.0);
	}
	Intersection intersection = Scene::intersect(ray);
	Material *m = intersection.m;
	Object *hitObject = intersection.obj;
	Vector3f hitColor = Vector3f(1.0f);
	//Vector3f hitColor = this->backgroundColor;
	//    float tnear = kInfinity;
	Vector2f uv;
	uint32_t index = 0;
	if (intersection.happened) {

		Vector3f hitPoint = intersection.coords;
		Vector3f N = intersection.normal; // normal
		Vector2f st; // st coordinates
		hitObject->getSurfaceProperties(hitPoint, ray.direction, index, uv, N, st);
		//        Vector3f tmp = hitPoint;
		switch (m->getType()) {
		/*case REFLECTION_AND_REFRACTION:
		{
			Vector3f reflectionDirection = normalize(reflect(ray.direction, N));
			Vector3f refractionDirection = normalize(refract(ray.direction, N, m->ior));
			Vector3f reflectionRayOrig = (dotProduct(reflectionDirection, N) < 0) ?
				hitPoint - N * EPSILON :
				hitPoint + N * EPSILON;
			Vector3f refractionRayOrig = (dotProduct(refractionDirection, N) < 0) ?
				hitPoint - N * EPSILON :
				hitPoint + N * EPSILON;
			Vector3f reflectionColor = castRay(Ray(reflectionRayOrig, reflectionDirection), depth + 1);
			Vector3f refractionColor = castRay(Ray(refractionRayOrig, refractionDirection), depth + 1);
			float kr;
			fresnel(ray.direction, N, m->ior, kr);
			hitColor = reflectionColor * kr + refractionColor * (1 - kr);
			break;
		}
		case REFLECTION:
		{
			float kr;
			fresnel(ray.direction, N, m->ior, kr);
			Vector3f reflectionDirection = reflect(ray.direction, N);
			Vector3f reflectionRayOrig = (dotProduct(reflectionDirection, N) < 0) ?
				hitPoint + N * EPSILON :
				hitPoint - N * EPSILON;
			hitColor = castRay(Ray(reflectionRayOrig, reflectionDirection), depth + 1) * kr;
			break;
		}*/
		default:	//diffuse and glossy
		{
			// [comment]
			// We use the Phong illumation model int the default case. The phong model
			// is composed of a diffuse and a specular reflection component.
			// [/comment]
			Vector3f direct_lightAmt = 0, specularColor = 0;
			Vector3f shadowPointOrig = (dotProduct(ray.direction, N) < 0) ?
				hitPoint + N * EPSILON :
				hitPoint - N * EPSILON;
			// [comment]
			// Loop over all lights in the scene and sum their contribution up
			// We also apply the lambert cosine law
			// [/comment]
			//direct lighting
			for (uint32_t i = 0; i < get_lights().size(); ++i)
			{
				auto area_ptr = dynamic_cast<AreaLight*>(this->get_lights()[i].get());
				if (area_ptr)
				{
					// handling area light
				}
				else
				{
					auto directionToLight = -1 * get_lights()[0]->position;
					auto lightDir = normalize(directionToLight);


					//Vector3f lightDir = get_lights()[i]->position - hitPoint;
					// square of the distance between hitPoint and the light
					float lightDistance2 = kInfinity;
					//lightDir = normalize(lightDir);
					float LdotN = std::max(0.f, dotProduct(lightDir, N));
					Object* shadowHitObject = nullptr;
					float tNearShadow = kInfinity;
					// is the point in shadow, and is the nearest occluding object closer to the object than the light itself?
					bool inShadow = bvh->Intersect(Ray(shadowPointOrig, lightDir)).happened;
					direct_lightAmt += (1 - inShadow) * get_lights()[i]->intensity * LdotN;
					Vector3f reflectionDirection = reflect(-lightDir, N);
					//specularColor += powf(std::max(0.f, -dotProduct(reflectionDirection, ray.direction)),
					//	m->specularExponent) * get_lights()[i]->intensity;
					
					//direct_lightAmt = direct_lightAmt /(float)No_of_samps;
				}
			}

			// [comment]
			// Compute indirect ligthing
			// [/comment]
			Vector3f indirectLigthing(0.0f);
//#ifdef GI
			uint32_t No_of_samps = 16;// / (depth + 1);
			Vector3f Nt, Nb;
			createCoordinateSystem(N, Nt, Nb);
			float pdf = 1 / (2 * M_PI);
			for (uint32_t n = 0; n < No_of_samps; ++n) {
				float r1 = distribution(generator);
				float r2 = distribution(generator);
				Vector3f sample = uniformSampleHemisphere(r1, r2);
				Vector3f sampleWorld(
					sample.x * Nb.x + sample.y * N.x + sample.z * Nt.x,
					sample.x * Nb.y + sample.y * N.y + sample.z * Nt.y,
					sample.x * Nb.z + sample.y * N.z + sample.z * Nt.z);
				// don't forget to divide by PDF and multiply by cos(theta)
				auto rL = castRay(Ray(hitPoint + sampleWorld * 0.0001,
					sampleWorld), depth + 1);
				indirectLigthing += rL;
				//indirectLigthing += castRay(hitPoint,
				//	sampleWorld, objects, lights, options, depth + 1);
			}
			// divide by N
			indirectLigthing = indirectLigthing / (float)No_of_samps;
//#endif
			
			//hitColor = (direct_lightAmt) * (hitObject->evalDiffuseColor(st) * m->Kd);
			hitColor = (direct_lightAmt / M_PI + indirectLigthing) * hitObject->evalDiffuseColor(st);// + specularColor * m->Ks;
			break;
			}
		}
	}

	//if (hitColor.x >= 1.0 || hitColor.y >= 1.0 || hitColor.z >= 1.0)
	//	std::cout << hitColor << "\n";

	hitColor.x = clamp(0.0f, 1.0f, hitColor.x);
	hitColor.y = clamp(0.0f, 1.0f, hitColor.y);
	hitColor.z = clamp(0.0f, 1.0f, hitColor.z);
	return hitColor;
}

typedef unsigned uint;
typedef Vector3f Vec3f;
typedef Vector2f Vec2f;
#include "BSDF.h"
#include "RNG.h"
Vector3f Scene::castRayIterative(Ray& ray, int depth) const
{
	//ray tracer constructor
	int maxPathLength = 8;
	int iter = 0;
	RNG randGen = int(rand() % 1000);
	//buffer.Setup(aScene.cam.res);

	Vec3f pathWeight(1.f);
	Vec3f color(0.f);
	//uint pathLength = 1;

	for (int pathLength=1;pathLength<= maxPathLength; ++pathLength)
	{
		Intersection intersection = Scene::intersect(ray);
		if (intersection.happened)
		{
			Vector3f N = intersection.normal; // normal
			Vec3f hitPoint = intersection.coords;
			BSDF bsdf(ray, intersection);

			if (pathLength >= maxPathLength || bsdf.ContinuationProb() < EPSILON)
				break;

			//direction light
			Vec3f directionToLight = -1 * get_lights()[0]->position;
			directionToLight = normalize(directionToLight);
			Vector3f lightDir = directionToLight;
			float lightDistance2 = kInfinity;

			float distance = 1e36f;
			Vec3f radiance = get_lights()[0]->intensity;
			float bsdfPdfW, cosThetaOut;

			Vec3f factor = bsdf.Evaluate(directionToLight, cosThetaOut, &bsdfPdfW);
			float LdotN = std::max(0.f, dotProduct(lightDir, N));
			Object* shadowHitObject = nullptr;
			float tNearShadow = kInfinity;

			Vector3f shadowPointOrig = (dotProduct(ray.direction, N) < 0) ?
				hitPoint + N * EPSILON :
				hitPoint - N * EPSILON;
			// is the point in shadow, and is the nearest occluding object closer to the object than the light itself?
			auto sIntersect = bvh->Intersect(Ray(shadowPointOrig, lightDir));
			bool inShadow = sIntersect.happened; //&& sIntersect.distance* sIntersect.distance < lightDistance2;
			if (!factor.IsZero() && !inShadow)
			{
				Vec3f contrib = cosThetaOut * radiance * factor;
				color += pathWeight * contrib;
			}

			Vec3f rn3 = randGen.GetVec3f();
			float pdf;

			Vec3f factorSample = bsdf.Sample(rn3, ray.direction, pdf, cosThetaOut);
			if (factorSample.IsZero())
				break;

			float continueProb = bsdf.ContinuationProb();
			pdf *= std::min(1.f, continueProb);
			pathWeight = pathWeight * (factorSample * (cosThetaOut / pdf));
			ray.origin = hitPoint + ray.direction * EPS_RAY;
			intersection.distance = 1e36f;
		}
	}
	color.x = clamp(0.0f, 1.0f, color.x);
	color.y = clamp(0.0f, 1.0f, color.y);
	color.z = clamp(0.0f, 1.0f, color.z);
	return color;
}
