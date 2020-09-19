#pragma once
#include "Ray.hpp"
#include "Intersection.hpp"
#include "Vector.hpp"
#include "Scene.hpp"
#include "Frame.h"
#include <string>
//#include <algorithm>

typedef Vector3f Vec3f;
typedef Vector2f Vec2f;

template<typename T>
T Sqr(const T& a) { return a * a; }

#define EPS_PHONG 1e-3f
#define PI_F     3.14159265358979f
#define INV_PI_F (1.f / PI_F)
#define EPS_COSINE 1e-6f
#define EPS_RAY    1e-4f

//////////////////////////////////////////////////////////////////////////
// Cosine lobe hemisphere sampling
Vec3f SamplePowerCosHemisphereW(
	const Vec2f& aSamples,
	const float  aPower,
	float* oPdfW)
{
	const float term1 = 2.f * PI_F * aSamples.x;
	const float term2 = std::pow(aSamples.y, 1.f / (aPower + 1.f));
	const float term3 = std::sqrt(1.f - term2 * term2);

	if (oPdfW)
	{
		*oPdfW = (aPower + 1.f) * std::pow(term2, aPower) * (0.5f * INV_PI_F);
	}

	return Vec3f(
		std::cos(term1) * term3,
		std::sin(term1) * term3,
		term2);
}

float PowerCosHemispherePdfW(
	const Vec3f& aNormal,
	const Vec3f& aDirection,
	const float  aPower)
{
	const float cosTheta = std::max(0.f, dotProduct(aNormal, aDirection));

	return (aPower + 1.f) * std::pow(cosTheta, aPower) * (INV_PI_F * 0.5f);
}

//////////////////////////////////////////////////////////////////////////
/// Sample direction in the upper hemisphere with cosine-proportional pdf
/** The returned PDF is with respect to solid angle measure */
Vec3f SampleCosHemisphereW(
	const Vec2f& aSamples,
	float* oPdfW)
{
	const float term1 = 2.f * PI_F * aSamples.x;
	const float term2 = std::sqrt(1.f - aSamples.y);

	const Vec3f ret(
		std::cos(term1) * term2,
		std::sin(term1) * term2,
		std::sqrt(aSamples.y));

	if (oPdfW)
	{
		*oPdfW = ret.z * INV_PI_F;
	}

	return ret;
}

float FresnelDielectric(
	float aCosInc,
	float IOR)
{
	if (IOR < 0)
		return 1.f;

	float eta;

	if (aCosInc < 0.f)
	{
		aCosInc = -aCosInc;
		eta = IOR;
	}
	else
	{
		eta = 1.f / IOR;
	}

	const float sinTrans2 = Sqr(eta) * (1.f - Sqr(aCosInc));
	const float cosTrans = std::sqrt(std::max(0.f, 1.f - sinTrans2));

	const float term1 = eta * cosTrans;
	const float rParallel =
		(aCosInc - term1) / (aCosInc + term1);

	const float term2 = eta * aCosInc;
	const float rPerpendicular =
		(term2 - cosTrans) / (term2 + cosTrans);

	return 0.5f * (Sqr(rParallel) + Sqr(rPerpendicular));
}

float Luminance(const Vec3f& aRGB)
{
	return 0.212671f * aRGB.x +
		0.715160f * aRGB.y +
		0.072169f * aRGB.z;
}

Vec3f ReflectLocal(const Vec3f& aVector)
{
	return Vec3f(-aVector.x, -aVector.y, aVector.z);
}

enum Events
{
	diffuse,
	Phong,
	reflect,
	refract,
};

typedef Vector3f Vec3f;
#define EPS_COSINE 1e-6f
class BSDF
{
	struct ComponentProbabilities
	{
		float diffProb;
		float phongProb;
	};
public:
	BSDF(const Ray& ray, const Intersection& intersectRes)
	{
		matType = "";
		frame.SetFromZ(intersectRes.normal);
		directionLFixed = frame.ToLocal(-ray.direction);

		if (std::abs(directionLFixed.z) < EPS_COSINE)
		{
			return;
		}

		matType = intersectRes.m->getType();
		mat = intersectRes.m;
		GetProbabilities(mat, probs);
	}

	Vec3f Evaluate(
		const Vec3f& directionW,
		float& theta,
		float* pdf) const
	{
		Vec3f result(0);
		*pdf = 0;

		const Vec3f directionL = frame.ToLocal(directionW);

		if (directionL.z * directionLFixed.z < 0)
			return result;

		theta = std::abs(directionL.z);

		//const Material& mat = scn.matVec[matID];

		result += EvaluateDiffuse(mat, directionL, pdf);
		//result += EvaluatePhong(mat, directionL, pdf);
		return result;
	}

	Vec3f Sample(
		const Vec3f& rn3,
		Vec3f& directionW,
		float& pdf,
		float& theta) const
	{
		Events sampledEvent;
		pdf = 0;
		Vec3f result(0);
		Vec3f directionL;

		if (rn3.z < probs.diffProb)
			sampledEvent = diffuse;
		//else if (rn3.z < probs.diffProb + probs.phongProb)
		//	sampledEvent = Phong;
		//else if (rn3.z < probs.diffProb + probs.phongProb + probs.reflectProb)
		//	sampledEvent = reflect;
		else
			sampledEvent = diffuse;


		//const Material& mat = scn.matVec[matID];

		switch (sampledEvent)
		{
		case diffuse:
		{
			result += SampleDiffuse(mat, rn3.GetXY(), directionL, pdf);

			if (result.IsZero())
				return Vec3f(0);

			//result += EvaluatePhong(mat, directionL, &pdf);
			break;
		}
		case Phong:
		{
			result += SamplePhong(mat, rn3.GetXY(), directionL, pdf);

			if (result.IsZero())
				return Vec3f(0);

			result += EvaluateDiffuse(mat, directionL, &pdf);
			break;
		}
		default:
			break;
		}

		theta = std::abs(directionL.z);

		if (theta < EPS_COSINE)
			return Vec3f(0.f);
		directionW = frame.ToWorld(directionL);
		return result;
	}

	float ContinuationProb() const
	{
		return probContinue;
	}

private:
	void GetProbabilities(const Material* mat, ComponentProbabilities& probs)
	{
		reflectCoeff = FresnelDielectric(directionLFixed.z, mat->ior);

		const float albedoDiffuse = Luminance(mat->m_color);
		const float albedoPhong = Luminance(mat->Ks);
		//const float albedoReflect = reflectCoeff * Luminance(mat->ior);//perfect specular
		//const float albedoRefract = (1.f - reflectCoeff) * (mat->ior > 0.f ? 1.f : 0.f);
		const float totalAlbedo = albedoDiffuse + albedoPhong;// + albedoReflect + albedoRefract;

		probs.diffProb = albedoDiffuse / totalAlbedo;
		probs.phongProb = albedoPhong / totalAlbedo;
		//probs.reflectProb = albedoReflect / totalAlbedo;
		//probs.refractProb = albedoRefract / totalAlbedo;

		probContinue = Vector3f::MaxVec((Vector3f(mat->m_color) + Vector3f(mat->Ks))); //+ reflectCoeff * Vector3f(0.5))) + (1.f - reflectCoeff);
		probContinue = std::min(1.f, std::max(0.f, probContinue));
	}

	Vec3f SampleDiffuse(
		const Material* mat,
		const Vec2f& rn3,
		Vec3f& directionL,
		float& pdf) const
	{
		if (directionLFixed.z < EPS_COSINE)
			return Vec3f(0);

		float weight;
		directionL = SampleCosHemisphereW(rn3, &weight);
		pdf += weight * probs.diffProb;

		return mat->m_color * INV_PI_F;
	}

	Vec3f SamplePhong(
		const Material* mat,
		const Vec2f& rn3,
		Vec3f& directionL,
		float& pdf) const
	{
		if (probs.phongProb == 0)
			return Vec3f(0.f);

		directionL = SamplePowerCosHemisphereW(rn3, mat->specularExponent, NULL);

		const Vec3f reflLocalDirFixed = ReflectLocal(directionLFixed);

		Frame frame;
		frame.SetFromZ(reflLocalDirFixed);
		directionL = frame.ToWorld(directionL);

		float dot_R_Wi = dotProduct(reflLocalDirFixed, directionL);

		if (dot_R_Wi <= EPS_PHONG)
			return Vec3f(0.f);

		pdf += PowerCosHemispherePdfW(reflLocalDirFixed, directionL, mat->specularExponent) * probs.phongProb;
		const Vec3f rho = mat->Ks * (mat->specularExponent + 2.f) * 0.5f * INV_PI_F;

		return rho * std::pow(dot_R_Wi, mat->specularExponent);
	}

	Vec3f EvaluateDiffuse(
		const Material* mat,
		const Vec3f& directionL,
		float* pdf) const
	{
		if (probs.diffProb == 0)
			return Vec3f(0);

		if (directionLFixed.z < EPS_COSINE || directionL.z < EPS_COSINE)
			return Vec3f(0);

		*pdf += probs.diffProb * std::max(0.f, directionL.z * INV_PI_F);

		return mat->m_color * INV_PI_F;
	}

	Vec3f EvaluatePhong(
		const Material* mat,
		const Vec3f& directionL,
		float* pdf) const
	{
		if (probs.phongProb == 0)
			return Vec3f(0);

		if (directionLFixed.z < EPS_COSINE || directionL.z < EPS_COSINE)
			return Vec3f(0);

		const Vec3f reflLocalDirIn = ReflectLocal(directionLFixed);
		const float dot_R_Wi = dotProduct(reflLocalDirIn, directionL);

		if (dot_R_Wi <= EPS_PHONG)
			return Vec3f(0.f);

		*pdf += probs.phongProb * PowerCosHemispherePdfW(reflLocalDirIn, directionL, mat->specularExponent);

		const Vec3f rho = mat->Ks * (mat->specularExponent + 2.f) * 0.5f * INV_PI_F;

		return rho * std::pow(dot_R_Wi, mat->specularExponent);
	}

private:
	Vec3f directionLFixed;
	ComponentProbabilities probs;
	float probContinue;
	float reflectCoeff;
	//int   matID;
	std::string matType;
	const Material* mat;
	Frame frame;
};