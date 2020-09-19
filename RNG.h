#pragma once
#include <random>
#include "Vector.hpp"

typedef Vector3f Vec3f;
typedef Vector2f Vec2f;
//////////////////////////////////////////////////////////////////////////
// Random Number Generator = RNG
class RNG
{
public:
	RNG(int aSeed = 1234) : randGen(aSeed) {}

	float GetFloat()
	{
		return mDistFloat(randGen);
	}

	Vec2f GetVec2f()
	{
		return Vec2f(GetFloat(), GetFloat());
	}

	Vec3f GetVec3f()
	{
		return Vec3f(GetFloat(), GetFloat(), GetFloat());
	}

private:
	std::mt19937_64 randGen;
	std::uniform_real_distribution<float> mDistFloat;
};