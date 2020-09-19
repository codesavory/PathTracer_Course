#pragma once
#include "Vector.hpp"
typedef Vector3f Vec3f;

//////////////////////////////////////////////////////////////////////////
// BMP TOOL
// References: http://www.smallvcm.com
#pragma region BMPUTIL
class Frame
{
public:

	Frame()
	{
		mX = Vec3f(1, 0, 0);
		mY = Vec3f(0, 1, 0);
		mZ = Vec3f(0, 0, 1);
	};

	Frame(const Vec3f& x, const Vec3f& y, const Vec3f& z) : mX(x), mY(y), mZ(z) {}

	void SetFromZ(const Vec3f& z)
	{
		Vec3f tmpZ = mZ = normalize(z);
		Vec3f tmpX = (std::abs(tmpZ.x) > 0.99f) ? Vec3f(0, 1, 0) : Vec3f(1, 0, 0);
		mY = normalize(crossProduct(tmpZ, tmpX));
		mX = crossProduct(mY, tmpZ);
	}

	Vec3f ToWorld(const Vec3f& a) const
	{
		return mX * a.x + mY * a.y + mZ * a.z;
	}

	Vec3f ToLocal(const Vec3f& a) const
	{
		return Vec3f(dotProduct(a, mX), dotProduct(a, mY), dotProduct(a, mZ));
	}

public:
	Vec3f mX, mY, mZ;
};
