/**
 * @file Math.h
 * @author wertrain
 */
#pragma once
#define _USE_MATH_DEFINES
#include <cmath>

struct V
{
    double x, y, z;
    V(){}
    V(const double v) : x(v), y(v), z(v) {}
    V(const double x, const double y, const double z) : x(x), y(y), z(z) {}
    double &operator[](int i) const { (&x)[i]; }
};

V operator+(V a, V b)
{
    return V(a.x + b.x, a.y + b.y, a.z + b.z);
}
V operator-(V a, V b)
{
    return V(a.x - b.x, a.y - b.y, a.z - b.z);
}
V operator*(V a, V b)
{
    return V(a.x * b.x, a.y * b.y, a.z * b.z);
}
V operator/(V a, V b)
{
    return V(a.x / b.x, a.y / b.y, a.z / b.z);
}
V operator-(V v)
{
    return V(-v.x, -v.y, -v.z);
}

double dot(const V& a, const V& b)
{
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

V cross(const V& a, const V& b)
{
    return V(
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x
    );
}

V normalize(const V &v)
{
    return v / sqrt(dot(v, v));
}

double DegreeToRadian(const double degree)
{
    return degree * M_PI / 180.0f;
}

double RadianToDegree(const double radian)
{
    return radian * 180.0f / M_PI;
}

std::tuple<V, V> tangentSpace(const V& n)
{
    const double s = std::copysign(1, n.z);
    const double a = -1 / (s + n.z);
    const double b = n.x*n.y*a;
    return
    {
        V(1 + s * n.x*n.x*a,s*b,-s * n.x),
        V(b,s + n.y*n.y*a,-n.y)
    };
}
