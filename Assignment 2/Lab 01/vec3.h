#pragma once

#include <cmath>

#define M_PI 3.141592654

class vec3 {

public:
    double x, y, z;

    inline vec3() :x(0.0f), y(0.0f), z(0.0f) {}
    ~vec3() {};

    inline vec3(double x, double y, double z) :x(x), y(y), z(z) {};

    inline vec3(double* arr) :x(arr[0]), y(arr[1]), z(arr[2]) {};

    /* Operator overloads */

    vec3& operator+=(const vec3& rhs) {
        x += rhs.x;
        y += rhs.y;
        z += rhs.z;
        return *this;
    }

    vec3& operator-=(const vec3& rhs) {
        x -= rhs.x;
        y -= rhs.y;
        z -= rhs.z;
        return *this;
    }

    vec3& operator*=(double scale) {
        x *= scale;
        y *= scale;
        z *= scale;
        return *this;
    }

    vec3& operator/=(double scale) {
        x /= scale;
        y /= scale;
        z /= scale;
        return *this;
    }

    vec3 operator-() const {
        return vec3(-x, -y, -z);
    }

    /* ------------------ */

    /* Utility */


    double magnitudeSq() const {
        return x * x + y * y + z * z;
    }

    double magnitude() const {
        return sqrt(magnitudeSq());
    }

    void normalize() {
        double mag = magnitude();
        x /= mag;
        y /= mag;
        z /= mag;
    }

    void copyToArray(double arr[3]) {
        arr[0] = x;
        arr[1] = y;
        arr[2] = z;
    }
};

inline vec3 operator+(vec3 lhs, const vec3 rhs) {
    lhs += rhs;
    return lhs;
}

inline vec3 operator-(vec3 lhs, const vec3 rhs) {
    lhs -= rhs;
    return lhs;
}

inline vec3 operator*(vec3 lhs, double rhs) {
    lhs *= rhs;
    return lhs;
}

inline vec3 operator/(vec3 lhs, double rhs) {
    lhs /= rhs;
    return lhs;
}

inline double dot(const vec3& a, const vec3& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

inline vec3 cross(const vec3& a, const vec3& b) {
    return vec3(
        a.y * b.z - b.y * a.z,
        -(a.x * b.z - b.x * a.z),
        a.x * b.y - b.x * a.y
    );
}

inline double rad2Deg(double angleRad) {
    return angleRad * 180.0f / M_PI;
}

inline double deg2Rad(double angleDeg) {
    return angleDeg * M_PI / 180.0f;
}

inline double magnitudeSq(const vec3& v) {
    return v.x * v.x + v.y * v.y + v.z * v.z;
}

inline double magnitude(const vec3& v) {
    return sqrt(magnitudeSq(v));
}

inline vec3 normalize(const vec3& v) {
    double mag = magnitude(v);
    return v / mag;
}

