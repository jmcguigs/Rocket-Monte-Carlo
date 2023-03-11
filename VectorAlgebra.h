#pragma once

#include <math.h>

namespace VectorAlgebra
{
    // BEGIN POINT3D
    class Point3D
    {
    public:
        float x, y, z;
        Point3D(float a, float b, float c)
        {
            x = a;
            y = b;
            z = c;
        }

        float distanceTo(Point3D other)
        {
            float dx = other.x - x;
            float dy = other.y - y;
            float dz = other.z - z;
            return sqrt(dx * dx + dy * dy + dz * dz);
        }

        Vector3D getDirectionTo(Point3D other)
        {
            float dx = other.x - x;
            float dy = other.y - y;
            float dz = other.z - z;

            float distance = sqrt(dx * dx + dy * dy + dz * dz);

            return Vector3D(dx / distance, dy / distance, dz / distance);
        }
    };

    // END POINT3D


    // BEGIN VECTOR3D

    class Vector3D
    {
    public:
        float x, y, z;
        Vector3D(float a, float b, float c)
        {
            x = a;
            y = b;
            z = c;
        }

        // take the dot product of this and another vector
        float dotProduct(Vector3D other)
        {
            return x * other.x + y * other.y + z * other.z;
        }

        // take the cross product of this and another vector
        Vector3D crossProduct(Vector3D other)
        {
            float i = y * other.z - z * other.y;
            float j = -1 * (x * other.z - other.x * z);
            float k = x * other.y - other.x * y;
            return Vector3D(i, j, k);
        }

        // return unit vector parallel to the 3D vector
        Vector3D direction(void)
        {
            float length = sqrt(x * x + y * y + z * z);
            return Vector3D(x / length, y / length, z / length);
        }

        friend Vector3D operator / (Vector3D self, float scalar);

        friend Vector3D operator * (Vector3D self, float scalar);

        friend Vector3D operator + (Vector3D a, Vector3D b);

        friend Vector3D operator - (Vector3D a, Vector3D b);
    };


    Vector3D operator * (Vector3D self, float scalar)
    {
        return Vector3D(self.x * scalar, self.y * scalar, self.z * scalar);
    }

    Vector3D operator / (Vector3D self, float scalar)
    {
        return Vector3D(self.x / scalar, self.y / scalar, self.z / scalar);
    }

    Vector3D operator + (Vector3D a, Vector3D b)
    {
        return Vector3D(a.x + b.x, a.y + b.y, a.z + b.z);
    }

    Vector3D operator - (Vector3D a, Vector3D b)
    {
        return Vector3D(a.x - b.x, a.y - b.y, a.z - b.z);
    }

    // END VECTOR3D

}
