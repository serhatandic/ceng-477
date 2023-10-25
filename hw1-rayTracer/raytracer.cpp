#include <iostream>
#include "parser.h"
#include "ppm.h"

typedef unsigned char RGB[3];

namespace parser {
    Vec3f::Vec3f() : x(0), y(0), z(0) {}

    Vec3f::Vec3f(float x, float y, float z) : x(x), y(y), z(z) {}

    Vec3f Vec3f::operator+(const Vec3f& other) const {
        return {x + other.x, y + other.y, z + other.z};
    }

    Vec3f Vec3f::operator-(const Vec3f& other) const {
        return {x - other.x, y- other.y, z - other.z};
    }

    Vec3f Vec3f::operator*(float scalar) const {
        return {x * scalar, y * scalar, z * scalar};
    }

    float Vec3f::dot(const Vec3f& other) const {
        return x * other.x + y * other.y + z * other.z;
    }

    Vec3f Vec3f::cross(const Vec3f& other) const {
        return {y * other.z - z * other.y, z * other.x - x * other.z, x * other.y - y * other.x};
    }

    Vec3f Vec3f::normalized() const {
        float length = sqrt(x * x + y * y + z * z);
        if (length > 0) {
            return *this * (1 / length);
        }
        return {};
    }

    Ray Scene::generateRay(int i, int j) { // ray goes through i,j th pixel
        float distFromLeft, distFromTop;
        Vec3f imageCenter, topLeftCorner, tipOfRay;
        Vec3f u, v, w;
        Camera cam = this -> cameras[0];

        w = cam.gaze * -1;
        v = cam.up;
        u = v.cross(w);

        float left = cam.near_plane.x;
        float right = cam.near_plane.y;
        float bottom = cam.near_plane.z;
        float top = cam.near_plane.w;

        int imageWidth = cam.image_width;
        int imageHeight = cam.image_height;

        float pixelWidth = (right - left) / (float) imageWidth;
        float pixelHeight = (top - bottom) / (float)imageHeight;

        distFromLeft = (i + 0.5) * pixelWidth;
        distFromTop = (j + 0.5) * pixelHeight;

        Vec3f gaze = cam.gaze;
        Vec3f camPosition = cam.position;

        float camImageDistance = cam.near_distance;

        imageCenter = camPosition + gaze * camImageDistance;
        topLeftCorner = imageCenter + u * left + v * top;
        tipOfRay = topLeftCorner + u*distFromLeft + v * (-distFromTop);

        return {camPosition, (tipOfRay + camPosition * -1).normalized()};
    }
}

std::ostream& operator<<(std::ostream& os, const parser::Vec3f& v) {
    return os << "(" << v.x << ", " << v.y << ", " << v.z << ")";
}

int main(int argc, char* argv[])
{
    // Sample usage for reading an XML scene file
    parser::Scene scene;

    scene.loadFromXml(argv[1]);

    // The code below creates a test pattern and writes
    // it to a PPM file to demonstrate the usage of the
    // ppm_write function.
    //
    // Normally, you would be running your ray tracing
    // code here to produce the desired image.

    const RGB BAR_COLOR[8] =
    {
        { 255, 255, 255 },  // 100% White
        { 255, 255,   0 },  // Yellow
        {   0, 255, 255 },  // Cyan
        {   0, 255,   0 },  // Green
        { 255,   0, 255 },  // Magenta
        { 255,   0,   0 },  // Red
        {   0,   0, 255 },  // Blue
        {   0,   0,   0 },  // Black
    };

    int width = 640, height = 480;
    int columnWidth = width / 8;

    unsigned char* image = new unsigned char [width * height * 3];

    int i = 0;
    for (int y = 0; y < height; ++y)
    {
        for (int x = 0; x < width; ++x)
        {
            int colIdx = x / columnWidth;
            image[i++] = BAR_COLOR[colIdx][0];
            image[i++] = BAR_COLOR[colIdx][1];
            image[i++] = BAR_COLOR[colIdx][2];
        }
    }

    write_ppm("test.ppm", image, width, height);

}
