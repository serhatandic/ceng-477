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

    Vec3f Vec3f::operator*(const Vec3f& other) const {
        return {x * other.x, y * other.y, z * other.z};
    }

    Vec3f Vec3f::operator/(const float other) const {
        return {x / other, y / other, z / other};
    }

    float Vec3f::dot(const Vec3f& other) const {
        return x * other.x + y * other.y + z * other.z;
    }

    Vec3f Vec3f::cross(const Vec3f& other) const {
        return {y * other.z - z * other.y, z * other.x - x * other.z, x * other.y - y * other.x};
    }
    Vec3f Vec3f::clamp() const {
        Vec3f color = {this->x, this->y, this->z};
        return {
                std::max(0.0f, std::min((float)255, color.x)),
                std::max(0.0f, std::min((float)255, color.y)),
                std::max(0.0f, std::min((float)255, color.z))
        };
    }
    Vec3f Vec3f::normalized() const {
        float length = sqrt(x * x + y * y + z * z);
        if (length > 0) {
            return *this * (1 / length);
        }
        return {};
    }

    Ray Scene::generateRay(int i, int j, Camera &cam) { // ray goes through i,j th pixel
        float distFromLeft, distFromTop;
        Vec3f imageCenter, topLeftCorner, tipOfRay;
        Vec3f u, v, w;

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

        distFromLeft = (float)(i + 0.5) * pixelWidth;
        distFromTop = (float)(j + 0.5) * pixelHeight;

        Vec3f gaze = cam.gaze;
        Vec3f camPosition = cam.position;

        float camImageDistance = cam.near_distance;

        imageCenter = camPosition + gaze * camImageDistance;
        topLeftCorner = imageCenter + u * left + v * top;
        tipOfRay = topLeftCorner + u*distFromLeft + v * (-distFromTop);

        return {camPosition, (tipOfRay + camPosition * -1).normalized()};
    }

    float Scene::intersect(parser::Sphere s, parser::Ray ray) const {
        float A, B, C; // constants for the quadratic equation
        float delta;
        float singleRoot, doubleRootFirst, doubleRootSecond;

        std::vector<Vec3f> vertexData = this->vertex_data;

        Vec3f sphereCenter = vertexData[s.center_vertex_id - 1];

        A = ray.direction.dot(ray.direction);
        B = 2*ray.direction.dot(ray.origin - sphereCenter);
        C = (ray.origin - sphereCenter).dot(ray.origin - sphereCenter) - s.radius * s.radius;

        delta = B*B - 4*A*C;
        if (delta == 0){
            singleRoot = -B / (2*A);
            return singleRoot;
        }else if (delta > 0){
            doubleRootFirst = (-B - sqrt(delta)) / (2*A);
            doubleRootSecond = (-B + sqrt(delta)) / (2*A);
            return std::min(doubleRootFirst, doubleRootSecond);
        }
        return -1;
    }

    Vec3f Scene::computeColor(Sphere s, PointLight pointLight, Vec3f ambientLight, Ray ray, Camera &cam) {
        float t = intersect(s, ray); // I'm assuming intersect returns a float and is properly defined elsewhere
        if (t < 0) return {0, 0, 0}; // No intersection, return black TODO: should return background color

        std::vector<Vec3f> vertexData = this->vertex_data;
        Vec3f sphereCenter = vertexData[s.center_vertex_id - 1];

        Vec3f intersectionPoint = ray.origin + ray.direction * t;
        Vec3f normal = (intersectionPoint - sphereCenter).normalized();
        Vec3f lightDir = (pointLight.position - intersectionPoint).normalized();
        // Vec3f viewDir = (cam.position - intersectionPoint).normalized(); // Assuming camera is in the scene
        Material material = materials[s.material_id - 1];

        // Ambient
        Vec3f ambient = ambientLight * material.ambient;

        // Diffuse
        float diff = std::max(normal.dot(lightDir), 0.0f);
        Vec3f diffuse = pointLight.intensity * material.diffuse * diff;

/*
        // Specular
        Vec3f reflectDir = reflect(-lightDir, normal);
        float spec = std::pow(std::max(viewDir.dot(reflectDir), 0.0f), material.phong_exponent);
        Vec3f specular = pointLight.intensity * material.specular * spec;*/

        // Sum up all components
        Vec3f result = ambient + diffuse;
        return result.clamp();
    }

    void Scene::renderScene(unsigned char* image) {
        Camera cam = this->cameras[0];

        int imageWidth = cam.image_width;
        int imageHeight = cam.image_height;
        int k = 0;
        for (int j = 0; j < imageHeight; j++) {
            for (int i = 0; i < imageWidth; i++) {
                Ray myRay = generateRay(i, j, cam);
                Vec3f rayColor = computeColor(this->spheres[0], this->point_lights[0], this->ambient_light, myRay, cam);

                image[k++] = (unsigned char) rayColor.x;
                image[k++] = (unsigned char) rayColor.y;
                image[k++] = (unsigned char) rayColor.z;

            }
        }
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
    parser::Camera cam = scene.cameras[0];
    int width = cam.image_width, height = cam.image_height;

    unsigned char* image = new unsigned char [width * height * 3];

    scene.renderScene(image);

    write_ppm("test.ppm", image, width, height);

}
