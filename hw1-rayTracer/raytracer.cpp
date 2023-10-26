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

    float Vec3f::determinant(const Vec3f& a, const Vec3f& b, const Vec3f& c) {
        return a.x * (b.y * c.z - c.y * b.z)
               - a.y * (b.x * c.z - c.x * b.z)
               + a.z * (b.x * c.y - c.x * b.y);
    }

    float Vec3f::distance(Vec3f &other) const {
        return sqrt((x - other.x)*(x - other.x) + (y - other.y)*(y - other.y) + (z - other.z)*(z - other.z));
    }
    Vec3f Vec3i::operator*(float scalar) const {
        return {x * scalar, y * scalar, z * scalar};
    }

    HitPoint::HitPoint(const Vec3f &normal, const Vec3f &hitPoint, const float &t,int materialId) : normal(normal), hitPoint(hitPoint), t(t),
                                                                                     materialId(materialId) {}

    HitPoint::HitPoint() {}

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

    float Scene::intersect(parser::Triangle triangle, parser::Ray ray) const {
        float Beta, Gamma, t;
        float detOfA, detForBeta, detForGamma, detFort;

        std::vector<Vec3f> vertexData = this->vertex_data;
        Vec3f triangleA = vertexData[triangle.indices.v0_id - 1];
        Vec3f triangleB = vertexData[triangle.indices.v1_id - 1];
        Vec3f triangleC = vertexData[triangle.indices.v2_id - 1];

        detOfA = parser::Vec3f::determinant(triangleA - triangleB, triangleA - triangleC, ray.direction);
        if (detOfA == 0){
            return -1;
        }

        detForBeta = parser::Vec3f::determinant(triangleA - ray.origin, triangleA - triangleC, ray.direction);
        detForGamma = parser::Vec3f::determinant(triangleA - triangleB, triangleA - ray.origin, ray.direction);
        detFort = parser::Vec3f::determinant(triangleA - triangleB, triangleA - triangleC, triangleA - ray.origin);

        Beta = detForBeta / detOfA;
        Gamma = detForGamma / detOfA;
        t = detFort / detOfA;
        if (Beta + Gamma <= 1 && Beta >= 0 && Gamma >= 0 && t > 0){
            return t;
        }
        return -1;
    }

    float Scene::intersect(parser::Face face, parser::Ray ray) const {
        float Beta, Gamma, t;
        float detOfA, detForBeta, detForGamma, detFort;

        std::vector<Vec3f> vertexData = this->vertex_data;
        Vec3f triangleA = vertexData[face.v0_id - 1];
        Vec3f triangleB = vertexData[face.v1_id - 1];
        Vec3f triangleC = vertexData[face.v2_id - 1];

        detOfA = parser::Vec3f::determinant(triangleA - triangleB, triangleA - triangleC, ray.direction);
        if (detOfA == 0){
            return -1;
        }

        detForBeta = parser::Vec3f::determinant(triangleA - ray.origin, triangleA - triangleC, ray.direction);
        detForGamma = parser::Vec3f::determinant(triangleA - triangleB, triangleA - ray.origin, ray.direction);
        detFort = parser::Vec3f::determinant(triangleA - triangleB, triangleA - triangleC, triangleA - ray.origin);

        Beta = detForBeta / detOfA;
        Gamma = detForGamma / detOfA;
        t = detFort / detOfA;
        if (Beta + Gamma <= 1 && Beta >= 0 && Gamma >= 0 && t > 0){
            return t;
        }
        return -1;
    }

    HitPoint Scene::closestIntersection(Ray ray){
        HitPoint hitPoint;
        float minT = 90000;
        hitPoint.t = -1;

        for (Sphere sphere: spheres){
            float sphereIntersection = intersect(sphere, ray);
            if (sphereIntersection < minT && sphereIntersection > 0){
                hitPoint.t = sphereIntersection;
                hitPoint.materialId = sphere.material_id;
                hitPoint.hitPoint = ray.origin + ray.direction * sphereIntersection;
                hitPoint.normal = (hitPoint.hitPoint - vertex_data[sphere.center_vertex_id - 1]).normalized();
                minT = sphereIntersection;
            }
        }

        for (Triangle triangle: triangles){
            float triangleIntersection = intersect(triangle, ray);
            if (triangleIntersection < minT && triangleIntersection > 0){
                hitPoint.t = triangleIntersection;
                hitPoint.materialId = triangle.material_id;
                hitPoint.hitPoint = ray.origin + ray.direction * triangleIntersection;
                Vec3f triangleA = vertex_data[triangle.indices.v0_id - 1];
                Vec3f triangleB = vertex_data[triangle.indices.v1_id - 1];
                Vec3f triangleC = vertex_data[triangle.indices.v2_id - 1];

                hitPoint.normal = (triangleA - triangleB).cross(triangleA - triangleC);
                minT = triangleIntersection;
            }
        }

        for (const Mesh& mesh: meshes){
            for (Face face: mesh.faces){
                float faceIntersection = intersect(face, ray);
                if (faceIntersection < minT && faceIntersection > 0){
                    hitPoint.t = faceIntersection;
                    hitPoint.materialId = mesh.material_id;
                    hitPoint.hitPoint = ray.origin + ray.direction * faceIntersection;

                    Vec3f faceA = vertex_data[face.v0_id - 1];
                    Vec3f faceB = vertex_data[face.v1_id - 1];
                    Vec3f faceC = vertex_data[face.v2_id - 1];

                    hitPoint.normal = (faceA - faceB).cross(faceA - faceC);
                    minT = faceIntersection;
                }
            }
        }

        return hitPoint;
    }

    Vec3f Scene::computeColor(HitPoint hitPoint, const std::vector<PointLight>& pointLights, Vec3f ambientLight, Ray ray, Camera &cam) {
        float t = hitPoint.t;
        Vec3f normal = hitPoint.normal.normalized();
        Vec3f intersectionPoint = hitPoint.hitPoint;

        std::cout << t << std::endl;
        if (t < 0) return this->background_color * 1.0f; // No intersection, return bgcolor

        Vec3f totalSpecular, totalDiffuse, ambient;
        Material material = materials[hitPoint.materialId- 1];
        ambient = (ambientLight * material.ambient);

        for (PointLight pointLight: pointLights){
            Vec3f lightDir = (pointLight.position - intersectionPoint).normalized();
            Vec3f viewDir = (cam.position - intersectionPoint).normalized(); // Assuming camera is in the scene

            // Diffuse
            float diff = std::max(normal.dot(lightDir), 0.0f);
            Vec3f diffuse = (pointLight.intensity)  * material.diffuse * diff;

            // Specular
            Vec3f halfVector = (lightDir + viewDir).normalized();
            float spec = std::pow(std::max(halfVector.dot(normal), 0.0f), material.phong_exponent);
            Vec3f specular = (pointLight.intensity  / pow(pointLight.position.distance(intersectionPoint),2)) * material.specular * spec;

            totalSpecular = specular + totalSpecular;
            totalDiffuse = diffuse + totalDiffuse;
        }

        // Sum up all components
        Vec3f result = totalDiffuse + totalSpecular + ambient;
        std::cout << result.x << " " << result.y << " " << result.z << std::endl;
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
                HitPoint hitPoint = closestIntersection(myRay);

                Vec3f rayColor = computeColor(hitPoint,this->point_lights, this->ambient_light, myRay, cam);

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
