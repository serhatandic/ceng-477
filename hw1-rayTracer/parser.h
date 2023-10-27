#ifndef __HW1__PARSER__
#define __HW1__PARSER__

#include <string>
#include <vector>

namespace parser
{
    //Notice that all the structures are as simple as possible
    //so that you are not enforced to adopt any style or design.
    struct Vec3f
    {
        float x, y, z;
        Vec3f();
        Vec3f(float x, float y, float z);
        Vec3f operator+(const Vec3f& other) const;
        Vec3f operator-(const Vec3f& other) const;
        Vec3f operator*(float scalar) const;
        Vec3f operator*(const Vec3f& other) const;
        Vec3f operator/(const float other) const;
        [[nodiscard]] float dot(const Vec3f& other) const;
        [[nodiscard]] Vec3f cross(const Vec3f& other) const;
        [[nodiscard]] Vec3f normalized() const;
        static float determinant(const Vec3f& a, const Vec3f& b, const Vec3f& c) ;
        Vec3f clamp() const;
        float distance(Vec3f &other) const;
        friend std::ostream& operator<<(std::ostream& os, const Vec3f& v);
    };

    struct Vec3i
    {
        int x, y, z;
        Vec3f operator*(float scalar) const;
    };

    struct Vec4f
    {
        float x, y, z, w;
    };

    struct Camera
    {
        Vec3f position;
        Vec3f gaze;
        Vec3f up;
        Vec4f near_plane;
        float near_distance;
        int image_width, image_height;
        std::string image_name;
    };

    struct PointLight
    {
        Vec3f position;
        Vec3f intensity;
    };

    struct Material
    {
        bool is_mirror;
        Vec3f ambient;
        Vec3f diffuse;
        Vec3f specular;
        Vec3f mirror;
        float phong_exponent;
    };

    struct Face
    {
        int v0_id;
        int v1_id;
        int v2_id;
    };

    struct Mesh
    {
        int material_id;
        std::vector<Face> faces;
    };

    struct Triangle
    {
        int material_id;
        Face indices;
    };

    struct Sphere
    {
        int material_id;
        int center_vertex_id;
        float radius;
    };

    struct Ray {
        Vec3f origin;      // The starting point of the ray
        Vec3f direction;   // The direction the ray is traveling
        int depth;

        Ray();

        Ray(Vec3f origin, Vec3f direction);
    };

    struct HitPoint {
        Vec3f normal;
        Vec3f hitPoint;
        float t;
        int materialId;

        HitPoint(const Vec3f &normal, const Vec3f &hitPoint, const float &t, int materialId);

        HitPoint();
    };

    struct Scene
    {
        //Data
        Vec3i background_color;
        float shadow_ray_epsilon;
        int max_recursion_depth;
        std::vector<Camera> cameras;
        Vec3f ambient_light;
        std::vector<PointLight> point_lights;
        std::vector<Material> materials;
        std::vector<Vec3f> vertex_data;
        std::vector<Mesh> meshes;
        std::vector<Triangle> triangles;
        std::vector<Sphere> spheres;
        //Functions
        void loadFromXml(const std::string &filepath);
        Ray generateRay(int i, int j, Camera &cam); // ray goes through i,j th pixel
        [[nodiscard]] float intersect(Sphere s, parser::Ray ray) const;
        [[nodiscard]] float intersect(Triangle t, parser::Ray ray) const;
        [[nodiscard]] float intersect(parser::Face face, parser::Ray ray) const;
        Vec3f computeColor(const std::vector<PointLight>& pointLights, Vec3f ambientLight, Ray ray, Camera &cam);
        Vec3f applyShading(HitPoint hitPoint, const std::vector<PointLight>& pointLights, Camera &cam, Ray ray);
        HitPoint closestIntersection(Ray ray);
        void renderScene(unsigned char* image);
        bool isShadow(Vec3f intersectionPoint, Vec3f lightSourceLocation);
    };
}

#endif
