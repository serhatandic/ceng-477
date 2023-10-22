#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>

class Vec3f {
public:
    double x, y, z;

    Vec3f(): x(0), y(0), z(0) {}

    Vec3f(double x, double y, double z) : x(x), y(y), z(z) {}

    Vec3f operator+(const Vec3f& other) const {
        return {x + other.x, y + other.y, z + other.z};
    }

    Vec3f operator-(const Vec3f& other) const {
        return {x - other.x, y- other.y, z - other.z};
    }

    Vec3f operator*(double scalar) const {
        return {x * scalar, y * scalar, z * scalar};
    }

    [[nodiscard]] double dot(const Vec3f& other) const {
        return x * other.x + y * other.y + z * other.z;
    }

    [[nodiscard]] Vec3f cross(const Vec3f& other) const {
        return {y * other.z - z * other.y, z * other.x - x * other.z, x * other.y - y * other.x};
    }

    [[nodiscard]] Vec3f normalized() const {
        double length = sqrt(x * x + y * y + z * z);
        if (length > 0) {
            return *this * (1 / length);
        }
        return {};
    }

    friend std::ostream& operator<<(std::ostream& os, const Vec3f& v) {
        return os << "(" << v.x << ", " << v.y << ", " << v.z << ")";
    }
};

class Sphere{
public:
    Vec3f center;
    Vec3f color;
    double radius;

    Sphere() {}

    Sphere(const Vec3f &center, const Vec3f &color, double radius) : center(center), color(color), radius(radius) {}

};

class Ray{
public:
    Vec3f origin, direction;

    Ray() = default;

    Ray(const Vec3f &origin, const Vec3f &direction) : origin(origin), direction(direction) {}

    [[nodiscard]] double intersectSphere(Sphere s) const{
        double A, B, C; // constants for the eqn
        double delta;
        double t, t1, t2;
        Vec3f c;
        c = s.center;

        A = direction.dot(direction);
        B = 2*direction.dot(origin - c);
        C = (origin - c).dot(origin - c) - s.radius*s.radius;

        delta = B*B - 4*A*C;

        if (delta == 0){
            t = -B / (2*A);
            return t;
        }else if (delta > 0){
            t1 = (-B - sqrt(delta)) / (2*A);
            t2 = (-B + sqrt(delta)) / (2*A);
            return std::min(t1, t2);
        }
        return -1;
    }

    [[nodiscard]] Vec3f computeColor(const std::vector<Sphere>& spheres) const {
        Vec3f c;
        double minT = 90000;
        double t;

        for(auto & sphere : spheres){
            t = intersectSphere(sphere);
            if (t < minT && t >= 1){
                c = sphere.color;
                minT = t;
            }
        }
        return c;
    }

    friend std::ostream& operator<<(std::ostream& os, const Ray& r) {
        return os << "Origin: " << r.origin << ", Direction: " << r.direction;
    }
};

typedef struct {
    int x,y,z;
} vec3i;

class Scene{
public:
    int nx, ny; // nx columns and ny rows
    double dist, left, right, top, bottom;
    Vec3f e; // eye location
    Vec3f gaze;
    Vec3f u, v, w;
    std::vector<std::vector<vec3i>> image;

    Scene(int width, int height): nx(width), ny(height), image(height, std::vector<vec3i>(width)){}

    void initializeScene() {
        w = gaze * -1;
        u = v.cross(w).normalized();

        for (int i = 0; i < ny; ++i) {
            for (int j = 0; j < nx; ++j) {
                image[i][j].x = image[i][j].y = image[i][j].z = 0;
            }
        }
    }

    [[nodiscard]] Ray generateRay(int i, int j) const {
        double su, sv;
        Vec3f m, q, s;

        double pixelWidth = (right - left) / nx;
        double pixelHeight = (top - bottom) / ny;

        su = (i + 0.5) * pixelWidth;
        sv = (j + 0.5) * pixelHeight;

        m = e + gaze.normalized() * dist; // center of the image window
        q = m + u * left + v * top; // top left corner of the image window
        s = q + u * su + v * (-1 * sv);

        return {e, s + e * -1};
    }

    void renderScene(const std::vector<Sphere>& spheres) {
        for (int j = 0; j < ny; j++) {
            for (int i = 0; i < nx; i++) {
                Ray myRay = generateRay(i, j);
                Vec3f pixel = myRay.origin + myRay.direction;
                Vec3f rayColor = myRay.computeColor(spheres);
                image[j][i].x = round(rayColor.x*255);
                image[j][i].y = round(rayColor.y*255);
                image[j][i].z = round(rayColor.z*255);

            }
        }
    }

    void writePPM() const {
        std::ofstream file("output.ppm");

        if (!file) {
            std::cerr << "Error: Could not open the file for writing.\n";
            return;
        }

        file << "P3\n";
        file << "#output.ppm\n";
        file << nx << " " << ny << "\n";
        file << "255\n";

        for (int i = 0; i < ny; i++) {
            for (int j = 0; j < nx; j++) {
                file << image[i][j].x << " " << image[i][j].y << " " << image[i][j].z << "\t";
            }
            file << "\n";
        }
    }

};

int main() {

    int nx, ny;
    std::cin >> nx >> ny;

    Scene scene(nx, ny);

    std::vector<Sphere> spheres(2);

    std::cin >> scene.dist >> scene.left >> scene.right >> scene.top >> scene.bottom;
    std::cin >> scene.e.x >> scene.e.y >> scene.e.z;
    std::cin >> scene.gaze.x >> scene.gaze.y >> scene.gaze.z;
    std::cin >> scene.v.x >> scene.v.y >> scene.v.z;
    std::cin >> scene.u.x >> scene.u.y >> scene.u.z;

    for (auto & sphere : spheres){
        std::cin >> sphere.center.x >> sphere.center.y >> sphere.center.z >> sphere.radius >> sphere.color.x >> sphere.color.y >> sphere.color.z ;
    }

    scene.initializeScene();
    scene.renderScene(spheres);
    scene.writePPM();

    return 0;
}
