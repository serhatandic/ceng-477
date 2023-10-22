#include <iostream>
#include <cstdio>
#include <cmath>

typedef struct {
    int x,y,z;
} vec3i;

typedef struct {
    double x,y,z;
} vec3f;

typedef struct {
    vec3f o, d;
} ray;

double dist, left, right, top, bottom;
vec3i ** image;
int nx, ny;

vec3f e; // eye location
vec3f gaze;
vec3f u, v, w;

vec3f multS(vec3f a, double s){
    vec3f result;
    result.x = a.x*s;
    result.y = a.y*s;
    result.z = a.z*s;
    return result;
}

vec3f add(vec3f a, vec3f b){
    vec3f result;
    result.x = a.x + b.x;
    result.y = a.y + b.y;
    result.z = a.z + b.z;
    return result;
}

double dotProduct(const vec3f& a, const vec3f& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

vec3f crossProduct(const vec3f& a, const vec3f& b) {
    vec3f result;

    result.x = a.y * b.z - a.z * b.y;
    result.y = a.z * b.x - a.x * b.z;
    result.z = a.x * b.y - a.y * b.x;

    return result;
}

vec3f normalize(vec3f vLocal) {
    double magnitude = sqrt(vLocal.x * vLocal.x + vLocal.y * vLocal.y + vLocal.z * vLocal.z);
    if (magnitude == 0) {
        std::cerr << "Warning: Trying to normalize a zero vector." << std::endl;
        return {0, 0, 0}; // Return a zero vector to avoid division by zero
    }
    return {vLocal.x / magnitude, vLocal.y / magnitude, vLocal.z / magnitude};
}

ray generateRay(int i, int j){
    ray result;
    double su, sv;
    vec3f m, q, s;

    double pixelWidth = (right - left) / nx;
    double pixelHeight = (top - bottom) / ny;

    su = (i + 0.5)*pixelWidth;
    sv = (j + 0.5)*pixelHeight;

    m = add(e, multS(normalize(gaze), dist)); // center of the image window
    q = add(m, add(multS(u, left), multS(v, top)));
    s = add(q, add(multS(u, su), multS(v, -1*sv)));

    result.o = e;
    result.d = add(s, multS(e, -1));

    return result;
}

void writePPM(vec3i** img){
    int i, j;
    FILE *fp;

    fp = fopen("output.ppm", "w");

    fprintf(fp, "P3\n");
    fprintf(fp, "#output.ppm\n");
    fprintf(fp, "%d %d\n", nx, ny);
    fprintf(fp, "255\n");

    for (i = 0; i < ny; i++){
        for (j = 0; j < nx; j++){
            fprintf(fp, "%d %d %d\t", img[i][j].x, img[i][j].y, img[i][j].z);
        }
        fprintf(fp, "\n");
    }
}

int main() {

    int i, j;


/*    scanf("%d %d", &nx, &ny);
    scanf("%f %f %f %f %f %f", &dist, &left, &right, &top, &bottom);*/

    std::cin >> nx >> ny;
    std::cin >> dist >> left >> right >> top >> bottom;
    std::cin >> e.x >> e.y >> e.z;
    std::cin >> gaze.x >> gaze.y >> gaze.z;
    std::cin >> v.x >> v.y >> v.z;

    w = multS(gaze, -1);
    u = crossProduct(v, w);
    image = (vec3i **) malloc(sizeof(vec3i*)*ny);

    for (i = 0; i < ny; i++){
        image[i] = (vec3i *) malloc(sizeof(vec3i) * nx);
    }

    for (i = 0; i < ny; i++){
        for (j = 0; j < nx; j++){
            image[i][j].x = image[i][j].y = image[i][j].z = 105;
        }
    }

    // main ray tracing loop
    for (i = 0; i < ny; i++) {
        for (j = 0; j < nx; j++) {
            ray myRay = generateRay(j,i);
            vec3f pixel;
            pixel = add(myRay.o, myRay.d);
            // std::cout << myRay.d.x << " " << myRay.d.y << " " << myRay.d.z << std::endl;
            std::cout << pixel.x << " " << pixel.y << " " << pixel.z << std::endl;
        }
    }


    writePPM(image);
    return 0;
}
