// FluidCube.h

#ifndef FLUID_H
#define FLUID_H

class Fluid {
public:
    // Constructor
    Fluid(int size, float dt, float diff, float visc);

    // Destructor
    ~Fluid();

    void printDensity();
    void printVelocity();
    void AddDensity(int x, int y, float amount);
    void AddVelocity(int x, int y, float amountX, float amountY);
    void Fade();

    static int IX(int x, int y, int N);
    static void Diffuse(int b, float* x, float* x0, float diff, float dt, int iter, int N);
    static void LinSolve(int b, float* x, float* x0, float a, float c, int iter, int N);
    static void SetBnd(int b, float* x, int N);

    static void project(float *velocX, float *velocY, float *p, float *div, int iter, int N);
    static void advect(int b, float *d, float *d0, float *velocX, float *velocY, float dt, int N);


    int size;
    float dt;
    int diff;
    int visc;

    float* s;
    float* density;
    float* Vx;
    float* Vy;
    float* Vx0;
    float* Vy0;

};

#endif // FLUID_CUBE_H
