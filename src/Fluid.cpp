#include <iostream>
#include <string>
#include "Fluid.h"

#include <cmath>

// Constructor
Fluid::Fluid(int size, float dt, float diff, float visc) : size(size), dt(dt), diff(diff), visc(visc) {
    // Allocate memory for arrays
    s = new float[size * size]();      // Using () initializes with zeros
    density = new float[size * size]();
    Vx = new float[size * size]();
    Vy = new float[size * size]();
    Vx0 = new float[size * size]();
    Vy0 = new float[size * size]();
}

// Destructor
Fluid::~Fluid() {
    // Deallocate memory
    delete[] s;
    delete[] density;
    delete[] Vx;
    delete[] Vy;
    delete[] Vx0;
    delete[] Vy0;
}
    //Goes form 3d to our 1d array
int Fluid::IX(int x, int y, int N) {
    x = (x < 0) ? 0 : (x >= N) ? N - 1 : x;
    y = (y < 0) ? 0 : (y >= N) ? N - 1 : y;
    
    return x + y * N;
}

void Fluid::AddDensity(int x, int y, float amount) {
    int N = size;
    density[Fluid::IX(x, y, size)] += amount;
}

void Fluid::Fade() {    
    for(int i = 0; i < this->size * this->size; i++){
        float d = density[i]-0.001;
        if(d < 0)
            d = 0;
        this->density[i] = d;

    }
}

void Fluid::AddVelocity(int x, int y, float amountX, float amountY) {    
    Vx[Fluid::IX(x, y, size)] += amountX;
    Vy[Fluid::IX(x, y, size)] += amountY;
}

void Fluid::Diffuse(int b, float* x, float* x0, float diff, float dt, int iter, int N) {
    float a = dt * diff * (N - 2) * (N - 2);
    Fluid::LinSolve(b, x, x0, a, 1 + 6 * a, iter, N);
}
bool isValidIndex(int i, int j, int N) {
    return i >= 0 && i < N && j >= 0 && j < N;
}


void Fluid::LinSolve(int b, float* x, float* x0, float a, float c, int iter, int N) {
    float cRecip = 1.0 / c;
    for (int k = 0; k < iter; k++) {
        for (int j = 1; j < N - 1; j++) {
            for (int i = 1; i < N - 1; i++) {
                x[Fluid::IX(i, j, N)] = (x0[Fluid::IX(i, j, N)] + a * (x[Fluid::IX(i+1, j,N)]
                                                                +x[Fluid::IX(i-1, j  , N)]
                                                                +x[Fluid::IX(i  , j+1, N)]
                                                                +x[Fluid::IX(i  , j-1, N)] )) * cRecip;
            }
        }
    
        Fluid::SetBnd(b, x, N);
    }
}

void Fluid::SetBnd(int b, float* x, int N) {
    for (int i = 1; i < N - 1; i++) {
        x[Fluid::IX(i, 0, N)] = b == 2 ? -x[Fluid::IX(i, 1, N)] : x[Fluid::IX(i, 1, N)];
        x[Fluid::IX(i, N - 1, N)] = b == 2 ? -x[Fluid::IX(i, N - 2, N)] : x[Fluid::IX(i, N - 2, N)];
    }

    for (int j = 1; j < N - 1; j++) {
        x[Fluid::IX(0, j, N)] = b == 1 ? -x[Fluid::IX(1, j, N)] : x[Fluid::IX(1, j, N)];
        x[Fluid::IX(N - 1, j, N)] = b == 1 ? -x[Fluid::IX(N - 2, j, N)] : x[Fluid::IX(N - 2, j, N)];
    }


    x[Fluid::IX(0, 0, N)] = 0.5f * (x[Fluid::IX(1, 0, N)] + x[Fluid::IX(0, 1, N)]);
    x[Fluid::IX(0, N - 1, N)] = 0.5f * (x[Fluid::IX(1, N - 1, N)] + x[Fluid::IX(0, N - 2, N)]);
    x[Fluid::IX(N - 1, 0, N)] = 0.5f * (x[Fluid::IX(N - 2, 0, N)] + x[Fluid::IX(N - 1, 1, N)]);
    x[Fluid::IX(N - 1, N - 1, N)] = 0.5f * (x[Fluid::IX(N - 2, N - 1, N)] + x[Fluid::IX(N - 1, N - 2,N)]);
    
}

void Fluid::project(float *velocX, float *velocY, float *p, float *div, int iter, int N)
{
    for (int j = 1; j < N - 1; j++) {
        for (int i = 1; i < N - 1; i++) {
            div[IX(i, j, N)] = -0.5f*(velocX[Fluid::IX(i+1, j  ,N)]
                                    -velocX[Fluid::IX(i-1, j  , N)]
                                    +velocY[Fluid::IX(i  , j+1, N)]
                                    -velocY[Fluid::IX(i  , j-1, N)])/N;
            p[IX(i, j, N)] = 0;
        }
    }

    SetBnd(0, div, N); 
    SetBnd(0, p, N);
    LinSolve(0, p, div, 1, 6, iter, N);
    
    for (int j = 1; j < N - 1; j++) {
        for (int i = 1; i < N - 1; i++) {
            velocX[Fluid::IX(i, j, N)] -= 0.5f * (  p[Fluid::IX(i+1, j, N)] -p[Fluid::IX(i-1, j, N)]) * N;
            velocY[Fluid::IX(i, j, N)] -= 0.5f * (  p[Fluid::IX(i, j+1, N)] -p[Fluid::IX(i, j-1, N)]) * N;
        }
    }
    
    SetBnd(1, velocX, N);
    SetBnd(2, velocY, N);
}

void Fluid::advect(int b, float *d, float *d0, float *velocX, float *velocY, float dt, int N)
{
    float i0, i1, j0, j1;
    
    float dtx = dt * (N - 2);
    float dty = dt * (N - 2);
    
    float s0, s1, t0, t1, u0, u1;
    float tmp1, tmp2, x, y;
    
    float Nfloat = N;
    float ifloat, jfloat;
    int i, j;
    
    
    for (j = 1, jfloat = 1; j < N - 1; j++, jfloat++) {
        for (i = 1, ifloat = 1; i < N - 1; i++, ifloat++) {
            tmp1 = dtx * velocX[Fluid::IX(i, j, N)];
            tmp2 = dty * velocY[Fluid::IX(i, j, N)];
            
            x = ifloat - tmp1;
            y = jfloat - tmp2;

            
            if (x < 0.5f) x = 0.5f;
            if (x > Nfloat + 0.5f) x = Nfloat + 0.5f;
            i0 = floorf(x);
            i1 = i0 + 1.0f;
            if (y < 0.5f) y = 0.5f;
            if (y > Nfloat + 0.5f) y = Nfloat + 0.5f;
            j0 = floorf(y);
            j1 = j0 + 1.0f;

            
            s1 = x - i0;
            s0 = 1.0f - s1;
            t1 = y - j0;
            t0 = 1.0f - t1;

            int i0i =i0;
            int i1i =i1;
            int j0i =j0;
            int j1i =j1;

            
            d[Fluid::IX(i, j, N)] =
                s0 * (t0 * d0[Fluid::IX(i0i, j0i, N)] + t1 * d0[Fluid::IX(i0i, j1i, N)])
                + s1 * (t0 * d0[Fluid::IX(i1i, j0i, N)] + t1 * d0[IX(i1i, j1i,  N)]);
        }
    }

    SetBnd(b, d, N);
}

void Fluid::printDensity() {
    std::cout << "Density values:" << std::endl;
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            std::cout << density[IX(i, j, size)] << " ";
        }
        std::cout << std::endl;
    }
}


void Fluid::printVelocity() {
    std::cout << "Density values:" << std::endl;
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            std::cout << Vx[IX(i, j, size)] << " ";
        }
        std::cout << std::endl;
    }
}
