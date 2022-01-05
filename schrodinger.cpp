/*
FTCS-method solver of the Schrodinger equation on a 2D grid. Initial condition
is a Gaussian wavepacket and the scalar potential is a hard wall on all 4 sides (a
large Gaussian spike).
*/

#include <iostream>
#include <fstream>
#include <cmath>

#define IX(i, j) ((i) + (N+2)*(j))

typedef struct complex  {
    float real;
    float im;
} complex;

// SIMULATION PARAMETERS. CHANGE AS YOU WISH
const int N = 99;
const int size = (N+2) * (N+2);
float dx = 0.01;
float dt = 2.5e-6;
int frames = 400;

// WAVEFUNCTION AND POTENTIAL
complex psi[size];
float potential[size];

// COMPLEX FUNCTIONS
float modSquared(complex z) {
    return z.real * z.real + z.im * z.im;
}

float modulus(complex z)    {
    return sqrt(modSquared(z));
}

complex reciprocal(complex z)   {
    float mod = modulus(z);
    z.real /= mod; z.im /= -mod;
    return z;
}

complex operator + (complex a, complex b)   {
    a.real += b.real; a.im += b.im;
    return a;
}

complex operator - (complex a, complex b)   {
    a.real -= b.real; a.im -= b.im;
    return a;
}

complex operator * (complex a, complex b)   {
    complex z;
    z.real = (a.real * b.real) - (a.im * b.im); z.im = (a.im * b.real) + (a.real * b.im);
    return z;
}

complex operator * (float k, complex z) {
    z.real *= k; z.im *= k;
    return z;
}

complex operator / (complex a, complex b)   {
    return a * reciprocal(b);
}

float scalarProduct2D(float * x, float * y) {
    return x[0] * y[0] + x[1] * y[1];
}

// ADDS YOUNG DOUBLE SLITS / HARD WALL. PROGRAMMED AS GAUSSIAN SPIKES WITH SHORT WIDTH
void addHorizontalSlits(float * potential, int slitWidth, int slitDistance, int y) {
    float constant = 0.5 / (dx * dx);
    int middle = (N+1) / 2;

    for (int i = 0; i < N+2; i++)  {
        if ((i <= middle - slitWidth - slitDistance / 2) || ((i >= middle - slitDistance / 2) && (i <= middle + slitDistance / 2)) || (i >= middle + slitWidth + slitDistance / 2)) {
            for (int j = 0; j < N+2; j++)   {
                potential[IX(i, j)] += 10000 * exp(-constant * (j-y) * (j-y));
            }
        }
    }
}

// SETS EVERYTHING TO 0. MUST BE INCLUDED IN MAIN.
void initialise(complex * psi, float * potential)   {
    complex zero; zero.real = 0; zero.im = 0;
    for (int i = 0; i < N+2; i++)   {
        for (int j = 0; j < N+2; j++)   {
            psi[IX(i, j)] = zero;
            potential[IX(i, j)] = 0;
        }
    }
}

// INSERTS GAUSSIAN WAVEPACKET AT POSITION r0 WITH WAVEVECTOR k AND WIDTH sigma. 
void setup(complex * psi, float * potential, float * r0, float * k, float sigma) {
    float constant1 = 0.5 / (sigma * sigma);
    float constant2 = 0.5 / (dx * dx);
    complex factor;
    float pos[2];
    for (int i = 0; i < N+2; i++)  {
        pos[0] = dx * i;
        for (int j = 0; j < N+2; j++)   {
            pos[1] = dx * j;
            float relPos[2] = {pos[0] - r0[0], pos[1] - r0[1]};
            factor.real = cos(scalarProduct2D(k, pos)); factor.im = sin(scalarProduct2D(k, pos));
            psi[IX(i, j)] = psi[IX(i, j)] + exp(-constant1 * scalarProduct2D(relPos, relPos)) * factor;
            potential[IX(i, j)] += 10000 * exp(-constant2 * (i * i + (j - N - 1) * (i - N - 1) + j * j + (j - N - 1) * (j - N - 1)));
            }
        }

    // THE FOLLOWING CODE MAY NOT BE NEEDED. NORMALISES WAVEFUNCTION SO INTEGRAL OVER GRID OF |psi|^2 = 1.
    float integral = 0;
    for (int i = 0; i < N+2; i++)  {
        for (int j = 0; j < N+2; j++)   {
            integral += dx * dx * modSquared(psi[IX(i, j)]);
        }
    }
    float a = 1/integral;
    for (int i = 0; i < N+2; i++)  {
        for (int j = 0; j < N+2; j++)   {
            psi[IX(i, j)] = a * psi[IX(i, j)];
        }
    }
}

// RETURNS i * [laplacian(psi) - potential * psi]. FOR USE IN RUNGE-KUTTA STEP.
complex function(complex * psi, complex RKConstant, float * potential, int i, int j)  {
    complex constant; constant.real = 0; constant.im = 0.5 / (dx * dx);

    return constant * (psi[IX(i+1, j)] + psi[IX(i-1, j)] + psi[IX(i, j+1)] + psi[IX(i, j-1)] - 4 * (psi[IX(i, j)] + RKConstant)) - potential[IX(i, j)] * (psi[IX(i, j)] + RKConstant);
}

// THE GOOD OLD EULER FORWARD STEP THAT DID NOT WORK VERY WELL. SHOULD PROBABLY BE REMOVED.
void eulerStep(complex * psi)   {
    complex zero; zero.real = 0; zero.im = 0;
    complex * store = psi;
    for (int i = 0; i < N+1; i++)   {
        for (int j = 0; j < N+1; j++)   {
            psi[IX(i, j)] = psi[IX(i, j)] + dt * function(psi, zero, potential, i, j);
        }
    }
}

// THE MUCH BETTER RUNGA-KUTTA 4th ORDER STEP.
void RK4Step(complex * psi) {
    complex k1, k2, k3, k4, zero;
    zero.real = 0; zero.im = 0; // 0 represented as 0 + 0i
    complex * store = psi; // temporary storage of psi
    for (int i = 1; i < N+1; i++)   {
        for (int j = 1; j < N+1; j++)   {
            k1 = function(psi, zero, potential, i, j);
            k2 = function(psi, 0.5 * dt * k1, potential, i, j);
            k3 = function(psi, 0.5 * dt * k2, potential, i, j);
            k4 = function(psi, dt * k3,       potential, i, j);
            psi[IX(i, j)] = store[IX(i, j)] + (dt / 6) * (k1 + 2 * (k2 + k3) + k4);
        }
    }
}

int main()  {
    float r0[2] = {dx * N/2, dx * N / 4};
    float k[2]  = {0, 60};
    float sigma = 0.1;

    initialise(psi, potential);
    setup(psi, potential, r0, k, sigma);
    addHorizontalSlits(potential, 0, 10, (N+1)/2);  // CAN BE DELETED, CHANGED, ETC.
    
    std::ofstream file;
    file.open("schrodinger_sim.dat");
    file << N << "," << frames << "," << "\n";

    for (int i = 0; i < frames; i++)    {
        std::cout << i << "\n";
        for (int k = 0; k < 16; k++)    {   // DONE 16 TIMES EVERY FRAME SO THE VIDEO ISN'T TOO SLOW
            RK4Step(psi);
        }
        
        for (int x = 0; x < N+2; x++)   {
            for (int y = 0; y < N+2; y++)   {
                file << modSquared(psi[IX(x, y)]) << ",";   // WE OUTPUT |psi|^2 TO THE FILE
            }
        }
        file << "\n";
    }

    file.close();
    return 0;
}
