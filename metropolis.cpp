#include <iostream>
#include <string>
#include <fstream>
#include <stdlib.h>
#include <cmath>

double const kB = 1.38065e-23; // Boltzmann constant
double const m = 39.95*1.6605e-27; // Argon atom mass
double const A = 1e-10; // Angstrom 

void initial_configuration(int N, double l, double m, double T, double** r0, double** v0);
double maxwell_boltzmann(double m, double T, double v);
double energy(double eps, double sigma, double m, int n, double l, double** r0, double** v0, double& T);
double potential_energy(double eps, double sigma, double m, int n, double l, double** r0);
double delta_potential_energy(double eps, double sigma, double m, int n, double l, int idx, double** r0, double** r1);
double drand(); // generates random numbers between 0 and 1 following uniform distribution

int main(int argc, const char* argv[])
{
    double rho, l, eps, sigma; // density, simulation cell length, Lennard-Jones parameters
    int N; // number of atoms
    long int Ntr; // number of relaxation time steps
    long int Nt; // number of simulation time steps
    double *r0mat, *r1mat, *v0mat, *v1mat;
    double **r0, **r1, **v0, **v1; // 2D arrays for positions and velocities
    double T; // temperature
    double smax; // maximum displacement (< 0.5*l)
    double ene, dene; // potential energy and delta potential energy
    std::ofstream fene, flast;
    
    eps = 120*kB;
    sigma = 3.4*A;
    std::cout << "Number of atoms: \n";
    std::cin >> N;
    std::cout << "Density: (kg/m^3) \n";
    std::cin >> rho;
    std::cout << "Temperature (K): \n";
    std::cin >> T;
    std::cout << "Maximum displacement (Angstrom): \n";
    std::cin >> smax;
    smax *= A;
    std::cout << "Relaxation time steps: \n";
    std::cin >> Ntr;
    std::cout << "Simulation time steps: \n";
    std::cin >> Nt;
    l = pow(m*N/rho, 1./3.);
    std::cout << "Simulation cell length (m): " << l << '\n';

    // dynamic allocation
    r0mat = new double[N*3];
    r1mat = new double[N*3];
    v0mat = new double[N*3];
    v1mat = new double[N*3];
    r0 = new double*[N];
    r1 = new double*[N];
    v0 = new double*[N];
    v1 = new double*[N];
    for (int i = 0; i < N; i++) {
        r0[i] = &r0mat[i*3];
        r1[i] = &r1mat[i*3];
        v0[i] = &v0mat[i*3];
        v1[i] = &v1mat[i*3];
    }

    std::cout << "Creating initial configuration... \n";
    initial_configuration(N, l, m, T, r0, v0);

    fene.open("energy.dat", std::ios::out);
    ene = potential_energy(eps, sigma, m, N, l, r0);

    // Metropolis algorithm implementation
    for (int it = 0; it < Ntr + Nt; it++) {
        int idx = (int) (N-1)*drand();
        for (int i = 0; i < N;  i++) {
            for (int j = 0; j < 3; j++) {
                r1[i][j] = r0[i][j];
            }
        }
        for (int j = 0; j < 3; j++) {
            r1[idx][j] += drand()*2*smax - smax;
        }
        // verify periodic boundary conditions
        for (int j = 0; j < 3; j++) {
            if (r1[idx][j] > l) r1[idx][j] -= l;
            else if (r1[idx][j] < 0.) r1[idx][j] += l;
        }
        dene = delta_potential_energy(eps, sigma, m, N, l, idx, r0, r1);
        // first acceptance attempt
        if (dene < 0) { 
            for (int i = 0; i < N; i++) {
                for (int j = 0; j < 3; j++) {
                    r0[i][j] = r1[i][j];
                }
            }
            ene += dene;
        }
        // second acceptance attempt
        else {
            double r = drand();
            if (exp(-dene/(kB*T)) > r) {
                for (int i = 0; i < N; i++) {
                    for (int j = 0; j < 3; j++) {
                        r0[i][j] = r1[i][j];
                    }
                }
                ene += dene;
            }
            else {
            }
        }
        fene << it << '\t' << ene << '\n';
    }

    // save final configuration 
    flast.open("argonlast.txt", std::ios::out);
    flast << N << '\t' << l/A << '\n';
    for (int i = 0; i < N; i++) {
        flast << r0[i][0]/A << '\t' << r0[i][1]/A << '\t' << r0[i][2]/A << '\n';
    }
    flast.close();

    //deallocation
    delete [] r0mat;
    delete [] r1mat;
    delete [] v0mat;
    delete [] v1mat;
    delete [] r0;
    delete [] r1;
    delete [] v0;
    delete [] v1;
}


void initial_configuration(int N, double l, double m, double T, double** r0, double** v0) {
    int i = 0;
    double f, c, v, d;
    int ni;
    ni = (int) pow((double) N, 1./3.);
    if (ni*ni*ni < N) ni++;
    double vp = sqrt(2.*kB*T/m); // most likely velocity
    double fmax = maxwell_boltzmann(m, T, vp);
    for (int ix = 0; ix < ni; ix++) {
        for (int iy = 0; iy < ni; iy++) {
            for (int iz = 0; iz < ni; iz++) {
                if (i < N) {
                    r0[i][0] = l / ((double) ni) * ix;
                    r0[i][1] = l / ((double) ni) * iy;
                    r0[i][2] = l / ((double) ni) * iz;
                    do {
                        v = 4.*vp*drand();
                        c = fmax*drand();
                        f = maxwell_boltzmann(m, T, v);
                    } while (c > f);
                    d = 0.;
                    for (int j = 0; j < 3; j++){
                        v0[i][j] = drand()-0.5;
                        d += pow(v0[i][j], 2.);
                    }
                    d = sqrt(d);
                    for (int j = 0; j < 3; j++) v0[i][j] *= v/d;
                    i++;
                }        
            } 
        }
    }
    // set center of mass velocity to 0
    double vcm[3] = {0.,0.,0.};
    for (int i = 0; i < N; i++) {
        for (int k = 0; k < 3; k++) {
            vcm[k] += v0[i][k];
        }
    }
    for (int k = 0; k < 3; k++) vcm[k] /= ((double) N);
    for (int i = 0; i < N; i++) {
        for (int k = 0; k < 3; k++) {
            v0[i][k] -= vcm[k];
        }
    }  
}


double maxwell_boltzmann(double m, double T, double v) {
    double f = sqrt(2./M_PI)*pow(m/(kB*T), 3./2.)*pow(v, 2.)*exp(-m*pow(v, 2.)/(2.*kB*T));
    return f;
}


double energy(double eps, double sigma, double m, int n, double l, double** r0, double** v0, double& T) {
    double ekin = 0.;
    double epot = 0.;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < 3; j++) ekin += 0.5*m*pow(v0[i][j], 2.);
    }
    T = 2.*ekin/(double (3*n-3))/kB;
    double dist_min, dist;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < i; j++) {
            dist_min = 1e10;
            for (int kx = -1; kx < 2; kx++) {
                for (int ky = -1; ky < 2; ky++) {
                    for (int kz = -1; kz < 2; kz++) {
                        dist = pow(r0[i][0]-r0[j][0]+l*kx, 2.)+pow(r0[i][1]-r0[j][1]+l*ky, 2.)+pow(r0[i][2]-r0[j][2]+l*kz, 2.);
                        dist = sqrt(dist);
                        if (dist < dist_min) dist_min = dist;
                    }
                }
            }
            epot += 4.*eps*(pow(sigma/dist_min, 12.)-pow(sigma/dist_min, 6.));
        }
    }
    double ene = ekin+epot;
    return ene;
}


double potential_energy(double eps, double sigma, double m, int n, double l, double** r0) {
    double epot = 0.;
    double dist_min, dist;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < i; j++) {
            dist_min=1e10;
            for (int kx = -1; kx < 2; kx++) {
                for (int ky = -1; ky < 2; ky++) {
                    for (int kz = -1; kz < 2; kz++) {
                        dist = pow(r0[i][0]-r0[j][0]+l*kx, 2.)+pow(r0[i][1]-r0[j][1]+l*ky, 2.)+pow(r0[i][2]-r0[j][2]+l*kz, 2.);
                        dist = sqrt(dist);
                        if (dist < dist_min) dist_min = dist;
                    }
                }
            }
            epot += 4.*eps*(pow(sigma/dist_min, 12.)-pow(sigma/dist_min, 6.));
        }
    }
    return epot;
}


double delta_potential_energy(double eps, double sigma, double m, int n, double l, int idx, double** r0, double** r1) {
    double depot=0.;
    double dist_min, dist;
    for (int j = 0; j < n; j++) {
        if (j != idx) {
            dist_min=1e10;
            for (int kx = -1; kx < 2; kx++) {
                for (int ky = -1; ky < 2; ky++) {
                    for (int kz = -1; kz < 2; kz++) {
                        dist = pow(r0[idx][0]-r0[j][0]+l*kx, 2.)+pow(r0[idx][1]-r0[j][1]+l*ky, 2.)+pow(r0[idx][2]-r0[j][2]+l*kz, 2.);
                        dist = sqrt(dist);
                        if (dist < dist_min) dist_min = dist;
                    }
                }
            }
            depot -= 2.*eps*(pow(sigma/dist_min, 12.)-pow(sigma/dist_min, 6.));
        }          
    }
    for (int j = 0; j < n; j++) {
        if (j != idx) {
            dist_min=1e10;
            for (int kx = -1; kx < 2; kx++) {
                for (int ky = -1; ky < 2; ky++) {
                    for (int kz = -1; kz < 2; kz++) {
                        dist = pow(r1[idx][0]-r1[j][0]+l*kx, 2.)+pow(r1[idx][1]-r1[j][1]+l*ky, 2.)+pow(r1[idx][2]-r1[j][2]+l*kz, 2.);
                        dist = sqrt(dist);
                        if (dist < dist_min) dist_min = dist;
                    }
                }
            }
            depot += 2.*eps*(pow(sigma/dist_min, 12.)-pow(sigma/dist_min, 6.));
        }          
    }
    return depot;
}


double drand() {
    double c;
    c = ((double) rand())/((double) RAND_MAX);
    return c;   
}