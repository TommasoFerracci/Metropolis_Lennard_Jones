// used to compute the pair correlation function

#include <iostream>
#include <string>
#include <fstream>
#include <stdlib.h>
#include <cmath>

double gaussian(double x, double x0, double sigma);

int main(int argc, const char* argv[])
{
    int N; // number of atoms
    double l; // simulation cell length
    double *r0mat;
    double **r0;
    double rmax, sigma;
    double *gr, *xr;
    double dist, dist_min;
    int ns;
    
    std::ifstream flast;
    std::ofstream fgr;
    
    flast.open("argonlast.txt", std::ios::in);
    flast >> N >> l;
    std::cout << "Atomi: " << N << "  Lato (A): " << l << '\n';

    r0mat = new double[N*3];
    r0 = new double*[N];
    for (int i = 0; i < N; i++) {
        r0[i] = &r0mat[i*3];
    }
    for (int i = 0; i < N; i++) {
        for (int k = 0;k < 3; k++) flast >> r0[i][k];
    }
    flast.close();
    
    std::cout << "rmax (A): \n";
    std::cin >> rmax;
    std::cout << "Guassian broadening (A): \n";
    std::cin >> sigma;
    std::cout << "Number of intervals: \n";
    std::cin >> ns;
    
    gr = new double[ns];
    xr = new double[ns];
    for (int i =0 ; i < ns; i++) gr[i] = 0.;
    for (int i = 0; i < ns; i++) xr[i] = rmax / ((double) ns) * ((double) i);
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            if (j != i) {
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
                for (int k = 0; k < ns; k++) gr[k] += gaussian(xr[k], dist_min, sigma);
            }
        }
    }
    
    for (int i = 0; i < ns; i++) gr[i] *= (1./(4.*M_PI*pow(xr[i]*N, 2.)))*pow(l, 3.);
    
    fgr.open("correlation.txt", std::ios::out);
    for (int i = 0; i < ns; i++) {
        fgr << xr[i] << '\t'  << gr[i] << '\n';
    }
    fgr.close();
    
    //deallocation
    delete [] r0mat;
    delete [] r0;
    delete [] gr;
    delete [] xr;    
}

double gaussian(double x, double x0, double sigma) {
    double g = (1./(sigma*sqrt(2.*M_PI)))*exp(-pow((x-x0)/sigma, 2.)/2.);
    return g;
}