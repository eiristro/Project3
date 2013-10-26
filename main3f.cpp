//packages
#include <iostream>
#include <armadillo>
#include <cmath>
#include <lib.cpp>
#include <math.h>
#include <iomanip>
#include <string>

//namespace
using namespace std;
using namespace arma;

//global constants
double G = 6.67384e-11;         //m^3 * kg^-1 * s^-2
double pi = 3.1415926535;
double Msun = 332946;           //earth masses
double r0 = 1.49597871e11;      //1 AU in meters
double m0 = 5.97219e24;         //earth mass in kg
double t0 = 31536000;           // year in seconds
double G0 = G * pow(t0,2)*m0/pow(r0,3); // AU^3 * M_earth^-1 * year^-2
double timestep = 1e-3;
double maxT = 300;           //years
int N = maxT / timestep;   //number of steps
int RunID = 0;


//creating planets
class planets{
    public:
        double mass, xpos, ypos, vx, vy, fx, fy;
        int N;
        const char* name;
        ofstream oFile;
        bool FileOpen;

        planets(const char *iname,double imass,double iradius,double ispeed);
        void force(planets** aSystem, int iObject, int iObjects);
        void correct(planets** aSystem, int iObjects);
        void Init(planets** aSystem, int iObjects, int iSun);
        void output();
};

planets::planets(const char* iname, double imass, double iradius, double ispeed){
    mass = imass/m0;
    xpos = iradius/r0;
    ypos = 0;
    vx = 0;
    vy = ispeed*(t0/r0);
    name = iname;

    FileOpen = false;
}

//calculating force
void planets::force(planets** aSystem, int iObject, int iObjects){
    double r;
    int i;

    fx = 0.0;
    fy = 0.0;
    for (i = 0; i<iObjects; i++) {
        if (i != iObject) {
            r = sqrt(pow(xpos - aSystem[i]->xpos, 2) + pow(ypos - aSystem[i]->ypos, 2));
            fx -= G0 * aSystem[i]->mass * mass * (xpos - aSystem[i]->xpos)/pow(r, 3);
            fy -= G0 * aSystem[i]->mass * mass * (ypos - aSystem[i]->ypos)/pow(r, 3);
        }
    }
}

void planets::correct(planets** aSystem, int iObjects) {
    /* Correcting the positions of the planets so that the center
       of mass is at coordinates (0, 0) */
    int i;
    double CMx = 0.0, CMy = 0.0, TMass = 0.0;

    for (i = 0; i < iObjects; i++) {
        CMx += aSystem[i]->xpos*aSystem[i]->mass;
        CMy += aSystem[i]->ypos*aSystem[i]->mass;
        TMass += aSystem[i]->mass;
    }
    CMx /= TMass;
    CMy /= TMass;

    for (i = 0; i < iObjects; i++) {
        aSystem[i]->xpos -= CMx;
        aSystem[i]->ypos -= CMy;
    }
}

void planets::Init(planets** aSystem, int iObjects, int iSun) {
    int i;
    double Px = 0.0, Py = 0.0;

    for (i = 0; i < iObjects; i++) {
        if (i != iSun) {
            Px += aSystem[i]->vx * aSystem[i]->mass;
            Py += aSystem[i]->vy * aSystem[i]->mass;
        }
    }

    aSystem[iSun]->vx = -Px / aSystem[iSun]->mass;
    aSystem[iSun]->vy = -Py / aSystem[iSun]->mass;
}

void planets::output(){
    int iPrecision=15;

    if(!FileOpen) {
        ostringstream sFileName;
        sFileName << "../Project3/" << RunID << "_" << name << ".dat";
        oFile.open(sFileName.str().c_str());
        FileOpen = true;
    }

    oFile << setprecision(iPrecision);
    oFile << setw(iPrecision+7) << xpos;
    oFile << setw(iPrecision+7) << ypos;
    oFile << endl;

}

int main()
{
    int iObjects = 10;
    int i, j, p = 0;

    planets** aSystem = new planets*[iObjects];
    //                       Name          Mass      Radius     Speed
    aSystem[p] = new planets("Sun",       1.9891e30,   0,          0      ); p++;
    aSystem[p] = new planets("Mercury",   3.3022e23,   69.82e9,    38.86e3); p++;
    aSystem[p] = new planets("Venus",     4.8685e24,   108.94e9,   34.79e3); p++;
    aSystem[p] = new planets("Earth",     5.9736e24,   152.10e9,   29.29e3); p++;
    aSystem[p] = new planets("Mars",      6.4185e23,   249.23e9,   21.97e3); p++;
    aSystem[p] = new planets("Jupiter",   1.8986e27,   816.62e9,   12.44e3); p++;
    aSystem[p] = new planets("Saturn",    5.6846e26,   1514.50e9,   9.09e3); p++;
    aSystem[p] = new planets("Uranus",    8.6810e25,   3003.62e9,   6.49e3); p++;
    aSystem[p] = new planets("Neptune",   1.0243e26,   4545.67e9,   5.37e3); p++;
    aSystem[p] = new planets("Pluto",     1.305e22,    7375.93e9,   3.71e3); p++;

    aSystem[0]->Init(aSystem, iObjects, 0);
    for (i = 0; i<iObjects; i++) {
        aSystem[i]->force(aSystem, i, iObjects);
        aSystem[i]->output();
    }

    double *axpos = new double[iObjects];
    double *aypos = new double[iObjects];
    double *avx = new double[4*iObjects];
    double *avy = new double[4*iObjects];
    double *aax = new double[4*iObjects];
    double *aay = new double[4*iObjects];
    double Step;
    double Time = 0.0;

    while(Time < maxT) {
        Time += timestep;
        for (i = 0; i<iObjects; i++) {
            axpos[i] = aSystem[i]->xpos;
            aypos[i] = aSystem[i]->ypos;
            avx[4*i] = aSystem[i]->vx;
            avy[4*i] = aSystem[i]->vy;
            aax[4*i] = aSystem[i]->fx/aSystem[i]->mass;
            aay[4*i] = aSystem[i]->fy/aSystem[i]->mass;
        }
        // Runge Kutta steps 1-3
        for (j = 0; j<3; j++) {
            for (i = 0; i<iObjects; i++) {
                if (j<2) {Step = timestep/2.0;} else {Step = timestep;};
                avx[i*4+j+1] = avx[i*4] + aax[i*4+j]*Step;
                avy[i*4+j+1] = avy[i*4] + aay[i*4+j]*Step;
                aSystem[i]->xpos = axpos[i] + avx[i*4+j]*Step;
                aSystem[i]->ypos = aypos[i] + avy[i*4+j]*Step;
            }
            for(i=0; i<iObjects; i++) {
                aSystem[i]->force(aSystem, i, iObjects);
                aax[i*4+j+1] = aSystem[i]->fx/aSystem[i]->mass;
                aay[i*4+j+1] = aSystem[i]->fy/aSystem[i]->mass;
            }

        }
        // Runge Kutta step 4
        for (i=0; i<iObjects; i++) {
            aSystem[i]->vx   = avx[i*4] + timestep*(aax[i*4] + 2.0*aax[i*4+1] + 2.0*aax[i*4+2] + aax[i*4+3])/6.0;
            aSystem[i]->vy   = avy[i*4] + timestep*(aay[i*4] + 2.0*aay[i*4+1] + 2.0*aay[i*4+2] + aay[i*4+3])/6.0;
            aSystem[i]->xpos = axpos[i] + timestep*(avx[i*4] + 2.0*avx[i*4+1] + 2.0*avx[i*4+2] + avx[i*4+3])/6.0;
            aSystem[i]->ypos = aypos[i] + timestep*(avy[i*4] + 2.0*avy[i*4+1] + 2.0*avy[i*4+2] + avy[i*4+3])/6.0;
        }
        // Updating system
        aSystem[0]->correct(aSystem, iObjects);
        for(i=0; i<iObjects; i++) {
            aSystem[i]->force(aSystem, i, iObjects);
            aSystem[i]->output();
        }

    }

    return 0;
}

