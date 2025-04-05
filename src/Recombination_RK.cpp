#include "Recombination_RK.h"
#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include <fstream>
#include <ostream>

#ifndef pi
#define pi 3.14159265358979323846
#endif

using namespace std;

static void derivatives6(double r1, double r2, double r3, double r4,
    double r5, double r6, double r7,
    double A,  double Fv, double Af,
    double Sv, double As, double A2,
    double& dAdt,  double& dFvdt, double& dAfdt,
    double& dSvdt, double& dAsdt, double& dA2dt)
{
// Compute instantaneous reaction rates:
double R1 = r1 * A * Fv;
double R2 = r2 * Af;
double R3 = r3 * A * Sv;
double R4 = r4 * A * As;
double R5 = r5 * Af * Sv;
double R6 = r6 * Af * As;
double R7 = r7 * Af * Af;  // Reaction 7: 2 Af -> A2 + 2Fv

// ODEs for the species concentrations:
dAdt  = -R1 + R2 - R3 - R4;
dFvdt = -R1 + R2 + R5 + R6 + 2.0 * R7;
dAfdt =  R1 - R2 - R5 - R6 - 2.0 * R7;
dSvdt = -R3 + R4 - R5 + R6;
dAsdt =  R3 - R4 + R5 - R6;
dA2dt =  R4 + R6 + R7;
}

static void rk4Step6(double r1, double r2, double r3, double r4,
double r5, double r6, double r7,
double& A,  double& Fv, double& Af,
double& Sv, double& As, double& A2,
double& t,  double dt)
{
double dA1, dFv1, dAf1, dSv1, dAs1, dA21;
derivatives6(r1, r2, r3, r4, r5, r6, r7,
A, Fv, Af, Sv, As, A2,
dA1, dFv1, dAf1, dSv1, dAs1, dA21);

double dA2_, dFv2_, dAf2_, dSv2_, dAs2_, dA22_;
{
double A_temp = A + 0.5 * dt * dA1;
double Fv_temp = Fv + 0.5 * dt * dFv1;
double Af_temp = Af + 0.5 * dt * dAf1;
double Sv_temp = Sv + 0.5 * dt * dSv1;
double As_temp = As + 0.5 * dt * dAs1;
double A2_temp = A2 + 0.5 * dt * dA21;
derivatives6(r1, r2, r3, r4, r5, r6, r7,
A_temp, Fv_temp, Af_temp, Sv_temp, As_temp, A2_temp,
dA2_, dFv2_, dAf2_, dSv2_, dAs2_, dA22_);
}

double dA3, dFv3, dAf3, dSv3, dAs3, dA23;
{
double A_temp = A + 0.5 * dt * dA2_;
double Fv_temp = Fv + 0.5 * dt * dFv2_;
double Af_temp = Af + 0.5 * dt * dAf2_;
double Sv_temp = Sv + 0.5 * dt * dSv2_;
double As_temp = As + 0.5 * dt * dAs2_;
double A2_temp = A2 + 0.5 * dt * dA22_;
derivatives6(r1, r2, r3, r4, r5, r6, r7,
A_temp, Fv_temp, Af_temp, Sv_temp, As_temp, A2_temp,
dA3, dFv3, dAf3, dSv3, dAs3, dA23);
}

double dA4, dFv4, dAf4, dSv4, dAs4, dA24;
{
double A_temp = A + dt * dA3;
double Fv_temp = Fv + dt * dFv3;
double Af_temp = Af + dt * dAf3;
double Sv_temp = Sv + dt * dSv3;
double As_temp = As + dt * dAs3;
double A2_temp = A2 + dt * dA23;
derivatives6(r1, r2, r3, r4, r5, r6, r7,
A_temp, Fv_temp, Af_temp, Sv_temp, As_temp, A2_temp,
dA4, dFv4, dAf4, dSv4, dAs4, dA24);
}

A  += (dt / 6.0) * (dA1 + 2.0 * dA2_ + 2.0 * dA3 + dA4);
Fv += (dt / 6.0) * (dFv1 + 2.0 * dFv2_ + 2.0 * dFv3 + dFv4);
Af += (dt / 6.0) * (dAf1 + 2.0 * dAf2_ + 2.0 * dAf3 + dAf4);
Sv += (dt / 6.0) * (dSv1 + 2.0 * dSv2_ + 2.0 * dSv3 + dSv4);
As += (dt / 6.0) * (dAs1 + 2.0 * dAs2_ + 2.0 * dAs3 + dAs4);
A2 += (dt / 6.0) * (dA21 + 2.0 * dA22_ + 2.0 * dA23 + dA24);

t += dt;
}

void RungeKuttaRecombination(double A,  double Fv, double Af,
                double Sv, double As, double A2,
                double Tw, double Tg, double M,
                double k1, double vd, double Ed,
                double k3, double k4, double Er, double ELHF,
                double vD, double ED,
                double dt, double tMax,
                const std::string& outputFilename)
{
double kb = 1.380649e-23;
double Na = 6.023e23;

double v_med = std::sqrt((8 * kb * Tg * Na) / (pi * M));
double phi_O = 0.25 * v_med * A;  // Use initial A

double r1 = k1 * phi_O;
double r2 = vd * std::exp(-Ed / (Na * kb * Tw));

double r3 = k3 * phi_O;
double r4 = k4 * std::exp(-Er / (Na * kb * Tw)) * r3;

double tau_d = vD * std::exp(-ED / (Na * kb * Tw));
double r5 = 0.75 * tau_d;

double r6 = tau_d * (k4 * std::exp(-Er / (Na * kb * Tw)));
double r7 = tau_d * (k4 * std::exp(-ELHF / (Na * kb * Tw)));

cout << "Computed reaction rates:" << "\n";
cout << "r1 = " << r1 << "\n";
cout << "r2 = " << r2 << "\n";
cout << "r3 = " << r3 << "\n";
cout << "r4 = " << r4 << "\n";
cout << "r5 = " << r5 << "\n";
cout << "r6 = " << r6 << "\n";
cout << "r7 = " << r7 << "\n";

double t = 0.0;
ofstream outFile(outputFilename);
if (!outFile) {
cerr << "Error opening file: " << outputFilename << "\n";
return;
}

// Write header (populations and reaction rates)
outFile << "Time\tA\tFv\tAf\tSv\tAs\tA2\tR1\tR2\tR3\tR4\tR5\tR6\tR7\n";

outFile << t << "\t" << A << "\t" << Fv << "\t" << Af << "\t" << Sv << "\t" << As << "\t" << A2;


while (t < tMax) {
rk4Step6(r1, r2, r3, r4, r5, r6, r7, A, Fv, Af, Sv, As, A2, t, dt);

double curr_R1 = r1 * A * Fv;
double curr_R2 = r2 * Af;
double curr_R3 = r3 * A * Sv;
double curr_R4 = r4 * A * As;
double curr_R5 = r5 * Af * Sv;
double curr_R6 = r6 * Af * As;
double curr_R7 = (Af >= 2) ? r7 * Af * Af : 0.0;

outFile << t << "\t" << A << "\t" << Fv << "\t" << Af << "\t" 
<< Sv << "\t" << As << "\t" << A2;
outFile << "\t" << curr_R1 << "\t" << curr_R2 << "\t" << curr_R3 
<< "\t" << curr_R4 << "\t" << curr_R5 << "\t" << curr_R6 
<< "\t" << curr_R7 << "\n";
}

outFile.close();

}

