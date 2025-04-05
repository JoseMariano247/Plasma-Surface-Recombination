#include "Recombination_RK.h"
#include "Recombination_MC_real.h"

#include <fstream>
#include <ostream>
#include <iostream>

using namespace std;

int main()
{

    double O         = 1e5;
    double Fv        = 1.5e5;
    double Sv        = 3e3;
    double A2        = 0.0;
    double M         = 16e-3;
    double Tg        = 500;
    double Tw        = 200;
    double k1        = 1;
    double k3        = 1;
    double k4        = 1;
    double vd        = 1e15;
    double vD        = 1e13;
    double Ed        = 30e3;
    double ED        = 15e3;
    double Er        = 17.5e3;
    double ELHF      = 17.5e3;
    double tstop     = 1e-10;

    MonteCarloRecombinationReal(
        O, Fv, Sv, A2,
        M, Tg, Tw,
        k1, k3, k4, vd,
        vD, Ed, ED, Er, ELHF,
        tstop, "Real_Test_MC.txt");

    RungeKuttaRecombination(O, Fv, 0.0,
        Sv, 0.0,  A2,
        Tw, Tg, M,
        k1,  vd, Ed,
        k3, k4, Er, ELHF,
        vD, ED,
        tstop/10000.0, tstop,
        "Real_Test_RK.txt");

    return 0;
}
