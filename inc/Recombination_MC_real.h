#pragma once

#include <string>
#include <vector>

std::vector<double> MonteCarloRecombinationReal(double initial_A, double initial_Fv,
    double initial_Sv, double initial_A2, double M, double Tg, double Tw,
    double k1, double k3, double k4, double vd,
    double vD, double Ed, double ED, double Er, double ELHF,
    double t_stop, const std::string& outputFilename);