#pragma once

#include <string>

void RungeKuttaRecombination(double A,  double Fv, double Af,
    double Sv, double As, double A2,
    double Tw, double Tg, double M,
    double k1, double vd, double Ed,
    double k3, double k4, double Er, double ELHF,
    double vD, double ED,
    double dt, double tMax,
    const std::string& outputFilename);