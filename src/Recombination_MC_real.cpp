#include "Recombination_MC_real.h"
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

std::vector<double> MonteCarloRecombinationReal(double initial_A, double initial_Fv,
    double initial_Sv, double initial_A2, double M, double Tg, double Tw,
    double k1, double k3, double k4, double vd,
    double vD, double Ed, double ED, double Er, double ELHF,
    double t_stop, const std::string& outputFilename)
    
{
    double gamma_ER    = 0.0;
    double gamma_LHS   = 0.0;
    double gamma_LHF   = 0.0;
    double gamma_total = 0.0;

    double kb = 1.380649e-23;
    double Na = 6.023e23;

    double v_med = std::sqrt((8*kb*Tg*Na)/(pi * M));

    double phi_O = 0.25 * v_med * initial_A;
    double Pr = k4 * std::exp(-Er/(Na*kb*Tw));
    double Prlh = k4 * std::exp(-ELHF/(Na*kb*Tw));
    double tau_d_1 = vD * std::exp(-ED/(Na*kb*Tw));

    //calculate reaction coefficents

    double r1 = k1 * phi_O;
    double r2 = vd * std::exp(-Ed/(Na*kb*Tw));
    double r3 = k3 * phi_O;
    double r4 = Pr * r3;
    double r5 = 0.75*tau_d_1;
    double r6 = tau_d_1 * Pr;
    double r7 = tau_d_1 * Prlh;

    double t = 0.0;
    double A = initial_A;
    double Fv = initial_Fv;
    double Af = 0.0;
    double Sv = initial_Sv;
    double As = 0.0;
    double A2 = initial_A2;

    const double S = Sv;
    const double F = Fv;

    std::vector<double> times{ t };
    std::vector<double> populations_A{ A };
    std::vector<double> populations_Fv{ Fv };
    std::vector<double> populations_Af{ Af };
    std::vector<double> populations_Sv{ Sv };
    std::vector<double> populations_As{ As };
    std::vector<double> populations_A2{ A2 };

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.0, 1.0);

    std::ofstream outFile(outputFilename);
    if (!outFile)
    {
        std::cerr << "Error opening file: " << outputFilename << std::endl;
    }
    
    outFile << "Time\tA\tFv\tAf\tSv\tAs\tA2\tR1\tR2\tR3\tR4\tR5\tR6\tR7\n";

    while (t < t_stop) {
        double R1 = r1 * A * Fv;      // A + Fv -> Af
        double R2 = r2 * Af;          // Af -> A + Fv
        double R3 = r3 * A * Sv;      // A + Sv -> As
        double R4 = r4 * A * As;      // A + As -> A2 + Sv
        double R5 = r5 * Af * Sv;     // Af + Sv -> Fv + As
        double R6 = r6 * Af * As;     // Af + As -> A2 + Sv + Fv
        double R7 = (Af >= 2) ? r7 * Af * Af : 0;  // Af + Af -> A2 + 2 Fv

        double totalRate = R1 + R2 + R3 + R4 + R5 + R6 + R7;
        if (totalRate <= 0) break;

        double r_time = dis(gen);
        double dt = -std::log(r_time) / totalRate;
        t += dt;

        double r_choice = dis(gen) * totalRate;
        double cumulative = 0.0;
        int reaction = -1;

        if ((cumulative += R1) >= r_choice)
            reaction = 1;
        else if ((cumulative += R2) >= r_choice)
            reaction = 2;
        else if ((cumulative += R3) >= r_choice)
            reaction = 3;
        else if ((cumulative += R4) >= r_choice)
            reaction = 4;
        else if ((cumulative += R5) >= r_choice)
            reaction = 5;
        else if ((cumulative += R6) >= r_choice)
            reaction = 6;
        else if ((cumulative += R7) >= r_choice)
            reaction = 7;

        // Update species counts based on the chosen reaction:
        switch (reaction)
        {
            case 1: // A + Fv -> Af
                //if (A > 0 && Fv > 0) { A--; Fv--; Af++; }
                if (A > 0 && Fv > 0) {Fv--; Af++; }
                break;
            case 2: // Af -> A + Fv
                //if (Af > 0) { Af--; A++; Fv++; }
                if (Af > 0) { Af--; Fv++; }
                break;
            case 3: // A + Sv -> As
                //if (A > 0 && Sv > 0) { A--; Sv--; As++; }
                if (A > 0 && Sv > 0) { Sv--; As++; }
                break;
            case 4: // A + As -> A2 + Sv
                //if (A > 0 && As > 0) { A--; As--; A2++; Sv++; }
                if (A > 0 && As > 0) { As--; A2++; Sv++; }
                break;
            case 5: // Af + Sv -> Fv + As
                if (Af > 0 && Sv > 0) { Af--; Sv--; Fv++; As++; }
                break;
            case 6: // Af + As -> A2 + Sv + Fv
                if (Af > 0 && As > 0) { Af--; As--; A2++; Sv++; Fv++; }
                break;
            case 7: // Af + Af -> A2 + 2 Fv
                if (Af >= 2) { Af -= 2; A2++; Fv += 2; }
                break;
            default:
                break;
        }

        // Record the updated state:
        times.push_back(t);
        populations_A.push_back(A);
        populations_Fv.push_back(Fv);
        populations_Af.push_back(Af);
        populations_Sv.push_back(Sv);
        populations_As.push_back(As);
        populations_A2.push_back(A2);


        // Print current state to console
        std::cout << t << "\t" << A << "\t" << Fv << "\t" << Af << "\t"
                  << Sv << "\t" << As << "\t" << A2  << "\n";

        outFile << t << "\t" << A << "\t" << Fv << "\t" << Af << "\t"
        << Sv << "\t" << As << "\t" << A2 << "\t" << R1 << "\t" << R2 << "\t" << R3 << "\t"
        << R4 << "\t" << R5 << "\t" << R6 << "\t" << R7 << "\n";
        
    }


    outFile.close();

    gamma_ER    = 2 * r4 * As * S / (phi_O * (S + F));
    gamma_LHS   = 2 * r6 * As * Af * S / (phi_O * (S + F));
    gamma_LHF   = 2 * r7 * Af * Af * F / (phi_O * (S + F));
    gamma_total = gamma_ER + gamma_LHF + gamma_LHS;

    return {Tw, gamma_ER, gamma_LHS, gamma_LHF, gamma_total};

}
