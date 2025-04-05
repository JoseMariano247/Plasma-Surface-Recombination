#include "Plasma-Surface-Recombination.h"
#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include <fstream>
#include <string>
#include <cstdlib>
#include <functional>
#include <map>
#include <algorithm>
#include <set>
#include <stdio.h>

#ifndef pi
#define pi 3.14159265358979323846
#endif

using namespace std;

// Global variables (defined in the Plasma-Surface-Recombination.cpp file)
extern double initial_A;
extern double global_Tw;
extern double global_Tg;
extern double global_M;
extern bool general_params_extracted;
extern double vD;
extern double ED;
extern double Er;
extern double k_4;


vector<ReactionEvent> buildEventsForReaction(const string& reaction, 
                                             const vector<double>& rates,
                                             int& rateIndex, 
                                             const map<string,int>& speciesIndex,
                                             bool chemPresent,
                                             bool surfPresent)
{
    double kb = 1.380649e-23;  
    double Na = 6.023e23;      

    vector<ReactionEvent> result;
    auto idx = [&](const string& s){ return speciesIndex.at(s); };

    if (reaction == "Basic") {
        double kA = rates[rateIndex++];
        double kB = rates[rateIndex++];
        {
            ReactionEvent e;
            e.k = kA;
            e.propensity = prop_single(idx("A"));
            e.delta.resize(speciesIndex.size(), 0.0);
            e.delta[idx("A")] = -1.0;
            e.delta[idx("B")] = +1.0;
            result.push_back(e);
        }
        {
            ReactionEvent e;
            e.k = kB;
            e.propensity = prop_single(idx("B"));
            e.delta.resize(speciesIndex.size(), 0.0);
            e.delta[idx("B")] = -1.0;
            e.delta[idx("A")] = +1.0;
            result.push_back(e);
        }
    }
    else {
        // For non‑Basic reactions, extract general parameters only once.
        if (!general_params_extracted) {
            global_Tw = rates[rateIndex++];
            global_Tg = rates[rateIndex++];
            global_M  = rates[rateIndex++];
            general_params_extracted = true;
        }
        double v_med = std::sqrt((8 * kb * global_Tg * Na) / (pi * global_M));
        double phi_O = 0.25 * v_med * initial_A;

        if (reaction == "Physisorption") {
            // Extract 3 parameters: k_1, vd, Ed.
            double k_1 = rates[rateIndex++];
            double vd  = rates[rateIndex++];
            double Ed  = rates[rateIndex++];
            {
                ReactionEvent e;
                e.k = k_1 * phi_O;
                e.propensity = prop_bimolecular(idx("A"), idx("Fv"));
                e.delta.resize(speciesIndex.size(), 0.0);
                e.delta[idx("A")]  = -1.0;
                e.delta[idx("Fv")] = -1.0;
                e.delta[idx("Af")] = +1.0;
                result.push_back(e);
            }
            {
                ReactionEvent e;
                e.k = vd * std::exp(-Ed / (Na * kb * global_Tw));
                e.propensity = prop_single(idx("Af"));
                e.delta.resize(speciesIndex.size(), 0.0);
                e.delta[idx("Af")] = -1.0;
                e.delta[idx("A")]  = +1.0;
                e.delta[idx("Fv")] = +1.0;
                result.push_back(e);
            }
        }
        else if (reaction == "Chemisorption") {
            // Extract 3 parameters: k_3, k_4, Er.
            double k_3 = rates[rateIndex++];
            k_4 = rates[rateIndex++];
            Er  = rates[rateIndex++];
            double Pr  = k_4 * std::exp(-Er / (Na * kb * global_Tw));
            {
                ReactionEvent e;
                e.k = k_3 * phi_O;
                e.propensity = prop_bimolecular(idx("A"), idx("Sv"));
                e.delta.resize(speciesIndex.size(), 0.0);
                e.delta[idx("A")]  = -1.0;
                e.delta[idx("Sv")] = -1.0;
                e.delta[idx("As")] = +1.0;
                result.push_back(e);
            }
            {
                ReactionEvent e;
                e.k = Pr * k_3 * phi_O;
                e.propensity = prop_bimolecular(idx("A"), idx("As"));
                e.delta.resize(speciesIndex.size(), 0.0);
                e.delta[idx("A")]  = -1.0;
                e.delta[idx("As")] = -1.0;
                e.delta[idx("A2")] = +1.0;
                e.delta[idx("Sv")] = +1.0;
                result.push_back(e);
            }
        }
        else if (reaction == "Surface Diffusion") {
            // Extract 2 parameters: vD and ED.
            vD = rates[rateIndex++];
            ED = rates[rateIndex++];
            double tau_d_1 = vD * std::exp(-ED / (Na * kb * global_Tw));
            {
                ReactionEvent e;
                e.k = 0.75 * tau_d_1;
                e.propensity = prop_bimolecular(idx("Af"), idx("Sv"));
                e.delta.resize(speciesIndex.size(), 0.0);
                e.delta[idx("Af")] = -1.0;
                e.delta[idx("Sv")] = -1.0;
                e.delta[idx("Fv")] = +1.0;
                e.delta[idx("As")] = +1.0;
                result.push_back(e);
            }
        }
        else if (reaction == "Langmuir-Hinshelwood recombination") {
            // LH normally requires 5 parameters: vD, ED, k4, Er, ELHF.
            double vD_local = 0.0, ED_local = 0.0, k4_local = 1.0, Er_local = 0.0, ELHF_local = 0.0;
            double tau_d_1 = 1.0;
            if (!chemPresent && !surfPresent) {
                // LH alone: extract all 5.
                vD_local = rates[rateIndex++];
                ED_local = rates[rateIndex++];
                k4_local = rates[rateIndex++];
                Er_local = rates[rateIndex++];
                ELHF_local = rates[rateIndex++];
                tau_d_1 = vD_local * std::exp(-ED_local / (Na * kb * global_Tw));
            } else if (chemPresent && !surfPresent) {
                // Chemisorption present: LH extracts vD, ED, and ELHF.
                vD_local = rates[rateIndex++];
                ED_local = rates[rateIndex++];
                ELHF_local = rates[rateIndex++];
                tau_d_1 = vD_local * std::exp(-ED_local / (Na * kb * global_Tw));
                // k4_local and Er_local assumed provided by Chemisorption.
            } else if (!chemPresent && surfPresent) {
                // Surface Diffusion present: LH extracts k4, Er, and ELHF.
                k4_local = rates[rateIndex++];
                Er_local = rates[rateIndex++];
                ELHF_local = rates[rateIndex++];
                // tau_d_1 assumed provided by Surface Diffusion.
            } else if (chemPresent && surfPresent) {
                // Both present: LH extracts only ELHF.
                ELHF_local = rates[rateIndex++];
                tau_d_1 = vD * std::exp(-ED / (Na * kb * global_Tw));
                k4_local = k_4;
                Er_local = Er;
            }
            double Pr   = k4_local * std::exp(-Er_local / (Na * kb * global_Tw));
            double Prlh = k4_local * std::exp(-ELHF_local / (Na * kb * global_Tw));
            {
                ReactionEvent e;
                e.k = tau_d_1 * Pr;
                e.propensity = prop_bimolecular(idx("Af"), idx("As"));
                e.delta.resize(speciesIndex.size(), 0.0);
                e.delta[idx("Af")] = -1.0;
                e.delta[idx("As")] = -1.0;
                e.delta[idx("A2")] = +1.0;
                e.delta[idx("Sv")] = +1.0;
                e.delta[idx("Fv")] = +1.0;
                result.push_back(e);
            }
            {
                ReactionEvent e;
                e.k = tau_d_1 * Prlh;
                e.propensity = prop_square(idx("Af"));
                e.delta.resize(speciesIndex.size(), 0.0);
                e.delta[idx("Af")] = -2.0;
                e.delta[idx("A2")] = +1.0;
                e.delta[idx("Fv")] = +2.0;
                result.push_back(e);
            }
        }
    }
    return result;
}


// Function that generates a progress bar.
void printProgressBar(double progress, double total) 
{
   int barWidth = 50;
   float percent = static_cast<float>(progress) / total;
   int filled = percent * barWidth;
 
   std::cout << "\r[";
   for (int i = 0; i <= barWidth; i++) {
       if (i < filled)
           std::cout << "=";
       else if (i == filled)
           std::cout << ">";
       else
           std::cout << " ";
   }
   std::cout << "] " << int(percent * 100.0) << "%" << std::flush;
}


// Function that runs the Monte Carlo simulation, storing both the populations
// and the instantaneous propensities (Big R values) at each time step.
void simulateMultiReaction(double t_stop, 
                           const vector<ReactionEvent>& events, 
                           vector<double>& state,
                           const vector<string>& speciesList,
                           const string& outputFilename)
{
    double t = 0.0;
    vector<double> times { t };
    vector<vector<double>> states { state };
    vector<vector<double>> propHistory; // To store propensities for each event at each time step.
    
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> dis(0.0, 1.0);

    printProgressBar(0.0, t_stop);

    while (t < t_stop) {
        double total_rate = 0.0;
        vector<double> rvec;
        rvec.reserve(events.size());
        for (auto &evt : events) {
            double r = evt.propensity(state, evt.k);
            rvec.push_back(r);
            total_rate += r;
        }

        if (total_rate <= 1e-15)
            break;

        double r1 = dis(gen);
        double dt = -log(r1) / total_rate;
        t += dt;
        if (t > t_stop)
            break;

        double r2 = dis(gen) * total_rate;
        double cum = 0.0;
        int chosen = -1;
        for (int i = 0; i < static_cast<int>(events.size()); i++) {
            cum += rvec[i];
            if (cum >= r2) {
                chosen = i;
                break;
            }
        }
        if (chosen < 0)
            break;
        
        for (int i = 0; i < static_cast<int>(state.size()); i++) {
            state[i] += events[chosen].delta[i];
            if (state[i] < 0)
                state[i] = 0;
        }

        times.push_back(t);
        states.push_back(state);
        propHistory.push_back(rvec);  // Store the propensities computed at this time step.

        printProgressBar(t, t_stop);
    }
    cout << "\n";

    ofstream outFile(outputFilename);
    if (!outFile) {
        cerr << "Error opening file: " << outputFilename << "\n";
        return;
    }
    // Write header: time, populations, then propensities.
    outFile << "Time";
    for (auto &s : speciesList) {
        outFile << "\tPopulation " << s;
    }
    for (size_t i = 0; i < events.size(); i++) {
        outFile << "\tR" << i+1;
    }
    outFile << "\n";
    
    // Write out each time step.
    for (size_t i = 0; i < times.size(); i++) {
        outFile << times[i];
        // Write species populations.
        for (size_t j = 0; j < speciesList.size(); j++) {
            outFile << "\t" << states[i][j];
        }
        // For the initial time step, we have no propensity values.
        if (i == 0) {
            for (size_t k = 0; k < events.size(); k++)
                outFile << "\t0";
        } else {
            int propIndex = i - 1; // propHistory is one element shorter.
            for (size_t k = 0; k < events.size(); k++) {
                outFile << "\t" << propHistory[propIndex][k];
            }
        }
        outFile << "\n";
    }
    
    outFile.close();
    cout << "Simulation complete. Output written to " << outputFilename << "\n";
}

