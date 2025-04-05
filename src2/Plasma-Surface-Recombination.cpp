/*
    Welcome to the code to perform the Monte Carlo simulation
    of the Plasma Surface recombination. This file uses a
    auxiliary .cpp file with all the necessary function and
    a header file.
*/

// Necessary libraries/header files
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


using namespace std;

// Global variable for initial concentration of A.
double initial_A = 0.0;

// Global general parameters (for non‑Basic reactions).
double global_Tw = 0.0, global_Tg = 0.0, global_M = 0.0;
bool general_params_extracted = false;

double vD;
double ED;
double Er;
double k_4;
double vD_sd;
double ED_sd;
double Er_ch;
double k4_ch; 

// This map specifies what type of species each reaction needs
static const map<string, vector<string>> reactionSpecies = {
    {"Basic", {"A", "B"}},
    {"Physisorption", {"A", "Fv", "Af"}},
    {"Chemisorption", {"A", "Sv", "As", "A2"}},
    {"Surface Diffusion", {"Af", "Sv", "Fv", "As"}},
    {"Langmuir-Hinshelwood recombination", {"Af", "As", "Fv", "Sv", "A2"}}
};

// This map specifies how many parameters each reaction requires
static const map<string, int> reactionRateCount = {
    {"Basic", 2},
    {"Physisorption", 3},
    {"Chemisorption", 3},
    {"Surface Diffusion", 2},
    {"Langmuir-Hinshelwood recombination", 5}
};


#include "Plasma-Surface-Recombination.h"
#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <set>
#include <algorithm>
#include <cmath>

#ifndef pi
#define pi 3.14159265358979323846
#endif

using namespace std;

int main(int argc, char* argv[]) {
    if (argc < 2) {
        cerr << "No arguments provided.\n";
        return 1;
    }

    // Parse reaction names until a numeric token is encountered.
    vector<string> reactions;
    int argIndex = 1;
    while (argIndex < argc) {
        string token = argv[argIndex];
        try {
            stod(token);
            break;
        } catch (...) {
            reactions.push_back(token);
            argIndex++;
        }
    }
    if (reactions.empty()) {
        cerr << "No reactions specified.\n";
        return 1;
    }

    // Compute nominal needed rate constants.
    int neededRateConstants = 0;
    for (auto &r : reactions) {
        auto it = reactionRateCount.find(r);
        if (it == reactionRateCount.end()) {
            cerr << "Unknown reaction: " << r << "\n";
            return 1;
        }
        neededRateConstants += it->second;
    }
    // Adjust for Langmuir-Hinshelwood recombination if shared.
    if (find(reactions.begin(), reactions.end(), "Langmuir-Hinshelwood recombination") != reactions.end()) {
        if ((find(reactions.begin(), reactions.end(), "Chemisorption") != reactions.end()) &&
            (find(reactions.begin(), reactions.end(), "Surface Diffusion") != reactions.end())) {
            neededRateConstants = neededRateConstants - reactionRateCount.at("Langmuir-Hinshelwood recombination") + 1;
        } else if (find(reactions.begin(), reactions.end(), "Chemisorption") != reactions.end()) {
            neededRateConstants = neededRateConstants - reactionRateCount.at("Langmuir-Hinshelwood recombination") + 3;
        } else if (find(reactions.begin(), reactions.end(), "Surface Diffusion") != reactions.end()) {
            neededRateConstants = neededRateConstants - reactionRateCount.at("Langmuir-Hinshelwood recombination") + 3;
        }
    }

    // Build union of species.
    set<string> usedSpecies;
    for (auto &r : reactions) {
        auto it = reactionSpecies.find(r);
        if (it == reactionSpecies.end()) {
            cerr << "No species info for reaction: " << r << "\n";
            return 1;
        }
        for (auto &s : it->second)
            usedSpecies.insert(s);
    }
    vector<string> allSpecies;
    vector<string> fixedOrder = {"A", "B", "Af", "As", "Fv", "Sv", "A2"};
    for (const auto &s : fixedOrder) {
        if (usedSpecies.find(s) != usedSpecies.end())
            allSpecies.push_back(s);
    }

    int totalNumericNeeded = (reactions.size() == 1 && reactions[0] == "Basic") ? 
                              (2 + allSpecies.size() + 1) : (3 + neededRateConstants + allSpecies.size() + 1);
    int numericAvailable = argc - argIndex;
    if (numericAvailable < totalNumericNeeded) {
        cerr << "Not enough numeric parameters provided. Expected " << totalNumericNeeded
             << ", got " << numericAvailable << ".\n";
        return 1;
    }

    vector<double> rates;
    rates.reserve(neededRateConstants);
    if (reactions.size() == 1 && reactions[0] == "Basic") {
        for (int i = 0; i < 2; i++) {
            rates.push_back(stod(argv[argIndex++]));
        }
    } else {
        for (int i = 0; i < 3; i++) { // General parameters: Tw, Tg, M.
            rates.push_back(stod(argv[argIndex++]));
        
        }
        for (int i = 0; i < neededRateConstants; i++) {
            rates.push_back(stod(argv[argIndex++]));
        }
    }


    vector<double> initState(allSpecies.size(), 0.0);
    for (int i = 0; i < (int)allSpecies.size(); i++) {
        initState[i] = stod(argv[argIndex++]);
    }

    double t_stop = stod(argv[argIndex++]);

    // Build mapping from species name to index.
    map<string,int> speciesMap;
    for (int i = 0; i < (int)allSpecies.size(); i++) {
        speciesMap[allSpecies[i]] = i;
    }

    if (speciesMap.find("A") != speciesMap.end()) {
        initial_A = initState[speciesMap["A"]];
    } else {
        cerr << "Species A is not in the union; cannot define initial_A.\n";
        return 1;
    }

    bool chemPresent = (find(reactions.begin(), reactions.end(), "Chemisorption") != reactions.end());
    bool surfPresent = (find(reactions.begin(), reactions.end(), "Surface Diffusion") != reactions.end());

    vector<ReactionEvent> events_MC;
    int ratePos = 0;
    for (auto &r : reactions) {
        vector<ReactionEvent> these = buildEventsForReaction(r, rates, ratePos, speciesMap, chemPresent, surfPresent);
        events_MC.insert(events_MC.end(), these.begin(), these.end());
    }

    // Run Monte Carlo simulation.
    string outputFilename_MC = "output.txt";
    vector<double> initStatecopy = initState;
    simulateMultiReaction(t_stop, events_MC, initState, allSpecies, outputFilename_MC);

    return 0;
}

