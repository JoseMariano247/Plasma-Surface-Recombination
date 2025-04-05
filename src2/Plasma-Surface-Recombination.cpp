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


int main(int argc, char* argv[]) {

    // In case no arguments are provided
    if (argc < 2) {
        cerr << "No arguments provided.\n";
        return 1;
    }  
    
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

    // In case no reactions are selected
    if (reactions.empty()) {
        cerr << "No reactions specified.\n";
        return 1;
    }


    int neededRateConstants = 0;
    for (auto &r : reactions) {
        auto it = reactionRateCount.find(r);

        // In case some reaction given in unknown
        if (it == reactionRateCount.end()) {
            cerr << "Unknown reaction: " << r << "\n";
            return 1;
        }
        neededRateConstants += it->second;
    }

    if (find(reactions.begin(), reactions.end(), "Langmuir-Hinshelwood recombination") != reactions.end()) {

        // If LH is present:
        if ( (find(reactions.begin(), reactions.end(), "Chemisorption") != reactions.end()) &&
             (find(reactions.begin(), reactions.end(), "Surface Diffusion") != reactions.end()) ) {
            // Both present: LH effective parameters = 1.
            neededRateConstants = neededRateConstants - reactionRateCount.at("Langmuir-Hinshelwood recombination") + 1;~


        } else if (find(reactions.begin(), reactions.end(), "Chemisorption") != reactions.end()) {
            // Chemisorption present: LH effective parameters = 3.
            neededRateConstants = neededRateConstants - reactionRateCount.at("Langmuir-Hinshelwood recombination") + 3;

        } else if (find(reactions.begin(), reactions.end(), "Surface Diffusion") != reactions.end()) {
            // Surface Diffusion present: LH effective parameters = 3.
            neededRateConstants = neededRateConstants - reactionRateCount.at("Langmuir-Hinshelwood recombination") + 3;

        }
    }

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

    int totalNumericNeeded = 0;
    if (reactions.size() == 1 && reactions[0] == "Basic") {
        totalNumericNeeded = 2 + (int)allSpecies.size() + 1;
    } else {
        totalNumericNeeded = 3 + neededRateConstants + (int)allSpecies.size() + 1;
    }

    int numericAvailable = argc - argIndex;
    if (numericAvailable < totalNumericNeeded) {
        cerr << "Not enough numeric parameters provided.\n";
        cerr << "Expected " << totalNumericNeeded << ", got " << numericAvailable << ".\n";
        return 1;
    }

    vector<double> rates;
    rates.reserve(neededRateConstants);
    if (reactions.size() == 1 && reactions[0] == "Basic") {
        for (int i = 0; i < 2; i++) {
            double val = stod(argv[argIndex++]);
            rates.push_back(val);
        }

    } else {

        for (int i = 0; i < 3; i++) {
            double val = stod(argv[argIndex++]);
            rates.push_back(val);
        }

        for (int i = 0; i < neededRateConstants; i++) {
            double val = stod(argv[argIndex++]);
            rates.push_back(val);
        }
    }

    vector<double> initState(allSpecies.size(), 0.0);
    for (int i = 0; i < (int)allSpecies.size(); i++) {
        initState[i] = stod(argv[argIndex++]);
    }

    double t_stop = stod(argv[argIndex++]);

    map<string,int> speciesMap;
    for (int i = 0; i < (int)allSpecies.size(); i++) {
        speciesMap[allSpecies[i]] = i;
    }

    // Safeguard to force A to be chosen
    if (speciesMap.find("A") != speciesMap.end()) {
        initial_A = initState[speciesMap["A"]];
    } else {
        cerr << "Species A is not in the union; cannot define initial_A.\n";
    }

    // Determine if Chemisorption and/or Surface Diffusion are present.
    bool chemPresent = false, surfPresent = false;
    for (auto &r : reactions) {
        if (r == "Chemisorption") chemPresent = true;
        if (r == "Surface Diffusion") surfPresent = true;
    }

    vector<ReactionEvent> events;
    int ratePos = 0;
    for (auto &r : reactions) {
        vector<ReactionEvent> these = buildEventsForReaction(r, rates, ratePos, speciesMap, chemPresent, surfPresent);
        events.insert(events.end(), these.begin(), these.end());
    }

    // File to store the values
    string outputFilename = "output.txt";
    simulateMultiReaction(t_stop, events, initState, allSpecies, outputFilename);

    return 0;
}
