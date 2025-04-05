#pragma once

#include <vector>
#include <functional>
#include <map>
#include <string>

auto prop_single = [](int idx) {
    return [=](const std::vector<double>& state, double k) -> double {
        return k * state[idx];
    };
};

auto prop_bimolecular = [](int idx1, int idx2) {
    return [=](const std::vector<double>& state, double k) -> double {
        return k * state[idx1] * state[idx2];
    };
};

auto prop_square = [](int idx) {
    return [=](const std::vector<double>& state, double k) -> double {
        return k * state[idx] * state[idx];
    };
};

// Structure for a reaction event.
struct ReactionEvent {
    std::function<double(const std::vector<double>&, double)> propensity;
    std::vector<double> delta;
    double k;
};

// Function declarations

// Builds ReactionEvent objects for a given reaction.
// The parameter 'general_params_extracted' is passed by value here, though you might choose
// to manage that globally instead.
std::vector<ReactionEvent> buildEventsForReaction(const std::string& reaction, 
    const std::vector<double>& rates, 
    int& rateIndex, 
    const std::map<std::string,int>& speciesIndex,
    bool chemPresent,
    bool surfPresent);

// Prints a progress bar to the terminal that updates in-place.
// 'progress' is the current progress, 'total' is the total value (e.g. simulation stop time).
void printProgressBar(double progress, double total);

// Runs the Gillespie simulation until t_stop and writes a full trajectory to outputFilename.
void simulateMultiReaction(
    double t_stop,
    const std::vector<ReactionEvent>& events,
    std::vector<double>& state,
    const std::vector<std::string>& speciesList,
    const std::string& outputFilename);