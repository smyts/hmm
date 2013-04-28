#include <vector>
#include <stdexcept>
#include <iostream>

#include "hmm.h"

using HMM::Data::Model;
using HMM::Data::ExperimentData;

namespace
{
    /**
     * \brief Converts first string element to char emission symbol for HMM
     *
     * \note
     * String assumed to be non-empty and symbol[0] is supposed to be in a..z ascii.
     */
    char symbolToInd(const std::string& symbol)
    {
        return symbol[0] - 'a';
    }
};

void Model::ReadModel(std::istream& modelSource)
{
    // part: states reading
    size_t nstates;
    std::string stateName;

    modelSource >> nstates;

    for (size_t i = 0; i < nstates; ++i) {
        modelSource >> stateName;
        stateNameToIndex[stateName] = i;
    }

    // part: alphabet reading
    modelSource >> alphabetSize;

    // part: transitions reading
    size_t ntransitions;
    std::string targetStateName;

    transitionProb.assign(nstates, std::vector<double> (nstates, 0));
    modelSource >> ntransitions;

    for (size_t i = 0; i < ntransitions; ++i) {
        double prob;
        modelSource >> stateName >> targetStateName >> prob;

        size_t fromInd = stateNameToIndex[stateName];
        size_t toInd   = stateNameToIndex[targetStateName];

        if (fromInd + 1 == nstates) {
            throw std::domain_error("Transition from the ending state is forbidden");
        }

        if (toInd == 0) {
            throw std::domain_error("Transition to the starting state is forbidden");
        }

        transitionProb[fromInd][toInd] = prob;
    }

    // part: state-symbol emission probabilities reading
    size_t nemissions;
    std::string symbol; // supposed to be single character, string is used for simpler reading code

    stateSymbolProb.assign(nstates, std::vector<double> (alphabetSize, 0));
    modelSource >> nemissions;

    for (size_t i = 0; i < nemissions; ++i) {
        double prob;
        modelSource >> stateName >> symbol >> prob;

        size_t stateInd = stateNameToIndex[stateName];
        size_t symbolInd = symbolToInd(symbol);

        if (stateInd == 0 || stateInd + 1 == nstates) {
            throw std::domain_error("Symbol emission from the beginning or the ending states is forbidden");
        }

        stateSymbolProb[stateInd][symbolInd] = prob;
    }
}

void ExperimentData::ReadExperimentData(const Model& model, std::istream& dataSource)
{
    size_t nsteps;
    size_t stepNumber;
    std::string stateName;
    std::string symbol;// supposed to be single character, string is used for simpler reading code

    dataSource >> nsteps;

    for (size_t i = 0; i < nsteps; ++i) {
        dataSource >> stepNumber >> stateName >> symbol;

        size_t stateInd = model.stateNameToIndex.at(stateName);
        size_t symbolInd = symbolToInd(symbol);

        timeStateSymbol.emplace_back(stepNumber, stateInd, symbolInd);
    }
}
