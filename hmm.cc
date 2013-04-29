#include <algorithm>
#include <vector>
#include <stdexcept>
#include <iostream>
#include <cstddef>

#include "hmm.h"


//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Data namespace definitions >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
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
//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< end of Data namespace definitions <<<<<<<<<<<<<<<<<<<<<<<<<<


//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Algorithms namespace definitions >>>>>>>>>>>>>>>>>>>>>>>>>>
const double HMM_BEGIN_STATE_PROBABILITY = 1.0;
const size_t HMM_UNDEFINED_STATE = -1;

namespace
{
    double CalcNewStateProbability(size_t stepNumber, size_t prevState,
                                   size_t curState, char curSymbol, const Model& model,
                                   const std::vector<std::vector<double> >& sequenceProbability)
    {
        double prevProbability = 1.;

        if (stepNumber == 0 && prevState == 0) {
            prevProbability = 1.;
        } else if (stepNumber == 0 && prevState != 0) {
            prevProbability = 0.;
        } else {
            prevProbability = sequenceProbability[stepNumber - 1][prevState];
        }

        return (prevProbability *
                model.transitionProb[prevState][curState] *
                model.stateSymbolProb[curState][curSymbol]);
    }

    size_t FindBestTransitionSource(size_t stepNumber, size_t curState,
                                    size_t curSymbol, const Model& model,
                                    const std::vector<std::vector<double> >& sequenceProbability)
    {
        if (stepNumber == 0) {
            return 0;
        }

        size_t nstates = model.transitionProb.size();
        double bestProbValue = -1;
        size_t bestPrevState = HMM_UNDEFINED_STATE;

        for (size_t prevState = 0; prevState < nstates; ++prevState) {
            double curProb = CalcNewStateProbability(stepNumber, prevState, curState,curSymbol, model, sequenceProbability);

            if (curProb > bestProbValue) {
                bestProbValue = curProb;
                bestPrevState = prevState;
            }
        }

        // there must be at least two states => the result won't be undefined
        return bestPrevState;
    }
};

std::vector<size_t>
HMM::Algorithms::FindMostProbableStateSequence(const Model& model, const ExperimentData& data)
{
    // part: prepare and initialize data structures for calculations
    size_t nstates = model.transitionProb.size();
    size_t maxtime = data.timeStateSymbol.size();

    std::vector<std::vector<double> > sequenceProbability(maxtime,
                                                          std::vector<double> (nstates, 0));
    std::vector<std::vector<size_t> > prevSeqState(maxtime,
                                                   std::vector<size_t> (nstates,
                                                                        HMM_UNDEFINED_STATE));

    // part: calculate probabilities for Viterbi algorithm using dynamic programming approach
    for (size_t t = 0; t < maxtime; ++t) {
        for (size_t curState = 0; curState < nstates; ++curState) {
            size_t curSymbol = std::get<2> (data.timeStateSymbol[t]);
            size_t bestPrevState = FindBestTransitionSource(t, curState, curSymbol, model, sequenceProbability);
            double bestProbValue = CalcNewStateProbability(t, bestPrevState, curState, curSymbol, model, sequenceProbability);

            sequenceProbability[t][curState] = bestProbValue;
            prevSeqState[t][curState] = bestPrevState;
        }
    }

    // part: collect most probable sequence in the reverse order
    std::vector<size_t> mostProbableSeq;
    ptrdiff_t curStep = maxtime - 1;
    size_t curState = std::distance(std::begin(sequenceProbability[curStep]),
                                    std::max_element(std::begin(sequenceProbability[curStep]),
                                                std::end(sequenceProbability[curStep])));

    for (; curStep > 0; --curStep) {
        curState = prevSeqState[curStep][curState];
        mostProbableSeq.push_back(curState);
    }

    mostProbableSeq.push_back(curState);


    // part: restore correct order and return results
    std::reverse(std::begin(mostProbableSeq), std::end(mostProbableSeq));

    return std::move(mostProbableSeq);
}

std::vector<std::pair<double, double> >
HMM::Algorithms::CalcForwardBackwardProbabiliies(const Model& model, const ExperimentData& data)
{
    // TODO: implement
    // TODO: add inverse stateInd to stateName conversions
}


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< end of Algorithms namespace definitions <<<<<<<<<<<<<<<<<<<<
