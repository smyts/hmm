#include <algorithm>
#include <vector>
#include <stdexcept>
#include <iostream>
#include <cstddef>

#include "hmm.h"

using std::vector;
using std::string;
using std::pair;

using HMM::Data::Model;
using HMM::Data::ExperimentData;

//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Data namespace definitions >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

/**
 * \note
 * Auxiliary functions, for internal usage only.
 */
namespace
{
    /**
     * \brief Converts first string element to char emission symbol for HMM
     *
     * \note
     * String assumed to be non-empty and symbol[0] is supposed to be in a..z ascii.
     */
    char symbolToInd(const string& symbol)
    {
        return symbol[0] - 'a';
    }
};

void Model::ReadModel(std::istream& modelSource)
{
    // section: states reading
    size_t nstates;
    string stateName;

    modelSource >> nstates;

    if (nstates < 2) {
        throw std::domain_error("There must be at least two states: begin and end");
    }

    for (size_t i = 0; i < nstates; ++i) {
        modelSource >> stateName;
        stateNameToIndex[stateName] = i;
        stateIndexToName.push_back(stateName);
    }

    // section: alphabet reading
    modelSource >> alphabetSize;

    // section: transitions reading
    size_t ntransitions;
    string targetStateName;

    transitionProb.assign(nstates, vector<double> (nstates, 0));
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

    // section: state-symbol emission probabilities reading
    size_t nemissions;
    string symbol; // supposed to be single character, string is used for simpler reading code

    stateSymbolProb.assign(nstates, vector<double> (alphabetSize, 0));
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
    string stateName;
    string symbol;// supposed to be single character, string is used for simpler reading code

    dataSource >> nsteps;

    if (nsteps == 0) {
        throw std::domain_error("Empty experiment data");
    }

    for (size_t i = 0; i < nsteps; ++i) {
        dataSource >> stepNumber >> stateName >> symbol;

        size_t stateInd = model.stateNameToIndex.at(stateName);
        size_t symbolInd = symbolToInd(symbol);

        timeStateSymbol.emplace_back(stepNumber, stateInd, symbolInd);
    }
}
//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< end of Data namespace definitions <<<<<<<<<<<<<<<<<<<<<<<<<<


//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Algorithms namespace definitions >>>>>>>>>>>>>>>>>>>>>>>>>>
const size_t HMM_UNDEFINED_STATE = -1;

/**
 * \note
 * Auxiliary functions, for internal usage only.
 */
namespace
{
    /**
     * \brief Aux. function to calculate new state probability for the Viterbi algorithm step
     */
    double CalcNewStateProbability(size_t stepNumber, size_t prevState,
                                   size_t curState, size_t curSymbol, const Model& model,
                                   const vector<vector<double> >& sequenceProbability)
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

    /**
     * \brief Aux. function to find the best previous state during the Viterbi algorithm step
     */
    size_t FindBestTransitionSource(size_t stepNumber, size_t curState,
                                    size_t curSymbol, const Model& model,
                                    const vector<vector<double> >& sequenceProbability)
    {
        if (stepNumber == 0) {
            return 0;
        }

        size_t nstates = model.transitionProb.size();
        double bestProbValue = -1;
        size_t bestPrevState = HMM_UNDEFINED_STATE;

        for (size_t prevState = 0; prevState < nstates; ++prevState) {
            double curProb = CalcNewStateProbability(stepNumber, prevState,
                                                     curState, curSymbol,
                                                     model, sequenceProbability);

            if (curProb > bestProbValue) {
                bestProbValue = curProb;
                bestPrevState = prevState;
            }
        }

        // there must be at least two states => the result won't be undefined
        return bestPrevState;
    }

    /**
     * \brief Aux. function to get the cumulative forward step transition probability
     *
     * \details
     * This is used inside forward-backward algorithm at forward probabilities calculation.
     */
    double CalcForwardStepProbability(size_t stepNumber, size_t curState,
                                      const Model& model, const ExperimentData& data,
                                      const vector<vector<double> >& forwardStateProbability)
    {
        size_t nstates = model.transitionProb.size();
        size_t curSymbol = std::get<2> (data.timeStateSymbol[stepNumber]);

        if (stepNumber == 0) {
            return model.transitionProb[0][curState] * model.stateSymbolProb[curState][curSymbol];
        } else {
            double prevCumulativeProb = 0;

            for (size_t prevState = 0; prevState < nstates; ++prevState) {
                prevCumulativeProb += (forwardStateProbability[stepNumber - 1][prevState] *
                                       model.transitionProb[prevState][curState]);
            }

            return prevCumulativeProb * model.stateSymbolProb[curState][curSymbol];
        }
    }

    /**
     * \brief Aux. function to get the cumulative probability for the backward step
     *
     * \details
     * This is used inside forward-backward algorithm at backward probabilities calculation.
     */
    double CalcBackwardStepProbability(size_t stepNumber, size_t curState,
                                       const Model& model, const ExperimentData& data,
                                       const vector<vector<double> >& backwardStateProbability)
    {
        size_t nstates = model.transitionProb.size();
        size_t maxtime = data.timeStateSymbol.size();

        if (stepNumber + 1 == maxtime) {
            return 1.; // probability to describe empty sequence is 1.
        } else {
            size_t nextSymbol = std::get<2> (data.timeStateSymbol[stepNumber + 1]);
            double nextCumulativeProb = 0.;

            for (size_t nextState = 0; nextState < nstates; ++nextState) {
                nextCumulativeProb += (model.transitionProb[curState][nextState] *
                                       model.stateSymbolProb[nextState][nextSymbol] *
                                       backwardStateProbability[stepNumber + 1][nextState]);
            }

            return nextCumulativeProb;
        }
    }
};

vector<size_t>
HMM::Algorithms::FindMostProbableStateSequence(const Model& model, const ExperimentData& data)
{
    // section: prepare and initialize data structures for calculations
    size_t nstates = model.transitionProb.size();
    size_t maxtime = data.timeStateSymbol.size();

    /**
     * \note
     * sequenceProbability[i][j] is the probability of the most probable sequence of states
     * for 1..i observations for which the last state is j-th
     */
    vector<vector<double> > sequenceProbability(maxtime,
                                                vector<double> (nstates, 0));

    /**
     * \note
     * prevSeqState[i][j] is the previous state from which the most probable
     * sequence (with probability sequenceProbability[i][j]) for 1..i observations with the last
     * state at j has been formed.
     * This information will help to recover the whole sequence.
     */
    vector<vector<size_t> > prevSeqState(maxtime,
                                         vector<size_t> (nstates, HMM_UNDEFINED_STATE));

    // section: calculate probabilities for Viterbi algorithm using dynamic programming approach
    for (size_t t = 0; t < maxtime; ++t) {
        for (size_t curState = 0; curState < nstates; ++curState) {
            size_t curSymbol = std::get<2> (data.timeStateSymbol[t]);
            size_t bestPrevState = FindBestTransitionSource(t, curState,
                                                            curSymbol, model,
                                                            sequenceProbability);
            double bestProbValue = CalcNewStateProbability(t, bestPrevState,
                                                           curState, curSymbol,
                                                           model, sequenceProbability);

            sequenceProbability[t][curState] = bestProbValue;
            prevSeqState[t][curState] = bestPrevState;
        }
    }

    // section: collect most probable sequence in the reverse order
    vector<size_t> mostProbableSeq;
    ptrdiff_t curStep = maxtime - 1;

    // find the last state of the most probable sequence to start recovery from it
    size_t curState = std::distance(std::begin(sequenceProbability[curStep]),
                                    std::max_element(std::begin(sequenceProbability[curStep]),
                                                     std::end(sequenceProbability[curStep])));

    for (; curStep > 0; --curStep) {
        curState = prevSeqState[curStep][curState];
        mostProbableSeq.push_back(curState);
    }

    mostProbableSeq.push_back(curState);

    // section: restore correct order and return results
    std::reverse(std::begin(mostProbableSeq), std::end(mostProbableSeq));

    return std::move(mostProbableSeq);
}

vector<vector<pair<double, double> > >
HMM::Algorithms::CalcForwardBackwardProbabiliies(const Model& model, const ExperimentData& data)
{
    size_t nstates = model.transitionProb.size();
    size_t maxtime = data.timeStateSymbol.size();

    /**
     * \note
     * forwardStateProbability[i][j] is the probability that any hidden sequence (with
     * the hidden state at i-th step equal to j) describes first 1..i observations.
     */
    vector<vector<double> > forwardStateProbability(maxtime, vector<double> (nstates, 0));

    // section: calculate forward probabilities of the forward-backward algorithm
    for (size_t t = 0; t < maxtime; ++t) {
        for (size_t curState = 0; curState < nstates; ++curState) {
            double cumulativePrevProbability =
                CalcForwardStepProbability(t, curState, model, data, forwardStateProbability);

            forwardStateProbability[t][curState] = cumulativePrevProbability;
        }
    }

    /**
     * \note
     * backwardStateProbability[i][j] is the probability that any hidden sequence (with
     * the hidden state at i+1 step equal to j) describes last i+1..T observations.
     */
    vector<vector<double> > backwardStateProbability(maxtime, vector<double> (nstates, 0.));

    // section: calculate backward probabilities of the forward-backward algorithm
    for (ptrdiff_t t = maxtime - 1; t >= 0; --t) {
        for (size_t curState = 0; curState < nstates; ++curState) {
            double cumulativeNextProbability =
                CalcBackwardStepProbability(t, curState, model, data, backwardStateProbability);

            backwardStateProbability[t][curState] = cumulativeNextProbability;
        }
    }

    // section: return joined results
    vector<vector<pair<double, double> > > forwardBackwardProbability
        (maxtime, vector<pair<double, double> > (nstates, pair<double, double>()));

    for (size_t t = 0; t < maxtime; ++t) {
        for (size_t curState = 0; curState < nstates; ++curState) {
            forwardBackwardProbability[t][curState] =
                pair<double, double>(forwardStateProbability[t][curState],
                                     backwardStateProbability[t][curState]);
        }
    }

    return std::move(forwardBackwardProbability);
}
//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< end of Algorithms namespace definitions <<<<<<<<<<<<<<<<<<<<


//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Estimation namespace definitions >>>>>>>>>>>>>>>>>>>>>>>>>>>
vector<size_t> HMM::Estimation::GetMostProbableStates(
    const vector<vector<pair<double, double> > >& forwardBackwardProb)
{
    size_t maxtime = forwardBackwardProb.size();
    vector<size_t> mostProbableStates;

    for (size_t t = 0; t < maxtime; ++t) {
        size_t mostProbableState =
            std::distance(std::begin(forwardBackwardProb[t]),
                          std::max_element(std::begin(forwardBackwardProb[t]),
                                           std::end(forwardBackwardProb[t]),
                                           [](const pair<double, double>& prev,
                                              const pair<double, double>& next)
                                           {return (prev.first * prev.second <
                                                    next.first * next.second);}));
        mostProbableStates.push_back(mostProbableState);
    }

    return std::move(mostProbableStates);
}

vector<vector<size_t> >
HMM::Estimation::CombineConfusionMatrix(const ExperimentData& realData,
                                        const vector<size_t>& predictedStates,
                                        const Model& model)
{
    size_t maxtime = predictedStates.size();
    size_t nstates = model.transitionProb.size();
    vector<vector<size_t> > confusionMatrix(nstates, vector<size_t> (nstates, 0));

    for (size_t t = 0; t < maxtime; ++t) {
        size_t predictedInd = predictedStates[t];
        size_t realInd      = std::get<1>(realData.timeStateSymbol[t]);

        ++confusionMatrix[predictedInd][realInd];
    }

    return std::move(confusionMatrix);
}

vector<HMM::Data::PredictionEstimation>
HMM::Estimation::GetStatePredictionEstimations(const vector<vector<size_t> >& confusionMatrix)
{
    size_t nstates = confusionMatrix.size();
    vector<PredictionEstimation> estimations(nstates);
    vector<size_t> colSums(nstates, 0);
    vector<size_t> rowSums(nstates, 0);

    // section: prepare auxiliary columns sums and row sums for further usage
    for (size_t i = 0; i < nstates; ++i) {
        for (size_t j = 0; j < nstates; ++j) {
            rowSums[i] += confusionMatrix[i][j];
            colSums[i] += confusionMatrix[j][i];
        }
    }

    size_t totalObservations = std::accumulate(std::begin(rowSums),
                                               std::end(rowSums), 0UL);

    // section: calculate prediction estimations for each state
    for (size_t state = 0; state < nstates; ++state) {
        estimations[state].truePositives = confusionMatrix[state][state];
        estimations[state].falsePositives = rowSums[state] - confusionMatrix[state][state];

        // neither predicted to be current state nor its real state is the current one
        estimations[state].trueNegatives = totalObservations - rowSums[state] - colSums[state] + confusionMatrix[state][state]; 
        estimations[state].falseNegatives = colSums[state] - confusionMatrix[state][state];

        // calculate f-measure
        double precision = 0;
        double recall = 0;

        if (rowSums[state] != 0) {
            precision = static_cast<double> (confusionMatrix[state][state]) / static_cast<double> (rowSums[state]);
        }

        if (colSums[state] != 0) {
            recall = static_cast<double> (confusionMatrix[state][state]) / static_cast<double> (colSums[state]);
        }

        if (rowSums[state] == 0 && colSums[state] == 0) {
            estimations[state].fMeasure = 0;
        } else {
            estimations[state].fMeasure = 2. * (precision * recall) / (precision + recall);
        }

    }

    return std::move(estimations);
}
//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< end of Estimation namespace definitions <<<<<<<<<<<<<<<<<<<<
