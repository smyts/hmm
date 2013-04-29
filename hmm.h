#ifndef HMM_H
#define HMM_H

#include <map>
#include <tuple>
#include <utility>
#include <vector>
#include <iostream>


/**
 * \note
 * Data structures structures, algorithms and prediction estimation for hidden markov models.
 */
namespace HMM
{
    namespace Data
    {
        /**
         * \brief Represents hidden markov model description
         */
        struct Model
        {
            /**
             * \brief Read model description from the stream
             *
             * \details 
             * It reads model description according to the specification file.
             * \note
             * It is supposed that the source stream is correct and contains all necessary data.
             * In order to catch errors make sure to enable exceptions for the stream before passing it here.
             */
            void ReadModel(std::istream& modelSource);

            /// number of different emission symbols (first such from a..z range in ascii)
            size_t alphabetSize; 

            /// conversion of state name string to state index
            std::map<std::string, size_t> stateNameToIndex;

            /// inverse conversion
            std::vector<std::string> stateIndexToName;

            /// element[i][j] here is the probability of transition from state i to j
            /// very first state is the begin state, the last is the end state
            std::vector<std::vector<double> > transitionProb;

            /// element[i][j] here is the probability to emit symbol j from state i
            std::vector<std::vector<double> > stateSymbolProb;
        };

        /**
         * \brief Represents experiment data for some particular model
         */
        struct ExperimentData
        {
            /**
             * \brief Read experiment data from the stream
             *
             * \details 
             * It reads experiment data according to the specification file.
             * \note
             * It is supposed that the source stream is correct and contains all necessary data.
             * In order to catch errors make sure to enable exceptions for the stream before passing it here.
             */
            void ReadExperimentData(const Model& model, std::istream& dataSource);

            /// Data triples as (time, state, symbol_emitted)
            std::vector<std::tuple<size_t, size_t, size_t> > timeStateSymbol;
        };

        /**
         * \brief State prediction estimation results for different hidden markov models algorithms
         */
        struct PredictionEstimation
        {
            size_t truePositives;
            size_t falsePositives;
            size_t trueNegatives;
            size_t falseNegatives;
            double fMeasure;
        };
    };

    namespace Algorithms
    {
        using Data::Model;
        using Data::ExperimentData;

        /**
         * \brief Finds most probable sequence of hidden states
         *
         * \details
         * Implementation is based on the Viterbi algorithm.
         *
         * \returns vector with predicted hidden state indices
         */
        std::vector<size_t>
        FindMostProbableStateSequence(const Model& model, const ExperimentData& data);

        /**
         * \brief Calculates alpha-beta value pairs for each time moment
         *
         * \details
         * Implementation is based on the Forward-Backward algorithm.
         * This function calculates pairs of alpha and beta values
         * for each time moment.
         * alpha -> a(t, i) is the dependent probability of the i-th
         * hidden state based on the first 0..t emitted symbols
         * beta  -> b(t, i) is the dependent probability of the i-th
         * hidden state based on the last t+1..END symbols.
         *
         * \returns vector result[t][i], where result[t][i].first is a(t, i)
         *          and result[t][i].second is b(t, i)
         */
        std::vector<std::vector<std::pair<double, double> > >
        CalcForwardBackwardProbabiliies(const Model& model, const ExperimentData& data);
    };

    namespace Estimation
    {
        using Data::Model;
        using Data::ExperimentData;
        using Data::PredictionEstimation;

        using std::vector;
        using std::pair;

        /**
         * \brief Use forward-backward probabilities to get the most probable state at each step
         */
        vector<size_t> GetMostProbableStates(
            const vector<vector<pair<double, double> > >& forwardBackwardProb);

        /**
         * \note
         * Confusion matrix element[i][j] is the number of elements with the
         * predicted state i when their real state is j.
         * It is used as an auxiliary data structure for calculation of various estimations
         * (true positives and etc., f-measure).
         */
        vector<vector<size_t> > CombineConfusionMatrix(const ExperimentData& realData,
                                                       const vector<size_t>& predictedStates, const Model& model);

        /**
         * \brief Use confusion matrix to calculate estimations of the prediction results
         *
         * \returns vector of prediction estimations for each state
         */
        vector<PredictionEstimation>
            GetStatePredictionEstimations(const vector<vector<size_t> >& confusionMatrix);
    };
};

#endif // HMM_H
