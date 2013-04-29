#include <map>
#include <tuple>
#include <utility>
#include <vector>
#include <iostream>

namespace HMM
{
    namespace Data
    {
        struct Model
        {
            /**
             * \brief Read model description from the stream
             *
             * \details 
             * Constructor reads model description according to the specification file.
             * \note
             * It is supposed that stream is correct and contains all necessary data.
             * In order to catch errors make sure to enable exceptions for the stream before passing it here.
             */
            void ReadModel(std::istream& modelSource);

            /// number of different emission symbols (first such from a..z range in ascii)
            size_t alphabetSize; 

            /// conversion of state name string to state index
            std::map<std::string, size_t> stateNameToIndex;

            /// element[i][j] here is the probability of transition from state i to j
            /// very first state is the begin state, the last is the end state
            std::vector<std::vector<double> > transitionProb;

            /// element[i][j] here is the probability to emit symbol j from state i
            std::vector<std::vector<double> > stateSymbolProb;
        };

        struct ExperimentData
        {
            /**
             * \brief Read experiment data from the stream
             *
             * \details 
             * Constructor reads experiment data according to the specification file.
             * \note
             * It is supposed that stream is correct and contains all necessary data.
             * In order to catch errors make sure to enable exceptions for the stream before passing it here.
             */
            void ReadExperimentData(const Model& model, std::istream& dataSource);

            /// Data triples as (time, state, symbol_emitted)
            std::vector<std::tuple<size_t, size_t, size_t> > timeStateSymbol;
        };

        /**
         * \note
         * Estimation in the form:
         * True positives, False Positives, True Negatives, False Negatives, F-measure
         */
        typedef std::tuple<size_t, size_t, size_t, size_t, double> PredictionEstimation;
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
         * \returns vector result[t], where result[t].first is a(t, i) and result[t].second is b(t, i)
         */
        std::vector<std::pair<double, double> >
        CalcForwardBackwardProbabiliies(const Model& model, const ExperimentData& data);
    };

    namespace Estimations
    {
        using Data::Model;
        using Data::ExperimentData;
        using Data::PredictionEstimation;

        std::vector<PredictionEstimation>
        GetStatePredictionEstimations(const ExperimentData& realData,
                                      const std::vector<size_t>& predictedStates);
    };
};
