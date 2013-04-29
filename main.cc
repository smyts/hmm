#include <fstream>
#include <iostream>

#include "hmm.h"

void showUsage(std::string programName)
{
    std::cerr << "Usage: " << programName
              << " path_to_model path_to_data " << std::endl;
}

void printPredictionEstimation(size_t stateInd,
                               const HMM::Data::PredictionEstimation& estimation,
                               const HMM::Data::Model& model)
{
    std::cout << "State " << model.stateIndexToName[stateInd]
              << " => "
              << "True Positives=" << estimation.truePositives << ", "
              << "False Positives=" << estimation.falsePositives << ", "
              << "True Negatives=" << estimation.trueNegatives << ", "
              << "False Negatives=" << estimation.falseNegatives << ", "
              << "f-measure=" << estimation.fMeasure << '\n';
}

int main(int argc, char* argv[])
{
    // section: check arguments and prepare input streams
    if (argc < 3) {
        showUsage(argv[0]);
        return -1;
    }

    std::ifstream modelSource(argv[1]);
    std::ifstream dataSource(argv[2]);

    if (! modelSource.good()) {
        std::cerr << "ERROR: Failed to open model file properly." << std::endl;
        return -1;
    }
    if (! dataSource.good()) {
        std::cerr << "ERROR: Failed to open data file properly." << std::endl;
        return -1;
    }

    // enable exceptions to signal errors later while reading model and data
    std::ios_base::iostate ioExcept = (std::ifstream::failbit |
                                       std::ifstream::badbit  |
                                       std::ifstream::eofbit);
    modelSource.exceptions(ioExcept);
    dataSource.exceptions(ioExcept);


    // section: read model and data
    HMM::Data::Model model;
    HMM::Data::ExperimentData data;

    try
    {
        model.ReadModel(modelSource);
    } catch(std::exception& e) {
        std::cerr << "ERROR: fatal problem while reading model. Details: '" << e.what()
                  << "'" << std::endl;
        return -1;
    } catch (...) {
        std::cerr << "ERROR: unknown exception while reading model" << std::endl;
        return -1;
    }

    try
    {
        data.ReadExperimentData(model, dataSource);
    } catch(std::exception& e) {
        std::cerr << "ERROR: fatal problem while reading experiment data. Details: '" << e.what()
                  << "'" << std::endl;
        return -1;
    } catch (...) {
        std::cerr << "ERROR: unknown exception while reading experiment data " << std::endl;
        return -1;
    }

    // secton: run and estimate viterbi predictions
    std::vector<size_t> mostProbableSeq =
        HMM::Algorithms::FindMostProbableStateSequence(model, data);
    std::vector<std::vector<size_t> > confusionMatrix =
        HMM::Estimation::CombineConfusionMatrix(data, mostProbableSeq, model);
    std::vector<HMM::Data::PredictionEstimation> estimations =
        HMM::Estimation::GetStatePredictionEstimations(confusionMatrix);

    std::cout << "Viterbi algorithm state prediction estimations:\n";

    // skip first and last states (begin and end)
    for (size_t i = 1; i + 1 < estimations.size(); ++i) {
        printPredictionEstimation(i, estimations[i], model);
    }

    std::cout << "\n";

    // section: run and estimate forward-backward predictions
    std::vector<std::vector<std::pair<double, double> > > forwardBackwardProb =
        HMM::Algorithms::CalcForwardBackwardProbabiliies(model, data);
    std::vector<size_t> mostProbableStates =
        HMM::Estimation::GetMostProbableStates(forwardBackwardProb);
    confusionMatrix = HMM::Estimation::CombineConfusionMatrix(data, mostProbableStates, model);
    estimations = HMM::Estimation::GetStatePredictionEstimations(confusionMatrix);

    std::cout << "Forward-backward algorithm state prediction estimations:\n";

    // skip first and last states (begin and end)
    for (size_t i = 1; i + 1 < estimations.size(); ++i) {
        printPredictionEstimation(i, estimations[i], model);
    }

    std::cout << "\n";

    return 0;
}
