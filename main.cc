#include <fstream>
#include <iostream>

#include "hmm.h"

void showUsage(std::string programName)
{
    std::cerr << "Usage: " << programName
              << " path_to_model path_to_data " << std::endl;
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
        std::cerr << "ERROR: fatal problem while reading model. Details: '" << e.what() << "'" << std::endl;
        return -1;
    } catch (...) {
        std::cerr << "ERROR: unknown exception while reading model" << std::endl;
        return -1;
    }

    try
    {
        data.ReadExperimentData(model, dataSource);
    } catch(std::exception& e) {
        std::cerr << "ERROR: fatal problem while reading experiment data. Details: '" << e.what() << "'" << std::endl;
        return -1;
    } catch (...) {
        std::cerr << "ERROR: unknown exception while reading experiment data " << std::endl;
        return -1;
    }

    // secton: run and estimate viterbi predictions

    // section: run and estimate forward-backward predictions

    return 0;
}
