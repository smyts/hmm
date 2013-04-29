What is it
===========
Simple C++ program that implements data structures for
hidden markov model (hmm) and two algorithms (Viterbi for most likely
sequence of hidden states and forward-backward for caclulating
forward and backward state probabilities).
There is also functionality for estimation of the state prediction
results of both algorithms.


Implementation note
===================
External dependencies (logging, build system, libraries) are
evaded on purpose in order to provide pure implementation and make it
easier to build and run for demonstration purposes.


Description of files and directories
====================================
* README.md  - this readme file, uses basic markdown markup
* hmm.h      - header file with declarations of data structures,
               algorithms and estimation functionality
* hmm.cc     - source file with implemenation of the hmm.h delcrarations
* model.spec - description of the file and data format
               for the hmm model description
* data.spec  - description of the file and data format
               for the hmm experiment data with the corresponding model
* model/     - directory for the model description files,
               currently contains only default model and failure tests
* data/      - directory for experiment data, currently contains only default data
* main.cc    - contains code that reads model and experiment data from given files
               and then runs Virterbi and forward-backward algorithms to use them
               as the hidden state predictors. The results of this program are
               the state prediction estimations for both algorithms and all states.
               Estimation is printed to the standard output and contains
               True Positives, False Positives, True Negatives, False Negatives
               and f-measure for each particular state prediction.

Technical stuff
===============

Compilation
-----------
* Just do it from the project directory:
  g++ main.cc hmm.cc -o app -std=c++11 -Wall -Wextra

Run with default example data
-----------------------------
* Compile as above and run from the project directory as:
  ./app models/default.model data/default.data

Simple testing
--------------
* There are models inside 'model/' dir as test cases for some trivial model validation.
  All of them, except one (default), are supposed to fail with different errors, which correspond to their file names.
  It is possible to use the following command to test against those test cases:
  g++ main.cc hmm.cc -o app -std=c++11 -Wall -Wextra && ls -1 models/*.model | xargs -r -n 1 -d '\n' -I 'modelfile' sh -c "./app modelfile data/default.data || true"
