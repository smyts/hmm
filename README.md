What is it
===========
Simple single-file C++ program that implements data structures for
hidden markov model and two algorithms (Viterbi for most likely
sequence of hidden states and forward-backward for caclulating
forward and backward state probabilities).


Implementation note
===================
External dependencies (logging, build system, libraries) are
evaded on purpose in order to provide pure implementation and make it
easier to build and run for demonstration purposes.


Description of files and directories
====================================
* README.md  - this readme file, uses basic markdown markup
* hmm.cc     - source file of this program
* model.spec - description of the file format for the hmm model description
* data.spec  - description of the file format for the hmm experiment data with the corresponding model
* model/     - directory for the model description files, currently contains only default model and failure tests
* data/      - directory for experiment data, currently contains only default data


Technical stuff
===============

Compilation
-----------
* Just do it:
  g++ main.cc hmm.cc -o app -std=c++11 -Wall -Wextra

Simple testing
--------------
* There are models inside 'model/' dir as test cases for some trivial model validation.
  All of them, except one (default), are supposed to fail with different errors, which correspond to their file names.
  It is possible to use the following command to test against those test cases:
  g++ main.cc hmm.cc -o app -std=c++11 -Wall -Wextra && ls -1 models/*.model | xargs -r -n 1 -d '\n' -I 'modelfile' sh -c "./app modelfile data/default.data || true"