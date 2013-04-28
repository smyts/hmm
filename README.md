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
* model/     - directory for the model description files, currently contains only default model
* data/      - directory for experiment data, currently contains only default data
