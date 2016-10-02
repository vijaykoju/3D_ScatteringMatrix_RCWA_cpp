#!/bin/bash

#g++ -L ../lib/ -I Eigen -O2 -DNDEBUG -msse2 rcwa_test.cpp -lfftw3 -lm
#g++ -L ../lib/ -I Eigen -Ofast -DNDEBUG rcwa_test.cpp -lfftw3 -lm -fopenmp 
icc -L ../lib/ -I Eigen -Ofast -DNDEBUG -msse2 rcwa_test.cpp -lfftw3 -lm -openmp 
#g++ -Ofast -DNDEBUG rcwa_test.cpp -lfftw3 -lm -fopenmp

time -p ./a.out
