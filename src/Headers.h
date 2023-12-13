#if !defined(HEADERS_INCLUDED)
#define HEADERS_INCLUDED

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]

//C and C++
#include <unistd.h>
#include <ios>
#include <iostream>
#include <fstream>

#include <stdio.h>
#include <iostream>
#include <string>
#include <limits>
#include <tuple>
#include <unordered_map>
#include <cstdint>     //defines uint_32, uint64 etc...
#include <cstddef>     //defines std::size_t (unsigned) and std::ssize_t (signed)
#include <stdlib.h>    // to use malloc and free
#include <mm_malloc.h>
#include <cmath>

#include <complex>

//Armadillo
#include <RcppArmadillo.h>

//Name space
using namespace arma;
using namespace Rcpp;
using namespace std;

//store coordinates, and weight, for an edge
typedef std::tuple<double, double, double> EDGE;

#endif
