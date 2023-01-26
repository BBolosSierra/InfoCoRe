// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// test
void test(const arma::SpMat<double>& Adj);
RcppExport SEXP _infoCoRe_test(SEXP AdjSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::SpMat<double>& >::type Adj(AdjSEXP);
    test(Adj);
    return R_NilValue;
END_RCPP
}
// test2
void test2(const arma::SpMat<double>& Adj);
RcppExport SEXP _infoCoRe_test2(SEXP AdjSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::SpMat<double>& >::type Adj(AdjSEXP);
    test2(Adj);
    return R_NilValue;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_infoCoRe_test", (DL_FUNC) &_infoCoRe_test, 1},
    {"_infoCoRe_test2", (DL_FUNC) &_infoCoRe_test2, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_infoCoRe(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
