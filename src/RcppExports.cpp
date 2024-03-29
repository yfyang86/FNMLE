// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// prodbyrow_c
NumericMatrix prodbyrow_c(NumericMatrix mat, NumericVector vec);
RcppExport SEXP _FNMLE_prodbyrow_c(SEXP matSEXP, SEXP vecSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type mat(matSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type vec(vecSEXP);
    rcpp_result_gen = Rcpp::wrap(prodbyrow_c(mat, vec));
    return rcpp_result_gen;
END_RCPP
}
// permsign_c
NumericMatrix permsign_c(int n);
RcppExport SEXP _FNMLE_permsign_c(SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(permsign_c(n));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_FNMLE_prodbyrow_c", (DL_FUNC) &_FNMLE_prodbyrow_c, 2},
    {"_FNMLE_permsign_c", (DL_FUNC) &_FNMLE_permsign_c, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_FNMLE(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
