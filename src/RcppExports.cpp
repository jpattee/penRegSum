// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// genotypeMatrix
arma::mat genotypeMatrix(const std::string fileName, int N, int P, arma::Col<int> col_skip_pos, arma::Col<int> col_skip, arma::Col<int> keepbytes, arma::Col<int> keepoffset, const int fillmissing);
RcppExport SEXP _penRegSum_genotypeMatrix(SEXP fileNameSEXP, SEXP NSEXP, SEXP PSEXP, SEXP col_skip_posSEXP, SEXP col_skipSEXP, SEXP keepbytesSEXP, SEXP keepoffsetSEXP, SEXP fillmissingSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::string >::type fileName(fileNameSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< int >::type P(PSEXP);
    Rcpp::traits::input_parameter< arma::Col<int> >::type col_skip_pos(col_skip_posSEXP);
    Rcpp::traits::input_parameter< arma::Col<int> >::type col_skip(col_skipSEXP);
    Rcpp::traits::input_parameter< arma::Col<int> >::type keepbytes(keepbytesSEXP);
    Rcpp::traits::input_parameter< arma::Col<int> >::type keepoffset(keepoffsetSEXP);
    Rcpp::traits::input_parameter< const int >::type fillmissing(fillmissingSEXP);
    rcpp_result_gen = Rcpp::wrap(genotypeMatrix(fileName, N, P, col_skip_pos, col_skip, keepbytes, keepoffset, fillmissing));
    return rcpp_result_gen;
END_RCPP
}
// normalize2
arma::vec normalize2(arma::mat& genotypes);
RcppExport SEXP _penRegSum_normalize2(SEXP genotypesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type genotypes(genotypesSEXP);
    rcpp_result_gen = Rcpp::wrap(normalize2(genotypes));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_penRegSum_genotypeMatrix", (DL_FUNC) &_penRegSum_genotypeMatrix, 8},
    {"_penRegSum_normalize2", (DL_FUNC) &_penRegSum_normalize2, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_penRegSum(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
