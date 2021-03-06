// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// gaussian1
arma::vec gaussian1(arma::vec x);
RcppExport SEXP monfuncreg_gaussian1(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(gaussian1(x));
    return rcpp_result_gen;
END_RCPP
}
// estimate1
double estimate1(arma::vec TT1, arma::vec YYYY, arma::vec NN, double h1, arma::vec weight, double t);
RcppExport SEXP monfuncreg_estimate1(SEXP TT1SEXP, SEXP YYYYSEXP, SEXP NNSEXP, SEXP h1SEXP, SEXP weightSEXP, SEXP tSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type TT1(TT1SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type YYYY(YYYYSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type NN(NNSEXP);
    Rcpp::traits::input_parameter< double >::type h1(h1SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type weight(weightSEXP);
    Rcpp::traits::input_parameter< double >::type t(tSEXP);
    rcpp_result_gen = Rcpp::wrap(estimate1(TT1, YYYY, NN, h1, weight, t));
    return rcpp_result_gen;
END_RCPP
}
