// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(openmp)]]

#include <omp.h>
#include <RcppEigen.h>

// [[Rcpp::export]]
SEXP eigenMapMatMult2(const Eigen::Map<Eigen::MatrixXd> A,
                      Eigen::Map<Eigen::MatrixXd> B,
                      int n_cores){

  Eigen::setNbThreads(n_cores);
  Eigen::MatrixXd C = A * B;
  return Rcpp::wrap(C);

}
