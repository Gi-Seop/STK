#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix crossatan(NumericMatrix m1, NumericMatrix m2) {
  int nrow1 = m1.nrow();
  int nrow2 = m2.nrow();
  int ncol = m1.ncol();

  if (ncol != m2.ncol()) {
    throw std::runtime_error("Incompatible number of dimensions");
  }

  NumericMatrix out(nrow1, nrow2);

  for (int r1 = 0; r1 < nrow1; r1++) {
    for (int r2 = 0; r2 < nrow2; r2++) {
      double dx = 0;
      double dy = 0;
      dx = m1(r1, 0) - m2(r2, 0);
      dy = m1(r1, 1) - m2(r2, 1);
      if(dx == 0) {
        out(r1, r2) = 1.570796;
      } else {
        out(r1, r2) = atan(dy/dx);
      }
    }
  }

  return out;
}
