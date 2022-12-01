#include <Rcpp.h>
using namespace Rcpp;

//' @title A Gibbs sampler using Rcpp
//' @description A Gibbs sampler using Rcpp
//' @param N the number of samples
//' @return a random sample of size \code{N}
//' @examples
//' \dontrun{
//' rnC <- gibbsC(100)
//' par(mfrow=c(2,1));
//' plot(rnC[,1],type='l')
//' plot(rnC[,2],type='l')
//' }
//' @export
// [[Rcpp::export]]
NumericMatrix gibbsC(int N) {
  NumericMatrix mat(N, 2);
  double  rho=0.9, mu1=0, mu2=0, sigma1=1, sigma2=1, s1=0, s2=0;
  s1 = sqrt(1-rho*rho)*sigma1;
  s2 = sqrt(1-rho*rho)*sigma2;
  double x = 0, y = 0;
  for(int i = 0; i < N; i++) {
    x = rnorm(1, mu1 + rho * (y - mu2) * sigma1/sigma2, s1)[0];
    y = rnorm(1, mu2 + rho * (x - mu1) * sigma2/sigma1, s2)[0];
    mat(i, 0) = x;
    mat(i, 1) = y;
  }
  return(mat);
}
