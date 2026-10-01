// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::interfaces(r, cpp)]]

#define ARMA_64BIT_WORD 1
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <math.h>

using namespace Rcpp;

// [[Rcpp::export]]
double K2_Binary_SPA_alt(double x, arma::vec muhat, arma::vec G)
{
    double res = 0.0;
	double temp1 = 0.0;
	double temp2 = 1.0;
	
    const int n = muhat.size();

    for(int i = 0; i < n; i++)
    {
		temp1 = muhat(i) * (1 - muhat(i)) * pow(G(i),2.0) * exp(x * G(i));
		temp2 = muhat(i) * exp(x * G(i)) + (1 - muhat(i));
		
        res = res + temp1/pow(temp2,2.0);
    }
	
    return res;
}
