// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::interfaces(r, cpp)]]

#define ARMA_64BIT_WORD 1
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <math.h>

using namespace Rcpp;

static void K1_K2_Binary_SPA(double x, const arma::vec& muhat, const arma::vec& G, double q, double& first, double& second)
{
	first = 0.0;
	second = 0.0;

	for(arma::uword i = 0; i < muhat.n_elem; i++)
	{
		double exponential = exp(-x * G(i));
		first = first - muhat(i) * G(i);
		first = first + muhat(i) * G(i)/(muhat(i) + (1 - muhat(i))*exponential);
		double temp1 = muhat(i) * (1 - muhat(i)) * pow(G(i),2.0)*exponential;
		double temp2 = muhat(i) + (1 - muhat(i)) * exponential;
		second = second + temp1/pow(temp2,2.0);
	}

	first = first - q;
}



// declare K_Binary_SPA
double K_Binary_SPA(double x, arma::vec muhat, arma::vec G);
// declare K_Binary_SPA_alt (alternative way if not converge)
double K_Binary_SPA_alt(double x, arma::vec muhat, arma::vec G);
// declare K1_Binary_SPA (first derivative)
double K1_Binary_SPA(double x, arma::vec muhat, arma::vec G, double q);
// declare K1_Binary_SPA_alt (first derivative, alternative way if not converge)
double K1_Binary_SPA_alt(double x, arma::vec muhat, arma::vec G, double q);
// declare K2_Binary_SPA (second derivative)
double K2_Binary_SPA(double x, arma::vec muhat, arma::vec G);
// declare K2_Binary_SPA_alt (second derivative, alternative way if not converge)
double K2_Binary_SPA_alt(double x, arma::vec muhat, arma::vec G);
// check na
bool check_is_na(double x); 
// check same sign
bool haveSameSign(double a, double b);

// [[Rcpp::export]]
double NR_Binary_SPA(arma::vec muhat, arma::vec G, double q, double init, double tol, int max_iter)
{
	double xi = 0.0;
	double xi_update = 0.0;
	
	// first step
	xi = init;
	xi_update = init;
	
	if(fabs(K1_Binary_SPA(xi, muhat, G, q)) > tol)
	{
		double first = 0.0;
		double second = 0.0;
		K1_K2_Binary_SPA(xi,muhat,G,q,first,second);
		xi_update = xi - first/second;
	}
	
	// iteration number
	int no_iter = 0;
	double numerator = 0.0;
	double denominator = 1.0;
	
	while (!(check_is_na(xi_update))&&(R_finite(xi_update)==1)&&(fabs(xi_update - xi) > tol)&&(fabs(K1_Binary_SPA(xi_update, muhat, G, q)) > tol)&&(no_iter < max_iter))
    {
		no_iter = no_iter + 1;
		xi = xi_update;	
		
		// calculate numerator
		K1_K2_Binary_SPA(xi,muhat,G,q,numerator,denominator);
		if((R_finite(numerator)==0)||(check_is_na(numerator)))
		{
			numerator = K1_Binary_SPA_alt(xi, muhat, G, q);
		}
		
		// calculate denominator
		
		if((R_finite(denominator)==0)||(check_is_na(denominator)))
		{
			denominator = K2_Binary_SPA_alt(xi, muhat, G);
		}
		
		xi_update = xi - numerator/denominator;
    }

	// test finite
	if((R_finite(xi_update)==0)||(check_is_na(xi_update)))
	{
		xi_update = xi;
	}
	
	
	// test convergence
	double w = 0.0; 
	double xhat = 0.0; 
	double ki = 0.0;

	if(no_iter == max_iter)
	{
		w = sqrt(2*(xhat*q - K_Binary_SPA(xhat,muhat,G)));
		xhat = xi_update;
		
		if (xhat < 0){
			w = -w;
		}

		ki = xhat*sqrt(K2_Binary_SPA(xhat,muhat,G));
		if((R_finite(ki)==0)||(check_is_na(ki)))
		{
			ki = xhat*sqrt(K2_Binary_SPA_alt(xhat,muhat,G));
		}
	
		if(fabs(w+log(ki/w)/w) < 38)
		{
			xi_update = 0.0;	
		}
	}
	
	return xi_update;
}
