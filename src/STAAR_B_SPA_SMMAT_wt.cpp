// [[Rcpp::depends(RcppArmadillo)]]

#define ARMA_64BIT_WORD 1
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <math.h>

using namespace Rcpp;

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
// declare NR
double NR_Binary_SPA(arma::vec muhat, arma::vec G, double q, double init, double tol, int max_iter);
// declare golden section search method for selecting initial value of bisection algorothm
arma::vec goldenSectionSearchForSignChange(double a, double b, arma::vec muhat, arma::vec G, double q, double tol, int max_iter);
// declare Bisection method
double Bisection_Binary_SPA(arma::vec muhat, arma::vec G, double q, double xmin, double xmax, double tol);
// declare Saddle_Binary_SPA using NR
double Saddle_Binary_SPA(double q, arma::vec muhat, arma::vec G, double tol, int max_iter, bool lower);
// declare Saddle_Binary_SPA using Bisection
double Saddle_Binary_SPA_Bisection(double q, arma::vec muhat, arma::vec G, double tol, int max_iter, double xmin, double xmax, bool lower);
// check na
bool check_is_na(double x); 
// check same sign
bool haveSameSign(double a, double b);

// [[Rcpp::export]]
arma::vec STAAR_B_SPA_SMMAT_wt(arma::mat G, arma::vec residuals, arma::vec muhat, arma::mat weights_B, double tol, int max_iter,
double p_filter_cutoff, arma::mat P)
{
	int i,k;
	double sum0 = 0.0;
	
	bool lower1 = false;
	bool lower2 = true;

	int vn = G.n_rows;
	int un = G.n_cols;
	int wn = weights_B.n_cols;

	// p-values
	arma::vec res;
	res.ones(wn);
	
	// p-values components
	double respart1 = 0.0;
	double respart2 = 0.0;

	// Score statistics
	arma::rowvec x = trans(residuals)*G;
	int n = x.size();
	// cumulative genotype function
	arma::mat G_cumu;
	G_cumu.zeros(vn,wn);
	G_cumu = G*weights_B;
	
	// init value of bisection
	double xmin = -100.0;
	double xmax = 100.0;
	
	// parameter for normal approximation in p-value calculation
	arma::mat Cov;
	Cov.zeros(un,un);
	
	Cov = trans(P*G)*G;
	
	arma::mat Covw;
	Covw.zeros(un,un);
	
	arma::mat Wleft;
	Wleft.zeros(un,un);
	
	arma::mat Wright;
	Wright.zeros(un,un);
	
	double sumw = 0.0;
	double sumx = 0.0;
	
	bool logp = false;
	
	for(i = 0; i < wn; i++)
	{
		sum0 = 0.0;
		for (k = 0; k < n; k++)
		{
			sum0 = sum0 + x(k)*weights_B(k,i);
		}
		
		// p-value calculation using normal approximation
		Wright.each_row() = trans(weights_B.col(i));
		Wleft.each_col() = weights_B.col(i);
		
		Covw = Wleft%Cov%Wright;
		sumw = arma::accu(Covw);
		
		sumx = pow(sum0, 2) / sumw;
		
		res(i) = R::pchisq(sumx,1,lower1,logp);
		
		// p-value calculation using SPA approximation
		if(res(i) < p_filter_cutoff)
		{
			
			// calculate p-value
			respart1 = 1.0;
			respart2 = 1.0;
			respart1 = Saddle_Binary_SPA(fabs(sum0), muhat, G_cumu.col(i), tol, max_iter, lower1);
			
			if((check_is_na(respart1))||(respart1 == 1.0))
			{
				respart1 = Saddle_Binary_SPA_Bisection(fabs(sum0), muhat, G_cumu.col(i), tol, max_iter, xmin, xmax, lower1);
			}
			
			if((check_is_na(respart1))||(respart1 == 1.0))
			{
				res(i) = 1;
			}else
			{
				respart2 = Saddle_Binary_SPA(-fabs(sum0), muhat, G_cumu.col(i), tol, max_iter, lower2);
				
				if((check_is_na(respart2))||(respart2 == 1.0))
				{
					respart2 = Saddle_Binary_SPA_Bisection(-fabs(sum0), muhat, G_cumu.col(i), tol, max_iter, xmin, xmax, lower2);
				}
				
				if((check_is_na(respart2))||(respart2 == 1.0))
				{
					res(i) = 1;
				}else
				{
					res(i) = respart1 + respart2;
				}
			}
			
			// if saddle point approximation fails, set the p-value as 1.0
			if(res(i) > 1)
			{
				res(i) = 1;		
			}
		}
	}

	return res;
}
