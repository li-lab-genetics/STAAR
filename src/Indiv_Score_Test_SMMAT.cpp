// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List Indiv_Score_Test_SMMAT(arma::sp_mat G, const arma::mat& P, arma::vec residuals)
{
	int i;

	// number of markers
	int p = G.n_cols;

	// Uscore
	arma::rowvec Uscore = trans(residuals)*G;

	arma::vec pvalue;
	pvalue.zeros(p);

	arma::vec Uscore_se;
	Uscore_se.zeros(p);

	double test_stat = 0;

	arma::mat P_G;
	P_G = P*G;

	for(i = 0; i < p; i++)
	{
		double Cov_ii = 0.0;
		for(auto g = G.begin_col(i); g != G.end_col(i); ++g)
		{
			Cov_ii += P_G(g.row(),i)*(*g);
		}
		Uscore_se(i) = sqrt(Cov_ii);

		if (Cov_ii == 0)
		{
			pvalue(i) = 1;
		}
		else
		{
			test_stat = pow(Uscore(i),2)/Cov_ii;
			pvalue(i) = R::pchisq(test_stat,1,false,false);
		}

	}

	return List::create(Named("Uscore") = trans(Uscore), Named("Uscore_se") = Uscore_se, Named("pvalue") = pvalue);
}

