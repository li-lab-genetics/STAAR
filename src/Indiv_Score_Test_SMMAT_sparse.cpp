// [[Rcpp::depends(RcppArmadillo)]]

#define ARMA_64BIT_WORD 1
#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

static double sparse_col_dot(const arma::sp_mat& A, const arma::sp_mat& B, arma::uword col)
{
	arma::sp_mat::const_col_iterator a_it = A.begin_col(col);
	arma::sp_mat::const_col_iterator a_end = A.end_col(col);
	arma::sp_mat::const_col_iterator b_it = B.begin_col(col);
	arma::sp_mat::const_col_iterator b_end = B.end_col(col);

	double out = 0;

	while ((a_it != a_end) && (b_it != b_end))
	{
		if (a_it.row() == b_it.row())
		{
			out += (*a_it) * (*b_it);
			++a_it;
			++b_it;
		}
		else if (a_it.row() < b_it.row())
		{
			++a_it;
		}
		else
		{
			++b_it;
		}
	}

	return out;
}

// [[Rcpp::export]]
List Indiv_Score_Test_SMMAT_sparse(arma::sp_mat G, arma::sp_mat Sigma_i, arma::mat Sigma_iX, arma::mat cov, arma::vec residuals)
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

	int q = Sigma_iX.n_cols;

	arma::mat tSigma_iX_G;
	tSigma_iX_G.zeros(q,p);

	arma::sp_mat S_G;
	S_G = Sigma_i*G;
	tSigma_iX_G = trans(Sigma_iX)*G;

	for(i = 0; i < p; i++)
	{
		double Cov_ii = sparse_col_dot(S_G, G, i);
		if(q > 0)
		{
			Cov_ii -= arma::as_scalar(trans(tSigma_iX_G.col(i))*cov*tSigma_iX_G.col(i));
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

