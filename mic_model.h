#ifndef MIC_MODEL_H
#define MIC_MODEL_H

#include <vector>

class CMicModel
{
public:
	CMicModel(  double t_value );
	~CMicModel();

public:
	int parameter_estimation( std::vector<int> const & citation_time );

	void check_log_likelihood( std::vector<int> const & citation_time );

	void compute_log_likelihood( std::vector<int> const & citation_time, double lambda, double mu, double sigma );

public:
	double get_mu();
	double get_sigma();
	double get_lambda();
	double get_alpha();
	double get_beta();
	double get_gamma();

private:
	double gaussian_cumulative_distribution_function( double upper );

	inline double gaussian_density_distribtuion_fuction( double x );

	double compute_lambda();

	double compute_gradient_alpha();

	double compute_gradient_beta();

	double compute_gradient_gamma();

	double compute_alpha_alpha();

	double compute_alpha_beta();

	double compute_alpha_gamma();

	double compute_gamma_gamma();

	double compute_gamma_beta();

	double compute_beta_beta();

	double log_likelihood();

	double log_normal( double x, double mu, double sigma );

	void hessian_matrix(double & a_a, double & b_b, double & g_g, double & a_b, double & a_g, double & b_g);

	void numerical_hessian_matrix( double & mu_mu, double & mu_sigma, double & sigma_sigma );

	void transform_hessian_matrix(double & a_a, double & b_b, double & g_g, double & a_b, double & a_g, double & b_g);

	void line_search( double delta_alpha, double delta_beta, double delta_gamma, double &LL );

private:
	double m_m;
	double m_T;
	double m_lambda;
	double m_mu;
	double m_sigma;
	double m_gamma;
	double m_beta;
	double m_alpha;

	std::vector<int> m_citations;
};

#endif
