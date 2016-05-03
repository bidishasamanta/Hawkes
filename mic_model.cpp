#include "mic_model.h"
#include <iostream>
#include <fstream>
#include <ctime>
#include <cmath>
#include <cstdlib>
#include <iomanip>

using namespace std;

const double sqrt2pi = 2.506628274631000502415765284811;
//const double exp = 2.71828;
const double sqrt2 = 1.4142135623730950488016887242097;
const double pi = 3.14159265358979323846;
const double sigma_min = 0.1;
const double g_squre_resolution = 1e-15;
const int max_iter = 5000;

CMicModel::CMicModel( double t_value )
{
	m_T = t_value;
	//m_m = 10.0;
};

CMicModel::~CMicModel()
{
}

/*
 *	get the value of parameter alpha
 *
 */
double CMicModel::get_alpha()
{
	return m_alpha;
}

/*
 *	get the value of parameter beta
 *
 */
double CMicModel::get_beta()
{
	return m_beta;
}

/*
 *	get the value of parameter gamma
 *
 */
double CMicModel::get_gamma()
{
	return m_gamma;
}

int CMicModel::parameter_estimation( vector<int> const & citation_time )
{
	if ( citation_time.empty() )
	{
		//cout << "lambda:\t" << 0 << endl;
		//cout << "mu:\t any value" << endl;
		//cout << "sigma:\t any value" << endl;
		return 0;
	}

	m_citations = citation_time;

	double gradient_alpha;
	double gradient_beta;
	double gradient_gamma;
	double a_a;
	double a_b;
	double a_g;
	double b_b;
	double b_g;
	double g_g;

	// initializing seed for random function
	srand( (int)time(NULL) );

	// initializing parameters mu and lambda, computing the parameter lambda
	//m_mu = 6 + 1.0 * rand() / RAND_MAX;    // generate random value between 0.0 and 2.0;
	//m_sigma = 1.0 + 0.5 * rand() / RAND_MAX; // generate random value between 0.01 and 1.01;

	m_alpha = 0.01;
	m_beta = 0.01;
	m_gamma = 0.01;

	double LL = log_likelihood();

	double g_square = 0.0;
	double delta_alpha;
	double delta_beta;
	double delta_gamma;

	int round;
	for ( round = 1; round < max_iter; ++round )
	{
		gradient_alpha = compute_gradient_alpha();
		//cout << gradient_alpha;
		gradient_beta = compute_gradient_beta();
		//cout << gradient_beta;
		gradient_gamma = compute_gradient_gamma();
		//cout << gradient_gamma;
		hessian_matrix(a_a, b_b, g_g, a_b, a_g, b_g);
		transform_hessian_matrix(a_a, b_b, g_g, a_b, a_g, b_g);
		delta_alpha = a_a * gradient_alpha + a_b * gradient_beta + a_g * gradient_gamma;//mu_mu * gradient_mu + mu_sigma * gradient_sigma;
		delta_beta = a_b * gradient_alpha + b_b * gradient_beta + b_g * gradient_gamma;
		delta_gamma = a_g * gradient_alpha + b_g * gradient_beta + g_g * gradient_gamma;

		// update parameters mu, sigma and lambda
		line_search( delta_alpha, delta_beta, delta_gamma, LL );

		// checking convergence
		g_square = gradient_alpha * gradient_alpha + gradient_beta * gradient_beta + gradient_gamma * gradient_gamma;
		if ( g_square < g_squre_resolution )
			break;

		//cout << "log likelihood: " << LL << "\t" << "g_square: " << g_square << endl;
	}

	//cout << "\nround: " << round << endl;
	//cout << "lambda:\t" << m_lambda << endl;
	//cout << "mu:\t" << m_mu << endl;
	//cout << "sigma:\t" << m_sigma << endl;

	return 0;
}


/*
 * compute the log likelihood for the given parameters: alpha, gamma, beta
 */
double CMicModel::log_likelihood()
{
	if ( m_citations.empty() )
		return 0;

	double LL = 0.0;

	double A = (m_gamma / m_beta) *(exp(-m_beta * m_T ) - 1);
	double B = 0.0;
	for(int i = 0; i< m_citations.size(); i++){
		B +=  exp(-m_beta*(m_T-m_citations[i])) - 1;
	}
	B = B * (m_alpha / m_beta);

	double C = 0.0;
	for(int i = 0; i< m_citations.size(); i++){
		double logsum = 0.0;
		for(int j = 0; j< i; j++){
			logsum += exp(-m_beta*(m_citations[i] - m_citations[j]));
		}
		C += log(m_gamma * exp(-m_beta*m_citations[i]) + logsum);
	}
	LL = A + B + C;
	return LL;
}


/*
 * compute the gradient of mu
 *
 */
double CMicModel::compute_gradient_alpha( )
{
	double gradient_alpha = 0.0;

	double B = 0.0;
	double C = 0.0;
	double denominator = 0.0;
	double numerator = 0.0;
	for(int i = 0; i< m_citations.size(); i++){
		B +=  exp( -m_beta*(m_T-m_citations[i])) - 1;
	}
	B /= m_beta;

	for(int i = 0; i< m_citations.size(); i++){
		double logsum = 0.0;
		for(int j = 0; j< i; j++){
			logsum += exp(-m_beta*(m_citations[i] - m_citations[j]));
		}
		denominator = m_gamma * exp(-m_beta*m_citations[i]) + m_alpha * logsum;
		numerator = logsum;
		C += numerator / denominator;
		}
	gradient_alpha = B + C;
	return gradient_alpha;
}

double CMicModel::compute_gradient_gamma()
{
	double gradient_gamma = 0.0;

	double B = 0.0;
	double C = 0.0;
	double denominator = 0.0;
	double numerator = 0.0;

	B = (exp( -m_beta*m_T) -1 )/m_beta;

	for(int i = 0; i< m_citations.size(); i++){
		double logsum = 0.0;
		for(int j = 0; j< i; j++){
			logsum += exp(-m_beta*(m_citations[i] - m_citations[j]));
		}
		denominator = m_gamma * exp(-m_beta*m_citations[i]) + m_alpha * logsum;
		numerator = -exp(-m_beta * m_citations[i]);
		C += numerator / denominator;
	}
	gradient_gamma = B + C;
	return gradient_gamma;
}


/*
 *  compute lambda by letting derivative equal to 0.
 *
 */
double CMicModel::compute_gradient_beta( )
{
	double gradient_beta = 0.0;

	double B = 0.0;
	double C = 0.0;
	double A = 0.0;
	double denominator = 0.0;
	double numerator = 0.0;

	A = (m_gamma / m_beta) * (-m_T * exp(-m_beta * m_T)) - (m_gamma / (m_beta * m_beta))*(exp( -m_beta*m_T) - 1);
	for(int i = 0; i< m_citations.size(); i++){
		B +=  exp(-m_beta*(m_T-m_citations[i])) - 1;
	}
	B = B * (m_alpha/(m_beta * m_beta));
	double sum = 0.0;
	for(int i = 0; i< m_citations.size(); i++){
			sum +=  exp( -m_beta*(m_T-m_citations[i])) * (m_T -m_citations[i]);
	}
	B = B + (m_alpha/m_beta) * sum;
	for(int i = 0; i< m_citations.size(); i++){
		double logsumdenom = 0.0, logsum=0.0;
		for(int j = 0; j< m_citations[m_citations.size()-1]; j++){
			//logsumdenom += exp( -m_beta*(m_citations[i] - m_citations[j]));
			logsum += exp(m_beta * m_citations[j]);//(m_citations[i] - m_citations[j]) * exp( -m_beta*(m_citations[i] - m_citations[j]));

		}
		//denominator = m_gamma * exp(-m_beta*m_citations[i]) + m_alpha * logsumdenom;

		//numerator = -m_gamma * m_citations[i] * exp(-m_beta * m_citations[i])- m_alpha *logsum;
		C +=(-m_citations[i] + 1) - (m_gamma / (m_gamma + m_alpha * logsum));
	}
	gradient_beta = A - B + C;
	return gradient_beta;
}

double CMicModel::compute_alpha_alpha( )
{
	double alpha_alpha = 0.0;

	double sum = 0.0;

	for(int i = 0; i< m_citations.size(); i++){
		sum = 0.0;
		for(int j = 0; j< i; j++){
			sum +=  exp(m_beta*(m_citations[j]));
		}
		alpha_alpha = alpha_alpha + (sum * sum)/pow(m_gamma + m_alpha * sum, 2);
	}
	return -alpha_alpha;
}

double CMicModel::compute_gamma_gamma( )
{
	double gamma_gamma = 0.0;

	double sum = 0.0;

	for(int i = 0; i< m_citations.size(); i++){
		sum = 0.0;
		for(int j = 0; j< i; j++){
			sum +=  exp( m_beta*(m_citations[j]));
		}
		gamma_gamma = 1 /pow((m_gamma + m_alpha * sum), 2);
	}
	return gamma_gamma;
}

double CMicModel::compute_beta_beta()
{
	double beta_beta = 0.0;

	double sum = 0.0;
	double A = 0.0;
	double B1 = 0.0, B2 = 0.0, B3 = 0.0, B = 0.0;
	double C1 = 0.0, C2 = 0.0, C3 = 0.0, C = 0.0;
	double sumNumerator, sumDenominator;

	A = ((2 * m_gamma * m_T)/pow(m_beta, 2) + (pow(m_T,2) * m_gamma) / m_beta + (2 * m_gamma)/pow(m_beta, 3) ) * exp( -m_beta*m_T);
	A = A - 2 * m_gamma / pow(m_beta, 3) - 2 * m_alpha * m_citations.size()/pow(m_beta, 3);

	for(int i = 0; i < m_citations.size(); i++){
		C1 += exp( -m_beta *(m_T - m_citations[i]));
		C2 += (m_T - m_citations[i]) * exp(-m_beta *(m_T - m_citations[i]));
		C3 += pow((m_T - m_citations[i]), 2) * exp(-m_beta *(m_T - m_citations[i]));
	}
	C = (2 * m_alpha * C1) / pow(m_beta, 3) + (2 * m_alpha * C2) / pow(m_beta, 2) + (m_alpha * C3) / pow(m_beta, 1);

	for(int i = 0; i< m_citations.size(); i++){
		sumNumerator = 0.0;
		sumDenominator = 0.0;
		for(int j = 0; j< i; j++){
			sumDenominator +=  exp(m_beta*(m_citations[j]));
			sumNumerator += exp(m_beta*m_citations[j]) * m_citations[j];
		}
		B += (m_gamma * m_alpha * sumNumerator) / pow(m_gamma + m_alpha * sumDenominator , 2);
	}
	beta_beta = A + B + C;
	return beta_beta;
}


double CMicModel::compute_alpha_gamma( )
{
	double alpha_gamma = 0.0;

	double sum = 0.0;

	for(int i = 0; i< m_citations.size(); i++){
		sum = 0.0;
		for(int j = 0; j< i; j++){
			sum +=  exp( m_beta*(m_citations[j]));
		}
		alpha_gamma = alpha_gamma +  sum/pow(m_gamma + m_alpha * sum, 2);
	}
	return -alpha_gamma;
}

double CMicModel::compute_alpha_beta()
{
	double alpha_beta = 0.0;

	double B = 0.0;
	double C = 0.0;
	double A = 0.0;
	double denominator = 0.0;
	double numerator = 0.0;

	//A = (m_gamma / m_beta) * (-m_T * pow(exp, -m_beta * m_T)) - (m_gamma / (m_beta * m_beta))*(pow(exp, -m_beta*m_T) - 1);
	for(int i = 0; i< m_citations.size(); i++){
		A +=  exp(-m_beta*(m_T-m_citations[i])) - 1;
	}
	A = -A / (m_beta * m_beta);

	for(int i = 0; i< m_citations.size(); i++){
		B +=  exp(-m_beta*(m_T-m_citations[i])) * (m_T - m_citations[i]);
	}
	B = -B/ m_beta;

	double sum = 0.0;
	for(int i = 0; i< m_citations.size(); i++){
		double sumDenom = 0.0;
		for (int j = 0; j<i; j++){
			sumDenom +=  m_citations[j] * exp( m_beta * m_citations[j]);
		}
		denominator = pow(m_gamma + m_alpha * sumDenom, 2);
		numerator = m_gamma * m_alpha *sumDenom ;
		C = C + numerator / denominator ;
	}

	alpha_beta = A + B + C;
	return alpha_beta;
}

double CMicModel::compute_gamma_beta()
{
	double gamma_beta = 0.0;

	double B = 0.0;
	double C = 0.0;
	double A = 0.0;
	double denominator = 0.0;
	double numerator = 0.0;

	//A = (m_gamma / m_beta) * (-m_T * pow(exp, -m_beta * m_T)) - (m_gamma / (m_beta * m_beta))*(pow(exp, -m_beta*m_T) - 1);

	A = -(exp(-m_beta * m_T) - 1) / (m_beta * m_beta);

	B = -(exp( -m_beta * m_T) * m_T) / 	m_beta;

	double sum = 0.0;
	for(int i = 0; i< m_citations.size(); i++){
		double sumDenom = 0.0;
		double sumNumer = 0.0;
		for (int j = 0; j<i; j++){
			sumNumer += exp(m_beta * m_citations[j]) * m_citations[j];
			sumDenom +=  exp( m_beta * m_citations[j]);
		}
		denominator = pow(m_gamma + m_alpha * sumDenom, 2);
		numerator = m_citations[i] * m_alpha * sumNumer ;
		C = C + numerator / denominator ;
	}

	gamma_beta = A + B + C;
	return gamma_beta;
}
/*
 * compute the value of Gaussian density function at the place of x
 *
 */
inline double CMicModel::gaussian_density_distribtuion_fuction( double x )
{
	return exp(-x*x/2.0) / sqrt2pi;
}

/*
 * the error function used to compute the cumulative Gasussian distribution function
 *
 */
inline double erf (double x)
{
	// constants
	const double a1 =  0.254829592;
	const double a2 = -0.284496736;
	const double a3 =  1.421413741;
	const double a4 = -1.453152027;
	const double a5 =  1.061405429;
	const double p =   0.3275911;

	// A&S formula 7.1.26
	double t = 1.0/(1.0 + p*std::fabs(x));
	double y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*exp(-x*x);

	if (x < 0)
		return -y;
	return y;
}

/*
 *	compute the value of cumulative Gasussian distribution function at the place of x
 *
 */
inline double CMicModel::gaussian_cumulative_distribution_function( double x )
{
	return (1.0+erf(x/sqrt2))*0.5;
}


/*
 *	compute the value of log normal density function at the place of x, with the parameters mu and sigma
 *
 */
inline double CMicModel::log_normal( double x, double mu, double sigma )
{
	double t = (log(x)-mu) / sigma;
	double value = exp( -0.5*t*t ) / x / sqrt2pi / sigma;

	return value;
}

void CMicModel::hessian_matrix( double & a_a, double & b_b, double & g_g, double & a_b, double & a_g, double & b_g)
{
	// compute second derivate with respect to \mu
	a_a =  compute_alpha_alpha();
	a_b = compute_alpha_beta();
	a_g = compute_alpha_gamma();
	b_g = compute_gamma_beta();
	b_b = compute_beta_beta();
	g_g = compute_gamma_gamma();
}



/*
 *	transform the hessian matrix into the inverse of negative hessian matrix
 *
 */
void CMicModel::transform_hessian_matrix(double & a_a, double & b_b, double & g_g, double & a_b, double & a_g, double & b_g)
{
	double a = a_a;
	double b = a_b;
	double c = a_g;
	double d = b_b;
	double e = b_g;
	double f = g_g;
	double det = 0.0;
	double a_temp, b_temp, c_temp, d_temp, e_temp, f_temp;
	a_temp = d * f - e * e;
	b_temp = - b * f + c * e;
	c_temp = b * e - c * d;
	d_temp = a * f - c * c;
	e_temp = - a* e + b * c;
	f_temp = a * d - b * b;

	det = (a * d * f + 2 * b * c * e - (a * e * e + f * b * b + d * c * c));

	a_a = a_temp/det;
	a_b = b_temp/det;
	a_g = c_temp/det;
	b_b = d_temp/det;
	b_g = e_temp/det;
	g_g = f_temp/det;


}

/*
 *	search the new value of mu and sigma which can increase the log likelihood, along the direction of gradient
 *
 */
void CMicModel::line_search( double delta_alpha, double delta_beta, double delta_gamma, double & LL )
{
	double old_LL = LL;
	double old_alpha = delta_alpha;
	double old_beta = delta_beta;
	double old_gamma = delta_gamma;
	// check the boundary of sigma
	/*if ( m_alpha + delta_alpha < alpha_min )
	{
		//cout << "boundary check\n";
		delta_alpha *= sigma_min-m_sigma / delta_alpha;
		delta_alpha = alpha_min-m_alpha; // delta_sigma *= (sigma_min-m_sigma) / delta_sigma
	}*/

	//double old_delta_sigma = delta_sigma;
	//double old_delta_mu = delta_mu;

	// guarantee the increase of log likelihood
	for ( int ii = 0; ii < 100; ++ii )
	{
		m_alpha = old_alpha + delta_alpha;
		m_beta = old_beta + delta_beta;
		m_gamma = old_gamma + delta_gamma;
		LL = log_likelihood();
		if ( LL >= old_LL ) return;

		delta_alpha /= 2.0;
		delta_gamma /= 2.0;
		delta_beta /= 2.0;
	}

	m_alpha = old_alpha;
	m_beta = old_beta;
	m_gamma = old_gamma;


}

/*
 *	print the landscape of log likelihood with respect to parameter mu and sigma
 *
 */
void CMicModel::check_log_likelihood( vector<int> const & citation_time )
{
	m_citations = citation_time;

	ofstream outFile( "log_likelihood.txt");
	if ( outFile.fail() )
	{
		cout << "Can not open the file: 'log_likelihood.txt'" << endl;
		return;
	}

	double l_mu = 7.4;
	double u_mu = 7.5;
	double l_sigma = 1.6;
	double u_sigma = 1.7;
	for ( m_mu = l_mu; m_mu < u_mu; m_mu += 0.001 )
	{
		for ( m_sigma = l_sigma; m_sigma < u_sigma; m_sigma += 0.001 )
		{
			double LL = log_likelihood();

			outFile << m_mu << "\t" << m_sigma << "\t" << LL << "\n";
		}
	}

	outFile.close();
}

/*
 *	compute the hessian matrix numerically, i.e., using the gradient function
 *
 */
/*void CMicModel::numerical_hessian_matrix( double & mu_mu, double & mu_sigma, double & sigma_sigma )
{
	double delta = 0.00001;

	// compute mu_mu
	double old_mu = m_mu;
	m_mu -= delta;
	double left_gradient_mu_delta = compute_gradient_mu();
	m_mu = old_mu + delta;
	double right_gradient_mu_delta = compute_gradient_mu();
	mu_mu = (right_gradient_mu_delta - left_gradient_mu_delta) / (2*delta);
	m_mu = old_mu;

	// compute mu_sigma
	old_mu = m_mu;
	m_mu -= delta;
	double left_gradient_delta = compute_gradient_sigma();
	m_mu = old_mu + delta;
	double right_gradient_delta = compute_gradient_sigma();
	mu_sigma = (right_gradient_delta-left_gradient_delta) / (2*delta);
	m_mu = old_mu;

	// compute sigma_sigma
	double old_sigma = m_sigma;
	m_sigma -= delta;
	double left_gradient_sigma_delta = compute_gradient_sigma();
	m_sigma = old_sigma + delta;
	double right_gradient_sigma_delta = compute_gradient_sigma();
	sigma_sigma = (right_gradient_sigma_delta - left_gradient_sigma_delta) / (2*delta);
	m_sigma = old_sigma;
}*/

void CMicModel::compute_log_likelihood( vector<int> const & citation_time, double lambda, double mu, double sigma )
{
	m_lambda = lambda;
	m_mu = mu;
	m_sigma = sigma;
	m_citations = citation_time;

	cout << "LL: " << log_likelihood() << endl;
}
