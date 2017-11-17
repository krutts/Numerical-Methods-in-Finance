#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include <chrono>
#include <random>
#define E 2.718281828459045

using namespace std;

double lambda, mu;
int no_of_trials;
unsigned seed = (unsigned)std::chrono::system_clock::now().time_since_epoch().count();
std::default_random_engine generator(seed);

// u.i.i.d. generator
double get_uniform() {
	std::uniform_real_distribution <double> distribution(0.0, 1.0);
	double number = distribution(generator);
	return (number);
}

//Normal cdf
double N(const double& z) { 
	if (z > 6.0) { return 1.0; }; // this guards against overflow 
	if (z < -6.0) { return 0.0; };
	double b1 = 0.31938153;
	double b2 = -0.356563782;
	double b3 = 1.781477937;
	double b4 = -1.821255978;
	double b5 = 1.330274429;
	double p = 0.2316419;
	double c2 = 0.3989423;
	double a = fabs(z);
	double t = 1.0 / (1.0 + a*p);
	double b = c2*exp((-z)*(z / 2.0));
	double n = ((((b5*t + b4)*t + b3)*t + b2)*t + b1)*t;
	n = 1.0 - b*n;
	if (z < 0.0) n = 1.0 - n;
	return n;
}

// returns exponential value of the number in parameter
double exp(double n) {
	return pow(E, n);
}

//inverse of Normal cdf
double inverseNormal(double lambda, double mu, double X) {
	return (N(sqrt(lambda / X)*(X / mu - 1)) + exp(2 * lambda / mu)*N(-sqrt(lambda / X)*(X / mu - 1)));
}


double monte_carlo_V(int no_of_trials, double lambda, double mu) {
	double phiXt = 0;
	for (int i = 0; i < no_of_trials; i++)
	{
		double x = get_uniform();
		double Xt = inverseNormal(lambda, mu, x); //since Xt is to be normalized
		phiXt += (1 / (1 + exp(1-Xt)))/no_of_trials; //using the logistic function given
	}
	return phiXt;
}

double monte_carlo_standard_deviation(int no_of_trials, double lambda, double mu) {
	double phiXt_MCEstimator = 0;
	double *phiXt = new double[no_of_trials + 1];
	double sd_phi = 0;
	for (int i = 0; i < no_of_trials; i++)
	{
		double x = get_uniform();
		double Xt = inverseNormal(lambda, mu, x);  //since Xt is to be normalized
		phiXt[i] = 1 / (1 + exp(1-Xt)); //storing phi values to use to calculate variance
		phiXt_MCEstimator += phiXt[i] / no_of_trials; //calculating the sum to calculate the estimated mean
	}
	for (int i = 0; i < no_of_trials; i++)
		sd_phi += pow(phiXt[i] - phiXt_MCEstimator, 2) / (no_of_trials - 1); //calculating the estimated variance using the Monte Carlo variance estimator formula
	return sqrt(sd_phi/no_of_trials);
}

double antithetic_variates_V(int no_of_trials, double lambda, double mu) {
	double phiXt = 0;
	for (int i = 0; i < no_of_trials/2; i++) //runs up to 5000 trials, half of the case in MC since we are generating 2 normal variables in every simulation here
	{
		double x = get_uniform();
		double Xt = inverseNormal(lambda, mu, x);
		double Xtcomp = inverseNormal(lambda, mu, 1-x); //using the complement of the uniform random variable generated to obtain another standard normal variable
		phiXt += (1 / (1 + exp(1-Xt)) + 1 / (1 + exp(1-Xtcomp)))/no_of_trials; //using both the above generated normal variables in the logistic function payoff
	}
	return phiXt;
}

double antithetic_variates_standard_deviation(int no_of_trials, double lambda, double mu) {
	double phiXt_AVEstimator = 0;
	double *phiXt = new double[no_of_trials + 1];
	double sd_phi = 0;
	for (int i = 0; i < no_of_trials/2; i++) //runs up to 5000 trials, half of the case in MC since we are generating 2 normal variables in every simulation here
	{
		double x = get_uniform();
		double Xt = inverseNormal(lambda, mu, x);
		double Xtcomp = inverseNormal(lambda, mu, 1 - x); //using the complement of the uniform random variable generated to obtain another standard normal variable
		phiXt[i] = 1 / (1 + exp(1-Xt)); //storing phi values from Xt in the first half of the matrix to use to calculate variance
		phiXt[no_of_trials/2 + i] = 1 / (1 + exp(1-Xtcomp)); //storing phi values from 1-Xt in the second half of the matrix to use to calculate variance
		phiXt_AVEstimator += (phiXt[i] + phiXt[no_of_trials / 2 + i] )/ no_of_trials; //calculating the sum to calculate the estimated mean
	}
	for (int i = 0; i < no_of_trials/2; i++)
		sd_phi += pow(0.5*(phiXt[i] + phiXt[no_of_trials / 2 + i]) - phiXt_AVEstimator, 2) / (no_of_trials - 1); //calculating the estimated variance using the given formula
	return sqrt(sd_phi/0.5/no_of_trials);
}



int main(int argc, char* argv[]) {

	no_of_trials = 10000;
	lambda = 1;
	mu = 1; //lambda,mu > 0

	cout << "--------------------------------" << endl;
	cout << "Variance Reduction Techniques: Antithetic Variates" << endl;
	cout << "Number of Trials = " << no_of_trials << endl;
	cout << "--------------------------------" << endl;
	cout << "Value of the derivative using Monte Carlo Simulation = " << monte_carlo_V(no_of_trials, lambda, mu) << endl;
	cout << "Standard Deviation of the derivative using Monte Carlo Simulation = " << monte_carlo_standard_deviation(no_of_trials, lambda, mu) << endl;
	cout << "Value of the derivative using Antithetic Variates = " << antithetic_variates_V(no_of_trials, lambda, mu) << endl;
	cout << "Standard Deviation of the derivative using Antithetic Variates = " << antithetic_variates_standard_deviation(no_of_trials, lambda, mu) << endl;

	return 0;
}