//Asian option using Monte Carlo Simulation
#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <random>
#include <chrono>
using namespace std;

double risk_free_rate, strike_price, initial_stock_price, expiration_time, volatility, average_price = 0.0, variance = 0.0;
int no_of_trials;
double* simulated_prices = new double[no_of_trials + 1];


unsigned seed = (unsigned)std::chrono::system_clock::now().time_since_epoch().count();
std::default_random_engine generator(seed);

double max(double a, double b) {
	return (b < a) ? a : b;
}

// u.i.i.d. generator
double get_uniform(double mu)
{
	std::uniform_real_distribution <double> distribution(mu, 1.0);
	double number = distribution(generator);
	return (number);
}

double simulate_asian_call_monte_carlo(double initial_stock_price, double risk_free_rate, double volatility, double expiration_time, int no_of_trials) {
	double avg_price = 0.0;
	double R = (risk_free_rate - 0.5*pow(volatility, 2))*expiration_time;
	double SD = volatility*sqrt(expiration_time);
	for (int i = 0; i < no_of_trials; i++) {
		// generate unit-normals using Box-Muller Transform
		double x = get_uniform(0);
		double y = get_uniform(0);
		double a = sqrt(-2.0*log(x)) * cos(6.283185307999998*y);
		double b = sqrt(-2.0*log(x)) * sin(6.283185307999998*y);

		
		double S1 = initial_stock_price * exp(R + SD*a);
		double S2 = initial_stock_price * exp(R - SD*a);
		double S3 = initial_stock_price * exp(R + SD*b);
		double S4 = initial_stock_price * exp(R - SD*b);

		avg_price += (S1 + S2 + S3 + S4) / 4;
	}
	
	return avg_price;
}

double variance_monte_carlo_asian_call(double initial_stock_price, double volatility, double expiration_time, double strike_price, double risk_free_rate, int no_of_trials)
{
	double average_price = 0;
	double *prices = new double[no_of_trials + 1];
	double variance = 0;
	for (int i = 0; i < no_of_trials; i++) {
		prices[i] = exp(-risk_free_rate*expiration_time)*((simulate_asian_call_monte_carlo(initial_stock_price, risk_free_rate, volatility, expiration_time, no_of_trials) / (double)(no_of_trials + 1)) - strike_price); //storing simulated option prices in an array to use to calculate variance and the Monte Carlo estimated mean
		average_price += prices[i] / no_of_trials; //adding all the values to find the Monte Carlo estimated mean
	}
	for (int i = 0; i < no_of_trials; i++)
		variance += pow(prices[i] - average_price, 2) / (no_of_trials - 1); //using the Monte Carlo variance formula to find the estimated put option variance
	return variance;
}


int main(int argc, char* argv[])
{
	sscanf(argv[1], "%lf", &expiration_time);
	sscanf(argv[2], "%lf", &risk_free_rate);
	sscanf(argv[3], "%lf", &volatility);
	sscanf(argv[4], "%lf", &initial_stock_price);
	sscanf(argv[5], "%lf", &strike_price);
	sscanf(argv[6], "%d", &no_of_trials);



	cout << "--------------------------------" << endl;
	cout << "Asian Call Option Pricing via Monte Carlo Simulation" << endl;
	cout << "Expiration Time (Years) = " << expiration_time << endl;
	cout << "Risk Free Interest Rate = " << risk_free_rate << endl;
	cout << "Volatility (%age of stock value) = " << volatility * 100 << endl;
	cout << "Initial Stock Price = " << initial_stock_price << endl;
	cout << "Strike Price = " << strike_price << endl;
	cout << "Number of Trials = " << no_of_trials << endl;
	cout << "--------------------------------" << endl;


	cout << "Average Asian Call Option Price = " << exp(-risk_free_rate*expiration_time)*((simulate_asian_call_monte_carlo(initial_stock_price, risk_free_rate, volatility, expiration_time, no_of_trials) / (double)(no_of_trials + 1)) - strike_price) << endl;
	cout << "Standard error in Asian Call Option Price = " << sqrt(variance_monte_carlo_asian_call(initial_stock_price, volatility, expiration_time, strike_price, risk_free_rate, no_of_trials) / no_of_trials) << endl;

	cout << "--------------------------------" << endl;
}
*/


