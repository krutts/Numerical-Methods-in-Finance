#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include <chrono>
#include <random>
#define E 2.718281828459045

using namespace std;

double risk_free_rate, drift_rate, strike_price, initial_stock_price, expiration_time, volatility, delT;
int no_of_trials;

unsigned seed = (unsigned)std::chrono::system_clock::now().time_since_epoch().count();
std::default_random_engine generator(seed);

// u.i.i.d. generator
double get_uniform() {
	std::uniform_real_distribution <double> distribution(0.0, 1.0);
	double number = distribution(generator);
	return (number);
}

double max(double a, double b) {
	return (b < a) ? a : b;
}

// returns exponential value of the number in parameter
double exp(double n) {
	return pow(E, n);
}

//simulates GBM path
double simulatePath(double initial_stock_price, double drift_rate, double volatility, double delT, double expiration_time, int flag) {
	int n = round(expiration_time / delT); //to ensure n is even
	double* S = new double[n + 1]; //stores simulated stock price values
	S[0] = initial_stock_price;

	for (int i = 1; i <= n/2; i++) {
		double x = get_uniform();
		double y = get_uniform();
		double a = sqrt(-2.0*log(x)) * cos(6.283185307999998*y); 
		double b = sqrt(-2.0*log(x)) * sin(6.283185307999998*y); //using Box-Mueller method to generate random values
		S[i] = S[i - 1] * exp((drift_rate - 0.5*pow(volatility, 2))*delT + volatility*sqrt(delT)*a);
		S[i + 1] = S[i] * exp((drift_rate - 0.5*pow(volatility, 2))*delT + volatility*sqrt(delT)*b); //using the normally distributed random values generated above to simulate stock prices
	}
	if (flag == 1) //if Asian
	{
		double average_price = 0;
		for (int i = 0; i < n; i++)
			average_price += S[i] / n;
		return average_price;
	}
	return S[n];
}

double option_price_monte_carlo_european_call(double initial_stock_price, double drift_rate, double volatility, double delT, double expiration_time, double strike_price, double risk_free_rate, int no_of_trials)
{
	double option_price = 0;
	for (int i = 0; i < no_of_trials; i++)
		option_price += exp(-risk_free_rate*expiration_time)*max(0, simulatePath(initial_stock_price, drift_rate, volatility, delT, expiration_time, 0) - strike_price) / no_of_trials; //calculating call option prices from payoffs with simulated stock prices
	return option_price;
}

double option_price_monte_carlo_european_put(double initial_stock_price, double drift_rate, double volatility, double delT, double expiration_time, double strike_price, double risk_free_rate, int no_of_trials)
{
	double option_price = 0;
	for (int i = 0; i < no_of_trials; i++)
		option_price += exp(-risk_free_rate*expiration_time)*max(0, strike_price - simulatePath(initial_stock_price, drift_rate, volatility, delT, expiration_time, 0)) / no_of_trials;  //calculating Eput option prices from payoffs with simulated stock prices
	return option_price; 
}

double option_price_monte_carlo_asian_call(double initial_stock_price, double drift_rate, double volatility, double delT, double expiration_time, double strike_price, double risk_free_rate, int no_of_trials)
{
	double option_price = 0;
	for (int i = 0; i < no_of_trials; i++)
		option_price += exp(-risk_free_rate*expiration_time)*max(0, simulatePath(initial_stock_price, drift_rate, volatility, delT, expiration_time, 1) - strike_price) / no_of_trials;  //calculating Asian call option prices from payoffs with simulated stock prices - in this case the average of the simulated stock prices is considered
	return option_price;
}

double variance_monte_carlo_european_call(double initial_stock_price, double drift_rate, double volatility, double delT, double expiration_time, double strike_price, double risk_free_rate, int no_of_trials)
{
	double MC_estimator = 0;
	double *option_prices = new double[no_of_trials + 1];
	double variance = 0;
	for (int i = 0; i < no_of_trials; i++) {
		option_prices[i] = exp(-risk_free_rate*expiration_time)*max(0, simulatePath(initial_stock_price, drift_rate, volatility, delT, expiration_time, 0) - strike_price); //storing simulated option prices in an array to use to calculate variance and the Monte Carlo estimated mean
		MC_estimator += option_prices[i] / no_of_trials; //adding all the values to find the Monte Carlo estimated mean
	}
	for (int i = 0; i < no_of_trials; i++)
		variance += pow(option_prices[i] - MC_estimator, 2) / (no_of_trials - 1); //using the Monte Carlo variance formula to find the estimated call option variance
	return variance;
}

double variance_monte_carlo_european_put(double initial_stock_price, double drift_rate, double volatility, double delT, double expiration_time, double strike_price, double risk_free_rate, int no_of_trials)
{
	double MC_estimator = 0;
	double *option_prices = new double[no_of_trials + 1];
	double variance = 0;
	for (int i = 0; i < no_of_trials; i++) {
		option_prices[i] = exp(-risk_free_rate*expiration_time)*max(0, strike_price - simulatePath(initial_stock_price, drift_rate, volatility, delT, expiration_time, 0)); //storing simulated option prices in an array to use to calculate variance and the Monte Carlo estimated mean
		MC_estimator += option_prices[i] / no_of_trials; //adding all the values to find the Monte Carlo estimated mean
	}
	for (int i = 0; i < no_of_trials; i++)
		variance += pow(option_prices[i] - MC_estimator, 2) / (no_of_trials - 1); //using the Monte Carlo variance formula to find the estimated put option variance
	return variance;
}

double variance_monte_carlo_asian_call(double initial_stock_price, double drift_rate, double volatility, double delT, double expiration_time, double strike_price, double risk_free_rate, int no_of_trials)
{
	double MC_estimator = 0;
	double *G = new double[no_of_trials + 1];
	double variance = 0;
	for (int i = 0; i < no_of_trials; i++) {
		G[i] = exp(-risk_free_rate*expiration_time)*max(0, simulatePath(initial_stock_price, drift_rate, volatility, delT, expiration_time, 1) - strike_price); //storing simulated option prices in an array to use to calculate variance and the Monte Carlo estimated mean
		MC_estimator += G[i] / no_of_trials; //adding all the values to find the Monte Carlo estimated mean
	}
	for (int i = 0; i < no_of_trials; i++)
		variance += pow(G[i] - MC_estimator, 2) / (no_of_trials - 1); //using the Monte Carlo variance formula to find the estimated put option variance
	return variance;
}


int main()
{
	drift_rate = 0.1;
	volatility = 0.15;
	initial_stock_price = 15;
	strike_price = 16;
	expiration_time = 1;
	risk_free_rate = 0;
	no_of_trials = 10000;
	delT = 0.4167; //=1/24 = half month

	cout << "--------------------------------" << endl;
	cout << "Monte Carlo Pricing for European and Asian Calls" << endl;
	cout << "Expiration Time (Years) = " << expiration_time << endl;
	cout << "Drift Rate = " << drift_rate << endl;
	cout << "Risk Free Interest Rate = " << risk_free_rate << endl;
	cout << "Volatility (%age of stock value) = " << volatility * 100 << endl;
	cout << "Initial Stock Price = " << initial_stock_price << endl;
	cout << "Strike Price = " << strike_price << endl;
	cout << "Number of Trials = " << no_of_trials << endl;
	cout << "--------------------------------" << endl;
	cout << "European Call Option Price = " << option_price_monte_carlo_european_call(initial_stock_price, drift_rate, volatility, delT, expiration_time, strike_price, risk_free_rate, no_of_trials) << endl;
	cout << "European Put Option Price = " << option_price_monte_carlo_european_put(initial_stock_price, drift_rate, volatility, delT, expiration_time, strike_price, risk_free_rate, no_of_trials) << endl;
	cout << "European Call Variance = " << variance_monte_carlo_european_call(initial_stock_price, drift_rate, volatility, delT, expiration_time, strike_price, risk_free_rate, no_of_trials) << endl;
	cout << "European Put Variance = " << variance_monte_carlo_european_put(initial_stock_price, drift_rate, volatility, delT, expiration_time, strike_price, risk_free_rate, no_of_trials) << endl;
	cout << "Asian Call Option Price = " << option_price_monte_carlo_asian_call(initial_stock_price, drift_rate, volatility, delT, expiration_time, strike_price, risk_free_rate, no_of_trials) << endl;
	cout << "Asian Call Variance = " << variance_monte_carlo_asian_call(initial_stock_price, drift_rate, volatility, delT, expiration_time, strike_price, risk_free_rate, no_of_trials) << endl;

	return 0;
}

