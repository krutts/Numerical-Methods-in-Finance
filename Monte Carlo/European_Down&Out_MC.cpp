
//European Down-and-Out Discrete Barrier Options Pricing via Monte Carlo Simulation

#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <random>
#include <chrono>
using namespace std;

double risk_free_rate, strike_price, initial_stock_price, expiration_time, volatility, barrier_price;
int no_of_trials, no_of_barriers;

unsigned seed = (unsigned)std::chrono::system_clock::now().time_since_epoch().count();
std::default_random_engine generator(seed);

double max(double a, double b) {
	return (b < a) ? a : b;
}

// u.i.i.d. generator
double get_uniform()
{
	std::uniform_real_distribution <double> distribution(0.0, 1.0);
	double number = distribution(generator);
	return (number);
}

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

double option_price_call_black_scholes(const double& S,       // spot (underlying) price
	const double& K,       // strike (exercise) price,
	const double& r,       // interest rate
	const double& sigma,   // volatility 
	const double& time) {  // time to maturity 
	double time_sqrt = sqrt(time);
	double d1 = (log(S / K) + r*time) / (sigma*time_sqrt) + 0.5*sigma*time_sqrt;
	double d2 = d1 - (sigma*time_sqrt);
	return S*N(d1) - K*exp(-r*time)*N(d2);
}


int main(int argc, char* argv[])
{
	sscanf(argv[1], "%lf", &expiration_time);
	sscanf(argv[2], "%lf", &risk_free_rate);
	sscanf(argv[3], "%lf", &volatility);
	sscanf(argv[4], "%lf", &initial_stock_price);
	sscanf(argv[5], "%lf", &strike_price);
	sscanf(argv[6], "%d", &no_of_trials);
	sscanf(argv[7], "%d", &no_of_barriers);
	sscanf(argv[8], "%lf", &barrier_price);


	cout << "--------------------------------" << endl;
	cout << "European Down-and-Out Discrete Barrier Options Pricing via Monte Carlo Simulation" << endl;
	cout << "Expiration Time (Years) = " << expiration_time << endl;
	cout << "Risk Free Interest Rate = " << risk_free_rate << endl;
	cout << "Volatility (%age of stock value) = " << volatility * 100 << endl;
	cout << "Initial Stock Price = " << initial_stock_price << endl;
	cout << "Strike Price = " << strike_price << endl;
	cout << "Barrier Price = " << barrier_price << endl;
	cout << "Number of Trials = " << no_of_trials << endl;
	cout << "Number of Discrete Barriers = " << no_of_barriers << endl;
	cout << "--------------------------------" << endl;

	double *simulated_barrier_call_option_price = new double[no_of_trials + 1];
	double *squared_simulated_barrier_call_option_price = new double[no_of_trials + 1];
	double *simulated_vanilla_call_option_price = new double[no_of_trials + 1];
	double *squared_simulated_vanilla_call_option_price = new double[no_of_trials + 1];
	double *predicted_barrier_call_option_price = new double[no_of_trials + 1];
	double covariance = 0.0, price_simulated_barrier_call_option = 0.0, variance_simulated_barrier_call_option = 0.0, price_simulated_vanilla_call_option = 0.0, variance_simulated_vanilla_call_option = 0.0;
	double delta_T = expiration_time / ((double)no_of_barriers);
	double delta_R = (risk_free_rate - 0.5*pow(volatility, 2))*delta_T;
	double delta_SD = volatility*sqrt(delta_T);

	for (int k = 0; k < no_of_trials; k++)
	{
		// by sharing random variables we create 4 paths
		double current_stock_price1 = initial_stock_price;
		double current_stock_price2 = initial_stock_price;
		double current_stock_price3 = initial_stock_price;
		double current_stock_price4 = initial_stock_price;

		//these will store the barrier-breach-adjusted price
		double S1 = 0, S2 = 0, S3 = 0, S4 = 0;

		//these will track whether the price path ever hits the barrier during the many simulations
		bool current_stock_price1_hit_barrier = false;
		bool current_stock_price2_hit_barrier = false;
		bool current_stock_price3_hit_barrier = false;
		bool current_stock_price4_hit_barrier = false;


		for (int j = 0; j < no_of_barriers; j++) {
			// create the unit normal variates using the Box-Muller Transform

			double x = get_uniform();
			double y = get_uniform();
			double a = sqrt(-2.0*log(x)) * cos(6.283185307999998*y);
			double b = sqrt(-2.0*log(x)) * sin(6.283185307999998*y);

			//if these prices ever hit the barrier, the corresponding hit_barrier variable will be set to true
			current_stock_price1 = current_stock_price1*exp(delta_R + delta_SD*a);
			current_stock_price1_hit_barrier = (current_stock_price1 <= barrier_price ? true  : current_stock_price1_hit_barrier);

			current_stock_price2 = current_stock_price2*exp(delta_R - delta_SD*a);
			current_stock_price2_hit_barrier = (current_stock_price2 <= barrier_price ? true : current_stock_price2_hit_barrier);

			current_stock_price3 = current_stock_price3*exp(delta_R + delta_SD*b);
			current_stock_price3_hit_barrier = (current_stock_price3 <= barrier_price ? true : current_stock_price3_hit_barrier);

			current_stock_price4 = current_stock_price4*exp(delta_R - delta_SD*b);
			current_stock_price4_hit_barrier = (current_stock_price4 <= barrier_price ? true : current_stock_price4_hit_barrier);


		}

		simulated_vanilla_call_option_price[k]= (max(0.0, current_stock_price1 - strike_price) +
			max(0.0, current_stock_price2 - strike_price) +
			max(0.0, current_stock_price3 - strike_price) +
			max(0.0, current_stock_price4 - strike_price)) / 4.0;

		price_simulated_vanilla_call_option += simulated_vanilla_call_option_price[k]; //stores the vanilla option price

		squared_simulated_vanilla_call_option_price[k]= (pow(max(0.0, current_stock_price1 - strike_price), 2) +
			pow(max(0.0, current_stock_price2 - strike_price), 2) +
			pow(max(0.0, current_stock_price3 - strike_price), 2) +
			pow(max(0.0, current_stock_price4 - strike_price), 2)) / 4.0;

		variance_simulated_vanilla_call_option += squared_simulated_vanilla_call_option_price[k]; //stores the variance of vanilla option

		//if breached, invalidate the price path by setting it to 0
		S1 = current_stock_price1_hit_barrier == true ? 0 : current_stock_price1;
		S2 = current_stock_price2_hit_barrier == true ? 0 : current_stock_price2;
		S3 = current_stock_price3_hit_barrier == true ? 0 : current_stock_price3;
		S4 = current_stock_price4_hit_barrier == true ? 0 : current_stock_price4;

		//calculating the average payoff from the simulated price paths that didn't hit the barrier
		simulated_barrier_call_option_price[k]= (max(0.0, S1 - strike_price) +
			max(0.0, S2 - strike_price) +
			max(0.0, S3 - strike_price) +
			max(0.0, S4 - strike_price)) / 4.0;

		price_simulated_barrier_call_option += simulated_barrier_call_option_price[k]; //stores the barrier option price

		squared_simulated_barrier_call_option_price[k]= (pow(max(0.0, S1 - strike_price), 2) +
			pow(max(0.0, S2 - strike_price), 2) +
			pow(max(0.0, S3 - strike_price), 2) +
			pow(max(0.0, S4 - strike_price), 2)) / 4.0;

		variance_simulated_barrier_call_option += squared_simulated_barrier_call_option_price[k]; //stores the variance of barrier option

		covariance += (max(0.0, current_stock_price1 - strike_price)*max(0.0, S1 - strike_price) +
			max(0.0, current_stock_price2 - strike_price)*max(0.0, S2 - strike_price) +
			max(0.0, current_stock_price3 - strike_price)*max(0.0, S3 - strike_price) +
			max(0.0, current_stock_price4 - strike_price)*max(0.0, S4 - strike_price)) / 4.0;
	}

	price_simulated_barrier_call_option = exp(-risk_free_rate*expiration_time)*(price_simulated_barrier_call_option / ((double)no_of_trials)); //discounted payoff of barrier option
	price_simulated_vanilla_call_option = exp(-risk_free_rate*expiration_time)*(price_simulated_vanilla_call_option / ((double)no_of_trials)); //discounted payoff of vanilla option

	double barrier_standard_error = sqrt((variance_simulated_barrier_call_option - pow(price_simulated_barrier_call_option, 2)) / (no_of_trials - 1) / no_of_trials);
	double vanilla_standard_error = sqrt((variance_simulated_vanilla_call_option - pow(price_simulated_vanilla_call_option, 2)) / (no_of_trials - 1) / no_of_trials);

	double beta = covariance / variance_simulated_vanilla_call_option;
	double alpha = price_simulated_barrier_call_option - beta*price_simulated_barrier_call_option;
	double sum_of_squared_residuals = 0.0;

	for (int k = 0; k < no_of_trials; k++)
	{
		predicted_barrier_call_option_price[k] = alpha + beta*simulated_vanilla_call_option_price[k];
		sum_of_squared_residuals += pow(predicted_barrier_call_option_price[k] - simulated_barrier_call_option_price[k], 2);
	}



	cout << "The average Barrier Call Price via explicit simulation of price paths = " << price_simulated_barrier_call_option << endl;
	cout << "Standard error in Barrier Call Option Price = " << barrier_standard_error << endl;
	cout << "Actual error in Barrier Call Option Price = " << price_simulated_barrier_call_option - 4.647650 << endl;
	cout << "The average Vanilla Call Price via explicit simulation of price paths = " << price_simulated_vanilla_call_option << endl;
	cout << "Standard error in Vanilla Call Option Price = " << vanilla_standard_error << endl;
	cout << "Actual error in Vanilla Call Option Price = " << price_simulated_vanilla_call_option - option_price_call_black_scholes(initial_stock_price, strike_price, risk_free_rate, volatility, expiration_time) << endl;
	cout << "R squared = " << 1-sum_of_squared_residuals/variance_simulated_barrier_call_option << endl;
	cout << "--------------------------------" << endl;
}
/**/