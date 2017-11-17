// Pricing an European Option using Simulation

#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <random>
#include <chrono>
#define UNIT_STEP(x) ((x)>0?(1):(0))
using namespace std;

double risk_free_rate, strike_price1, strike_price2, strike_price3, initial_stock_price, expiration_time, volatility;
int no_of_trials, m1, m2, m3, m4;


double max(double a, double b) {
	return (b < a) ? a : b;
}

unsigned seed = (unsigned)std::chrono::system_clock::now().time_since_epoch().count();
std::default_random_engine generator(seed);

// u.i.i.d. generator
double get_uniform()
{
	std::uniform_real_distribution <double> distribution(0.0, 1.0);
	double number = distribution(generator);
	return (number);
}

// unit-normal i.i.d. generator
double get_gaussian()
{
	return (sqrt(-2.0*log(get_uniform()))*cos(6.283185307999998*get_uniform()));
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
};

double d(double initial_stock_price, double volatility, double expiration_time, double strike_price, double risk_free_rate) {
	return (log(initial_stock_price / strike_price) + risk_free_rate*expiration_time) / (volatility*sqrt(expiration_time)) + 0.5*volatility*sqrt(expiration_time);
}

double simulatePath(double initial_stock_price, double volatility, double expiration_time, int flag) {
	double total = 0;
	double S = 0;
	for (int i = 1; i <= no_of_trials; i++) {
		double x = get_uniform();
		double y = get_uniform();
		double a = sqrt(-2.0*log(x)) * cos(6.283185307999998*y);
		double b = sqrt(-2.0*log(x)) * sin(6.283185307999998*y); //using Box-Mueller method to generate random values
		double S1 = initial_stock_price*exp((risk_free_rate - 0.5*pow(volatility, 2))*expiration_time + volatility*sqrt(expiration_time)*a);
		double S2 = initial_stock_price*exp((risk_free_rate - 0.5*pow(volatility, 2))*expiration_time + volatility*sqrt(expiration_time)*b); //using the normally distributed random values generated above to simulate stock 
		total += (S1 + S2) / 2;
		S = (S1 + S2) / 2;
	}
	if (flag == 1) //if Asian
		return total/no_of_trials;
	return S;
}

double simulatePathOnMDates(double initial_stock_price, double volatility, double expiration_time, int m) {
	double total = 0;
	double S = 0;
	for (int i = 1; i <= no_of_trials; i++) {
		double x = get_uniform();
		double y = get_uniform();
		double a = sqrt(-2.0*log(x)) * cos(6.283185307999998*y);
		double b = sqrt(-2.0*log(x)) * sin(6.283185307999998*y); //using Box-Mueller method to generate random values
		double S1 = initial_stock_price*exp((risk_free_rate - 0.5*pow(volatility, 2))*expiration_time + volatility*sqrt(expiration_time)*a);
		double S2 = initial_stock_price*exp((risk_free_rate - 0.5*pow(volatility, 2))*expiration_time + volatility*sqrt(expiration_time)*b); //using the normally distributed random values generated above to simulate stock 
		if(i%m==0) 
			total += (S1 + S2) / 2;
	}
	return total / no_of_trials;
}

double delta_BSM_european_put(double initial_stock_price, double volatility, double expiration_time, double strike_price, double risk_free_rate, int no_of_trials)
{
	return exp(-risk_free_rate*expiration_time)*(N(d(initial_stock_price, volatility, expiration_time, strike_price, risk_free_rate))-1);
}

double gamma_BSM_european_put(double initial_stock_price, double volatility, double expiration_time, double strike_price, double risk_free_rate, int no_of_trials)
{
	return exp(-risk_free_rate*expiration_time)*exp(-pow(d(initial_stock_price, volatility, expiration_time, strike_price, risk_free_rate),2)/2)/(initial_stock_price*volatility*sqrt(expiration_time)*2.5066);
}

double option_price_monte_carlo_european_put(double initial_stock_price, double volatility, double expiration_time, double strike_price, double risk_free_rate, int no_of_trials)
{
	double option_price = 0;
	for (int i = 0; i < no_of_trials; i++)
		option_price += exp(-risk_free_rate*expiration_time)*max(0, strike_price - simulatePath(initial_stock_price, volatility, expiration_time, 0)) / no_of_trials;  //calculating Eput option prices from payoffs with simulated stock prices
	return option_price;
}

double option_price_monte_carlo_asian_call(double initial_stock_price, double volatility, double expiration_time, double strike_price, double risk_free_rate, int no_of_trials)
{
	double option_price = 0;
	for (int i = 0; i < no_of_trials; i++)
		option_price += exp(-risk_free_rate*expiration_time)*max(0, simulatePath(initial_stock_price, volatility, expiration_time, 1) - strike_price) / no_of_trials;  //calculating Asian call option prices from payoffs with simulated stock prices - in this case the average of the simulated stock prices is considered
	return option_price;
}

double delta_monte_carlo_european_put(double initial_stock_price, double volatility, double expiration_time, double strike_price, double risk_free_rate, int no_of_trials)
{
	double St = 0;
	double delta = 0;
	for (int i = 0; i < no_of_trials; i++)
	{
		St = simulatePath(initial_stock_price, volatility, expiration_time, 0);
		if(strike_price>St)
			delta += exp(-risk_free_rate*expiration_time)*St/initial_stock_price/no_of_trials;  //calculating Eput option prices from payoffs with simulated stock prices
	}
	return delta;
}

double delta_pathwise_european_put(double initial_stock_price, double volatility, double expiration_time, double strike_price, double risk_free_rate, int no_of_trials)
{
	double St = 0;
	double delta = 0;
	for (int i = 0; i < no_of_trials; i++)
	{
		St = simulatePath(initial_stock_price, volatility, expiration_time, 0);
		if (strike_price>St)
			delta += exp(-risk_free_rate*expiration_time)*St / initial_stock_price / no_of_trials;  //calculating Eput option prices from payoffs with simulated stock prices
	}
	return delta-1;
}

double delta_LRM_european_put(double initial_stock_price, double volatility, double expiration_time, double strike_price, double risk_free_rate, int no_of_trials)
{
	double St = 0;
	double delta = 0;
	double X = 0;
	for (int i = 0; i < no_of_trials; i++)
	{
		St = simulatePath(initial_stock_price, volatility, expiration_time, 0);
		double Z = get_gaussian();
		X+= max(strike_price - St, 0)*Z / initial_stock_price / volatility / sqrt(expiration_time);
	}
	delta = exp(-risk_free_rate*expiration_time)*X / no_of_trials;  //calculating Eput option prices from payoffs with simulated stock prices
	return delta;
}

double gamma_LRM_european_put(double initial_stock_price, double volatility, double expiration_time, double strike_price, double risk_free_rate, int no_of_trials)
{
	double St = 0;
	double gamma = 0;
	double X = 0;
	for (int i = 0; i < no_of_trials; i++)
	{
		St = simulatePath(initial_stock_price, volatility, expiration_time, 0);
		double Z = get_gaussian();
		X += max(strike_price - St, 0)*(((pow(Z,2) -1) / pow(initial_stock_price*volatility,2) / sqrt(expiration_time)) - (Z / pow(initial_stock_price, 2) / volatility / sqrt(expiration_time)));
	}
	gamma += exp(-risk_free_rate*expiration_time)*X / no_of_trials;  //calculating Eput option prices from payoffs with simulated stock prices
	return gamma;
}

double delta_LRM_asian_call(double initial_stock_price, double volatility, double expiration_time, double strike_price, double risk_free_rate, int no_of_trials, int m)
{
	double Sbar = 0;
	double delta = 0;
	double X = 0;
	for (int i = 0; i < no_of_trials; i++)
	{
		Sbar = simulatePathOnMDates(initial_stock_price, volatility, expiration_time, m);
		double Z = get_gaussian();
		X += max(Sbar - strike_price, 0) + Z / initial_stock_price / volatility / sqrt(expiration_time);
	}
	delta += exp(-risk_free_rate*expiration_time)*X / no_of_trials;  //calculating Eput option prices from payoffs with simulated stock prices
	return delta;
}

double delta_pathwise_asian_call(double initial_stock_price, double volatility, double expiration_time, double strike_price, double risk_free_rate, int no_of_trials, int m)
{
	double option_price = 0;
	double Sbar = 0;
	double delta = 0;
	for (int i = 0; i < no_of_trials; i++)
	{
		Sbar = simulatePathOnMDates(initial_stock_price, volatility, expiration_time, m);
		if (strike_price>Sbar)
			delta += exp(-risk_free_rate*expiration_time)*Sbar / initial_stock_price / no_of_trials;  //calculating Eput option prices from payoffs with simulated stock prices
	}
	return delta;
}

void variance_delta_asian_call(double initial_stock_price, double volatility, double expiration_time, double strike_price, double risk_free_rate, int no_of_trials, int m) {
	double* D_LRM = new double[no_of_trials + 1];
	double* D_PM = new double[no_of_trials + 1];
	double LRM_Estimator = 0, PM_Estimator = 0, Var_LRM = 0, Var_PM = 0;
	for (int i = 0; i < no_of_trials; i++) {
		D_LRM[i] = delta_LRM_asian_call(initial_stock_price, volatility, expiration_time, strike_price, risk_free_rate, no_of_trials, m);
		D_PM[i] = delta_pathwise_asian_call(initial_stock_price, volatility, expiration_time, strike_price, risk_free_rate, no_of_trials, m);
		LRM_Estimator += D_LRM[i]/no_of_trials;
		PM_Estimator += D_PM[i]/no_of_trials;
	}
	for (int i = 0; i < no_of_trials; i++) {
		Var_LRM += pow((D_LRM[i] - LRM_Estimator), 2) / (no_of_trials - 1);
		Var_PM += pow((D_PM[i] - PM_Estimator), 2) / (no_of_trials - 1);
	}
	cout << "Variance in Delta using LRM with " << m << " dates = " << Var_LRM << endl;
	cout << "Variance in Delta using Pathwise Method with " << m << " dates = " << Var_PM << endl;
}

int main(int argc, char* argv[])
{
	initial_stock_price = 100;
	volatility = 0.3;
	risk_free_rate = 0.05;
	expiration_time = 0.1;
	strike_price1 = 90;
	strike_price2 = 100;
	strike_price3 = 110;
	no_of_trials = 250;
	m1 = 32;
	m2 = 64;
	m3 = 96;
	m4 = 128;

	double R = (risk_free_rate - 0.5*pow(volatility, 2))*expiration_time;
	double SD = volatility*sqrt(expiration_time);

	cout << "--------------------------------" << endl;
	cout << "Greeks Estimation via Monte Carlo Simulation" << endl;
	cout << "Expiration Time (Years) = " << expiration_time << endl;
	cout << "Risk Free Interest Rate = " << risk_free_rate << endl;
	cout << "Volatility (%age of stock value) = " << volatility * 100 << endl;
	cout << "Initial Stock Price = " << initial_stock_price << endl;
	cout << "Number of Trials = " << no_of_trials << endl;
	cout << "--------------------------------" << endl;

	/*cout << "European Put Option with K = 90:" << endl;
	cout << "Price = " << option_price_monte_carlo_european_put(initial_stock_price, volatility, expiration_time, strike_price1, risk_free_rate, no_of_trials) << endl;
	cout << "Delta from BSM = " << delta_BSM_european_put(initial_stock_price, volatility, expiration_time, strike_price1, risk_free_rate, no_of_trials) << endl;
	cout << "Delta using Pathwise Estimators = " << delta_pathwise_european_put(initial_stock_price, volatility, expiration_time, strike_price1, risk_free_rate, no_of_trials) << endl;
	cout << "Delta using LRM = " << delta_LRM_european_put(initial_stock_price, volatility, expiration_time, strike_price1, risk_free_rate, no_of_trials) << endl;
	cout << "Gamma from BSM = " << gamma_BSM_european_put(initial_stock_price, volatility, expiration_time, strike_price1, risk_free_rate, no_of_trials) << endl;
	cout << "Gamma using LRM = " << gamma_LRM_european_put(initial_stock_price, volatility, expiration_time, strike_price1, risk_free_rate, no_of_trials) << endl;
	cout << "European Put Option with K = 100:" << endl;
	cout << "Price = " << option_price_monte_carlo_european_put(initial_stock_price, volatility, expiration_time, strike_price2, risk_free_rate, no_of_trials) << endl;
	cout << "Delta from BSM = " << delta_BSM_european_put(initial_stock_price, volatility, expiration_time, strike_price2, risk_free_rate, no_of_trials) << endl;
	cout << "Delta using Pathwise Estimators = " << delta_pathwise_european_put(initial_stock_price, volatility, expiration_time, strike_price2, risk_free_rate, no_of_trials) << endl;
	cout << "Delta using LRM = " << delta_LRM_european_put(initial_stock_price, volatility, expiration_time, strike_price2, risk_free_rate, no_of_trials) << endl;
	cout << "Gamma from BSM = " << gamma_BSM_european_put(initial_stock_price, volatility, expiration_time, strike_price2, risk_free_rate, no_of_trials) << endl;
	cout << "Gamma using LRM = " << gamma_LRM_european_put(initial_stock_price, volatility, expiration_time, strike_price2, risk_free_rate, no_of_trials) << endl;
	cout << "European Put Option with K = 110:" << endl;
	cout << "Price = " << option_price_monte_carlo_european_put(initial_stock_price, volatility, expiration_time, strike_price3, risk_free_rate, no_of_trials) << endl;
	cout << "Delta from BSM = " << delta_BSM_european_put(initial_stock_price, volatility, expiration_time, strike_price3, risk_free_rate, no_of_trials) << endl;
	cout << "Delta using Pathwise Estimators = " << delta_pathwise_european_put(initial_stock_price, volatility, expiration_time, strike_price3, risk_free_rate, no_of_trials) << endl;
	cout << "Delta using LRM = " << delta_LRM_european_put(initial_stock_price, volatility, expiration_time, strike_price3, risk_free_rate, no_of_trials) << endl;
	cout << "Gamma from BSM = " << gamma_BSM_european_put(initial_stock_price, volatility, expiration_time, strike_price3, risk_free_rate, no_of_trials) << endl;
	cout << "Gamma using LRM = " << gamma_LRM_european_put(initial_stock_price, volatility, expiration_time, strike_price3, risk_free_rate, no_of_trials) << endl;
	cout << "Asian Call Option Price with K = 90:" << endl;
	cout << "Price = " << option_price_monte_carlo_asian_call(initial_stock_price, volatility, expiration_time, strike_price1, risk_free_rate, no_of_trials) << endl;
	cout << "Delta using Pathwise Estimators  with 32 dates = " << delta_pathwise_asian_call(initial_stock_price, volatility, expiration_time, strike_price1, risk_free_rate, no_of_trials, m1) << endl;
	cout << "Delta using Pathwise Estimators  with 64 dates = " << delta_pathwise_asian_call(initial_stock_price, volatility, expiration_time, strike_price1, risk_free_rate, no_of_trials, m2) << endl;
	cout << "Delta using Pathwise Estimators  with 96 dates = " << delta_pathwise_asian_call(initial_stock_price, volatility, expiration_time, strike_price1, risk_free_rate, no_of_trials, m3) << endl;
	cout << "Delta using Pathwise Estimators  with 128 dates = " << delta_pathwise_asian_call(initial_stock_price, volatility, expiration_time, strike_price1, risk_free_rate, no_of_trials, m4) << endl;
	cout << "Delta using LRM with 32 dates = " << delta_LRM_asian_call(initial_stock_price, volatility, expiration_time, strike_price1, risk_free_rate, no_of_trials, m1) << endl;
	cout << "Delta using LRM with 64 dates = " << delta_LRM_asian_call(initial_stock_price, volatility, expiration_time, strike_price1, risk_free_rate, no_of_trials, m2) << endl;
	cout << "Delta using LRM with 96 dates = " << delta_LRM_asian_call(initial_stock_price, volatility, expiration_time, strike_price1, risk_free_rate, no_of_trials, m3) << endl;
	cout << "Delta using LRM with 128 dates = " << delta_LRM_asian_call(initial_stock_price, volatility, expiration_time, strike_price1, risk_free_rate, no_of_trials, m4) << endl;
	cout << "Asian Call Option Price with K = 100:" << endl;
	cout << "Price = " << option_price_monte_carlo_asian_call(initial_stock_price, volatility, expiration_time, strike_price2, risk_free_rate, no_of_trials) << endl;
	cout << "Delta using Pathwise Estimators  with 32 dates = " << delta_pathwise_asian_call(initial_stock_price, volatility, expiration_time, strike_price2, risk_free_rate, no_of_trials, m1) << endl;
	cout << "Delta using Pathwise Estimators  with 64 dates = " << delta_pathwise_asian_call(initial_stock_price, volatility, expiration_time, strike_price2, risk_free_rate, no_of_trials, m2) << endl;
	cout << "Delta using Pathwise Estimators  with 96 dates = " << delta_pathwise_asian_call(initial_stock_price, volatility, expiration_time, strike_price2, risk_free_rate, no_of_trials, m3) << endl;
	cout << "Delta using Pathwise Estimators  with 128 dates = " << delta_pathwise_asian_call(initial_stock_price, volatility, expiration_time, strike_price2, risk_free_rate, no_of_trials, m4) << endl;
	cout << "Delta using LRM with 32 dates = " << delta_LRM_asian_call(initial_stock_price, volatility, expiration_time, strike_price2, risk_free_rate, no_of_trials, m1) << endl;
	cout << "Delta using LRM with 64 dates = " << delta_LRM_asian_call(initial_stock_price, volatility, expiration_time, strike_price2, risk_free_rate, no_of_trials, m2) << endl;
	cout << "Delta using LRM with 96 dates = " << delta_LRM_asian_call(initial_stock_price, volatility, expiration_time, strike_price2, risk_free_rate, no_of_trials, m3) << endl;
	cout << "Delta using LRM with 128 dates = " << delta_LRM_asian_call(initial_stock_price, volatility, expiration_time, strike_price2, risk_free_rate, no_of_trials, m4) << endl;
	cout << "Asian Call Option Price with K = 110:" << endl;
	cout << "Price = " << option_price_monte_carlo_asian_call(initial_stock_price, volatility, expiration_time, strike_price3, risk_free_rate, no_of_trials) << endl;
	cout << "Delta using Pathwise Estimators  with 32 dates = " << delta_pathwise_asian_call(initial_stock_price, volatility, expiration_time, strike_price3, risk_free_rate, no_of_trials, m1) << endl;
	cout << "Delta using Pathwise Estimators  with 64 dates = " << delta_pathwise_asian_call(initial_stock_price, volatility, expiration_time, strike_price3, risk_free_rate, no_of_trials, m2) << endl;
	cout << "Delta using Pathwise Estimators  with 96 dates = " << delta_pathwise_asian_call(initial_stock_price, volatility, expiration_time, strike_price3, risk_free_rate, no_of_trials, m3) << endl;
	cout << "Delta using Pathwise Estimators  with 128 dates = " << delta_pathwise_asian_call(initial_stock_price, volatility, expiration_time, strike_price3, risk_free_rate, no_of_trials, m4) << endl;
	cout << "Delta using LRM with 32 dates = " << delta_LRM_asian_call(initial_stock_price, volatility, expiration_time, strike_price3, risk_free_rate, no_of_trials, m1) << endl;
	cout << "Delta using LRM with 64 dates = " << delta_LRM_asian_call(initial_stock_price, volatility, expiration_time, strike_price3, risk_free_rate, no_of_trials, m2) << endl;
	cout << "Delta using LRM with 96 dates = " << delta_LRM_asian_call(initial_stock_price, volatility, expiration_time, strike_price3, risk_free_rate, no_of_trials, m3) << endl;
	cout << "Delta using LRM with 128 dates = " << delta_LRM_asian_call(initial_stock_price, volatility, expiration_time, strike_price3, risk_free_rate, no_of_trials, m4) << endl;*/
	variance_delta_asian_call(initial_stock_price, volatility, expiration_time, strike_price3, risk_free_rate, no_of_trials, m1);
	variance_delta_asian_call(initial_stock_price, volatility, expiration_time, strike_price3, risk_free_rate, no_of_trials, m2);
	variance_delta_asian_call(initial_stock_price, volatility, expiration_time, strike_price3, risk_free_rate, no_of_trials, m3);
	variance_delta_asian_call(initial_stock_price, volatility, expiration_time, strike_price3, risk_free_rate, no_of_trials, m4);

	cout << "--------------------------------" << endl;
}



