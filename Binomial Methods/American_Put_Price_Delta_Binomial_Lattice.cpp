
#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <vector>
using namespace std;

float up_factor, uptick_prob, risk_free_rate, strike_price;
float initial_stock_price, expiration_time, volatility, R;
vector <float> Svals, Pvals;
int no_of_steps;

float max(float a, float b) {
	return (b < a) ? a : b;
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

int main(int argc, char* argv[])
{
	sscanf(argv[1], "%f", &expiration_time);
	sscanf(argv[2], "%d", &no_of_steps);
	sscanf(argv[3], "%f", &risk_free_rate);
	sscanf(argv[4], "%f", &volatility);
	sscanf(argv[5], "%f", &initial_stock_price);
	sscanf(argv[6], "%f", &strike_price);

	double dt = expiration_time / ((float)no_of_steps);
	up_factor = exp(volatility*sqrt(dt));
	R = exp(risk_free_rate*dt);
	uptick_prob = (R - (1 / up_factor)) / (up_factor - (1 / up_factor));
	
	//to calculate delta
	double d = (log(initial_stock_price / strike_price) + (risk_free_rate + pow(volatility, 2)*expiration_time*0.5)) / (volatility*pow(expiration_time, 0.5));
	double delta = N(d) - 1;

	for (int i = 0; i <= 2 * no_of_steps + 1; i++)
	{
		Svals.push_back(0); //initializing the vector that will store stock prices
		Pvals.push_back(0); //initializing the vector that will store option values at various times
	}

	Svals[no_of_steps + 1] = initial_stock_price;

	for (int i = 1; i <= no_of_steps; i++)
	{
		Svals[no_of_steps + 1 + i] = up_factor*Svals[no_of_steps + i];
		Svals[no_of_steps + 1 - i] = (1/up_factor)*Svals[no_of_steps + 2 - i];
	}

	for (int i = 1; i <= 2 * no_of_steps + 1; i += 2)
		Pvals[i] = max(strike_price - Svals[i], 0);


	for (int t = 1; t <= no_of_steps; t++)
	{
		for (int i = t + 1; i <= 2 * no_of_steps + 1 - t; i += 2)
			Pvals[i] = max((uptick_prob*Pvals[i + 1] + ((1 - uptick_prob)*Pvals[i - 1])) / R, strike_price - Svals[i]); // //since an American option can be exercised at any time, we check what has the higher value at the time of exercising

	}

	cout << "American Put Option Pricing using a Binomial Lattice" << endl;
	cout << "Expiration Time (Years) = " << expiration_time << endl;
	cout << "Number of Steps = " << no_of_steps << endl;
	cout << "Risk Free Interest Rate = " << risk_free_rate << endl;
	cout << "Volatility (%age of stock value) = " << volatility * 100 << endl;
	cout << "Initial Stock Price = " << initial_stock_price << endl;
	cout << "Strike Price = " << strike_price << endl;
	cout << "--------------------------------------" << endl;
	cout << "Up Factor = " << up_factor << endl;
	cout << "Uptick Probability = " << uptick_prob << endl;
	cout << "--------------------------------------" << endl;
	cout << "Price of an American Put Option = " << Pvals[no_of_steps+1] << endl;
	cout << "Delta of an American Put Option = " << delta << endl;
	cout << "--------------------------------------" << endl;
}
