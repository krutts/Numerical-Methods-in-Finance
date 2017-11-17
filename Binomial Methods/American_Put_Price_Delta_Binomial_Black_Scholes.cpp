
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <vector>
using namespace std;


double up_factor, risk_free_rate, strike_price, uptick_prob;
double initial_stock_price, expiration_time, volatility, R;
double **memoized_array, bsm, richardson;
int no_of_steps;


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

double option_price_put_black_scholes(const double& S,      // spot price
	const double& K,      // Strike (exercise) price,
	const double& r,      // interest rate
	const double& sigma,  // volatility
	const double& time) {
	double time_sqrt = sqrt(time);
	double d1 = (log(S / K) + r*time) / (sigma*time_sqrt) + 0.5*sigma*time_sqrt;
	double d2 = d1 - (sigma*time_sqrt);
	return K*exp(-r*time)*N(-d2) - S*N(-d1);
};



double max(double a, double b) {
	return (b < a) ? a : b;
}


void reinitialize(double **memoized_array) //clearing/iniitializing the array
{
	for(int i = 0; i <= no_of_steps; i++)
		for (int j = 0; j <= 2 * no_of_steps + 1; j++)  
			memoized_array[i][j] = -1;
}

double get_American_option_price_BBS(int k, int i, double current_stock_price, int no_of_steps) {
	if (memoized_array[k][i] != -1)
		return memoized_array[k][i];
	else if (k == no_of_steps - 1)
		memoized_array[k][i] = max((strike_price - current_stock_price), (option_price_put_black_scholes(current_stock_price, strike_price, risk_free_rate, volatility, expiration_time / ((double)no_of_steps))));
	else
		memoized_array[k][i] = max((strike_price - current_stock_price), (uptick_prob*get_American_option_price_BBS(k + 1, i + 1, current_stock_price*up_factor, no_of_steps) + (1 - uptick_prob)*get_American_option_price_BBS(k + 1, i - 1, current_stock_price / up_factor, no_of_steps))/R);
	return memoized_array[k][i];
}


int main(int argc, char* argv[])
{

	sscanf(argv[1], "%lf", &expiration_time);
	sscanf(argv[2], "%d", &no_of_steps);
	sscanf(argv[3], "%lf", &risk_free_rate);
	sscanf(argv[4], "%lf", &volatility);
	sscanf(argv[5], "%lf", &initial_stock_price);
	sscanf(argv[6], "%lf", &strike_price);

	memoized_array = new double *[no_of_steps + 1]; //creating a multi-dimensional array that will save values in order to avoid unnecessary repetitive calculations and save computational time
	for (int i = 0; i <= no_of_steps; i++) //
		memoized_array[i] = new double[2 * no_of_steps + 1];

	cout << "Binomial American Option Pricing" << endl;
	cout << "Expiration Time (Years) = " << expiration_time << endl;
	cout << "Number of Steps = " << no_of_steps << endl;
	cout << "Risk Free Interest Rate = " << risk_free_rate << endl;
	cout << "Volatility (%age of stock value) = " << volatility * 100 << endl;
	cout << "Initial Stock Price = " << initial_stock_price << endl;
	cout << "Strike Price = " << strike_price << endl;
	cout << "--------------------------------------" << endl;

	cout << "Converging values from the Binomial Black Scholes Method :" << endl;
	for (int i = 2; i <= no_of_steps; i+=2)
	{
		reinitialize(memoized_array);
		double dt = expiration_time / ((float)i);
		up_factor = exp(volatility*sqrt(dt));
		R = exp(risk_free_rate*dt);
		uptick_prob = (R - (1 / up_factor)) / (up_factor - (1 / up_factor));
		cout << "Price from BSM = " << get_American_option_price_BBS(0, no_of_steps, initial_stock_price, i) << endl;
		bsm = get_American_option_price_BBS(0, no_of_steps, initial_stock_price, no_of_steps - 1);
	}

	cout << "----------------------------------------------------------------------------------------" << endl;
	cout << "Converging values from the Richardson extrapolation Method :" << endl;
	for (int i = 2; i <= no_of_steps; i += 2)
	{
		reinitialize(memoized_array);
		double dt = expiration_time / ((float)i);
		up_factor = exp(volatility*sqrt(dt));
		R = exp(risk_free_rate*dt);
		uptick_prob = (R - (1 / up_factor)) / (up_factor - (1 / up_factor));
		cout << "Price from Richardson Extrapolation = " << 2 * get_American_option_price_BBS(0, i, initial_stock_price, i) - get_American_option_price_BBS(0, i, initial_stock_price, i / 2) << endl;
		richardson = 2 * get_American_option_price_BBS(0, i, initial_stock_price, (no_of_steps - 1)) - get_American_option_price_BBS(0, i, initial_stock_price, (no_of_steps - 1) / 2);
	}

	
	cout << "Price of an American Put Option from BSM = " << bsm << endl;
	cout << "Price of an American Put Option from Richardson Extrapolation = " << richardson << endl;

}
