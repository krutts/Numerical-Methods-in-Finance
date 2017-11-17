
#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <vector>
#include <ctime>
#include <chrono>
using namespace std;
using namespace chrono;

float up_factor1, up_factor2, uptick_prob, risk_free_rate, strike_price, rho, div1 = 0.0, div2 = 0.0, nu1, nu2;
float initial_stock_price1, initial_stock_price2, expiration_time, volatility1, volatility2, R, p_uu, p_ud, p_du, p_dd;
vector <float> S1vals, S2vals, temp;
vector <vector<float>> Cvals;
int no_of_steps;

float max(float a, float b) {
	return (b < a) ? a : b;
}

float get_European_spread_option_price(float expiration_time, int no_of_steps, float risk_free_rate, float volatility1, float volatility2, float initial_stock_price1, float initial_stock_price2, float strike_price, float rho)
{

	up_factor1 = exp(volatility1*sqrt(expiration_time / ((float)no_of_steps)));
	up_factor2 = exp(volatility2*sqrt(expiration_time / ((float)no_of_steps)));
	nu1 = risk_free_rate - div1 - 0.5*pow(volatility1, 2);
	nu2 = risk_free_rate - div2 - 0.5*pow(volatility2, 2);
	R = exp(risk_free_rate*expiration_time / ((float)no_of_steps));

	p_uu = (1 / R)*0.25*(1 + sqrt(expiration_time / ((float)no_of_steps))*(nu1 / volatility1 + nu2 / volatility2) + rho);
	p_ud = (1 / R)*0.25*(1 + sqrt(expiration_time / ((float)no_of_steps))*(nu1 / volatility1 - nu2 / volatility2) - rho);
	p_du = (1 / R)*0.25*(1 + sqrt(expiration_time / ((float)no_of_steps))*(-nu1 / volatility1 + nu2 / volatility2) - rho);
	p_dd = (1 / R)*0.25*(1 + sqrt(expiration_time / ((float)no_of_steps))*(-nu1 / volatility1 - nu2 / volatility2) + rho);

	for (int i = 0; i <= 2 * no_of_steps + 1; i++)
	{
		S1vals.push_back(0);
		S2vals.push_back(0);
		temp.push_back(0);
	}

	S1vals[1] = initial_stock_price1 / pow(up_factor1, no_of_steps);
	S2vals[1] = initial_stock_price2 / pow(up_factor2, no_of_steps);

	for (int i = 2; i <= 2 * no_of_steps + 1; i++)
	{
		S1vals[i] = up_factor1*S1vals[i - 1];
		S2vals[i] = up_factor2*S2vals[i - 1];
	}

	for (int i = 0; i <= 2 * no_of_steps + 1; i++)
		Cvals.push_back(temp);


	for (int i = 1; i <= 2 * no_of_steps + 1; i += 2)
	{
		for (int j = 1; j <= 2 * no_of_steps + 1; j += 2)
			Cvals[i][j] = max(S1vals[i] - S2vals[j] - strike_price, 0);
	}

	for (int t = 1; t <= no_of_steps; t++)
	{
		for (int i = t + 1; i <= 2 * no_of_steps + 1 - t; i += 2)
		{
			for (int j = t + 1; j <= 2 * no_of_steps + 1 - t; j += 2)
				Cvals[i][j] = p_uu*Cvals[i + 1][j + 1] + p_ud*Cvals[i + 1][j - 1] + p_du*Cvals[i - 1][j + 1] + p_dd*Cvals[i - 1][j - 1];
		}

	}

	return Cvals[no_of_steps + 1][no_of_steps + 1];
}

int main(int argc, char* argv[])
{

	sscanf(argv[1], "%f", &expiration_time);
	sscanf(argv[2], "%d", &no_of_steps);
	sscanf(argv[3], "%f", &risk_free_rate);
	sscanf(argv[4], "%f", &volatility1);
	sscanf(argv[5], "%f", &volatility2);
	sscanf(argv[6], "%f", &initial_stock_price1);
	sscanf(argv[7], "%f", &initial_stock_price2);
	sscanf(argv[8], "%f", &strike_price);
	sscanf(argv[9], "%f", &rho);
	

	cout << "European Spread Option Pricing" << endl;
	cout << "Expiration Time (Years) = " << expiration_time << endl;
	cout << "Number of Steps = " << no_of_steps << endl;
	cout << "Risk Free Interest Rate = " << risk_free_rate << endl;
	cout << "Volatility (%age of stock value) 1 = " << volatility1 * 100 << endl;
	cout << "Volatility (%age of stock value) 2 = " << volatility2 * 100 << endl;
	cout << "Initial Stock Price 1 = " << initial_stock_price1 << endl;
	cout << "Initial Stock Price 2 = " << initial_stock_price2 << endl;
	cout << "Strike Price = " << strike_price << endl;
	cout << "Correlation co-efficient = " << rho << endl;
	cout << "--------------------------------------" << endl;
	auto repeatedTime1 = system_clock::now(); //stores the current CPU time, i.e. in this case, the time just before the start of the execution of the option pricing method
	cout << "Price of an European-type Spread Option = " << get_European_spread_option_price(expiration_time, no_of_steps, risk_free_rate, volatility1, volatility2, initial_stock_price1, initial_stock_price2, strike_price, rho) << endl;
	auto repeatedTime2 = system_clock::now(); //stores the current CPU time, i.e. in this case, the time just after the end of the execution of the option pricing method
	cout << "It took " << duration_cast<chrono::nanoseconds>(repeatedTime2 - repeatedTime1).count() / pow(10, 9) << " seconds to complete." << endl; //calculating the difference in the CPU times calculated above
	cout << "--------------------------------------" << endl;
}
