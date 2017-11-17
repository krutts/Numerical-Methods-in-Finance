
#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <vector>
using namespace std;

float up_factor1, up_factor2, up_factor3, uptick_prob, risk_free_rate, strike_price, correlation12, correlation23, correlation13, nu1, nu2, nu3;
float initial_stock_price1, initial_stock_price2, initial_stock_price3, expiration_time, volatility1, volatility2, volatility3, R, p_uu, p_ud, p_du, p_dd, p_uuu, p_uud, p_udu, p_udd, p_duu, p_dud, p_ddu, p_ddd;
vector <float> S1vals, S2vals, S3vals, temp;
vector <vector<float>> C1vals, temp2;
vector <vector <vector<float>>> C2vals;
int no_of_steps;

float max(float a, float b, float c) {
return ((a > b) ? ((a > c) ? a : c) : ((b > c) ? b : c));
}

int main(int argc, char* argv[])
{

sscanf(argv[1], "%f", &expiration_time);
sscanf(argv[2], "%d", &no_of_steps);
sscanf(argv[3], "%f", &risk_free_rate);
sscanf(argv[4], "%f", &volatility1);
sscanf(argv[5], "%f", &volatility2);
sscanf(argv[6], "%f", &volatility3); //input volatility3 as 0 if no. of assets = 2
sscanf(argv[7], "%f", &initial_stock_price1);
sscanf(argv[8], "%f", &initial_stock_price2);
sscanf(argv[9], "%f", &initial_stock_price3);
sscanf(argv[10], "%f", &strike_price);
sscanf(argv[11], "%f", &correlation12); //only this value is required if the no. of assets = 2
sscanf(argv[12], "%f", &correlation23);
sscanf(argv[13], "%f", &correlation13);

double dt = expiration_time / ((float)no_of_steps);
up_factor1 = exp(volatility1*sqrt(dt));
up_factor2 = exp(volatility2*sqrt(dt));
up_factor3 = exp(volatility3*sqrt(dt));
nu1 = risk_free_rate - 0.5*pow(volatility1, 2);
nu2 = risk_free_rate - 0.5*pow(volatility2, 2);
nu3 = risk_free_rate - 0.5*pow(volatility3, 2);
R = exp(risk_free_rate*dt);

if (volatility3 == 0) //i.e. no. of assets = 2
{
	p_uu = (1 / R)*0.25*(1 + sqrt(dt)*(nu1 / volatility1 + nu2 / volatility2) + correlation12);
	p_ud = (1 / R)*0.25*(1 + sqrt(dt)*(nu1 / volatility1 - nu2 / volatility2) - correlation12);
	p_du = (1 / R)*0.25*(1 + sqrt(dt)*(-nu1 / volatility1 + nu2 / volatility2) - correlation12);
	p_dd = (1 / R)*0.25*(1 + sqrt(dt)*(-nu1 / volatility1 - nu2 / volatility2) + correlation12);
}

else //risk-neutral probabilities for no. of assets = 3
{
	p_uuu = (1 / R)*0.125*(1 + correlation12 + correlation13 + correlation23 + sqrt(dt)*(nu1 / volatility1 + nu2 / volatility2 + nu3 / volatility3));
	p_uud = (1 / R)*0.125*(1 + correlation12 - correlation13 - correlation23 + sqrt(dt)*(nu1 / volatility1 + nu2 / volatility2 - nu3 / volatility3));
	p_udu = (1 / R)*0.125*(1 - correlation12 + correlation13 - correlation23 + sqrt(dt)*(nu1 / volatility1 - nu2 / volatility2 + nu3 / volatility3));
	p_udd = (1 / R)*0.125*(1 - correlation12 - correlation13 + correlation23 + sqrt(dt)*(nu1 / volatility1 - nu2 / volatility2 - nu3 / volatility3));
	p_duu = (1 / R)*0.125*(1 - correlation12 - correlation13 + correlation23 - sqrt(dt)*(nu1 / volatility1 - nu2 / volatility2 - nu3 / volatility3));
	p_dud = (1 / R)*0.125*(1 - correlation12 + correlation13 - correlation23 - sqrt(dt)*(nu1 / volatility1 - nu2 / volatility2 + nu3 / volatility3));
	p_ddu = (1 / R)*0.125*(1 + correlation12 - correlation13 - correlation23 - sqrt(dt)*(nu1 / volatility1 + nu2 / volatility2 - nu3 / volatility3));
	p_ddd = (1 / R)*0.125*(1 + correlation12 + correlation13 + correlation23 - sqrt(dt)*(nu1 / volatility1 + nu2 / volatility2 + nu3 / volatility3));

}

for (int i = 0; i <= 2 * no_of_steps + 1; i++)
{
S1vals.push_back(0); //initializing vector for stock values of stock 1
S2vals.push_back(0); //initializing vector for stock values of stock 2
S3vals.push_back(0); //initializing vector for stock values of stock 3
temp.push_back(0); //initializing a temporary vector to initialize the vector-of-vectors further ahead
}

for (int i = 0; i <= 2 * no_of_steps + 1; i++)
temp2.push_back(temp); //initializing a temporary vector-of-vectors to initialize the vector-of-vector-of-vectors further ahead


S1vals[1] = initial_stock_price1 / pow(up_factor1, no_of_steps);
S2vals[1] = initial_stock_price2 / pow(up_factor2, no_of_steps);
S3vals[1] = initial_stock_price3 / pow(up_factor3, no_of_steps);

for (int i = 2; i <= 2 * no_of_steps + 1; i++)
{
S1vals[i] = up_factor1*S1vals[i - 1];
S2vals[i] = up_factor2*S2vals[i - 1];
S3vals[i] = up_factor3*S3vals[i - 1];
}

for (int i = 0; i <= 2 * no_of_steps + 1; i++)
{
C1vals.push_back(temp); //initializing the vector-of-vectors
C2vals.push_back(temp2); //initializing the vetcor-of-vector-of-vectors
}

if (volatility3 == 0) //execute this only if no. of assets = 2
{
	for (int i = 1; i <= 2 * no_of_steps + 1; i += 2)
	{
		for (int j = 1; j <= 2 * no_of_steps + 1; j += 2)
		{
			C1vals[i][j] = max(max(S1vals[i], S2vals[j], 0) - strike_price, 0, 0);
		}
	}

	for (int t = 1; t <= no_of_steps; t++)
	{
		for (int i = t + 1; i <= 2 * no_of_steps + 1 - t; i += 2)
		{
			for (int j = t + 1; j <= 2 * no_of_steps + 1 - t; j += 2)
			{
				double hold = p_uu*C1vals[i + 1][j + 1] + p_ud*C1vals[i + 1][j - 1] + p_du*C1vals[i - 1][j + 1] + p_dd*C1vals[i - 1][j - 1];
				C1vals[i][j] = max(hold, max(S1vals[i], S2vals[j], 0) - strike_price, 0); //since an American option can be exercised at any time, we check what has the higher value at the time of exercising
			}
		}

	}
}
else //if no. of assets = 3
{
	for (int i = 1; i <= 2 * no_of_steps + 1; i += 2)
	{
		for (int j = 1; j <= 2 * no_of_steps + 1; j += 2)
		{
			for (int k = 1; k <= 2 * no_of_steps + 1; k += 2)
				C2vals[i][j][k] = max(max(S1vals[i], S2vals[j], S3vals[k]) - strike_price, 0, 0);
		}
	}

	for (int t = 1; t <= no_of_steps; t++)
	{
		for (int i = t + 1; i <= 2 * no_of_steps + 1 - t; i += 2)
		{
			for (int j = t + 1; j <= 2 * no_of_steps + 1 - t; j += 2)
			{
				for (int k = t + 1; k <= 2 * no_of_steps + 1 - t; k += 2)
				{
					double hold = p_uuu*C2vals[i + 1][j + 1][k + 1] + p_uud*C2vals[i + 1][j + 1][k - 1] + p_udu*C2vals[i + 1][j - 1][k + 1] + p_udd*C2vals[i + 1][j - 1][k - 1] + p_duu*C2vals[i - 1][j + 1][k + 1] + p_dud*C2vals[i - 1][j + 1][k - 1] + p_ddu*C2vals[i - 1][j - 1][k + 1] + p_ddd*C2vals[i - 1][j - 1][k - 1];
					C2vals[i][j][k] = max(hold, max(S1vals[i], S2vals[j], S3vals[k]) - strike_price, 0); //since an American option can be exercised at any time, we check what has the higher value at the time of exercising
				}
			}
		}

	}
}
cout << "American Max Option Pricing" << endl;
cout << "Expiration Time (Years) = " << expiration_time << endl;
cout << "Number of Steps = " << no_of_steps << endl;
cout << "Risk Free Interest Rate = " << risk_free_rate << endl;
volatility3 == 0 ? cout << "Number of assets = 2" << endl : cout << "Number of assets = 3" << endl;
cout << "Volatility (%age of stock value) 1 = " << volatility1 * 100 << endl;
cout << "Volatility (%age of stock value) 2 = " << volatility2 * 100;
volatility3 == 0 ? cout << endl : cout << endl << "Volatility (%age of stock value) 3 = " << volatility3 * 100 << endl;
cout << "Initial Stock Price 1 = " << initial_stock_price1 << endl;
cout << "Initial Stock Price 2 = " << initial_stock_price2;
volatility3 == 0 ? cout << endl : cout << endl << "Initial Stock Price 3 = " << initial_stock_price3 << endl;
cout << "Strike Price = " << strike_price << endl;
volatility3 == 0? cout << "Probabilities = " << p_uu << ", " << p_ud << ", " << p_du << ", " << p_dd << endl : cout << endl << "Probabilities = " << p_uuu << ", " << p_uud << ", " << p_udu << ", " << p_udd << ", " << p_duu << ", " << p_dud << ", " << p_ddu << ", " << p_ddd << endl;
volatility3 == 0 ? cout << "Correlation co-efficient = " << correlation12 << endl : cout << endl << "Correlation co-efficients = " << correlation12 << ", " << correlation13 << ", " << correlation23 << endl;
cout << "--------------------------------------" << endl;
volatility3 == 0 ? cout << "Price of an American Max Option = " << C1vals[no_of_steps + 1][no_of_steps + 1] << endl: cout << endl << "Price of an American Max Option = " << C2vals[no_of_steps + 1][no_of_steps + 1][no_of_steps + 1] << endl;
cout << "--------------------------------------" << endl;
}

