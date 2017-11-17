
#include<iostream>
#include<vector>
#include <cmath>
#include "newmat.h"
#include "newmatap.h"
#include "newmatio.h"
using namespace std;

double initial_stock_price, strike_price, barrier_price, max_stock_price;
double dS, dt, sol, rate, volatility, time, price;

double max(double a, double b) {
	return (a > b) ? a : b;
}

//cdf of Gaussian distribution to compute values for Black-Scholes formula
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

//function to compute option price using Black-Scholes formula
float closed_form_down_and_out_european_put_option(double initial_stock_price, double strike_price, double rate, double time, double volatility, double barrier_price)
{
	double time_sqrt = sqrt(time);
	double d1 = (log(initial_stock_price / strike_price) + rate*time) / (volatility*time_sqrt) + 0.5*volatility*time_sqrt;
	double d2 = d1 - (volatility*time_sqrt);
	float vanilla_put = strike_price*exp(-rate*time)*N(-d2) - initial_stock_price*N(-d1);
	float lambda = (rate + ((volatility*volatility) / 2)) / (volatility*volatility);
	float temp = 2 * lambda - 2.0;
	float x1 = (log(initial_stock_price / barrier_price) / (volatility*sqrt(time))) + (lambda*volatility*sqrt(time));
	float y = (log(barrier_price*barrier_price / (initial_stock_price*strike_price)) / (volatility*sqrt(time))) + (lambda*volatility*sqrt(time));
	float y1 = (log(barrier_price / initial_stock_price) / (volatility*sqrt(time))) + (lambda*volatility*sqrt(time));
	return (vanilla_put - (-initial_stock_price*N(-x1) + strike_price*exp(-rate*time)*N(-x1 + volatility*sqrt(time)) +
		initial_stock_price*pow(barrier_price / initial_stock_price, 2 * lambda)*(N(y) - N(y1)) -
		strike_price*exp(-rate*time)*pow(barrier_price / initial_stock_price, temp)*(N(y - volatility*sqrt(time)) - N(y1 - volatility*sqrt(time)))));
}

//function to compute option price using Crank-Nicolson FDS
void getSolution(double initial_stock_price, double strike_price, double rate, double time, double volatility, double barrier_price, double max_stock_price, double dS, double dt) {
	int M = round((max_stock_price - barrier_price) / dS);
	dS = (max_stock_price - barrier_price) / M;
	int N = round(time / dt);
	dt = time / N;

	////mesh to store u(x,t) values
	Matrix solution(M+1, N+1);
	solution = 0;
	ColumnVector vetS(M + 1), veti(M+1);
	int index = 1;
	for (double i = barrier_price; index <= M+1; i+=dS)
	{
		vetS(index) = i;
		veti(index) = i / dS;
		index++;
	}

	//initializing mesh
	for (int i = 1; i <= M + 1; i++)
		solution(i, N + 1) = max(strike_price - vetS(i), 0);

	ColumnVector alpha(M + 1), beta(M + 1), gamma(M + 1);

	//calculating alpha, beta, gamma values as derived in the previous question
	for (int i = 1; i <= M+1; i++)
	{
		alpha(i) = (0.25*dt*(pow(volatility, 2)*pow(veti(i), 2) - rate*veti(i)));
		beta(i) = (-dt*0.5*(pow(volatility, 2)*pow(veti(i), 2) + rate));
		gamma(i) = (0.25*dt*(pow(volatility, 2)*pow(veti(i), 2) + rate*veti(i)));
	}

	Matrix M1(M - 1, M - 1);
	Matrix M2(M - 1, M - 1);

	M1 = 0;
	M2 = 0;

	//initializing the matrices as specified in the precious question
	for (int i = 1; i <= M - 1; i++)
	{
		M1(i, i) = 1 - beta(i + 1);
		M2(i, i) = 1 + beta(i + 1);
	}

	for (int i = 1; i <= M - 2; i++)
	{
		M1(i + 1, i) = -alpha(i + 2);
		M1(i, i + 1) = -gamma(i + 1);
		M2(i + 1, i) = alpha(i + 2);
		M2(i, i + 1) = gamma(i + 1);
	}

	//using submatrices of M1 and M2 to calculate mesh values
	for (int j = N; j >= 1; j--)
		solution.SubMatrix(2, M, j, j) = M1.i()*(M2*solution.SubMatrix(2, M, j + 1, j + 1));

	//to find option price
	for (int i = 1; i <= M + 1; i++)
	{
		if (vetS(i) == initial_stock_price)
		{
			price = solution(i, 1);
			break;
		}
		else if (vetS(i) > initial_stock_price)
		{
			price = (solution(i, 1) - solution(i - 1, 2)) / 2;
			break;
		}
	}
}

int main(int argc, char* argv[])
{
	sscanf(argv[1], "%lf", &initial_stock_price);
	sscanf(argv[2], "%lf", &strike_price);
	sscanf(argv[3], "%lf", &rate);
	sscanf(argv[4], "%lf", &time);
	sscanf(argv[5], "%lf", &volatility);
	sscanf(argv[6], "%lf", &barrier_price);
	sscanf(argv[7], "%lf", &max_stock_price);
	sscanf(argv[8], "%lf", &dS);
	sscanf(argv[9], "%lf", &dt);

	getSolution(initial_stock_price, strike_price, rate, time, volatility, barrier_price, max_stock_price, dS, dt);

	cout << "European Down-and-Out Barrier Option Pricing" << endl;
	cout << "Expiration time (Years) = " << time << endl;
	cout << "Risk Free Interest Rate = " << rate << endl;
	cout << "Volatility (%age of stock value) = " << volatility * 100 << endl;
	cout << "Initial Stock Price = " << initial_stock_price << endl;
	cout << "Maximum Stock Price = " << max_stock_price << endl;
	cout << "Strike Price = " << strike_price << endl;
	cout << "Barrier Price = " << barrier_price << endl;
	cout << "--------------------------------------" << endl;
	cout << "Price of an European Down and Out Put Option using Crank-Nicolson FDS method= " << price << endl;
	cout << "Price of an European Down and Out Put Option from Theory = " <<
		closed_form_down_and_out_european_put_option(initial_stock_price, strike_price, rate, time, volatility, barrier_price) << endl;

	return 0;
}
