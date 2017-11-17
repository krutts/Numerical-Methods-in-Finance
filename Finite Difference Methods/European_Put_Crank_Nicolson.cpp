
#include<iostream>
#include<vector>
#include <cmath>
#include "newmat.h"
#include "newmatap.h"
#include "newmatio.h"
#include <iomanip>
using namespace std;

#define e 2.71828

double initial_stock_price, strike_price, barrier_price, max_stock_price;
double dS, dt, sol, rate, volatility, time_of_maturity, price, delta, gamma, theta;

double max(double a, double b) {
	return (a > b) ? a : b;
}


//function to compute option price using Crank-Nicolson FDS
double getSolution(double initial_stock_price, double strike_price, double rate, double time_of_maturity, double volatility, double max_stock_price, double dS, double dt) {
	int M = round(max_stock_price / dS);
	int N = round(time_of_maturity / dt);

	////mesh to store u(x,t) values
	Matrix solution(M + 1, N + 1);
	solution = 0;
	ColumnVector vetS(M + 1), veti(M + 1), vetj(N+1);
	int index = 1;
	//vector to store x values
	for (double i = 0; i <= max_stock_price; i += dS)
		vetS(index++) = i;
	index = 1;
	//vector to store i values
	for (int i = 0; i <= M; i += 1)
		veti(index++) = i;
	index = 1;
	//vector to store t values
	for (int i = 0; i <= N; i += 1)
		vetj(index++) = i;

	//initializing mesh, setting boundary conditions for European put
	for (int i = 1; i <= M + 1; i++)
		solution(i, N + 1) = max(strike_price - vetS(i), 0);

	for (int j = 1; j <= N + 1; j++)
		solution(1, j) = (strike_price)*pow(e, -rate*dt*(N-vetj(j)));


	ColumnVector alpha(M + 1), beta(M + 1), gamma(M + 1);

	//calculating alpha, beta, gamma values as derived in the previous question
	for (int i = 1; i <= M + 1; i++)
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

	ColumnVector aux(M - 1);
	aux = 0;
	//using submatrices of M1 and M2 to calculate mesh values
	for (int j = N; j >= 1; j--)
	{
		aux(1) = alpha(2)*(solution(1, j)+solution(1,j+1)); //to include the boundary conditions in the submatrix
		aux(M - 1) = gamma(M + 1)*(solution(M+1, j) + solution(M + 1, j + 1));
		solution.SubMatrix(2, M, j, j) = M1.i()*(M2*solution.SubMatrix(2, M, j + 1, j + 1) + aux);
	}

	
	//to find option price
	for (int i = 1; i <= M + 1; i++)
	{
		if (vetS(i) == initial_stock_price)
		{
			price = solution(i, 1);
			delta = 0.5*(solution(i + 1, 1) - solution(i - 1, 1))/dS;
			theta = 0.5*(solution(i, 3) - solution(i, 1))/dt;
			break;
		}
		else if (vetS(i) > initial_stock_price)
		{
			price = (solution(i, 1) + solution(i + 1, 2)) / 2;
			delta = 0.5*(solution(i + 1, 1) - solution(i - 1, 1)) / dS;
			theta = 0.5*(solution(i, 3) - solution(i, 1)) / dt;
			break;
		}
	}
	return price;
	
}

int main(int argc, char* argv[])
{
	sscanf(argv[1], "%lf", &initial_stock_price);
	sscanf(argv[2], "%lf", &strike_price);
	sscanf(argv[3], "%lf", &rate);
	sscanf(argv[4], "%lf", &time_of_maturity);
	sscanf(argv[5], "%lf", &volatility);
	sscanf(argv[6], "%lf", &max_stock_price);
	sscanf(argv[7], "%lf", &dS);
	sscanf(argv[8], "%lf", &dt);

	cout << "European Put Option Pricing using Crank-Nicolson FDS method " << endl;
	cout << "Expiration time (Years) = " << time_of_maturity << endl;
	cout << "Risk Free Interest Rate = " << rate << endl;
	cout << "Volatility (%age of stock value) = " << volatility * 100 << endl;
	cout << "Initial Stock Price = " << initial_stock_price << endl;
	cout << "Maximum Stock Price = " << max_stock_price << endl;
	cout << "Strike Price = " << strike_price << endl;
	cout << "--------------------------------------" << endl;

	double price0 = getSolution(initial_stock_price, strike_price, rate, time_of_maturity, volatility, max_stock_price, dS, dt); //p(S,t)

	cout << "Price of European Put Option = " << price0 << endl;
	cout << "--------------------------------------" << endl;
	cout << "Approximations of sensitivities: " << endl;
	cout << "Delta = " << delta << endl;

	initial_stock_price += dS;
	double price1 = getSolution(initial_stock_price, strike_price, rate, time_of_maturity, volatility, max_stock_price, dS, dt); //p(S+dS,t)
	initial_stock_price -= 2*dS;
	double price2 = getSolution(initial_stock_price, strike_price, rate, time_of_maturity, volatility, max_stock_price, dS, dt); //p(S-dS,t)

	gamma = (price2 - 2 * price0 + price1) / pow(dS, 2);

	cout << "Gamma = " << gamma << endl;
	cout << "Theta = " << theta << endl;

	return 0;
}
