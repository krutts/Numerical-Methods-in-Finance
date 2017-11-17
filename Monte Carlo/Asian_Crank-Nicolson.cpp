
#include<iostream>
#include<vector>
#include <cmath>
#include "newmat.h"
#include "newmatap.h"
#include "newmatio.h"
#include <iomanip>
#include <chrono>
using namespace std;
using namespace chrono;

#define e 2.71828

double initial_stock_price, strike_price, barrier_price, max_stock_price;
double dS, dt, sol, rate, volatility, time_of_maturity, price, delta, gamma, theta;

double max(double a, double b) {
	return (a > b) ? a : b;
}


//function to compute average strike option price using Crank-Nicolson FDS
double getSolution(double initial_stock_price, double rate, double time_of_maturity, double volatility, double dS, double dt) {
	int M = round(initial_stock_price*2/dS);
	int N = round(time_of_maturity /dt);

	max_stock_price = M*dS;

	//mesh to store u(x,t) values
	Matrix solution(M + 1, N + 1);
	solution = 0;
	ColumnVector vetS(M + 1), veti(M + 1), vetj(N + 1), vetR(M + 1);
	int index = 1;

	vetR = 0;
	//vector to store x values
	for (double i = 0; i <= max_stock_price; i += dS)
		vetS(index++) = i;
	for (int i = 1; i <= M + 1; i++) {
		for (int j = 1; j <= i; j++)
			if(vetS(i)!=0)
			vetR(i) += vetS(j) / vetS(i);
	}
	vetR(1)=0.5;
	double dR = vetR(M+1) - vetR(M);
	index = 1;
	//vector to store i values
	for (int i = 0; i <= M; i += 1)
		veti(index++) = i;

	index = 1;


	//initializing mesh, setting boundary conditions for Asian average strike option
	for (int i = 1; i <= M + 1; i++)
		solution(i, N + 1) = initial_stock_price*max(1 - vetR(i)/ time_of_maturity, 0); //payoff


	ColumnVector D(M + 1), E(M + 1), F(M + 1), G(M + 1), K(M + 1), J(M + 1);

	//calculating alpha, beta, gamma values as derived in the previous question
	for (int i = 1; i <= M + 1; i++)
	{
		double d = 1 - rate*vetR(i);
		double c = 0.5*pow(volatility, 2)*pow(vetR(i), 2);
		D(i) = 0.5*c / pow(dR, 2) + 0.25*d / dR;
		E(i) = -c / pow(dR, 2) +1/dt;
		F(i) = 0.5*c / pow(dR, 2) - 0.25*d / dR;
		G(i) = -0.5*c / pow(dR, 2) - 0.25*d / dR;
		K(i) = c / pow(dR, 2) + 1 / dt;
		J(i) = -0.5*c / pow(dR, 2) + 0.25*d / dR;
	}

	Matrix M1(M - 1, M - 1);
	Matrix M2(M - 1, M - 1);

	M1 = 0;
	M2 = 0;

	//initializing the matrices as specified in the precious question
	for (int i = 1; i <= M - 1; i++)
	{
		M1(i, i) = E(i + 1);
		M2(i, i) = K(i + 1);
	}

	for (int i = 1; i <= M - 2; i++)
	{
		M1(i + 1, i) = F(i + 2);
		M1(i, i + 1) = D(i + 1);
		M2(i + 1, i) = J(i + 2);
		M2(i, i + 1) = G(i + 1);
	}

	ColumnVector aux(M - 1);
	aux = 0;

	//using submatrices of M1 and M2 to calculate mesh values
	for (int j = N; j >= 1; j--)
	{
		aux(1) = F(2)*solution(1, j + 1) - J(2)*solution(1, j); //to include the boundary conditions in the submatrix
		aux(M - 1) = D(M)*solution(M + 1, j + 1) - G(M)*solution(M + 1, j);
		solution.SubMatrix(2, M, j, j) = M2.i()*(M1*solution.SubMatrix(2, M, j + 1, j + 1) + aux);
	}

	
	//to find option price
	for (int i = 1; i <= M + 1; i++)
	{
		if (vetS(i) == initial_stock_price)
		{
			cout << "i=" << i << endl;
			price = initial_stock_price*solution(i, 1);
			delta = 0.5*initial_stock_price*(solution(i + 1, 1) - solution(i - 1, 1)) / dS;
			theta = 0.5*initial_stock_price*(solution(i, 3) - solution(i, 1)) / dt;
			break;
		}
		else if (vetS(i) > initial_stock_price)
		{
			price = initial_stock_price*(solution(i, 1) + solution(i + 1, 2)) / 2;
			delta = 0.5*initial_stock_price*(solution(i + 1, 1) - solution(i - 1, 1)) / dS;
			theta = 0.5*initial_stock_price*(solution(i, 3) - solution(i, 1)) / dt;
			break;
		}
	}
	return price;

}

int main(int argc, char* argv[])
{
	initial_stock_price = 50;
	rate = 0.06;
	time_of_maturity = 10;
	volatility = 0.05;
	dS = 2;
	dt = 0.1;

	cout << "Asian average strike Option Pricing using Crank-Nicolson FDS method " << endl;
	cout << "Expiration time (Years) = " << time_of_maturity << endl;
	cout << "Risk Free Interest Rate = " << rate << endl;
	cout << "Volatility (%age of stock value) = " << volatility * 100 << endl;
	cout << "Initial Stock Price = " << initial_stock_price << endl;
	cout << "--------------------------------------" << endl;
	auto repeatedTime1 = system_clock::now();
	double price0 = getSolution(initial_stock_price, rate, time_of_maturity, volatility, dS, dt); //p(S,t)
	auto repeatedTime2 = system_clock::now();
	cout << "Price of Asian Option = " << price0 << endl;
	cout << "It took " << duration_cast<chrono::nanoseconds>(repeatedTime2 - repeatedTime1).count() / pow(10, 9) << " seconds to complete." << endl; //calculating the difference in the CPU times calculated above
	cout << "--------------------------------------" << endl;
	cout << "Approximations of sensitivities: " << endl;
	cout << "Delta = " << delta << endl;

	initial_stock_price += dS;
	double price1 = getSolution(initial_stock_price, rate, time_of_maturity, volatility, dS, dt); //p(S+dS,t)
	initial_stock_price -= 2 * dS;
	double price2 = getSolution(initial_stock_price, rate, time_of_maturity, volatility, dS, dt); //p(S-dS,t)

	gamma = (price2 - 2 * price0 + price1) / pow(dS, 2);

	cout << "Gamma = " << gamma << endl;
	cout << "Theta = " << theta << endl;
	
	return 0;
}
