
#include "stdafx.h"
#include<iostream>
#include<vector>
#include <fstream>
using namespace std;
int xmin, xmax, tmax, c;
double dx, dt, sol;

void getSolution(int xmin, double dx, int xmax, double dt, int tmax, int c) {
	int N = ceil((xmax - xmin) / dx);
	int M = ceil(tmax / dt);
	double rho = dt / dx;

	//mesh to store u(x,t) values
	double** solution;
	solution = new double*[N+1];

	//initializing mesh
	for (int i = 0; i <= N; i++)
		solution[i] = new double[M+1];
	for (int i = 0; i <= N; i++)
		for (int j = 0; j <= M; j++)
			solution[i][j] = 0;

	//inclusing boundary conditions
	vector<double> vetx;
	for (double i = xmin; i <= xmax; i += dx)
		vetx.push_back(i<=0?((i<-1)?0:i+1):1);

	for (int i = 0; i <= N; i++)
		solution[i][0] = vetx[i];

	//calculating values using the FDS
	for (int j = 1; j <= M; j++)
		for (int i = 1; i <= N; i++)
			solution[i][j] = (1 - c*rho)*solution[i][j-1] + (c*rho)*solution[i - 1][j-1];
	
	ofstream of0("out0"), of1("out20"), of2("out40"), of3("out60"), of4("out80");
	
	//writing mesh values at t=0,20,40,60,80
	for (int i = 0; i <= N; i++)
		for (int j = 0; j <= M; j++)
		{
				of0 << solution[i][0] << endl;
				of1 << solution[i][20] << endl;
				of2 << solution[i][40] << endl;
				of3 << solution[i][60] << endl;
				of4 << solution[i][80] << endl;
		}
}

int main(int argc, char* argv[])
{
	sscanf(argv[1], "%d", &xmin);
	sscanf(argv[2], "%d", &xmax);
	sscanf(argv[3], "%lf", &dx);
	sscanf(argv[4], "%d", &tmax);
	sscanf(argv[5], "%lf", &dt);
	sscanf(argv[6], "%d", &c);

	getSolution(xmin, dx, xmax, dt, tmax, c);

	cout << "Transport equation FDS" << endl;
	cout << "Writing output to files..." << endl;

    return 0;
}

