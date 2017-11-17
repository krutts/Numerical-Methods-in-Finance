
#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <random>
#include <chrono>
using namespace std;
using namespace chrono;
std::default_random_engine generator;

double r, K;
double S0, T, sigma, mu, delta_t;
int N;

double max(double a, double b) {
	return (b < a) ? a : b;
}

double get_uniform()
{
	std::uniform_real_distribution <double> distribution(0.0, 1.0);
	double number = distribution(generator);
	return (number);
}

double Arithmetic_Asian_SimPath(double S0, double T, double sigma, double delta_t, double mu) {
	int n = round(T / delta_t); // n needs to be even
	double *t = new double[n + 1];
	double *S = new double[n + 1];
	double sum = 0;
	t[0] = 0;
	S[0] = S0;
	for (int i = 1; i <= n / 2; i++) {
		t[2 * i - 1] = delta_t + t[2 * i - 2];
		t[2 * i] = delta_t + t[2 * i - 1];
		// create the unit normal variates using the Box-Muller Transform
		double x = get_uniform();
		double y = get_uniform();
		double a = sqrt(-2.0*log(x)) * cos(6.283185307999998*y);
		double b = sqrt(-2.0*log(x)) * sin(6.283185307999998*y);
		S[2 * i - 1] = S[2 * i - 2] * exp((mu - pow(sigma, 2) / 2)*delta_t + sigma*sqrt(delta_t)*a);
		S[2 * i] = S[2 * i - 1] * exp((mu - pow(sigma, 2) / 2)*delta_t + sigma*sqrt(delta_t)*b);
	}
	for (int i = 0; i <= n; i++)
		sum = sum + S[i] / (n + 1);
	return sum;
}

double Arithmetic_asian_call(double S0, double T, double sigma, double delta_t, double mu, double r, double K, int N) {
	double MC_call = 0;
	for (int i = 1; i <= N; i++) {
		MC_call = MC_call + exp(-r*T)*max(Arithmetic_Asian_SimPath(S0, T, sigma, delta_t, mu) - K, 0) / N;
	}
	return MC_call;
}

//Path-Dependent method to compute delta
double Delta_Arithmetic_asian_call(double S0, double T, double sigma, double delta_t, double mu, double r, double K, int N) {
	double MC_call = 0;
	for (int i = 1; i <= N; i++) {
		double alpha = Arithmetic_Asian_SimPath(S0, T, sigma, delta_t, mu);
		if (alpha - K > 0) {
			MC_call = MC_call + exp(-r*T)*alpha / S0 / N;
		}
	}
	return MC_call;
}

//Finite Difference method to compute Gamma
double Gamma_Arithmetic_asian_call(double S0, double T, double sigma, double delta_t, double mu, double r, double K, int N) {
	double MC_call = 0;
	double h = S0 / N;
	for (int i = 1; i <= N; i++) {
		MC_call = MC_call + (exp(-r*T)*max(Arithmetic_Asian_SimPath(S0 + h, T, sigma, delta_t, mu) - K, 0)
			- 2 * exp(-r*T)*max(Arithmetic_Asian_SimPath(S0, T, sigma, delta_t, mu) - K, 0)
			+ exp(-r*T)*max(Arithmetic_Asian_SimPath(S0 - h, T, sigma, delta_t, mu) - K, 0)) / (N * h * h);
	}
	return MC_call;
}

double Arithmetic_asian_var_call(double S0, double T, double sigma, double delta_t, double mu, double r, double K, int N) {
	double MC_call = 0;
	double *varas_call = new double[N + 1];
	varas_call[0] = 0;
	for (int i = 1; i <= N; i++) {
		varas_call[i] = exp(-r*T)*max(Arithmetic_Asian_SimPath(S0, T, sigma, delta_t, mu) - K, 0);
		MC_call = MC_call + varas_call[i] / N;
	}
	//MC_call=mean(call)
	for (int i = 1; i <= N; i++) {
		varas_call[0] = varas_call[0] + pow(MC_call - varas_call[i], 2) / (N - 1);
	}
	return varas_call[0] / N;
}

double Strike_Asian_SimPath(double S0, double T, double sigma, double delta_t, double mu) {
	int n = round(T / delta_t); // n needs to be even
	double *t = new double[n + 1];
	double *S = new double[n + 1];
	double sum = 0;
	t[0] = 0;
	S[0] = S0;
	for (int i = 1; i <= n / 2; i++) {
		t[2 * i - 1] = delta_t + t[2 * i - 2];
		t[2 * i] = delta_t + t[2 * i - 1];
		// create the unit normal variates using the Box-Muller Transform
		double x = get_uniform();
		double y = get_uniform();
		double a = sqrt(-2.0*log(x)) * cos(6.283185307999998*y);
		double b = sqrt(-2.0*log(x)) * sin(6.283185307999998*y);
		S[2 * i - 1] = S[2 * i - 2] * exp((mu - pow(sigma, 2) / 2)*delta_t + sigma*sqrt(delta_t)*a);
		S[2 * i] = S[2 * i - 1] * exp((mu - pow(sigma, 2) / 2)*delta_t + sigma*sqrt(delta_t)*b);
	}
	for (int i = 0; i <= n; i++)
		sum = sum + S[i] / (n + 1);
	return S[n] - sum;
}

double Strike_asian_call(double S0, double T, double sigma, double delta_t, double mu, double r, int N) {
	double MC_call = 0;
	for (int i = 1; i <= N; i++) {
		MC_call = MC_call + exp(-r*T)*max(Strike_Asian_SimPath(S0, T, sigma, delta_t, mu), 0) / N;
	}
	return MC_call;
}

//Path-Dependent method to compute delta
double Delta_Strike_asian_call(double S0, double T, double sigma, double delta_t, double mu, double r, int N) {
	double MC_call = 0;
	for (int i = 1; i <= N; i++) {
		double alpha = Strike_Asian_SimPath(S0, T, sigma, delta_t, mu);
		if (alpha > 0) {
			MC_call = MC_call + exp(-r*T)*alpha / S0 / N;
		}
	}
	return MC_call;
}

//Finite Difference method to compute Gamma
double Gamma_Strike_asian_call(double S0, double T, double sigma, double delta_t, double mu, double r, int N) {
	double MC_call = 0;
	double h = S0 / N;
	for (int i = 1; i <= N; i++) {
		MC_call = MC_call + (exp(-r*T)*max(Strike_Asian_SimPath(S0 + h, T, sigma, delta_t, mu), 0)
			- 2 * exp(-r*T)*max(Strike_Asian_SimPath(S0, T, sigma, delta_t, mu), 0)
			+ exp(-r*T)*max(Strike_Asian_SimPath(S0 - h, T, sigma, delta_t, mu), 0)) / (N * h * h);
	}
	return MC_call;
}

double Strike_asian_var_call(double S0, double T, double sigma, double delta_t, double mu, double r, int N) {
	double MC_call = 0;
	double *varas_call = new double[N + 1];
	varas_call[0] = 0;
	for (int i = 1; i <= N; i++) {
		varas_call[i] = exp(-r*T)*max(Strike_Asian_SimPath(S0, T, sigma, delta_t, mu), 0);
		MC_call = MC_call + varas_call[i] / N;
	}
	//MC_call=mean(call)
	for (int i = 1; i <= N; i++) {
		varas_call[0] = varas_call[0] + pow(MC_call - varas_call[i], 2) / (N - 1);
	}
	return varas_call[0] / N;
}

double Norm(const double& z) {
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

double Geometric_Closed_Form_Asian_Call(double S0, double T, double sigma, double delta_t, double r) {
	double n = round(T / delta_t);
	double adj_sigma = sigma*sqrt((2 * n + 1) / (6 * (n + 1)));
	double rho = 0.5*(r - 0.5*pow(sigma, 2) + pow(adj_sigma, 2));
	double d1 = (log(S0 / K) + (rho + 0.5*pow(adj_sigma, 2))*T) / (adj_sigma*sqrt(T));
	double d2 = (log(S0 / K) + (rho - 0.5*pow(adj_sigma, 2))*T) / (adj_sigma*sqrt(T));
	return exp(-r*T)*(S0*exp(rho*T)*Norm(d1) - K*Norm(d2));
}

double Geometric_Asian_SimPath(double S0, double T, double sigma, double delta_t, double mu) {
	int n = round(T / delta_t); // n needs to be even
	double *t = new double[n + 1];
	double *S = new double[n + 1];
	double prod = 1;
	t[0] = 0;
	S[0] = S0;
	for (int i = 1; i <= n / 2; i++) {
		t[2 * i - 1] = delta_t + t[2 * i - 2];
		t[2 * i] = delta_t + t[2 * i - 1];
		// create the unit normal variates using the Box-Muller Transform
		double x = get_uniform();
		double y = get_uniform();
		double a = sqrt(-2.0*log(x)) * cos(6.283185307999998*y);
		double b = sqrt(-2.0*log(x)) * sin(6.283185307999998*y);
		S[2 * i - 1] = S[2 * i - 2] * exp((mu - pow(sigma, 2) / 2)*delta_t + sigma*sqrt(delta_t)*a);
		S[2 * i] = S[2 * i - 1] * exp((mu - pow(sigma, 2) / 2)*delta_t + sigma*sqrt(delta_t)*b);
	}
	for (int i = 0; i <= n; i++)
		prod = prod*S[i];
	prod = log(prod);
	prod = prod / (n + 1);
	prod = exp(prod);
	return prod;
}

double Geometric_asian_call(double S0, double T, double sigma, double delta_t, double mu, double r, double K, int N) {
	double MC_call = 0;
	for (int i = 1; i <= N; i++) {
		MC_call = MC_call + exp(-r*T)*max(Geometric_Asian_SimPath(S0, T, sigma, delta_t, mu) - K, 0) / N;
	}
	return MC_call;
}

double Geometric_asian_var_call(double S0, double T, double sigma, double delta_t, double mu, double r, double K, int N) {
	double MC_call = 0;
	double *varas_call = new double[N + 1];
	varas_call[0] = 0;
	for (int i = 1; i <= N; i++) {
		varas_call[i] = exp(-r*T)*max(Geometric_Asian_SimPath(S0, T, sigma, delta_t, mu) - K, 0);
		MC_call = MC_call + varas_call[i] / N;
	}
	for (int i = 1; i <= N; i++) {
		varas_call[0] = varas_call[0] + pow(MC_call - varas_call[i], 2) / (N - 1);
	}
	return varas_call[0] / N;
}

double CV_Arithmetic_Asian_Call(double S0, double T, double sigma, double delta_t, double mu, double r, double K, int N) {
	int n = round(T / delta_t); // n needs to be even
	double *t = new double[n + 1];
	double *S = new double[n + 1];
	double *p = new double[N];
	double *q = new double[N];
	double expected_x = Geometric_Closed_Form_Asian_Call(S0, T, sigma, delta_t, r);
	double sum = 0;
	double MC_call = 0;
	t[0] = 0;
	S[0] = S0;
	double beta; //beta is what is called b in the lecture notes
	double var = 0; // variance of X, * (4*N*(4*N-1))
	double covar = 0; // covariance of X, Y, * (4*N*(4*N-1))
	double ave = Arithmetic_asian_call(S0, T, sigma, delta_t, mu, r, K, N);
	for (int j = 1; j <= N; j++) {
		double prod = 1;
		for (int i = 1; i <= n / 2; i++) {
			t[2 * i - 1] = delta_t + t[2 * i - 2];
			t[2 * i] = delta_t + t[2 * i - 1];
			// create the unit normal variates using the Box-Muller Transform
			double x = get_uniform();
			double y = get_uniform();
			double a = sqrt(-2.0*log(x)) * cos(6.283185307999998*y);
			double b = sqrt(-2.0*log(x)) * sin(6.283185307999998*y);
			S[2 * i - 1] = S[2 * i - 2] * exp((mu - pow(sigma, 2) / 2)*delta_t + sigma*sqrt(delta_t)*a);
			S[2 * i] = S[2 * i - 1] * exp((mu - pow(sigma, 2) / 2)*delta_t + sigma*sqrt(delta_t)*b);
		}
		for (int i = 0; i <= n; i++) {
			sum = sum + S[i] / (n + 1);
			prod = prod*S[i];
		}
		prod = log(prod);
		prod = prod / (n + 1);
		prod = exp(prod);
		p[j] = exp(-r*T)*max(sum - K, 0) / 4 / N;
		q[j] = exp(-r*T)*max(prod - K, 0);
		covar = covar + (q[j] - expected_x)*(p[j] - ave);
		var = var + pow(q[j] - expected_x, 2);
	}
	beta = covar / var;
	for (int j = 1; j <= N; j++) {
		MC_call = MC_call + (p[j] + beta*(q[j] - expected_x)) / N;
	}
	return MC_call;
}

double CV_Arithmetic_Asian_var_Call(double S0, double T, double sigma, double delta_t, double mu, double r, double K, int N) {
	int n = round(T / delta_t); // n needs to be even
	double *t = new double[n + 1];
	double *S = new double[n + 1];
	double *p = new double[N];
	double *q = new double[N];
	double expected_x = Geometric_Closed_Form_Asian_Call(S0, T, sigma, delta_t, r);
	double sum = 0;
	double MC_call = 0;
	t[0] = 0;
	S[0] = S0;
	double beta; //beta is what is called b in the lecture notes
	double var = 0; // variance of X, * (4*N*(4*N-1))
	double covar = 0; // covariance of X, Y, * (4*N*(4*N-1))
	double covar2 = 0;
	double var2 = 0;
	double cor;
	double se;
	double ave = Arithmetic_asian_call(S0, T, sigma, delta_t, mu, r, K, N);
	for (int j = 1; j <= N; j++) {
		double prod = 1;
		for (int i = 1; i <= n / 2; i++) {
			t[2 * i - 1] = delta_t + t[2 * i - 2];
			t[2 * i] = delta_t + t[2 * i - 1];
			// create the unit normal variates using the Box-Muller Transform
			double x = get_uniform();
			double y = get_uniform();
			double a = sqrt(-2.0*log(x)) * cos(6.283185307999998*y);
			double b = sqrt(-2.0*log(x)) * sin(6.283185307999998*y);
			S[2 * i - 1] = S[2 * i - 2] * exp((mu - pow(sigma, 2) / 2)*delta_t + sigma*sqrt(delta_t)*a);
			S[2 * i] = S[2 * i - 1] * exp((mu - pow(sigma, 2) / 2)*delta_t + sigma*sqrt(delta_t)*b);
		}
		for (int i = 0; i <= n; i++) {
			sum = sum + S[i] / (n + 1);
			prod = prod*S[i];
		}
		prod = log(prod);
		prod = prod / (n + 1);
		prod = exp(prod);
		p[j] = exp(-r*T)*max(sum - K, 0) / 4 / N;
		q[j] = exp(-r*T)*max(prod - K, 0);
		covar = covar + (q[j] - expected_x)*(p[j] - ave);
		var = var + pow(q[j] - expected_x, 2);
	}
	beta = covar / var;
	for (int j = 1; j <= N; j++) {
		MC_call = MC_call + (p[j] + beta*(q[j] - expected_x)) / N;
		covar2 = covar2 + (q[j] - expected_x)*(p[j] - beta * (q[j] - expected_x) - ave);
		var2 = var2 + pow(p[j] - beta * (q[j] - expected_x) - ave, 2) / (N*(N - 1));
	}
	cor = covar2 / sqrt(var * var2 / N / (N - 1));
	se = sqrt(var2 / 4 / N);
	return var2;
}

double CV_Arithmetic_Asian_cor_Call(double S0, double T, double sigma, double delta_t, double mu, double r, double K, int N) {
	int n = round(T / delta_t); // n needs to be even
	double *t = new double[n + 1];
	double *S = new double[n + 1];
	double *p = new double[N];
	double *q = new double[N];
	double expected_x = Geometric_Closed_Form_Asian_Call(S0, T, sigma, delta_t, r);
	double sum = 0;
	double MC_call = 0;
	t[0] = 0;
	S[0] = S0;
	double beta; //beta is what is called b in the lecture notes
	double var = 0; // variance of X, * (4*N*(4*N-1))
	double covar = 0; // covariance of X, Y, * (4*N*(4*N-1))
	double covar2 = 0;
	double var2 = 0;
	double cor;
	double se;
	double ave = Arithmetic_asian_call(S0, T, sigma, delta_t, mu, r, K, N);
	for (int j = 1; j <= N; j++) {
		double prod = 1;
		for (int i = 1; i <= n / 2; i++) {
			t[2 * i - 1] = delta_t + t[2 * i - 2];
			t[2 * i] = delta_t + t[2 * i - 1];
			// create the unit normal variates using the Box-Muller Transform
			double x = get_uniform();
			double y = get_uniform();
			double a = sqrt(-2.0*log(x)) * cos(6.283185307999998*y);
			double b = sqrt(-2.0*log(x)) * sin(6.283185307999998*y);
			S[2 * i - 1] = S[2 * i - 2] * exp((mu - pow(sigma, 2) / 2)*delta_t + sigma*sqrt(delta_t)*a);
			S[2 * i] = S[2 * i - 1] * exp((mu - pow(sigma, 2) / 2)*delta_t + sigma*sqrt(delta_t)*b);
		}
		for (int i = 0; i <= n; i++) {
			sum = sum + S[i] / (n + 1);
			prod = prod*S[i];
		}
		prod = log(prod);
		prod = prod / (n + 1);
		prod = exp(prod);
		p[j] = exp(-r*T)*max(sum - K, 0) / 4 / N;
		q[j] = exp(-r*T)*max(prod - K, 0);
		covar = covar + (q[j] - expected_x)*(p[j] - ave);
		var = var + pow(q[j] - expected_x, 2);
	}
	beta = covar / var;
	for (int j = 1; j <= N; j++) {
		MC_call = MC_call + (p[j] + beta*(q[j] - expected_x)) / N;
		covar2 = covar2 + (q[j] - expected_x)*(p[j] - beta * (q[j] - expected_x) - ave);
		var2 = var2 + pow(p[j] - beta * (q[j] - expected_x) - ave, 2) / (N*(N - 1));
	}
	cor = covar2 / sqrt(var * var2 / N / (N - 1));
	se = sqrt(var2 / 4 / N);
	return cor;
}

double CV_Arithmetic_Asian_se_Call(double S0, double T, double sigma, double delta_t, double mu, double r, double K, int N) {
	int n = round(T / delta_t); // n needs to be even
	double *t = new double[n + 1];
	double *S = new double[n + 1];
	double *p = new double[N];
	double *q = new double[N];
	double expected_x = Geometric_Closed_Form_Asian_Call(S0, T, sigma, delta_t, r);
	double sum = 0;
	double MC_call = 0;
	t[0] = 0;
	S[0] = S0;
	double beta; //beta is what is called b in the lecture notes
	double var = 0; // variance of X, * (4*N*(4*N-1))
	double covar = 0; // covariance of X, Y, * (4*N*(4*N-1))
	double covar2 = 0;
	double var2 = 0;
	double cor;
	double se;
	double ave = Arithmetic_asian_call(S0, T, sigma, delta_t, mu, r, K, N);
	for (int j = 1; j <= N; j++) {
		double prod = 1;
		for (int i = 1; i <= n / 2; i++) {
			t[2 * i - 1] = delta_t + t[2 * i - 2];
			t[2 * i] = delta_t + t[2 * i - 1];
			// create the unit normal variates using the Box-Muller Transform
			double x = get_uniform();
			double y = get_uniform();
			double a = sqrt(-2.0*log(x)) * cos(6.283185307999998*y);
			double b = sqrt(-2.0*log(x)) * sin(6.283185307999998*y);
			S[2 * i - 1] = S[2 * i - 2] * exp((mu - pow(sigma, 2) / 2)*delta_t + sigma*sqrt(delta_t)*a);
			S[2 * i] = S[2 * i - 1] * exp((mu - pow(sigma, 2) / 2)*delta_t + sigma*sqrt(delta_t)*b);
		}
		for (int i = 0; i <= n; i++) {
			sum = sum + S[i] / (n + 1);
			prod = prod*S[i];
		}
		prod = log(prod);
		prod = prod / (n + 1);
		prod = exp(prod);
		p[j] = exp(-r*T)*max(sum - K, 0) / 4 / N;
		q[j] = exp(-r*T)*max(prod - K, 0);
		covar = covar + (q[j] - expected_x)*(p[j] - ave);
		var = var + pow(q[j] - expected_x, 2);
	}
	beta = covar / var;
	for (int j = 1; j <= N; j++) {
		MC_call = MC_call + (p[j] + beta*(q[j] - expected_x)) / N;
		covar2 = covar2 + (q[j] - expected_x)*(p[j] - beta * (q[j] - expected_x) - ave);
		var2 = var2 + pow(p[j] - beta * (q[j] - expected_x) - ave, 2) / (N*(N - 1));
	}
	cor = covar2 / sqrt(var * var2 / N / (N - 1));
	se = sqrt(var2 / 4 / N);
	return se;
}

double CV_Strike_Asian_Call(double S0, double T, double sigma, double delta_t, double mu, double r, double K, int N) {
	int n = round(T / delta_t); // n needs to be even
	double *t = new double[n + 1];
	double *S = new double[n + 1];
	double *p = new double[N];
	double *q = new double[N];
	double expected_x = Geometric_Closed_Form_Asian_Call(S0, T, sigma, delta_t, r);
	double sum = 0;
	double MC_call = 0;
	t[0] = 0;
	S[0] = S0;
	double beta; //beta is what is called b in the lecture notes
	double var = 0; // variance of X, * (4*N*(4*N-1))
	double covar = 0; // covariance of X, Y, * (4*N*(4*N-1))
	double ave = Strike_asian_call(S0, T, sigma, delta_t, mu, r, N);
	for (int j = 1; j <= N; j++) {
		double prod = 1;
		for (int i = 1; i <= n / 2; i++) {
			t[2 * i - 1] = delta_t + t[2 * i - 2];
			t[2 * i] = delta_t + t[2 * i - 1];
			// create the unit normal variates using the Box-Muller Transform
			double x = get_uniform();
			double y = get_uniform();
			double a = sqrt(-2.0*log(x)) * cos(6.283185307999998*y);
			double b = sqrt(-2.0*log(x)) * sin(6.283185307999998*y);
			S[2 * i - 1] = S[2 * i - 2] * exp((mu - pow(sigma, 2) / 2)*delta_t + sigma*sqrt(delta_t)*a);
			S[2 * i] = S[2 * i - 1] * exp((mu - pow(sigma, 2) / 2)*delta_t + sigma*sqrt(delta_t)*b);
		}
		for (int i = 0; i <= n; i++) {
			sum = sum + S[i] / (n + 1);
			prod = prod*S[i];
		}
		prod = log(prod);
		prod = prod / (n + 1);
		prod = exp(prod);
		p[j] = exp(-r*T)*max(S[n] - sum, 0);
		q[j] = exp(-r*T)*max(prod - K, 0);
		covar = covar + (q[j] - expected_x)*(p[j] - ave);
		var = var + pow(q[j] - expected_x, 2);
	}
	beta = covar / var;
	cout << beta << endl;
	for (int j = 1; j <= N; j++) {
		MC_call = MC_call + (p[j] + beta*(q[j] - expected_x)) / N;
	}
	return MC_call;
}

double CV_Strike_Asian_var_Call(double S0, double T, double sigma, double delta_t, double mu, double r, double K, int N) {
	int n = round(T / delta_t); // n needs to be even
	double *t = new double[n + 1];
	double *S = new double[n + 1];
	double *p = new double[N];
	double *q = new double[N];
	double expected_x = Geometric_Closed_Form_Asian_Call(S0, T, sigma, delta_t, r);
	double sum = 0;
	double MC_call = 0;
	t[0] = 0;
	S[0] = S0;
	double beta; //beta is what is called b in the lecture notes
	double var = 0; // variance of X, * (4*N*(4*N-1))
	double covar = 0; // covariance of X, Y, * (4*N*(4*N-1))
	double covar2 = 0;
	double var2 = 0;
	double cor;
	double se;
	double ave = Strike_asian_call(S0, T, sigma, delta_t, mu, r, N);
	for (int j = 1; j <= N; j++) {
		double prod = 1;
		for (int i = 1; i <= n / 2; i++) {
			t[2 * i - 1] = delta_t + t[2 * i - 2];
			t[2 * i] = delta_t + t[2 * i - 1];
			// create the unit normal variates using the Box-Muller Transform
			double x = get_uniform();
			double y = get_uniform();
			double a = sqrt(-2.0*log(x)) * cos(6.283185307999998*y);
			double b = sqrt(-2.0*log(x)) * sin(6.283185307999998*y);
			S[2 * i - 1] = S[2 * i - 2] * exp((mu - pow(sigma, 2) / 2)*delta_t + sigma*sqrt(delta_t)*a);
			S[2 * i] = S[2 * i - 1] * exp((mu - pow(sigma, 2) / 2)*delta_t + sigma*sqrt(delta_t)*b);
		}
		for (int i = 0; i <= n; i++) {
			sum = sum + S[i] / (n + 1);
			prod = prod*S[i];
		}
		prod = log(prod);
		prod = prod / (n + 1);
		prod = exp(prod);
		p[j] = exp(-r*T)*max(S[n] - sum, 0);
		q[j] = exp(-r*T)*max(prod - K, 0);
		covar = covar + (q[j] - expected_x)*(p[j] - ave);
		var = var + pow(q[j] - expected_x, 2);
	}
	beta = covar / var;
	for (int j = 1; j <= N; j++) {
		MC_call = MC_call + (p[j] + beta*(q[j] - expected_x)) / N;
		covar2 = covar2 + (q[j] - expected_x)*(p[j] - beta * (q[j] - expected_x) - ave);
		var2 = var2 + pow(p[j] - beta * (q[j] - expected_x) - ave, 2) / (N*(N - 1));
	}
	cor = covar2 / sqrt(var * var2 / N / (N - 1));
	se = sqrt(var2 / 4 / N);
	return var2;
}

double CV_Strike_Asian_cor_Call(double S0, double T, double sigma, double delta_t, double mu, double r, double K, int N) {
	int n = round(T / delta_t); // n needs to be even
	double *t = new double[n + 1];
	double *S = new double[n + 1];
	double *p = new double[N];
	double *q = new double[N];
	double expected_x = Geometric_Closed_Form_Asian_Call(S0, T, sigma, delta_t, r);
	double sum = 0;
	double MC_call = 0;
	t[0] = 0;
	S[0] = S0;
	double beta; //beta is what is called b in the lecture notes
	double var = 0; // variance of X, * (4*N*(4*N-1))
	double covar = 0; // covariance of X, Y, * (4*N*(4*N-1))
	double covar2 = 0;
	double var2 = 0;
	double cor;
	double se;
	double ave = Strike_asian_call(S0, T, sigma, delta_t, mu, r, N);
	for (int j = 1; j <= N; j++) {
		double prod = 1;
		for (int i = 1; i <= n / 2; i++) {
			t[2 * i - 1] = delta_t + t[2 * i - 2];
			t[2 * i] = delta_t + t[2 * i - 1];
			// create the unit normal variates using the Box-Muller Transform
			double x = get_uniform();
			double y = get_uniform();
			double a = sqrt(-2.0*log(x)) * cos(6.283185307999998*y);
			double b = sqrt(-2.0*log(x)) * sin(6.283185307999998*y);
			S[2 * i - 1] = S[2 * i - 2] * exp((mu - pow(sigma, 2) / 2)*delta_t + sigma*sqrt(delta_t)*a);
			S[2 * i] = S[2 * i - 1] * exp((mu - pow(sigma, 2) / 2)*delta_t + sigma*sqrt(delta_t)*b);
		}
		for (int i = 0; i <= n; i++) {
			sum = sum + S[i] / (n + 1);
			prod = prod*S[i];
		}
		prod = log(prod);
		prod = prod / (n + 1);
		prod = exp(prod);
		p[j] = exp(-r*T)*max(S[n] - sum, 0);
		q[j] = exp(-r*T)*max(prod - K, 0);
		covar = covar + (q[j] - expected_x)*(p[j] - ave);
		var = var + pow(q[j] - expected_x, 2);
	}
	beta = covar / var;
	for (int j = 1; j <= N; j++) {
		MC_call = MC_call + (p[j] + beta*(q[j] - expected_x)) / N;
		covar2 = covar2 + (q[j] - expected_x)*(p[j] - beta * (q[j] - expected_x) - ave);
		var2 = var2 + pow(p[j] - beta * (q[j] - expected_x) - ave, 2) / (N*(N - 1));
	}
	cor = covar2 / sqrt(var * var2 / N / (N - 1));
	se = sqrt(var2 / 4 / N);
	return cor;
}

double CV_Strike_Asian_se_Call(double S0, double T, double sigma, double delta_t, double mu, double r, double K, int N) {
	int n = round(T / delta_t); // n needs to be even
	double *t = new double[n + 1];
	double *S = new double[n + 1];
	double *p = new double[N];
	double *q = new double[N];
	double expected_x = Geometric_Closed_Form_Asian_Call(S0, T, sigma, delta_t, r);
	double sum = 0;
	double MC_call = 0;
	t[0] = 0;
	S[0] = S0;
	double beta; //beta is what is called b in the lecture notes
	double var = 0; // variance of X, * (4*N*(4*N-1))
	double covar = 0; // covariance of X, Y, * (4*N*(4*N-1))
	double covar2 = 0;
	double var2 = 0;
	double cor;
	double se;
	double ave = Strike_asian_call(S0, T, sigma, delta_t, mu, r, N);
	for (int j = 1; j <= N; j++) {
		double prod = 1;
		for (int i = 1; i <= n / 2; i++) {
			t[2 * i - 1] = delta_t + t[2 * i - 2];
			t[2 * i] = delta_t + t[2 * i - 1];
			// create the unit normal variates using the Box-Muller Transform
			double x = get_uniform();
			double y = get_uniform();
			double a = sqrt(-2.0*log(x)) * cos(6.283185307999998*y);
			double b = sqrt(-2.0*log(x)) * sin(6.283185307999998*y);
			S[2 * i - 1] = S[2 * i - 2] * exp((mu - pow(sigma, 2) / 2)*delta_t + sigma*sqrt(delta_t)*a);
			S[2 * i] = S[2 * i - 1] * exp((mu - pow(sigma, 2) / 2)*delta_t + sigma*sqrt(delta_t)*b);
		}
		for (int i = 0; i <= n; i++) {
			sum = sum + S[i] / (n + 1);
			prod = prod*S[i];
		}
		prod = log(prod);
		prod = prod / (n + 1);
		prod = exp(prod);
		p[j] = exp(-r*T)*max(S[n] - sum, 0);
		q[j] = exp(-r*T)*max(prod - K, 0);
		covar = covar + (q[j] - expected_x)*(p[j] - ave);
		var = var + pow(q[j] - expected_x, 2);
	}
	beta = covar / var;
	for (int j = 1; j <= N; j++) {
		MC_call = MC_call + (p[j] + beta*(q[j] - expected_x)) / N;
		covar2 = covar2 + (q[j] - expected_x)*(p[j] - beta * (q[j] - expected_x) - ave);
		var2 = var2 + pow(p[j] - beta * (q[j] - expected_x) - ave, 2) / (N*(N - 1));
	}
	cor = covar2 / sqrt(var * var2 / N / (N - 1));
	se = sqrt(var2 / 4 / N);
	return se;
}

int main(int argc, char* argv[])
{
	T = 1;
	r = 0.03;
	sigma = 0.6;
	S0 = 99;
	K = 105;
	N = 200000;
	delta_t = 0.025;
	cout << N << endl;
	cout << "The Arithmetic Asian Call Price by Monte-Carlo = " << Arithmetic_asian_call(S0, T, sigma, delta_t, r, r, K, N) << endl;
	cout << "The Delta of the Arithmetic Asian Call Price by Monte-Carlo = " << Delta_Arithmetic_asian_call(S0, T, sigma, delta_t, r, r, K, N) << endl;
	cout << "The Gamma of the Arithmetic Asian Call Price by Monte-Carlo = " << Gamma_Arithmetic_asian_call(S0, T, sigma, delta_t, r, r, K, N) << endl;
	cout << "The Variance of the Arithmetic Asian Call Price by Monte-Carlo = " << Arithmetic_asian_var_call(S0, T, sigma, delta_t, r, r, K, N) << endl;
	cout << "The Strike Arithmetic Asian Call Price by Monte-Carlo = " << Strike_asian_call(S0, T, sigma, delta_t, r, r, N) << endl;
	cout << "The Delta of the Strike Arithmetic Asian Call Price by Monte-Carlo = " << Delta_Strike_asian_call(S0, T, sigma, delta_t, r, r, N) << endl;
	cout << "The Gamma of the Strike Arithmetic Asian Call Price by Monte-Carlo = " << Gamma_Strike_asian_call(S0, T, sigma, delta_t, r, r, N) << endl;
	cout << "The Variance of the Strike Arithmetic Asian Call Price by Monte-Carlo = " << Strike_asian_var_call(S0, T, sigma, delta_t, r, r, N) << endl;
	cout << "The Geometric Asian Call Price by Closed Form Solution = " << Geometric_Closed_Form_Asian_Call(S0, T, sigma, delta_t, r) << endl;
	cout << "The Geometric Asian Call Price by Monte-Carlo = " << Geometric_asian_call(S0, T, sigma, delta_t, r, r, K, N) << endl;
	cout << "The Variance of the Geometric Asian Call Price by Monte-Carlo = " << Geometric_asian_var_call(S0, T, sigma, delta_t, r, r, K, N) << endl;
	cout << "The Standard Error of the Geometric Asian Call Price by Monte-Carlo = " << sqrt(Geometric_asian_var_call(S0, T, sigma, delta_t, r, r, K, N) / N) << endl;
	cout << "The Arithmetic Asian Call Price by Monte-Carlo with Control Variates = " << CV_Arithmetic_Asian_Call(S0, T, sigma, delta_t, r, r, K, N) << endl;
	cout << "The Variance of the Arithmetic Asian Call Price by Monte-Carlo with Control Variates = " << CV_Arithmetic_Asian_var_Call(S0, T, sigma, delta_t, r, r, K, N) << endl;
	cout << "The Correlation of the Arithmetic Asian Call Price by Monte-Carlo with Control Variates = " << CV_Arithmetic_Asian_cor_Call(S0, T, sigma, delta_t, r, r, K, N) << endl;
	cout << "The Standard Error of the Arithmetic Asian Call Price by Monte-Carlo with Control Variates = " << CV_Arithmetic_Asian_se_Call(S0, T, sigma, delta_t, r, r, K, N) << endl;
	auto repeatedTime1 = system_clock::now();
	cout << "The Strike Arithmetic Asian Call Price by Monte-Carlo with Control Variates = " << CV_Strike_Asian_Call(S0, T, sigma, delta_t, r, r, K, N) << endl;
	auto repeatedTime2 = system_clock::now();
	cout << "It took " << duration_cast<chrono::nanoseconds>(repeatedTime2 - repeatedTime1).count() / pow(10, 9) << " seconds to complete." << endl; //calculating the difference in the CPU times calculated above
	cout << "The Variance of the Strike Arithmetic Asian Call Price by Monte-Carlo with Control Variates = " << CV_Strike_Asian_var_Call(S0, T, sigma, delta_t, r, r, K, N) << endl;
	cout << "The Correlation of the Strike Arithmetic Asian Call Price by Monte-Carlo with Control Variates = " << CV_Strike_Asian_cor_Call(S0, T, sigma, delta_t, r, r, K, N) << endl;
	cout << "The Standard Error of the Strike Arithmetic Asian Call Price by Monte-Carlo with Control Variates = " << CV_Strike_Asian_se_Call(S0, T, sigma, delta_t, r, r, K, N) << endl;

}