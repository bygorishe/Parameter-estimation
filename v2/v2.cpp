#include <iostream>
#include <cmath>
#include <vector>

using namespace std;

# define pi 3.14159265358979323846

vector<double> operator * (vector<double> vector, const double& w) {
	size_t size = vector.size();
	for (size_t i = 0; i < size; ++i)
		vector[i] *= w;
	return vector;
}

class MNK {
private:
	double delta = 10., delta0 = 0.01, I = 1;
	vector<int> A = { 0, 0, 0 }, B = { 100, 0, 0 };
	vector<vector<int>>	M = {{ 200, 0, 0 }, { 500, 0, 0 }, { 1000, 0, 0 }},
		N = { { 300, 0, 0 }, { 600, 0, 0 }, { 1100, 0, 0 } };

	vector<double> V_sint, w, Vt, V;

	double r(vector<int> a, vector<int> b) {
		return sqrt((a[0] - b[0]) * (a[0] - b[0]) + (a[1] - b[1]) * (a[1] - b[1]) + (a[2] - b[2]) * (a[2] - b[2]));
	}

	vector<double> potential_V(double delta) { //обуз с  - delta * delta
		vector<double> temp(3);
		for (int i = 0; i < 3; i++)
			temp[i] = (I * (1 / r(B, M[i]) - 1 / r(A, M[i])) - (1 / r(B, N[i]) - 1 / r(A, N[i]))) / (2 * pi * delta);
		return temp;
	}

	double Func(vector<double> Vt) {
		double sum = 0.;
		for (int i = 0; i < 3; i++)
			sum += w[i] * (Vt[i] - V_sint[i]) * w[i] * (Vt[i] - V_sint[i]);
		return sum;
	}

public:
	MNK() {
		V_sint.resize(3); w.resize(3), Vt.resize(3), V.resize(3);

		cout << "sigma sint = " << delta << " sigma0 = " << delta0 << endl << endl;

		for (int k = 0; k < 2; k++) {
			if (k == 0) {
				V_sint = potential_V(delta);
				cout << "test 1, without noise" << endl;
			}
			else {
				V_sint = potential_V(delta) * 1.1;
				cout << "test 2, with noise" << endl;
			}

			double sigma = delta0;

			for (int i = 0; i < 3; i++)
				w[i] = 1 / V_sint[i];

			Vt = potential_V(sigma);

			//cout << "delta sint = " << delta << " delta0 = " << delta0 << endl;

			double F = 1;
			int iter = 0;

			while (F > 0.0000001) {
				double A11 = 0, b1 = 0, delta_d;
				V = potential_V(-sigma * sigma); //obyz func
				for (int i = 0; i < 3; i++) {
					A11 += (w[i] * V[i]) * (w[i] * V[i]);
					b1 -= w[i] * w[i] * V[i] * (Vt[i] - V_sint[i]);
				}

				delta_d = b1 / A11;

				sigma += delta_d;

				Vt = potential_V(sigma);

				F = Func(Vt);

				cout << "iter " << iter++ << " |delta " << delta_d << " |sigma " << sigma << " |f " << F << endl;
			}
			cout << endl << endl;
		}
	}
};

int main() {
	MNK M;
}