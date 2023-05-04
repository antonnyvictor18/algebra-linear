#include <iostream>
#include <cmath>
#include <vector>

using namespace std;

const double EPS = 1e-6; // valor de tolerância para o critério de parada

// função que realiza o produto matriz-vetor
vector<double> prodMatVec(vector<vector<double>> A, vector<double> x) {
    int n = A.size();
    vector<double> y(n);
    for (int i = 0; i < n; i++) {
        y[i] = 0;
        for (int j = 0; j < n; j++) {
            y[i] += A[i][j] * x[j];
        }
    }
    return y;
}

// função que calcula a norma euclidiana de um vetor
double norma(vector<double> v) {
    double s = 0;
    for (double x : v) {
        s += x * x;
    }
    return sqrt(s);
}

// função que divide um vetor por um escalar
vector<double> divVec(double k, vector<double> v) {
    int n = v.size();
    vector<double> w(n);
    for (int i = 0; i < n; i++) {
        w[i] = v[i] / k;
    }
    return w;
}

vector<vector<double>> calcMatRot(vector<vector<double>> A, int p, int q) {
    int n = A.size();
    vector<vector<double>> J(n, vector<double>(n, 0));
    double theta = atan2(A[q][p], A[p][p]);
    double c = cos(theta);
    double s = sin(theta);
    for (int i = 0; i < n; i++) {
        J[i][i] = 1;
    }
    J[p][p] = c;
    J[q][q] = c;
    J[p][q] = -s;
    J[q][p] = s;
    return J;
}



int main() {
    int n;
    cout << "Digite a ordem da matriz A: ";
    cin >> n;

    vector<vector<double>> A(n, vector<double>(n));
    cout << "Digite os elementos da matriz A:" << endl;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            cin >> A[i][j];
        }
    }

    vector<double> x(n, 1); // vetor inicial x
    double lambda, lambda_ant = 0; // autovalor atual e autovalor anterior
    do {
        vector<double> y = prodMatVec(A, x); // y = A*x
        lambda_ant = lambda;
        lambda = y[0] / x[0]; // autovalor
        x = divVec(norma(y), y); // autovetor
    } while (abs(lambda - lambda_ant) > EPS); // critério de parada

    cout << "Autovalor: " << lambda << endl;
    cout << "Autovetor: ";
    for (double v : x) {
        cout << v << " ";
    }
    cout << endl;

    return 0;
}
