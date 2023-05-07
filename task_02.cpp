#include <iostream>
#include <cmath>
#include <vector>

using namespace std;

double residuo = 0.001; // valor de tolerância para o critério de parada

// função que realiza o produto matriz-vetor
vector<double> prodMatVec(vector<vector<double>> &A, vector<double> &x) {
    int n = A.size();
    vector<double> y(n,0);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            y[i] += A[i][j] * x[j];
        }
    }
    return y;
}

// função que calcula a norma euclidiana de um vetor
double normaEuclidiana(vector<double> v) {
    double s = 0;
    for (double x : v) {
        s += x * x;
    }
    return sqrt(s);
}

double normaMaiorValor(vector<double> &y){
    double maior = 0;
    for (int i = 0; i < y.size(); i++){
        if (y[i] > maior){
            maior = y[i];
        }
    }
    return maior;
}

// função que divide um vetor por um escalar
vector<double> divVec(double &lambda, vector<double>&y) {
    int n = y.size();
    vector<double> x(n,0);
    for (int i = 0; i < y.size(); i++) {
        x[i] = y[i]/lambda;
    }
    return x;
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
    int iter = 0;
    vector<double> x(n, 1);
    vector<double> y(n,0); // vetor inicial x
    double lambda, lambda_ant = 0; // autovalor atual e autovalor anterior
    do {
        lambda_ant = lambda;
        y = prodMatVec(A, x);
        lambda = normaMaiorValor(y);
        x = divVec(lambda,y);
        iter++;
    } while (abs((lambda - lambda_ant)/lambda) > residuo); // critério de parada
    
    cout << "Número de iterações: " << iter << endl;
    cout << "Autovalor: " << lambda << endl;
    cout << "Autovetor: ";
    for (double v : x) {
        cout << v << " ";
    }
    cout << endl;

    return 0;
}
