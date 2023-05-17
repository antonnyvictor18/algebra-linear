#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include "utils.h"
using namespace std;


int main() {
    int ICOD, iter;
    int n = 10;
    double tol;
    string arquivo = "Matriz_A.dat";
    vector<vector<double>> A(n, vector<double>(n));

    cout << "Lendo a Matriz A:" << endl;
    lerMatriz(A,n,arquivo);
    cout << "Matriz A Lida:" << endl;
    imprimirMatriz(A);

    cout << "Escolha o método de resolução (Metodo da Potência -> 1 ou Método de Jacobi -> 2): ";
    cin >> ICOD;

    cout << "Escolha uma tolerância (Entre com um valor de ponto flutuante entre 0 e 1): ";
    cin >> tol;

    if (ICOD == 1){
        iter = 0;
        vector<double> x(n, 1.0);
        vector<double> y(n,0.0); 
        double lambda = 1;
        double lambda_ant; 
        do {
            lambda_ant = lambda;
            y = prodMatVec(A, x);
            lambda = normaMaiorValor(y);
            x = divVec(lambda,y);
            iter++;
        } while (abs((lambda - lambda_ant)/lambda) > tol); 
            
        cout << "Número de iterações: " << iter << endl;
        cout << "Autovalor: " << lambda << endl;
        cout << "Autovetor: ";
        imprimirVetor(x);    
    }

    else if (ICOD == 2){
        iter = 0;
        int p, q;
        double max, determinante;
        vector<vector<double>> P(n, vector<double>(n));
        vector<vector<double>> P_trans(n, vector<double>(n));
        vector<vector<double>> X(n, vector<double>(n));
        identidade(P,n);
        identidade(X,n);
        max = maxElemento(A, p, q);
        while (max > tol) { // Critério de parada
        
            cout << "iteração: " << iter << ", valor máximo fora da diagonal: "<< max <<endl;
            cout << "Matriz P: " <<endl;
            imprimirMatriz(P);

            cout << "Matriz A: " <<endl;
            imprimirMatriz(A);

            cout << "Matriz X: " << endl;
            imprimirMatriz(X);

            cout << "\n";
            identidade(P,n);
            rotacaoJacobi(A, P, P_trans, X, p, q, n);
            max = maxElemento(A, p, q);
            iter++;
        }
        cout << "Matriz com Autovalores: "<<endl;
        imprimirMatriz(A);
        cout << "Matriz com Autovetores: "<<endl;
        imprimirMatriz(X);
        cout << "Determinante de A: " << det(A,n) << endl;

        
    }
    return 0;
}
