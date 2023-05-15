#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include "utils.h"
using namespace std;

double residuo = 0.0001; // valor de tolerância para o critério de parada


int main() {
    int ICOD, iter;
    int n = 10;
    int sair = 1;
    string arquivo = "Matriz_A.dat";
    vector<vector<double>> A(n, vector<double>(n));

    cout << "Lendo a Matriz A:" << endl;
    lerMatriz(A,n,arquivo);
    cout << "Matriz A Lida:" << endl;
    imprimirMatriz(A);

    while (sair){
        cout << "Escolha o método de resolução (Metodo da Potência -> 1 ou\n Método de Jacobi -> 2):\n";
        cin >> ICOD;

        if (ICOD == 1){
            iter = 0;
            vector<double> x(n, 1);
            vector<double> y(n,0); 
            double lambda, lambda_ant = 1; // autovalor atual e autovalor anterior
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
            imprimirVetor(x);    
            }

        else if (ICOD == 2){
            iter = 0;
            vector<vector<double>> V(n, vector<double>(n));
            jacobi(A,V);
            imprimeAutovaloreeAutovetores(A, V, n);
        }

        cout << "Digite 0 para encerrar o programa ou 1 para escolher outro método: ";
        cin >> sair;
    }
    
    return 0;
}
