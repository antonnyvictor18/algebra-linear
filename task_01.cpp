#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include "utils.h"
using namespace std;


// Função principal
int main() {
    int n, ICOD, maxIter;
    double tol = 0.0001;
    int sair  = 1;
    string arquivo = "Matriz_A.dat";
    string id;
    n = 10;

    // Cria as matrizes A, L e U e o vetor b
    vector<vector<double>> A(n, vector<double>(n, 0.0));
    vector<vector<double>> L(n, vector<double>(n, 0.0));
    vector<vector<double>> U(n, vector<double>(n, 0.0));
    vector<double> b(n, 0.0);
    vector<double> x;

    // Lê os valores de A e b
    cout << "Lendo a matriz A... " << endl;
    lerMatriz(A,n,arquivo);
    cout << "Matriz A Lida: " << endl;
    imprimirMatriz(A);

    
    while (sair){
        cout << "Escolha o método de resolução (Decomposição LU (ICOD =1) ou Decomposição de Cholesky (ICOD =2)): ";
        cin >> ICOD;
        if (ICOD == 1){
            // Executa a decomposição LU
            decomposicaoLU(A, L, n);
            for (int i = 1; i<= 3; i++){
                id = to_string(i);
                cout << "Lendo o vetor B"+id+"...\n";
                arquivo = "Vetor_B_0"+id+".dat";
                lerVetor(b,n,arquivo);
                cout << "Vetor B_0"+id+" Lido : " << endl;
                imprimirVetor(b);

                // Resolve o sistema linear Ax = b usando a decomposição LU
                x = resolverSistemaLU(A, L, b, n);
                cout << "A solução do sistema é:\n";
                imprimirVetor(x);
            }
                    
        }   
            
            
        else if (ICOD == 2){
            decomposicaoCholesky(A,L,n);
            for (int i = 1; i<= 3; i++){
                id = to_string(i);
                cout << "Lendo o vetor B_0"+id+"...\n";
                arquivo = "Vetor_B_0"+id+".dat";
                lerVetor(b,n,arquivo);
                cout << "Vetor B_0"+id+" Lido:" << endl;
                imprimirVetor(b);

                // Resolve o sistema linear Ax = b usando Cholesky
                x = resolverSistemaCholesky(A,L,b,n);
                cout << "A solução do sistema é:\n";
                imprimirVetor(x);
            }
        }

        else if (ICOD == 3){
            
            cout << "Qual a quantidade máxima de iterações desejada? ";
            cin >> maxIter;
            for (int i = 1; i<= 3; i++){
                id = to_string(i);
                cout << "Lendo o vetor B_0"+id+"...\n";
                arquivo = "Vetor_B_0"+id+".dat";
                lerVetor(b,n,arquivo);
                cout << "Vetor B_0"+id+" Lido:" << endl;
                imprimirVetor(b);

                // Resolve o sistema linear Ax = b usando Cholesky
                x = jacobi(A,b,tol,maxIter);
                cout << "A solução do sistema é:\n";
                imprimirVetor(x);
            }

        
        }
        else if (ICOD == 4){

            // cout << "Qual a tolerância ? ";
            // cin >> tol;
            for (int i = 1; i<= 3; i++){
                id = to_string(i);
                cout << "Lendo o vetor B_0"+id+"...\n";
                arquivo = "Vetor_B_0"+id+".dat";
                lerVetor(b,n,arquivo);
                cout << "Vetor B_0"+id+" Lido:" << endl;
                imprimirVetor(b);

                // Resolve o sistema linear Ax = b usando Cholesky
                x = gauss_seidel(A,b,tol);
                cout << "A solução do sistema é:\n";
                imprimirVetor(x);
            }
        }
        else {
            cerr << "ICOD não definido !"<<endl;
            return EXIT_FAILURE;
        }
        cout << "Deseja escolher outro método? Se sim, escola 1. Do contrário, entre com qualquer digito. ";
        cin >> sair;
    }
    return 0;
}
