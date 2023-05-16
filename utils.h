#ifndef UTILS_H
#define UTILS_H
#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <cmath>
using namespace std;

void identidade(vector<vector<double>> &P, int &n){
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (i == j) {
                P[i][j] = 1.0;
            } 
            else {
                P[i][j] = 0.0;
            }
        }
    }   
}


void imprimirVetor(vector<double>&x){
    cout << "[";
    for (int i = 0; i < x.size(); i++) {
        cout << x[i] << ", ";
    }
    cout << "]" << endl;
}

void imprimirMatriz(vector<vector<double>>& matriz) {
    int n = matriz.size();
    int m = matriz[0].size();
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            cout << matriz[i][j] << " ";
        }
        cout << "\n";
    }
}

void lerVetor(vector<double> &vetor, int &n, string &arquivo){
    ifstream fin(arquivo); 
    char ch;
    string last_ch = " ";
    bool negativo = false;
    int contador= 0;

    while(fin.get(ch)){ 
        if (ch == ' '){
            if(last_ch == " "){
                continue;
            }
            else if (last_ch != " "){
                if (negativo){
                    last_ch = '-' + last_ch;
                }
                
                vetor[contador] = stod(last_ch);
                last_ch = ch;
                negativo = false;
                contador++;
                continue;
            }
            
        }

        else if (ch == '-'){
            negativo =true;
            continue;
        }

        else if (ch != ' '){
            if (last_ch == " "){
                last_ch = ch;
                continue;
            }

            else if (last_ch != " "){
                last_ch = last_ch + ch;
                continue;
            }
        }
    }
    vetor[contador] = stod(last_ch);    
}

void lerMatriz(vector<vector<double>> &A, int &n, string &arquivo){
 ifstream fin(arquivo); 
 char ch;
 string last_ch = " ";
 bool negativo = false;
 int contador = 0;
 int m = n*n;
 vector<double> vetor(m,0.0);

 while(fin.get(ch)){ 
    if (ch == ' '){
        if( last_ch == " "){
            continue;
        }
        else if (last_ch != " "){
            if (negativo){
                last_ch = '-' + last_ch + ".0";
            }
            
            //cout << last_ch << endl;
            vetor[contador] = stod(last_ch);
            last_ch = ch;
            negativo = false;
            contador++;
            continue;
        }
        
    }

    else if (ch == '-'){
        negativo =true;
        continue;
    }

    else if (ch != ' '){
        if (last_ch == " "){
            last_ch = ch;
            continue;
        }

        else if (last_ch != " "){
            last_ch = last_ch + ch;
            continue;
        }
        

    }

    }
    vetor[contador] = stod(last_ch);
    contador = 0;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            A[i][j] = vetor[contador];
            contador++;
        }
    }
}


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

void imprimeAutovaloreeAutovetores(vector<vector<double>>&A, vector<vector<double>>&V, int &n){
    for (int i = 0; i < n; i++) {
    cout << "Autovalor " << i+1 << ": " << A[i][i] << endl;
    cout << "Autovetor " << i+1 << ": ";
        for (int j = 0; j < n; j++) {
            cout << V[j][i] << " ";
        }
    cout << endl << endl;
    }
}



vector<vector<double>> decomposicaoCholesky(vector<vector<double>>& A, vector<vector<double>>& L, int n) {
    // Calcular a decomposição de Cholesky
    for (int i = 0; i < n; i++) {
        for (int j = 0; j <= i; j++) {
            double s = 0;
            for (int k = 0; k < j; k++) {
                s += L[i][k] * L[j][k];
            }
            if (i == j) {
                L[i][j] = sqrt(A[i][i] - s);
            } else {
                L[i][j] = (1.0 / L[j][j]) * (A[i][j] - s);
            }
        }
    }

    return L;
}

// Função para decompor a matriz A em L e U
void decomposicaoLU(vector<vector<double>>& A, vector<vector<double>>& L, int n) {
    // Executa o processo de eliminação de gauss para obter L e U
    for (int k = 0; k < n; k++) {
        L[k][k] = 1.0;
        for (int i = k + 1; i < n; i++) {
            L[i][k] = A[i][k] / A[k][k];
            for (int j = k; j < n; j++) {
                A[i][j] = A[i][j] - L[i][k] * A[k][j];
            }
        }
    }

}

bool converge(vector<vector<double>>&A){
    int n = A.size();
    double soma_linha, soma_coluna;
    for (int i = 0; i < n; i++){
        soma_linha = 0.0;
        soma_coluna = 0.0;

        for(int j = 0; j < n; j++ ){

            if(j == i){
                continue;
            }

            soma_linha += abs(A[i][j]);
            soma_coluna += abs(A[j][i]);
        }

        if(soma_linha > abs(A[i][i]) or soma_coluna > abs(A[i][i]) ){
            return false;
        }
    }
    return true;
}
// Função que calcula a solução do sistema linear Ax = B usando o método iterativo de Jacobi
vector<double> jacobi(vector<vector<double>>& A, vector<double>& B, int &n, double &tol, int &maxIter) {
    if(!converge(A)){
        cerr << "Matriz não converge !!" << endl;
        exit(1);
    }

    int k = 0;
    vector<double> Xold(n, 1.0);
    vector<double> Xnew(n, 0.0);
    double numerador,denominador;
    double residuo = tol + 1.0;

    while (residuo > tol && k < maxIter) {
        for (int i = 0; i < n; i++) {
            double soma = 0.0;
            for (int j = 0; j < n; j++) {
                if (j != i) {
                    soma += A[i][j] * Xold[j];
                }
            }
            Xnew[i] = (B[i] - soma) / A[i][i];
        }

        residuo = 0;
        numerador = 0;
        denominador = 0.0;
        for (int i = 0; i < n; i++) {
            numerador += pow(Xnew[i] - Xold[i], 2);
            denominador += pow(Xnew[i],2);
        }
        residuo = sqrt(numerador)/sqrt(denominador);
        Xold = Xnew;
        k++;
    }
    if (k == maxIter) {
        cout << "O método iterativo de Jacobi não convergiu em " << maxIter << " iterações!" << endl;
    } else {
        cout << "O método iterativo de Jacobi convergiu em " << k << " iterações!" << endl;
    }
    return Xnew;
}


vector<double> gauss_seidel(vector<vector<double>> &A, vector<double> &B, int &n, double &tol) {
    int iter = 0; // número de iterações
    vector<double> Xold(n, 1.0);// estimativa inicial da solução
    vector<double> Xnew(n, 0.0); 
    double residuo = tol + 1; // residuo inicial (qualquer valor maior que tol)
    double numerador,denominador;

    while (residuo > tol) {
        numerador = 0.0;
        denominador = 0.0;
        for (int i = 0; i < n; i++) {
            double soma = 0;
            for (int j = 0; j <= i-1; j++) {
                soma += A[i][j] * Xnew[j];
            }

            for (int j = i + 1; j < n; j++){
                //cout << "passou no 2 quando o i era " << i << " e o j era " << j <<endl;
                soma += A[i][j] * Xold[j];
            }


            Xnew[i] = (B[i] - soma)/A[i][i];  // nova estimativa da solução
            numerador += pow(Xnew[i] - Xold[i], 2);
            denominador += pow(Xnew[i],2);
            
            }
            residuo = sqrt(numerador)/sqrt(denominador);
            Xold = Xnew; // atualiza a estimativa da solução
            iter++; // incrementa o número de iterações
        }
    // verifica se o método convergiu ou não
    cout << "O método de Gauss-Seidel convergiu em " << iter << " iterações." << endl;

    return Xold;
}







// Função para resolver o sistema linear Ax = b usando a decomposição LU
vector<double> resolverSistemaLU(vector<vector<double>>& A, vector<vector<double>>& L, vector<double>& b, int n) {
    // Encontra a solução de Ly = b
    vector<double> y(n, 0.0);
    for (int i = 0; i < n; i++) {
        double soma = 0.0;
        for (int j = 0; j < i; j++) {
            soma += L[i][j] * y[j];
        }
        y[i] = b[i] - soma;
    }

    // Encontra a solução de Ux = y
    vector<double> x(n, 0.0);
    for (int i = n - 1; i >= 0; i--) {
        double soma = 0.0;
        for (int j = i + 1; j < n; j++) {
            soma += A[i][j] * x[j];
        }
        x[i] = (y[i] - soma) / A[i][i];
    }

    return x;
}

vector<double> resolverSistemaCholesky(vector<vector<double>>& A, vector<vector<double>>& L, vector<double>& B, int n) {
    vector<double> y(n, 0.0);
    vector<double> x(n, 0.0);

    // Resolver Ly = B usando substituição direta
    for (int i = 0; i < n; i++) {
        double s = 0.0;
        for (int j = 0; j < i; j++) {
            s += L[i][j] * y[j];
        }
        y[i] = (1.0 / L[i][i]) * (B[i] - s);
    }

    // Resolver L^T x = y usando substituição reversa
    for (int i = n - 1; i >= 0; i--) {
        double s = 0.0;
        for (int j = i + 1; j < n; j++) {
            s += L[j][i] * x[j];
        }
        x[i] = (1.0 / L[i][i]) * (y[i] - s);
    }

    return x;
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

// Função para calcular o maior elemento fora da diagonal
double maxElemento(vector<vector<double>>&A, int& p, int& q) {
    double max = 0.0;
    int n = A.size();
    for (int i = 0; i < n; i++) {
        for (int j = i+1; j < n; j++) {
            if (abs(A[i][j]) > max) {
                max = abs(A[i][j]);
                p = i;
                q = j;
            }
        }
    }
    return max;
}

// Função para calcular o ângulo de rotação
double anguloRotacao(vector<vector<double>>&A, int &p, int &q) {
    if (A[p][p] == A[q][q]) {
        return M_PI/4.0;
    } else {
        double tau = 2.0*A[p][q]/(A[p][p]-A[q][q]);
        return atan(tau)/2.0;
    }
}

void transporMatriz(vector<vector<double>>& matriz, vector<vector<double>>& matrizTransposta, int &n) {
   
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            matrizTransposta[j][i] = matriz[i][j];
        }
    }
}

void prodMatMat(vector<vector<double>> &X, vector<vector<double>> &P, int &n){
    vector<vector<double>> multi(n, vector<double>(n));    
    for(int i =0; i<n; i++){
        for(int j =0; j<n;j++){
            multi[i][j] = 0.0;
            for(int k = 0; k<n; k++){
                multi[i][j] += X[i][k] * P[k][j]; 
            }
        }
    }
    X = multi;
}

void prodMatMatMat(vector<vector<double>> &P_trans, vector<vector<double>> &A,vector<vector<double>> &P, int &n){
    vector<vector<double>> multi(n, vector<double>(n));    
    for(int i =0; i<n; i++){
        for(int j =0; j<n;j++){
            multi[i][j] = 0.0;
            for(int k = 0; k<n; k++){
                multi[i][j] += P_trans[i][k] * A[k][j]; 
            }
            
        }
    }

    for(int i =0; i<n; i++){
        for(int j =0; j<n;j++){
            A[i][j] = 0.0;
            for(int k = 0; k<n; k++){
                A[i][j] += multi[i][k] * P[k][j]; 
            }
    
        }
    }
}



// Função para fazer a rotação de Jacobi
void rotacaoJacobi(vector<vector<double>>&A, vector<vector<double>>&P, vector<vector<double>>&P_trans, vector<vector<double>>&X, int &p, int &q, int &n) {
    double angulo = anguloRotacao(A,p,q);
    double c = cos(angulo);
    double s = sin(angulo);
    P[p][p] = c;
    P[q][q] = c;
    P[p][q] = -s;
    P[q][p] = s;
    transporMatriz(P,P_trans,n);
    prodMatMatMat(P_trans,A,P,n);
    prodMatMat(X,P,n);
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




#endif