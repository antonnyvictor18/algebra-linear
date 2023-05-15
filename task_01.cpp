#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
using namespace std;


void lerVetor(vector<double> &vetor, int &n, string &arquivo){
    ifstream fin(arquivo); 
    char ch;
    string last_ch = " ";
    bool negativo = false;
    int contador, contador2 = 0;

    while(fin.get(ch)){ 
        if (ch == ' '){
            if(last_ch == " " or contador2 == 0){
                contador2++;
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
 int contador, contador2 = 0;
 vector<double> vetor(n*n,0);

 while(fin.get(ch)){ 
    if (ch == ' '){
        if(last_ch == " " or contador2 == 0){
            contador2++;
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
    contador = 0;
    for (int i = 0; i < 10; i++) {
        for (int j = 0; j < 10; j++) {
            A[i][j] = vetor[contador];
            contador++;
        }
    } 
}


void imprimirVetor(vector<double>&x){
    cout << "A solução do sistema é:\n";
    for (int i = 0; i < x.size(); i++) {
        cout << "x" << i + 1 << " = " << x[i] << "\n";
    }
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
void decomposicaoLU(vector<vector<double>>& A, vector<vector<double>>& L, vector<vector<double>>& U, int n) {
    // Inicializa L com a matriz identidade e U com a matriz A
    for (int i = 0; i < n; i++) {
        L[i][i] = 1.0;
        for (int j = 0; j < n; j++) {
            U[i][j] = A[i][j];
        }
    }

    // Executa o processo de eliminação de gauss para obter L e U
    for (int k = 0; k < n; k++) {
        for (int i = k + 1; i < n; i++) {
            L[i][k] = U[i][k] / U[k][k];
            for (int j = k; j < n; j++) {
                U[i][j] = U[i][j] - L[i][k] * U[k][j];
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
vector<double> jacobi(vector<vector<double>>& A, vector<double>& B, double &tol, int &maxIter) {
    if(!converge(A)){
        cerr << "Matriz não converge !!" << endl;
        return;
    }

    int n = A.size();
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

        residuo, numerador, denominador = 0.0;
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


vector<double> gauss_seidel(vector<vector<double>> &A, vector<double> &B, double &tol, int &maxIter) {
    int n = A.size();
    int iter = 0; // número de iterações
    vector<double> Xold(n, 1);// estimativa inicial da solução
    vector<double> Xnew(n, 0); 
    double residuo = tol + 1; // residuo inicial (qualquer valor maior que tol)
    double numerador,denominador;

    while (residuo > tol && iter < maxIter) {
        residuo, numerador, denominador = 0.0;
        for (int i = 0; i < n; i++) {
            double soma = 0;
            for (int j = 0; j <= i-1; j++) {
                    soma += A[i][j] * Xnew[j];
            }

            for (int j = i + 1; j < n; j++){
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
    if (iter == maxIter && residuo > tol) {
        cerr << "O método de Gauss-Seidel não convergiu em " << maxIter << " iterações." << endl;
    }

    return Xnew;
}







// Função para resolver o sistema linear Ax = b usando a decomposição LU
vector<double> resolverSistemaLU(vector<vector<double>>& A, vector<vector<double>>& L, vector<vector<double>>& U, vector<double>& b, int n) {
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
            soma += U[i][j] * x[j];
        }
        x[i] = (y[i] - soma) / U[i][i];
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



// Função principal
int main() {
    int n, ICOD, sair;
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

    cout << "Escolha o método de resolução (Decomposição LU (ICOD =1) ou Decomposição de Cholesky (ICOD =2)): ";
    cin >> ICOD;

    if (ICOD == 1){
        // Executa a decomposição LU
        decomposicaoLU(A, L, U, n);
        for (int i = 1; i<= 3; i++){
            id = to_string(i);
            cout << "Lendo o vetor B"+id+"...\n";
            arquivo = "Vetor_B_0"+id+".dat";
            lerVetor(b,n,arquivo);
            cout << "Vetor B_0"+id+" Lido : " << endl;
            imprimirVetor(b);

            // Resolve o sistema linear Ax = b usando a decomposição LU
            x = resolverSistemaLU(A, L, U, b, n);
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
        double tol = 0.0001;
        int maxIter;

        cout << "Qual a quantidade máxima de iterações desejada? ";
        cin >> maxIter;

        x = jacobi(A,b,tol,maxIter);
        imprimirVetor(x);
    
    }
    else if (ICOD == 4){
        double tol;
        int maxIter;

        cout << "Qual a tolerância ? ";
        cin >> tol;

        cout << "Qual a quantidade máxima de iterações desejada? ";
        cin >> maxIter;

        x = gauss_seidel(A,b,tol,maxIter);
        imprimirVetor(x);
    }
    else {
        cerr << "ICOD não definido !"<<endl;
        return EXIT_FAILURE;
    }
    
    return 0;
}
