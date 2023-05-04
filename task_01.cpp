#include <iostream>
#include <vector>
#include <cmath>
using namespace std;

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
    
    cout << "Digite a ordem da matriz: ";
    cin >> n;

    // Cria as matrizes A, L e U e o vetor b
    vector<vector<double>> A(n, vector<double>(n, 0.0));
    vector<vector<double>> L(n, vector<double>(n, 0.0));
    vector<vector<double>> U(n, vector<double>(n, 0.0));
    vector<double> b(n, 0.0);
    vector<double> x;

    // Lê os valores de A e b
    cout << "Digite os valores da matriz A:\n";
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            cin >> A[i][j];
        }
    }

    cout << "Escolha o método de resolução (Decomposição LU (ICOD =1) ou Decomposição de Cholesky (ICOD =2)): ";
    cin >> ICOD;

    if (ICOD == 1){
        // Executa a decomposição LU
        decomposicaoLU(A, L, U, n);
        while (true){
            cout << "Digite os novoss valores do vetor b:\n";
            for (int i = 0; i < n; i++) {
                cin >> b[i];
            }

            // Resolve o sistema linear Ax = b usando a decomposição LU
            x = resolverSistemaLU(A, L, U, b, n);
            cout << "A solução do sistema é:\n";
            for (int i = 0; i < n; i++) {
                cout << "x" << i + 1 << " = " << x[i] << "\n";
            }

            cout << "desaja sair ? ('1' para sim ou '0' para não) ";
            cin >> sair;

            if (sair == 1){
                return 0;
            }
            
        }
                
    }   
        
        
    else if (ICOD == 2){
        decomposicaoCholesky(A,L,n);
        while (true){
            cout << "Digite os novoss valores do vetor b:\n";
            for (int i = 0; i < n; i++) {
                cin >> b[i];
            }
            // Resolve o sistema linear Ax = b usando a decomposição LU
            x = resolverSistemaCholesky(A,L,b,n);
            cout << "A solução do sistema é:\n";
            for (int i = 0; i < n; i++) {
                cout << "x" << i + 1 << " = " << x[i] << "\n";
            }

            cout << "desaja sair ? ('1' para sim ou '0' para não) ";
            cin >> sair;

            if (sair == 1){
                return 0;
            }
        }
    }
    else {
        cerr << "ICOD não definido"<<endl;
        return EXIT_FAILURE;
    }
    
    return 0;
}
