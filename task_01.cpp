#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

// Função para imprimir a matriz
void print_matrix(vector<vector<double>> &A, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            cout << A[i][j] << "\t";
        }
        cout << endl;
    }
    cout << endl;
}

// Função para imprimir o vetor
void print_vector(vector<double> &x, int n) {
    for (int i = 0; i < n; i++) {
        cout << x[i] << endl;
    }
    cout << endl;
}

// Função para decompor a matriz A em LU
void LU_decomposition(vector<vector<double>> &A, vector<vector<double>> &L, vector<vector<double>> &U, int n) {
    for (int i = 0; i < n; i++) {
        // Calcula os elementos da matriz L
        for (int j = 0; j <= i; j++) {
            double s = 0;
            for (int k = 0; k < j; k++) {
                s += L[i][k] * U[k][j];
            }
            L[i][j] = A[i][j] - s;
        }
        // Calcula os elementos da matriz U
        for (int j = i + 1; j < n; j++) {
            double s = 0;
            for (int k = 0; k < i; k++) {
                s += L[i][k] * U[k][j];
            }
            U[i][j] = (A[i][j] - s) / L[i][i];
        }
        L[i][i] = 1; // Define o elemento diagonal da matriz L como 1
    }
}

// Função para resolver o sistema linear de equações utilizando a decomposição LU
void LU_solver(vector<vector<double>> &L, vector<vector<double>> &U, vector<double> &B, vector<double> &x, int n) {
    vector<double> y(n);
    // Resolve o sistema Ly = B
    for (int i = 0; i < n; i++) {
        double s = 0;
        for (int j = 0; j < i; j++) {
            s += L[i][j] * y[j];
        }
        y[i] = (B[i] - s) / L[i][i];
    }
    // Resolve o sistema Ux = y
    for (int i = n - 1; i >= 0; i--) {
        double s = 0;
        for (int j = i + 1; j < n; j++) {
            s += U[i][j] * x[j];
        }
        x[i] = (y[i] - s) / U[i][i];
    }
}

// Função para decompor a matriz A em Cholesky
void Cholesky_decomposition(vector<vector<double>> &A, vector<vector<double>> &L, int n) {
    for (int i = 0; i < n; i++) {
    // Calcula os elementos da matriz L
    for (int j = 0; j <= i; j++) {
        double s = 0;
        for (int k = 0; k < j; k++) {
            s += L[i][k] * L[j][k];
        }
        if (i == j) {
            L[i][i] = sqrt(A[i][i] - s);
        } else {
            L[i][j] = (A[i][j] - s) / L[j][j];
        }
        }
    }
}
// Função para resolver o sistema linear de equações utilizando a decomposição de Cholesky
void Cholesky_solver(vector<vector<double>> &L, vector<double> &B, vector<double> &x, int n) {
    vector<double> y(n);
    // Resolve o sistema Ly = B
    for (int i = 0; i < n; i++) {
        double s = 0;
        for (int j = 0; j < i; j++) {
            s += L[i][j] * y[j];
        }

        y[i] = (B[i] - s) / L[i][i];
    }
    // Resolve o sistema L^Tx = y
    for (int i = n - 1; i >= 0; i--) {
        double s = 0;
        for (int j = i + 1; j < n; j++) {
            s += L[j][i] * x[j];
        }
        x[i] = (y[i] - s) / L[i][i];
    }
}

int main() {
    int n;
    cout << "Digite o tamanho da matriz A: ";
    cin >> n;
    // Aloca as matrizes A, L e U e o vetor B
vector<vector<double>> A(n, vector<double>(n));
vector<vector<double>> L(n, vector<double>(n));
vector<vector<double>> U(n, vector<double>(n));
vector<double> B(n);
vector<double> x(n);

// Preenche a matriz A e o vetor B
cout << "Digite os elementos da matriz A:" << endl;
for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
        cin >> A[i][j];
    }
}
cout << "Digite os elementos do vetor B:" << endl;
for (int i = 0; i < n; i++) {
    cin >> B[i];
}

int ICOD;
cout << "Digite o código do método a ser utilizado (1 para LU e 2 para Cholesky): ";
cin >> ICOD;

if (ICOD == 1) {
    // Realiza a decomposição LU da matriz A
    LU_decomposition(A, L, U, n);
} else if (ICOD == 2) {
    // Realiza a decomposição de Cholesky da matriz A
    Cholesky_decomposition(A, L, n);
} else {
    cout << "Código inválido." << endl;
    return 0;
}

// Resolve o sistema linear de equações para vários vetores B
while (true) {
    cout << "Digite os elementos do vetor B ou -1 para sair:" << endl;
    bool flag = false;
    for (int i = 0; i < n; i++) {
        cin >> B[i];
        if (B[i] == -1) {
            flag = true;
            break;
        }
    }
    if (flag) {
        break;
    }
    if (ICOD == 1) {
        // Resolve o sistema linear de equações utilizando a decomposição LU
        LU_solver(L, U, B, x, n);
    } else if (ICOD == 2) {
        // Resolve o sistema linear de equações utilizando a decomposição de Cholesky
        Cholesky_solver(L, B, x, n);
    }
    // Imprime o resultado
    cout << "A solução do sistema é: ";
    for (int i = 0; i < n; i++) {
        cout << x[i] << " ";
    }
    cout << endl;
}

return 0;
}