#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>

using namespace std;

vector<vector<double>> getSubmatrix(const vector<vector<double>>& A, int p, int q, int n) {
    vector<vector<double>> temp(n - 1, vector<double>(n - 1, 0.0));
    int i = 0, j = 0;
    
    for (int row = 0; row < n; row++) {
        for (int col = 0; col < n; col++) {
            if (row != p && col != q) {
                temp[i][j] = A[row][col];
                j++;
                if (j == n - 1) {
                    j = 0;
                    i++;
                }
            }
        }
    }
    return temp;
}

double getDeterminant(const vector<vector<double>>& A, int n) {
    if (n == 1) return A[0][0];
    if (n == 2) return A[0][0] * A[1][1] - A[0][1] * A[1][0];

    double det = 0;
    int sign = 1;
    for (int f = 0; f < n; f++) {
        vector<vector<double>> sub = getSubmatrix(A, 0, f, n);
        det += sign * A[0][f] * getDeterminant(sub, n - 1);
        sign = -sign;
    }
    return det;
}

bool inverseByDeterminant(const vector<vector<double>>& A, vector<vector<double>>& inv, int n) {
    double det = getDeterminant(A, n);
    
    if (det == 0) {
        cout << "[오류] 행렬식이 0이므로 역행렬이 존재하지 않습니다." << endl;
        return false;
    }

    if (n == 1) {
        inv[0][0] = 1.0 / A[0][0];
        return true;
    }

    vector<vector<double>> adj(n, vector<double>(n, 0.0));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            vector<vector<double>> sub = getSubmatrix(A, i, j, n);
            int sign = ((i + j) % 2 == 0) ? 1 : -1;
            adj[j][i] = sign * getDeterminant(sub, n - 1); 
        }
    }

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            inv[i][j] = adj[i][j] / det;
        }
    }
    return true;
}

bool inverseByGaussJordan(const vector<vector<double>>& A, vector<vector<double>>& inv, int n) {
    vector<vector<double>> aug(n, vector<double>(2 * n, 0.0));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            aug[i][j] = A[i][j];
        }
        aug[i][i + n] = 1.0; 
    }

    for (int i = 0; i < n; i++) {
        if (aug[i][i] == 0) {
            cout << "[오류] 피벗이 0이 되어 가우스-조던 소거법을 진행할 수 없습니다." << endl;
            return false;
        }

        double pivot = aug[i][i];
        for (int j = 0; j < 2 * n; j++) {
            aug[i][j] /= pivot;
        }

        for (int k = 0; k < n; k++) {
            if (k != i) {
                double factor = aug[k][i];
                for (int j = 0; j < 2 * n; j++) {
                    aug[k][j] -= factor * aug[i][j];
                }
            }
        }
    }

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            inv[i][j] = aug[i][j + n];
        }
    }
    return true;
}

void printMatrix(const vector<vector<double>>& mat, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            cout << setw(10) << fixed << setprecision(4) << mat[i][j] << " ";
        }
        cout << endl;
    }
}

bool compareMatrices(const vector<vector<double>>& M1, const vector<vector<double>>& M2, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (abs(M1[i][j] - M2[i][j]) > 0.0001) {
                return false;
            }
        }
    }
    return true;
}

int main() {
    int n;

    cout << "정방행렬의 차수를 입력하세요: ";
    if (!(cin >> n) || n <= 0) return 0;

    vector<vector<double>> A(n, vector<double>(n, 0.0));
    for (int i = 0; i < n; i++) {
        cout << i + 1 << "행: ";
        for (int j = 0; j < n; j++) {
            cin >> A[i][j];
        }
    }

    cout << endl;

    vector<vector<double>> invDet(n, vector<double>(n, 0.0));
    vector<vector<double>> invGJ(n, vector<double>(n, 0.0));

    cout << "    행렬식  으로 구한 역행렬:" << endl;
    bool successDet = inverseByDeterminant(A, invDet, n);
    if (successDet) {
        printMatrix(invDet, n);
    }

    cout << endl;

    cout << "가우스-조던 소거법으로 구한 역행렬:" << endl;
    bool successGJ = inverseByGaussJordan(A, invGJ, n);
    if (successGJ) {
        printMatrix(invGJ, n);
    }

    cout << endl;

    if (successDet && successGJ) {
        if (compareMatrices(invDet, invGJ, n)) {
            cout << "두 방법의 결과가 동일합니다." << endl;
        } else {
            cout << "두 방법의 결과가 다릅니다." << endl;
        }
    }

    return 0;
}
