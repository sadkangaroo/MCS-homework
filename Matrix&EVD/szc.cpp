#include "evd.h"
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <iomanip>
using namespace std;
using namespace splab;

const int MaxN = 20 + 5;
const int N = 20;

double dis[MaxN][MaxN];
string name[MaxN];

#define sq(x) ((x) * (x))

double xx[MaxN], yy[MaxN];

int main() {

    freopen("test.in", "r", stdin);
    freopen("test.out", "w", stdout);

    for (int i = 0; i < N; ++i) for (int j = 0; j < N; ++j) dis[i][j] = 0;
    for (int i = 0; i < N; ++i) for (int j = 0; j < i; ++j) {
        cin >> dis[i][j];
        dis[j][i] = dis[i][j];
    }
    for (int i = 0; i < N; ++i) cin >> name[i];
    cout << setiosflags(ios::fixed) << setprecision(3);

    Matrix<double> A(N, N); 
    for (int i = 0; i < N; ++i) for (int j = 0; j < N; ++j) {
        A[i][j] = sq(dis[i][j]);
        double tmp = 0;
        for (int k = 0; k < N; ++k) tmp += sq(dis[i][k]);
        tmp /= N;
        A[i][j] -= tmp;
        tmp = 0;
        for (int k = 0; k < N; ++k) tmp += sq(dis[k][j]);
        tmp /= N;
        A[i][j] -= tmp;
        tmp = 0;
        for (int k1 = 0; k1 < N; ++k1) for (int k2 = 0; k2 < N; ++k2)
           tmp += sq(dis[k1][k2]); 
        tmp /= N; tmp /= N;
        A[i][j] += tmp;
        A[i][j] *= -0.5;
    }
    EVD<double> eig;
    eig.dec(A);
    Matrix<double> D, V;
    V=eig.getV();
    D=diag(eig.getD());
    for (int i = 0; i < N; ++i) for (int j = 0; j < N; ++j)
        if (i != j || i >= 2) D[i][j] = 0;
        else D[i][j] = sqrt(D[i][j]);
    V = V * D;
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < 1; ++j) 
            printf("%.3lf ", V[i][j]);
        puts("");
    }
    puts("");
    for (int i = 0; i < N; ++i) {
        for (int j = 1; j < 2; ++j) 
            printf("%.3lf ", V[i][j]);
        puts("");
    }
    puts("");
    for (int i = 0; i < N; ++i) cout << name[i] << endl;

    return 0;

}
