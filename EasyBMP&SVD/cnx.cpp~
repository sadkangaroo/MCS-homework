#include "svd.h"
#include "svd.cpp"
#include "stdio.h"
#include <cstring>
#include <cstdlib>
#include <cmath>
#include "EasyBMP.h"
#include "EasyBMP.cpp"
#include <algorithm>

using namespace std;

const int MaxN = 3000 + 5, MaxM = 3000 + 5;
const double eps=1e-7;
double s[MaxN * MaxM], u[MaxN * MaxM], v[MaxN * MaxM], c[MaxN * MaxM];

int N, M, keep;

void solve(BMP &input, int typ) {
	for(int i=0;i<N;++i)
		for(int j=0;j<M;++j) {
            if (typ == 0) s[i*M+j]=input(j, i)->Red;
            if (typ == 1) s[i*M+j]=input(j, i)->Blue;
            if (typ == 2) s[i*M+j]=input(j, i)->Green;
        }
    double t1 = 0;
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < M; ++j)
            t1 += s[i * M + j] * s[i * M + j];
    t1 = sqrt(t1);
	memset(u,0,sizeof(u));
	memset(v,0,sizeof(v));
	if (dluav(s,N,M,u,v,eps) > 0) puts("OK");
    else puts("Oops");

	int x=(int)(min(N, M) / 100.0 * keep + 0.5);
	for(int i=0;i<N;++i)
		for(int j=0;j<M;++j)
			if(i!=j||i>=x)
				s[i*M+j]=0;
				
    
	damul(u,s,N, x, x, N, N, M,c);
	damul(c,v,N, x, M, N, M, M,s);

    double t2 = 0;
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < M; ++j)
            t2 += s[i * M + j] * s[i * M + j];
    t2 = sqrt(t2);
    printf("%.6lf\n", t2 / t1);

	for(int i=0;i<N;++i)
		for(int j=0;j<M;++j) {
			if (typ == 0) input(j, i)->Red=(int)(s[i*M+j]+0.5);
			if (typ == 1) input(j, i)->Blue=(int)(s[i*M+j]+0.5);
			if (typ == 2) input(j, i)->Green=(int)(s[i*M+j]+0.5);
        }
}

int main(int argc, char **arg) {

	BMP input;
    keep = atoi(arg[1]);
	input.ReadFromFile(arg[2]);
    N=input.TellHeight();
    M=input.TellWidth();
	printf("%d %d\n",N,M);
    solve(input, 0); solve(input, 1); solve(input, 2);
	input.WriteToFile(arg[3]);

	return 0;

}
