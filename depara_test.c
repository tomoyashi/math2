/*
新たに作り直す．
このプログラムでは，諸々の行列計算(行列積など)を行う関数をヘッダファイルに入れて呼び出せるようにする．


また，並列計算の見直しを行う．

できれば，行列をファイルから読み込んで計算できるようにしたい．
*/

#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <mpi.h>
#include <stdbool.h>
#include "matrix1.h"
#include "mytimer.h"

#define SIZE 1000
double eps = 1e-15;

//この関数の場合は行列を返しているわけではなく、第4引数目のtmpに計算結果が入っていることに注意する
void func1(double x, double a2[], int n, double tmp[]){
    //double tmp[n*n];
    int i, c = 0;
    for(i = 0; i < n*n; i++){
        tmp[i] = a2[i];
        if(i == c*n + c){
            tmp[i] += x*x;
            c++;
        }
    }
    return;
}

double func2(double x){
    return exp(M_PI/2.0*sinh(x));
}

double func3(double x){
    return M_PI/2.0*exp(M_PI /2.0*sinh(x)) * cosh(x);
}

void deformula(double h, int size, double A[], double a2[], int N, int nprocs, int myrank, double tmp[], int points[]){    
    //double eps1 = eps*size*size;    
    int size2 = size*size, i;
    double tmp2[size2], tmp3[size2];
    int info, ipiv[size];
    
    //プロセス0のときのみx=0のときの値をtmpに入るようにする．
    if(myrank == 0){
        //memcpy(tmp, A, sizeof(tmp2));
        matcpy(A, tmp, size);
        func1(func2(0.0), a2, size, tmp3);
        dgesv_(&size, &size, tmp3, &size, ipiv, tmp, &size, &info);
        mat_multidouble(tmp, func3(0.0), size2);
    }
    int n_loop, istart;
    if(N%nprocs == 0){
        n_loop = N/nprocs;
        istart = myrank*n_loop + 1;
    }else{
        n_loop = N/nprocs + 1;
        istart = myrank*n_loop + 1;
        int marge = N-(N/nprocs)*nprocs;
        if(myrank == marge){
            //istart = myrank*n_loop + 1;
            n_loop = N/nprocs;
        }else if(myrank >= marge+1){
            istart = myrank*(N/nprocs) + marge + 1;
            n_loop = N/nprocs;
        }
    }
    //
    for(i = istart; (i < istart + n_loop) && (i <= N) ; i++){
        double x = points[i]*h;
        matcpy(A, tmp2, size);
        func1(func2(x), a2, size, tmp3);
        dgesv_(&size, &size, tmp3, &size, ipiv, tmp2, &size, &info);
        mat_multidouble(tmp2, func3(x), size2);
        mat_sum(tmp, tmp2, size2, 1);
    }
    //printf("\n");
}


int main(int argc, char *argv[])
{
    double emin = 0.1, emax=10, A[SIZE*SIZE], A2[SIZE*SIZE], etime;
    //double I[SIZE*SIZE];    setmatrix2(A, I, SIZE, SIZE, emax, emin, 4);
    //int cond = atoi(argv[1]);
    
    //setmatrix_X(A, SIZE, emax, emin, 1, cond, 'Y');
    setmatrix_sym(A, SIZE, emax, emin, 1);

    double h = 0.1;    
    int nmax = 80, nn = SIZE*SIZE;//こちらの場合はnmaxは正負の分点数の合計
    int points[nmax+1];
    for(int i = 1; i < nmax+1; i+=2){
        points[i] = i/2+1;
        points[i+1] = -(i/2+1);
    }
    MPI_Init(&argc, &argv);
    int myrank, nprocs;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    //int npm[2]={0,0}, npmsum[2]={0,0};
    
    //tmpは初期化しておく
    double tmp[nn], tmp2[nn];
    for(int j = 0; j < nn; j++){
        tmp[j] = 0;
    }
    
    MPI_Barrier(MPI_COMM_WORLD);
    timer_start();
    matrixA(A, A2, SIZE);//A^2の計算
    deformula(h, SIZE, A, A2, nmax, nprocs, myrank, tmp, points);
    MPI_Reduce(tmp, tmp2, nn, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if(myrank == 0){
        mat_multidouble(tmp2, h/M_PI*2.0, SIZE*SIZE);
        etime = timer_end();
        //double rerr = calc_rerr_id(tmp2, SIZE);
        double rerr = calc_err_id(tmp2, SIZE);
        //print_mat(tmp2, SIZE, SIZE);
        //printf h:刻み幅, 分点数, 絶対誤差, log10(絶対誤差), 行列サイズ(行), プロセス数，計算時間    
        printf("%f %d %.15f %f %d %d %0.8f\n", h, nmax+1, rerr, log10(rerr), SIZE, nprocs, etime);
        //printf("%f %d %.13f %f %d %d %0.8f %d\n", h, nmax*2+1, rerr, log10(rerr), SIZE, nprocs, etime, cond);
    }


    MPI_Finalize();
    //fclose(fp);
    return 0;
}