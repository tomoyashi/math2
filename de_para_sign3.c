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

void deformula(double h, int size, double A[], double a2[], int N, int nprocs, int myrank, double tmp[]){    
    //double eps1 = eps*size*size;    
    int size2 = size*size, i;
    double tmp2[size2], tmp3[size2];
    int info, ipiv[size];
    
    //プロセス0のときのみx=0のときの値をtmpに入るようにする．
    if(myrank == 0){
        memcpy(tmp, A, sizeof(tmp2));
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
    //changeがtrueのとき前半分のプロセスでこれで正の部分を計算
    for(i = istart; (i < istart + n_loop) && (i <= N) ; i++){
        //printf("i: %d ", i);
        double x = i*h;
        //printf("x:%f i:%d ",x, i);
        memcpy(tmp2, A, sizeof(tmp2));
        func1(func2(x), a2, size, tmp3);
        dgesv_(&size, &size, tmp3, &size, ipiv, tmp2, &size, &info);
        mat_multidouble(tmp2, func3(x), size2);
        int da = 1;
        mat_sum(tmp, tmp2, size2, da);
        //print_mat(tmp, size, size);
        //int incx = 1, nn = size2;
        //double normf = dnrm2_(&nn, tmp2, &incx);
        //printf("%g %d\n", normf, i);
    }
    for(i = istart; (i < istart + n_loop) && (i <= N); i++){
        double x = -i*h;
        //printf("x:%f i:%d ",x, i);
        memcpy(tmp2, A, sizeof(tmp2));
        func1(func2(x), a2, size, tmp3);
        dgesv_(&size, &size, tmp3, &size, ipiv, tmp2, &size, &info);
        mat_multidouble(tmp2, func3(x), size2);
        int da = 1;
        mat_sum(tmp, tmp2, size2, da);
        //print_mat(tmp, size, size);
        //int incx = 1, nn = size2;
        //double normf = dnrm2_(&nn, tmp2, &incx);
        //printf("%g %d\n", normf, i);
    }
    //printf("\n");
}


int main(int argc,char *argv[])
{
    double emin = 0.1, emax=10, A[SIZE*SIZE], A2[SIZE*SIZE], I[SIZE*SIZE], etime;
    setmatrix2(A, I, SIZE, SIZE, emax, emin, 4);
    double h = 0.1;    
    int nmax = 45, nn = SIZE*SIZE;
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
    deformula(h, SIZE, A, A2, nmax, nprocs, myrank, tmp);
    MPI_Reduce(tmp, tmp2, nn, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if(myrank == 0){
        mat_multidouble(tmp2, h/M_PI*2.0, SIZE*SIZE);
        etime = timer_end();
        //double rerr = calc_rerr_id(tmp2, SIZE);
        double err = calc_err_id(tmp2, SIZE);
        //printf h:刻み幅, 分点数, 絶対誤差, log10(絶対誤差), 行列サイズ(行), プロセス数，計算時間    
        printf("%f %d %.10f %f %d %d %0.8f\n", h, nmax*2+1, err, log10(err), SIZE, nprocs, etime);
    }


    MPI_Finalize();
    //fclose(fp);
    return 0;
}
