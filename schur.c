#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include "mytimer.h"
#include "matrix1.h"
double eps = 1e-14;

/*
void dgemm_(char *transA, char *transB, int *m, int *n, int *k, double *alpha, double *A, int *lda, double *B, int *ldB, double *beta, double *C, int *ldc);

jobVS: if ='N' シューアベクトルは計算されない． else = 'V' 計算される



void dgees_(char *jobvs, char *sort, int *select, int *n, double *A, int *lda, int *sdim, double *wr, double *wi, double *vs, int *ldvs, double *work, int *lwork, bool *bwork, int *info);


static struct timespec time_start, time_end;

void timer_start(){
    clock_gettime(CLOCK_REALTIME, &time_start);
}

double timer_end(){
    clock_gettime(CLOCK_REALTIME, &time_end);
    return (time_end.tv_sec - time_start.tv_sec) + (time_end.tv_nsec - time_start.tv_nsec)/1e+9;
}
*/
double signs(double s){
    if(s > 0){
        return 1.0;
    }else if(s < 0){
        return -1.0;
    }
    else{
        printf("error %lf\n", s);
        return 0;
    }
}

//行列の転置 アルゴリズムでは共役転置をしているが実行列を対象としているので純粋な転置を行う．
//今後の課題は共役転置
void transm(double *a, double *at){
    return;
}


//matrix signum function via schur method 
void msf_schur(double *a, double *q, int size, double *ans){
    //まずschur分解する
    schur(a, q, size);
    double u[size*size];
    //uを0で初期化
    for(int i=0; i < size*size; i++){
        u[i] = 0;
    }
    for(int i=0; i < size; i++){
        //printf("%lf\n", a[i*size+i]);
        u[i*size + i] = signs(a[i*size+i]);
    }
    
    for(int j = 1; j < size; j++){
        for(int i=j-1; i > 0; i--){
            //printf("i = %d j = %d\n", i, j);
            if(u[i*size + i] +u[j*size + j] != 0){
                
                double tmp = 0.0;
                for(int k = i+1; k <= j-1; k++){
                    tmp += u[i*size + k]*u[k*size + j];
                }
                u[i*size + j] = -tmp/(u[i*size+i]+u[j*size+j]);
            }else{
                double tmp = 0.0;
                for(int k = i+1; k <= j-1; k++){
                    tmp += u[i*size + k]*a[k*size + j] - a[i*size + k]*u[k*size+j];
                }
                u[i*size + j] = a[i*size + j]*(u[i*size+i] - u[j*size + j])/(a[i*size + i] - a[j*size + j]) + tmp/(a[i*size + i] - a[j*size + j]);
            }
        }
    }
    double tmp[size*size], qt[size*size];
    mat_T(q, qt, size);
    //print_mat(u, size, size);
    mat_multi(q, u, tmp, size);
    mat_multi(tmp, qt, ans, size);
    return ;
}


int main(int argc, char *argv[]){

    int size= 1000, nn = size*size, incx = 1;
    double emin = 0.1, emax = 10;

    double a[size*size], I[size*size], q[size*size], ans[size*size], dif[size*size];

    setmatrix2(a, I, size, size, emax, emin, 3); 
    //print_mat(a, size, size);
    //print_mat(I, size, size);
    //clock_gettime(CLOCK_REALTIME, &ts);
    //printf("tv_sec= %ld, tv_nsec= %ld\n", ts.tv_sec, ts.tv_nsec);
    timer_start();
    msf_schur(a, q, size, ans);

    double time = timer_end();

    //print_mat(ans, size, size);
    mat_diff(dif, ans, I, size*size);
    double normf = dnrm2_(&nn, dif, &incx);
    double rel_err = calc_rerr_id(ans, size);
    //printf("size of matrix, minimum eigenvalue, max eigenvalue, log(err)\n");
    printf("%d %lf %lf %lf %lf %lf\n", size, emin, emax, log10(normf), log10(rel_err), time);
    return 0;

}
