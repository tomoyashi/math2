#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
//#include <mpi.h>

//#define SIZE 10

void dgemm_(char *transA, char *transB, int *m, int *n, int *k, double *alpha, double *A, int *lda, double *B, int *ldB, double *beta, double *C, int *ldc);
void dgetrf_(int *m, int *n, double *a, int *lda, int *ipiv, int *info);
void dgetri_(int *n, double *a, int *lda, int *ipiv, double *work, int *lwork, int *info);
double dnrm2_(int *n, double *x, int *incx);
void dgesv_(int *n, int *nrhs, double *a, int *lda, int *ipiv, double *B, int *ldb, int *info);
void dscal_(int *n, double *da, double *dx, int *incx);
void daxpy_(int *n, double *da, double *dx, int *incx, double *dy, int *incy);
void dgees_(char *jobvs, char *sort, bool (*SELECT)(double, double), int *n, double *A, int *lda, int *sdim, double *wr, double *wi, double *vs, int *ldvs, double *work, int *lwork, bool *bwork, int *info);

bool SELECT(const double wr, const double wi){
    return true;
}

//引数の行列aを逆行列にする関数
void matrix_inv(double *a, int size){
    int n = size; //行のサイズ
    int m = size; //列のサイズ
    int lda = size; //mと同じ値
    int info; //計算が成功すれば0を返す
    int ipiv[size]; //要素数はm, nのうち小さいほうとする
    int lwork = size; //nと同じ値
    double work[size]; //要素数はlworkと同じ値
    dgetrf_(&n, &m, a, &lda, ipiv, &info);
    dgetri_(&m, a, &lda, ipiv, work, &lwork, &info);
}

//行列A^2を求める, a2に結果が入る
//sizeには正方行列R^(n*n)のnが入る
void matrixA(double a[], double a2[], int size){
    double alpha = 1.0, beta = 0.0;
    int n = size;
    dgemm_("n", "n", &n, &n, &n, &alpha, a, &n, a, &n, &beta, a2, &n);

    return;
}

//行列の積を計算 計算結果は3引数目の変数に入る
//sizeには行列の行数
void mat_multi(double a[], double b[], double a2[], int size){
    double alpha = 1.0, beta = 0.0;
    int n = size;
    dgemm_("n", "n", &n, &n, &n, &alpha, a, &n, b, &n, &beta, a2, &n);

    return;
}

//行列の和を計算する関数 n*n正方行列 aにbを加える
//lapackサブルーチンdaxpyを利用
//出力結果は5番目の引数、DYに入っている(ここではdaxpyの5番目の引数には配列aを入れてるのでaの結果が更新される)
//dy = dy + αdxを計算。α=1とすれば純粋な加算
//sizeのところには行列の要素数を入れる
void mat_sum(double a[], double b[], int size, int da){
    int n=size, incx = 1, incy = 1;
    double da1= da; 
    daxpy_(&n, &da1, b, &incx, a, &incy);
    return;
}
//ある行列に対称行列を足すための関数
//対称性を利用する
//nは行数(正方行列)
//a = a + b
void mat_sum_sym(double a[], double b[], int n){
    for(int i=0; i < n; i++){
        for(int j = i; j < n; j++){
            if(i == j){
                a[i + j*n] += b[i + j*n];
            }else{
                a[i+j*n] += b[i+j*n];
                a[j + i*n] += b[i + j*n];
            }
        }
    }
}

void mat_sum_osoi(double a[], double b[], int n){
    for(int i=0; i<n*n; i++){
        a[i] += b[i];
    }
}


//行列を出力する関数
void print_mat(double a[], int n, int m){
    for(int i = 0; i < n; i++){
        for(int j = 0; j < m; j++){
            printf("%f ", a[i + j*m]);
        }
        printf("\n");
    }
    printf("\n");
    return ;
}


//サイズsize*sizeの正方行列aに実数xを掛ける関数
//行列のときsize = n*nが入るようにする。
void mat_multidouble(double a[], double x, int size){
    int s = size, incx = 1;
    dscal_(&s, &x, a, &incx);
    return;
}


//n*m行列, Iは単位行列
void setmatrix(double A[], double I[], int n, int m, double emax, double emin){
    double X[n*m]; //まず適当な行列を作成する
    double Xinv[n*m];
    srand(2);//計算結果比較のためシード固定
    //srand((unsigned)time(NULL));
    for(int i = 0; i < n*m; ++i){
        X[i] = ((double)rand()/ ((double)RAND_MAX + 1));
        Xinv[i] = X[i];
    }
    //print_mat(X, n, n);
    double a[n*m];//ここ初期化すればよくないか？後で確認
    for(int i = 0; i < n*m; i++){
        a[i] = 0;
        I[i] = 0;
    }
    //ここ上のループでできるかも
    for(int i = 0; i < n; i++){
        I[i*n+i] = 1;
        a[i*n + i] = ((double)rand()/ ((double)RAND_MAX))*(emax-emin) + emin;
    }
    a[0] = emax; a[1*n + 1] = emin;//対角成分はi*n + i  ここでは固有値を好みの値にセット
    
    matrix_inv(Xinv, n);
    //print_mat(Xinv, n, n);
    double tmp[n*m];
    mat_multi(X, a, tmp, n);
    mat_multi(tmp, Xinv, A, n);   

    return ;
}

//行列生成シードを引数で指定できるようにする
void setmatrix2(double A[], double I[], int n, int m, double emax, double emin, int seed){
    double X[n*m]; //まず適当な行列を作成する
    double Xinv[n*m];
    srand(seed);//計算結果比較のためシード固定
    //srand((unsigned)time(NULL));
    for(int i = 0; i < n*m; ++i){
        X[i] = ((double)rand()/ ((double)RAND_MAX + 1));
        Xinv[i] = X[i];
    }
    //print_mat(X, n, n);
    double a[n*m];//ここ初期化すればよくないか？後で確認
    for(int i = 0; i < n*m; i++){
        a[i] = 0;
        I[i] = 0;
    }
    //ここ上のループでできるかも
    for(int i = 0; i < n; i++){
        I[i*n+i] = 1;
        a[i*n + i] = ((double)rand()/ ((double)RAND_MAX))*(emax-emin) + emin;
    }
    a[0] = emax; a[1*n + 1] = emin;//対角成分はi*n + i  ここでは固有値を好みの値にセット
    
    matrix_inv(Xinv, n);
    //print_mat(Xinv, n, n);
    double tmp[n*m];
    mat_multi(X, a, tmp, n);
    mat_multi(tmp, Xinv, A, n);   

    return ;
}
//シューア分解を行う関数
//aをシューア分解し，計算後シューア標準形が入る．vsにはシューアベクトル行列が入る
void schur(double *a, double *vs, int size){
    double wr[size], wi[size], work[3*size];
    char jobvs = 'V', sort = 'N';
    int n = size, lda = size, ldvs = size, lwork = 3*size, sdim, info;
    bool bwork[size];
    dgees_(&jobvs, &sort, SELECT, &n, a, &lda, &sdim, wr, wi, vs, &ldvs, work, &lwork, bwork, &info);
    return;
}

//行列aの転置をatに入れる sizeは行の長さ
void mat_T(double *a, double *at, int size){
    for(int i = 0; i < size; i++){
        at[i*size + i] = a[i*size + i];
    }
    //print_mat(a, size, size);
    for(int j = 1; j < size; j++){
        for(int i = j-1; i >= 0; i--){
            at[i*size + j]= a[j*size + i];
            at[j*size + i] = a[i*size + j];
            //printf("%d %d ", i*size + j, j*size + i);
        }
    }
    return;
}

//行列の差を求める
void mat_diff(double dif[], double A[], double sol[], int size){
    for(int i=0; i < size; i++){
        dif[i] = A[i] - sol[i]; 
    }
    return;
}


//nは行のサイズ，aにはファイルの値がlapack用に入れられる．filenameはファイル名(拡張子含む)
void file_read(char *filename, double *a, int n){
    FILE *fp;
    if((fp = fopen(filename, "r")) == NULL){
        printf("file can not open\n");
        exit(1);
    }
    //無限ループ怖い．注意
    int mat_row_num = 0, mat_column_num = 0;
    for(;;){
        /*
        bufに行を読み込むバッファオーバーランするかもしれない．
        おかしくなったらここ確認
        エラー発生したらbufのサイズ大きくする方がよき
        buf
        1つの数字あたり18文字あると考えて，それがn個あり，空白がn個あり，一応10を足しておく．
        */
        char buf[22*n+n];

        if(fgets(buf, sizeof(buf), fp) == NULL){
            if(feof(fp)){
                break;
            }
            else{
                fputs("error happen!\n", stderr);
                exit(EXIT_FAILURE);
            }
        }

        //末尾が改行文字であれば，'\0'で上書きする
        char* p = strchr(buf, '\n');
        if(p != NULL){
            *p = '\0';
        }

        /*
        以下で行をdelimiterで区切られた数値をそれぞれ配列に入れていく．
        ここでは，delimeterは空白となっている．データファイルの種類によって変更する．
        stは各数値の始まりの位置を覚えておく．
        文字を走査していき，delimiterにあったらbufのstポインタからi-st+1(これが数値の長さ)をtmpに格納する．
        tmpをatofによってdouble型に変換．
        atofはエラー判定ができないため注意．strtofの方が安全そうである．
        */
        //puts(buf);
        int st = 0;
        for(int i = 0; buf[i] != '\0'; i++){
            if(buf[i+1] == '\0'){
                char tmp[32];
                strncpy(tmp, buf+st, i - st + 1);
                double num = atof(tmp);
                a[mat_row_num + n*mat_column_num] = num;
                //printf("%f ", num);
                mat_column_num = 0;
            }
            else if(buf[i] == ' '){
                char tmp[32];
                strncpy(tmp, buf+st, i - st);
                double num = atof(tmp);
                a[mat_row_num + n*mat_column_num] = num;
                //printf("%f ", a[mat_num]);
                st = i+1;
                mat_column_num++;
            }
        }
        mat_row_num++;
        //free(buf);
    }
    if(fclose(fp) == EOF){
        fputs("file cant close!!\n", stderr);
        exit(EXIT_FAILURE);
    }
}

//関数実行後，配列Iに単位行列が格納される．sizeは行の大きさ(正方行列なので列も同じ)
void identity(double *I, int n){
    for(int i = 0; i < n*n; i++){
        if(i/n == i%n){
            I[i] = 1.;
        }else{
            I[i] = 0.;
        }
    }
}

//単位行列との相対誤差を求める 正方行列とする．引数にはn*nのnを入れる(行か列のサイズ)．
//返り値として相対誤差を返す
double calc_rerr_id(double *A, int size){
    //incxは1にすることで全要素の計算ができる
    int nn = size*size, incx = 1;
    double dif[size*size], I[size*size], normdif, normsol;
    //対角成分のときは1でそれ以外は0なので，それぞれから引くようにする
    for(int i = 0; i < size*size; i++){
        if(i%size == (i/size)){//対角成分のとき
            dif[i] = 1.0 - A[i];
            I[i] = 1.0;
        }
        else{
            dif[i] = 0.0 - A[i];
            I[i] = 0.0;
        }
    }

    normdif = dnrm2_(&nn, dif, &incx);
    normsol = dnrm2_(&nn, I, &incx);

    return normdif/normsol;
}

//返り値として絶対誤差を求める
double calc_err_id(double *A, int size){
    //incxは1にすることで全要素の計算ができる
    int nn = size*size, incx = 1;
    double dif[size*size], normdif;
    //対角成分のときは1でそれ以外は0なので，それぞれから引くようにする
    for(int i = 0; i < size*size; i++){
        if(i%size == (i/size)){//対角成分のとき
            dif[i] = 1.0 - A[i];
        }
        else{
            dif[i] = 0.0 - A[i];
        }
    }

    normdif = dnrm2_(&nn, dif, &incx);

    return normdif;
}