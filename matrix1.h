#ifndef MATRIX1
#define MATRIX1


/*lapack関連プロトタイプ宣言*/
//-------------------------------------------------------------------------------------------//

//行列乗算
void dgemm_(char *transA, char *transB, int *m, int *n, int *k, double *alpha, double *A, int *lda, double *B, int *ldB, double *beta, double *C, int *ldc);
void dgetrf_(int *m, int *n, double *a, int *lda, int *ipiv, int *info);
void dgetri_(int *n, double *a, int *lda, int *ipiv, double *work, int *lwork, int *info);
//ユークリッドノルムを計算
double dnrm2_(int *n, double *x, int *incx);
//A*X = Bを計算
void dgesv_(int *n, int *nrhs, double *a, int *lda, int *ipiv, double *B, int *ldb, int *info);
//倍
void dscal_(int *n, double *da, double *dx, int *incx);
//足し算など
void daxpy_(int *n, double *da, double *dx, int *incx, double *dy, int *incy);
//シューア分解
void dgees_(char *jobvs, char *sort, bool (*SELECT)(double, double), int *n, double *A, int *lda, int *sdim, double *wr, double *wi, double *vs, int *ldvs, double *work, int *lwork, bool *bwork, int *info);
//qr分解用
void dgeqr2_(int *m, int *n, double *a, int *lda, double *tau, double *work, int *info);
void dgeqrf_(int *m, int *n, double *a, int *lda, double *tau, double *work, int *lwork, int *info); //こっちでいい

//行列normの計算 norm='f'ならdnrm2と(たぶん)同じ
double dlange_(char *norm, int *m, int *n, double *a, int *lda, double *work);

//行列コピー
void dlacpy_(char *uplo, int *m, int *n, double *A, int *lda, double *B, int *ldb);

void dorgqr_(int *m, int *n, int *k, double *A, int *lda, double *tau, double *work, int *lwork, int *info);
//------------------------------------------------------------------------------------------//

//逆行列計算
void matrix_inv(double *a, int size);

//Aの累乗(A^2), a2に結果が入る
void matrixA(double *a, double *a2, int size);

//行列の積を計算 a*b = a2
void mat_multi(double *a, double *b, double *a2, int size);

//行列の和を計算する関数 n*n正方行列 aにbを加える
//lapackサブルーチンdaxpyを利用
//出力結果は5番目の引数、DYに入っている
//dy = dy + αdxを計算。α=1とすれば純粋な加算
//α=-1とすれば減算ができる
//sizeには行列の要素数を入れる
void mat_sum(double *a, double *b, int size, int da);

//入力行列が対称な場合のときに2つの行列の和を計算するための関数
void mat_sum_sym(double *a, double *b, int n);

//mat_sumでlapackルーチンを使わない版．これはlapack版との性能を比較したいときに使うのがいいかも
void mat_sum_osoi(double *a, double *b, int n);


//行列を出力する関数
void print_mat(double *a, int n, int m);

//行列の差を求める
void mat_diff(double *dif, double *A, double *sol, int size);

//サイズsize*sizeの正方行列aに実数xを掛ける関数
//行列のときsize = n*nが入るようにする。
void mat_multidouble(double *a, double x, int size);

//n*m行列, Iは単位行列
void setmatrix(double *A, double *I, int n, int m, double emax, double emin);

//seedを引数によって変えられるようにする
void setmatrix2(double *A, double *I, int n, int m, double emax, double emin, int seed);

//単位行列を計算しない関数
void setmatrix3(double *A, int n, double emax, double emin, int seed);

//aをシューア分解し，計算後シューア標準形が入る．vsにはシューアベクトル行列が入る
void schur(double *a, double *vs, int size);


void mat_T(double *a, double *at, int size);

void mat_diff(double *dif, double *A, double *sol, int size);

void file_read(char *filename, double *a, int n);

void identity(double *I, int n);

double calc_rerr_id(double *A, int size);

double calc_err_id(double *A, int size);

void qr(double *q, double *r, int nsize);

double dif_norm(double *A, double *B, int n);

void matcpy(double *A, double *B, int n);

void setmatrix_X(double *A, int n, double emax, double emin, int seed, double cond, char qrs);

void setmatrix_sym(double A[], int n, double emax, double emin, int seed);

double normf(double *A, int n);

void sign_jordan(double *eig, double *evec, int size);

void setmatrix_sym2(double *A, int n, double emax, double emin, int seed, double *ans, char ch);

#endif
