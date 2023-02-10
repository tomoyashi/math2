#Last modified: 2023/01/30 

#通常コンパイル時にはCCを利用．MPIを使うときはMPICC
#lapackを使う計算のときにはLAPACK
#必要であればコンパイルオプションを追加する
#コピペするとタブ入力が空白に変ってしまい正しく動作しなくなってしまうらしいので注意


CC           = gcc
MPICC        = mpicc
CFLAGS       = -O3 -Wall
LAPACK       = -l lapack -l blas 
MATRIX_LIB   = matrix1.o
TIMER        = mytimer.o

all: schur 

depara3: de_para_sign3.o $(TIMER) $(MATRIX_LIB)
	$(MPICC) $(CFLAGS) -o $@ $^ $(LAPACK) -lm

schur: schur.o $(TIMER) $(MATRIX_LIB)
	$(CC) $(CFLAGS) -o $@ $^ $(LAPACK) -lm 

matrix1.o: matrix1.c
	$(CC) $(CFLAGS) -c $< $(LAPACK) -lm

schur.o: schur.c
	$(CC) $(CFLAGS) -c $< $(LAPACK) -lm

mytimer.o: mytimer.c 
	$(CC) $(CFLAGS) -c $< 

de_para_sign3.o: de_para_sign3.c
	$(MPICC) $(CFLAGS) -c $< $(LAPACK) -lm

#生成した実行ファイルなど消去
.PHONY: clean
clean:
	rm -f schur matrix1.o schur.o depara3 de_para_sign3.o
 