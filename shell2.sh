#!/bin/bash                                                                                 
#PJM -L "node=1"               # 1ノード                                                    
#PJM -L "rscgrp=small"         # リソースグループの指定                                     
#PJM -L "elapse=10:00"         # ジョブの経過時間制限値                                     
#PJM -g ra000006            # グループ指定                                                  
#PJM -x PJM_LLIO_GFSCACHE=/vol0003 # ジョブで使用するデータ領域のvolume                     
#PJM --mpi "max-proc-per-node=48" # 1ノードあたりに生成するMPIプロセス数の上限値             
#PJM -s                                                                                     

export PLE_MPI_STD_EMPTYFILE=off # 標準出力/標準エラー出力への出力がない場合はファイルを作成しない                                                                                    
export OMP_NUM_THREADS=8

# execute job                                                                               
for i in `seq 1 48`
do
    mpiexec -stdout-proc file_depara3 -n $i ./depara3
done

#depara3用file_depara3.4.0
#mpifccpx -Kfast,SVE,openmp matrix1.c de_para_sign3.c mytimer.c  -o depara3 -SSL2 -lm
# mpifccpx -Kfast,SVE,openmp matrix1.c schur.c mytimer.c  -o shcur -SSL2 -lm

#mpifccpx -Kfast,SVE,openmp,parallel matrix1.c de_para_sign3.c mytimer.c  -o depara3 -SSL2 -lm
# mpifccpx -Kfast,SVE,openmp,parallel matrix1.c schur.c mytimer.c  -o shcur -SSL2BLAMP -lm に変更 2023/2/28 
#1000 0.100000 10.000000 -12.302306 -13.802306 2.042561 変更後1.35倍速くなった