#include <stdio.h>
#include <stdlib.h>

#define N 256  // 可根據需要調整大小，太大 Cachegrind 會跑很久

int main() {
    int i, j, k;
    double sum_ijk = 0.0, sum_ikj = 0.0, sum_kji = 0.0;

    // 建立 3D 陣列
    double ***A = (double ***)malloc(N * sizeof(double **));
    for (i = 0; i < N; i++) {
        A[i] = (double **)malloc(N * sizeof(double *));
        for (j = 0; j < N; j++) {
            A[i][j] = (double *)malloc(N * sizeof(double));
        }
    }

    // 初始化
    for (i = 0; i < N; i++)
        for (j = 0; j < N; j++)
            for (k = 0; k < N; k++)
                A[i][j][k] = i + j + k;

    // ===== 巢狀迴圈 ijk =====
    for (i = 0; i < N; i++)
        for (j = 0; j < N; j++)
            for (k = 0; k < N; k++)
                sum_ijk += A[i][j][k];

    // ===== 巢狀迴圈 ikj =====
    for (i = 0; i < N; i++)
        for (k = 0; k < N; k++)
            for (j = 0; j < N; j++)
                sum_ikj += A[i][j][k];

    // ===== 巢狀迴圈 kji =====
    for (k = 0; k < N; k++)
        for (j = 0; j < N; j++)
            for (i = 0; i < N; i++)
                sum_kji += A[i][j][k];

    printf("sum_ijk = %f\n", sum_ijk);
    printf("sum_ikj = %f\n", sum_ikj);
    printf("sum_kji = %f\n", sum_kji);

    // 釋放記憶體
    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++)
            free(A[i][j]);
        free(A[i]);
    }
    free(A);

    return 0;
}