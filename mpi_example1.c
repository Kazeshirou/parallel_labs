#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

int all, rank;

void print_current_row(double* row, int n) {
    printf("[ ");
    for (int j = 0; j < n + 1; j++) {
        printf("%lf\t", row[j]);
    }
    printf("]");
}

// Внутри каждой итерации строки обновляются параллельно.
// i -  неизменяемая строка, используется для вычисления коэффициента.
// k - строка, которая будет обновлена.
// n - порядок СЛАУ.
// j - номер переменной, столбец которой обнуляется (по сути диагональный
// элемент, но так как в текущей реализации не требует, чтобы ненулевые элементы
// шли строго по диагонали, достаточно того, чтобы осталось по одному ненулевому
// элемент в каждой строке и в каждом столбце).
void update_row(double* i, double* k, int n, int j) {
    if (k[j] == 0) {
        return;
    }

    double t = k[j] / i[j];
    k[j]     = 0;
    for (int l = j + 1; l < n + 1; l++) {
        k[l] -= i[l] * t;
    }
}


int gaussian(double** src_matrix, double* x, int n) {
    int    n1 = n + 1;
    double recv_k[n1];
    int    iter  = n / all;
    int    iter1 = iter + 1;
    int    tail  = n % all;
    int    rcounts_row[all];
    int    rcounts[all];
    int    displs_row[all];
    int    displs[all];
    double recv_rows[iter1][n1];
    displs[0]     = 0;
    displs_row[0] = 0;
    for (int k = 0; k < all; k++) {
        rcounts_row[k] = iter;
        if (k && (k <= tail)) {
            rcounts_row[k] = iter1;
            rcounts[k]     = iter1 * n1;
        } else {
            rcounts_row[k] = iter;
            rcounts[k]     = iter * n1;
        }
        if (k) {
            displs_row[k] = displs_row[k - 1] + rcounts_row[k - 1];
            displs[k]     = displs[k - 1] + rcounts[k - 1];
        }
    }

    size_t  size   = sizeof(double) * n * (n + 1);
    double* matrix = (double*)malloc(sizeof(double) * n1 * n);
    double  a[n];
    if (!rank) {
        memcpy(matrix, src_matrix, size);
    }

    for (int i = 0; i < n; i++) {
        int offset = i * n1;
        int j      = 0;
        if (!rank) {
            while (matrix[offset + j] == 0) {
                j++;
                if (j > n) {
                    printf("Вырожденная матрица\n");
                    // print_current_row(matrix[i], n);
                    // printf("\n");
                    MPI_Abort(MPI_COMM_WORLD, 1);
                    return 1;
                }
            }
            a[i] = j;
        }
        // if (all == 1 && !rank) {
        //     for (int k = 0; k < n; k++) {
        //         if (k == i) {
        //             continue;
        //         }
        //         // printf("i = %d j = %d k = %d\ncurren_matrix = [\n", i, j,
        //         k); update_row(matrix[i], matrix[k], n, j);
        //         // for (int l = 0; l < n; l++) {
        //         //     printf("\t");
        //         //     print_current_row(matrix[l], n);
        //         //     printf("\n");
        //         // }
        //         // printf("]\n");
        //     }

        //     continue;
        // }

        MPI_Bcast(matrix + offset, n + 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(&j, 1, MPI_INT, 0, MPI_COMM_WORLD);

        MPI_Scatterv(matrix, rcounts, displs, MPI_DOUBLE, recv_rows[0],
                     rcounts[rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
        // printf("rank = %d i = %d rows = [\n", rank, i);

        for (int k = 0; k < rcounts_row[rank]; k++) {
            if (i == (displs_row[rank] + k)) {
                continue;
            }
            // printf("\t");
            // print_current_row(recv_rows[k], n);
            // printf("\n");
            update_row(matrix + offset, recv_rows[k], n, j);
        }

        MPI_Gatherv(recv_rows[0], rcounts[rank], MPI_DOUBLE, matrix, rcounts,
                    displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }

    if (!rank) {
        // printf("create answer\n");
        for (int i = 0; i < n; i++) {
            int j = (int)a[i];
            x[j]  = matrix[i * n1 + n] / matrix[i * n1 + j];
        }
    }

    free(matrix);
    return 0;
}

int main(int argc, char* argv[]) {
    // Ответ: -1, 1, 3.
    // int    n               = 3;
    // double matrix[][3 + 1] = {{2., 1., 1., 2.},
    //                           {1., -1., 0., -2.},
    //                           {3., -1., 2., 2.}};
    // double x[3]            = {0.};


    // Ответ: 4, -17, 13, -2.
    // int    n               = 4;
    // double matrix[][4 + 1] = {{3., 3., 3., 1., 1.},
    //                           {3., -2., -4., -2., 1.},
    //                           {3., -1., -2., 2., 2.},
    //                           {1., 2., 2., -2., 1.}};
    // double x[4]            = {0.};

    // Ответ: 1, 2, 3, 4, 10.
    // int    n               = 5;
    // double matrix[][5 + 1] = {{1., 1., 0., 0., 0., 3.},
    //                           {1., 1., 1., 0., 0., 6.},
    //                           {0., 1., 1., 1., 0., 9.},
    //                           {0., 0., 1., 1., 1., 17.},
    //                           {0., 0., 0., 0., 1., 10.}};
    // double x[5]            = {0.};

    // Ответ: 1, 2, 3, 4, 5, 6.
    // int    n               = 6;
    // double matrix[][6 + 1] = {
    //     {1., 1., 0., 0., 0., 0., 3.},  {1., 1., 1., 0., 0., 0., 6.},
    //     {0., 1., 1., 1., 0., 0., 9.},  {0., 0., 1., 1., 1., 0., 12.},
    //     {0., 0., 0., 1., 1., 1., 15.}, {0., 0., 0., 0., 1., 1., 11.}};
    // double x[6] = {0.};

    int     n = 1000;
    double* matrix;
    double* x;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &all);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (!rank) {
        matrix = (double*)malloc(sizeof(double) * n * (n + 1));
        x      = (double*)malloc(sizeof(double) * n);
        srand(time(NULL));

        int n1 = n + 1;
        for (int i = 0; i < n; i++) {
            int offset = i * n1;
            for (int j = 0; j < n + 1; j++) {
                matrix[offset + j] = rand();
            }
        }
        // printf("GAUSSIAN start with n = %d and \nmatrix = \n[\n", n);
        // for (int i = 0; i < n; i++) {
        //     printf("\t");
        //     print_current_row(matrix[i], n);
        //     printf("\n");
        // }
        // printf("]\n\n\n");
    }

    double start = MPI_Wtime();

    int res = gaussian((double**)matrix, x, n);

    double end = MPI_Wtime();

    if (!rank) {
        if (res) {
            printf("Что-то пошло не так\n");
        } else {
            printf("Время = %lf\n", end - start);
            // printf("Ответ= \n\t");
            // print_current_row(x, n - 1);
            // printf("\n");
        }
        free(matrix);
        free(x);
    }
    MPI_Finalize();


    return 0;
}
