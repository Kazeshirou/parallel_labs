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
    double recv_k[n + 1];

    int iter = n / all;
    int tail = n % all;
    int rcounts[all];
    int displs[all];
    for (int k = 0; k < all; k++) {
        if (k < tail) {
            rcounts[k] = n + 1;
            displs[k]  = k * rcounts[k];
        } else {
            rcounts[k] = 0;
            displs[k]  = 0;
        }
    }

    double matrix[n][n + 1];
    double a[n];
    if (!rank) {
        size_t size = sizeof(double) * n * (n + 1);
        memcpy(matrix, src_matrix, size);
    }

    for (int i = 0; i < n; i++) {
        int j = 0;
        if (!rank) {
            while (matrix[i][j] == 0) {
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
        if (all == 1 && !rank) {
            for (int k = 0; k < n; k++) {
                if (k == i) {
                    continue;
                }
                // printf("i = %d j = %d k = %d\ncurren_matrix = [\n", i, j, k);
                update_row(matrix[i], matrix[k], n, j);
                // for (int l = 0; l < n; l++) {
                //     printf("\t");
                //     print_current_row(matrix[l], n);
                //     printf("\n");
                // }
                // printf("]\n");
            }

            continue;
        }

        MPI_Bcast(matrix[i], n + 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(&j, 1, MPI_INT, 0, MPI_COMM_WORLD);

        for (int k = 0; k < iter; k++) {
            int current_begin = k * all;
            MPI_Scatter(matrix[current_begin], n + 1, MPI_DOUBLE, recv_k, n + 1,
                        MPI_DOUBLE, 0, MPI_COMM_WORLD);
            if ((current_begin + rank) != i) {
                update_row(matrix[i], recv_k, n, j);
            }
            MPI_Gather(recv_k, n + 1, MPI_DOUBLE, matrix[current_begin], n + 1,
                       MPI_DOUBLE, 0, MPI_COMM_WORLD);
        }

        if (tail && rcounts[rank]) {
            if (tail == 1) {
                // printf("tail %d\n", i);
                int current_begin = iter * all;
                if ((current_begin + rank) != i) {
                    update_row(matrix[i], matrix[current_begin + rank], n, j);
                }
            } else {
                int current_begin = iter * all;
                MPI_Scatterv(matrix[current_begin], rcounts, displs, MPI_DOUBLE,
                             recv_k, n + 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
                if ((current_begin + rank) != i) {
                    update_row(matrix[i], recv_k, n, j);
                }
                MPI_Gatherv(recv_k, n + 1, MPI_DOUBLE, matrix[current_begin],
                            rcounts, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            }
        }
    }

    if (!rank) {
        // printf("create answer\n");
        for (int i = 0; i < n; i++) {
            int j = (int)a[i];
            x[j]  = matrix[i][n] / matrix[i][j];
        }
    }
    return 0;
}

#define N 700

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

    int    n = N;
    double matrix[N][N + 1];
    double x[N];

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &all);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (!rank) {
        srand(time(NULL));

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n + 1; j++) {
                matrix[i][j] = rand();
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

    MPI_Finalize();
    if (!rank) {
        if (res) {
            printf("Что-то пошло не так\n");
        } else {
            printf("Время = %lf\n", end - start);
            // printf("Ответ= \n\t");
            // print_current_row(x, n - 1);
            // printf("\n");
        }
    }
    return 0;
}
