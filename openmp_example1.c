#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

int check_res(double* matrix, double* x, int n) {
    int n1 = n + 1;
    for (int i = 0; i < n; i++) {
        int    offset = i * n1;
        double lpart  = 0;
        for (int j = 0; j < n; j++) {
            lpart += matrix[offset + j] * x[j];
        }
        if (abs(lpart - matrix[offset + n]) > 1e-3) {
            return 0;
        }
    }
    return 1;
}

void print_current_row(double* row, int n) {
    printf("[ ");
    for (int j = 0; j < n + 1; j++) {
        printf("%lf\t", row[j]);
    }
    printf("]");
}


void update_row(double* i, double* k, int n, int j) {
    double t = k[j] / i[j];
    k[j]     = 0;
    for (int l = j + 1; l < n + 1; l++) {
        k[l] -= i[l] * t;
    }
}


int gaussian(double** src_matrix, double* x, int n) {
    size_t size = sizeof(double) * n * (n + 1);
    double matrix[n][n + 1];
    double a[n];
    memcpy(matrix, src_matrix, size);

    for (int i = 0; i < n; i++) {
        int j = 0;
        while (matrix[i][j] == 0) {
            j++;
            if (j > n) {
                return 1;
            }
        }
        a[i] = j;
#pragma omp parallel
        {
            int k;
#pragma omp for private(k)
            for (k = 0; k < n; k++) {
                if (k == i) {
                    continue;
                }
                update_row(matrix[i], matrix[k], n, j);
            }
        }
    }
#pragma omp parallel
    {
        int i;
#pragma omp for private(i)
        for (i = 0; i < n; i++) {
            int j = (int)a[i];
            x[j]  = matrix[i][n] / matrix[i][j];
        }
    }
    return 0;
}

int main(int argc, char** argv) {
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


    double start = omp_get_wtime();

    int rv = gaussian((double**)matrix, x, n);

    double end = omp_get_wtime();
    int    f   = check_res((double*)matrix, x, n);
    printf("rv = %d Корректность решения = %d Время = %lf\n", rv, f,
           end - start);
    // printf("Ответ= \n\t");
    // print_current_row(x, n - 1);
    // printf("\n");

    free(matrix);
    free(x);

    return 0;
}
