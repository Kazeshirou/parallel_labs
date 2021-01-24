#include "gaussian.h"

#include <string.h>

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

        for (int k = 0; k < n; k++) {
            if (k == i) {
                continue;
            }

            update_row(matrix[i], matrix[k], n, j);
        }
    }

    for (int i = 0; i < n; i++) {
        int j = (int)a[i];
        x[j]  = matrix[i][n] / matrix[i][j];
    }

    return 0;
}
