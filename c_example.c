#include <stdio.h>

#include "gaussian.h"

void print_current_row(double* row, int n) {
    printf("[ ");
    for (int j = 0; j < n + 1; j++) {
        printf("%lf\t", row[j]);
    }
    printf("]");
}

int main() {
    // int n = 3;
    // double matrix[][4] = {{2., 1., 1., 2.},
    //                       {1., -1., 0., -2.},
    //                       {3., -1., 2., 2.}};
    // double x[3]        = {0.};

    // int n = 4;
    // double matrix[][4 + 1] = {{3., 3., 3., 1., 1.},
    //                           {3., -2., -4., -2., 1.},
    //                           {3., -1., -2., 2., 2.},
    //                           {1., -2., 2., -2., 1.}};
    // double x[4]            = {0.};

    // Ответ: 1, 2, 3, 4, 10.
    int    n               = 5;
    double matrix[][5 + 1] = {{1., 1., 0., 0., 0., 3.},
                              {1., 1., 1., 0., 0., 6.},
                              {0., 1., 1., 1., 0., 9.},
                              {0., 0., 1., 1., 1., 17.},
                              {0., 0., 0., 0., 1., 10.}};
    double x[5]            = {0.};

    // int    n               = 6;
    // double matrix[][6 + 1] = {
    //     {1., 1., 0., 0., 0., 0., 3.},  {1., 1., 1., 0., 0., 0., 6.},
    //     {0., 1., 1., 1., 0., 0., 9.},  {0., 0., 1., 1., 1., 0., 12.},
    //     {0., 0., 0., 1., 1., 1., 15.}, {0., 0., 0., 0., 1., 1., 11.},
    // };
    // double x[6] = {0.};


    gaussian((double**)matrix, (double*)x, n);

    printf("GAUSSIAN end with answer = \n\t");
    print_current_row(x, n - 1);
    printf("\n");
    return 0;
}