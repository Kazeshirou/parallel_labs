#include <omp.h>
#include <stdio.h>

int main(int argc, char** argv) {
    int partial_sum, total_sum;

#pragma omp parallel private(partial_sum) shared(total_sum)
    {
        partial_sum = 0;
        total_sum   = 0;

#pragma omp for
        for (int i = 1; i <= 1000; i++) {
            partial_sum += i;
        }


#pragma omp critical
        { total_sum += partial_sum; }
    }
    printf("Total sum : % d\n", total_sum);

    return 0;
}
