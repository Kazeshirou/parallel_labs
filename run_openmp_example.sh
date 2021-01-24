OMP_NUM_THREADS=$1;
export OMP_NUM_THREADS

gcc openmp_example1.c -o openmp_example -fopenmp
./openmp_example
