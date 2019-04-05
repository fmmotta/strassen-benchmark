#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include <math.h>
#define sgemm sgemm_


void sgemm (char *transA, char *transB, 
            int *m, int *n, int *k, 
            float *alpha, float *A, int *ldA, float *B, int *ldB,
            float *beta, float *C, int *ldC ) ; 



void imp_matriz(int m, int n, int ldA, float *A)
{
int i,j ;

for (i=0; i<m; i++)
  {
  if (i==0) printf("[");
     else printf(" ") ;
  for (j=0; j<n; j++)
     {
     printf( " %+10.6f ", A [i + j * ldA]);
     if (j < n-1) printf(",") ;
        else if (i<m-1) printf(";\n") ;
             else printf("];\n") ;
     }
  }

}

void strassen( float *A, float *B, float *C,int n, int min_size, int start_a, int start_b, int start_c)
    {
        float alpha = 1.0;
        float beta = 0.0;
        int i;
        int odd = 0;        
        if (n < min_size) {
            sgemm("N","N", &n, &n, &n, &alpha, A + start_a, &n, B + start_b, &n, &beta, C + start_c, &n);
            return;
        }
        else if (n == 1){
            C[start_c] = A[start_a]*B[start_b];
            return;
        }
        else
        {
            int old_n = n;
            while (ceil(log2(n)) != floor(log2(n))){
                n += 1;
                for (i = (n*n - 1); i >= 0; --i){
                    if ((i > (n*n - n - 1)) || (i+1)%n == 0){
                         A[i] = 0;
                         B[i] = 0;
                    }
                    else {
                        A[i] = A[i - (i/n)];
                        B[i] = B[i - (i/n)];
                    }
                }            
            }
            int block_size = n/2;
            int sq_size = block_size*block_size;
            
            for (i = 0; i < sq_size; ++i) {
                C[start_c + n*n + i] = A[start_a + i+block_size*((2*sq_size + i)/block_size)] - A[start_a + i+block_size*((2*sq_size + i + block_size)/block_size)];
                C[start_c + n*n + sq_size + i] = B[start_b + i+block_size*((block_size + i)/block_size)] + B[start_b + i+block_size*((2*sq_size + i + block_size)/block_size)]; 
            }
            strassen(C, C, C, block_size, min_size, start_c + n*n, start_c + n*n + sq_size, start_c + n*n + 2 * sq_size);
            
            for (i = 0; i < sq_size ; ++i) {
                C[start_c + i+block_size*(i/block_size)] = C[start_c + n*n + 2 * sq_size + i];
            }

            for (i = 0; i < sq_size; ++i) {     
                C[start_c + n*n + i] = A[start_a + i+block_size*(i/block_size)] + A[start_a + i+block_size*((2*sq_size + i + block_size)/block_size)];  
                C[start_c + n*n + sq_size + i] = B[start_b + i+block_size*(i/block_size)] + B[start_b + i+block_size*((2*sq_size + i + block_size)/block_size)];
            }
            strassen(C, C, C, block_size, min_size, start_c + n*n, start_c + n*n + sq_size, start_c + n*n + 2 * sq_size);
            
            for (i = 0; i < sq_size ; ++i) {
                C[start_c + i+block_size*(i/block_size)] += C[start_c + n*n + 2 * sq_size + i];
                C[start_c + i+block_size*((2*sq_size + block_size + i)/block_size)] = C[start_c + n*n + 2 * sq_size + i]; 
            }

            for (i = 0; i < sq_size; ++i) {     
                C[start_c + n*n + i] = A[start_a + i+block_size*(i/block_size)] - A[start_a + i+block_size*((block_size + i)/block_size)];
                C[start_c + n*n + sq_size + i] = B[start_b + i+block_size*(i/block_size)] + B[start_b + i+block_size*((2*sq_size + i)/block_size)];
            }
            strassen(C, C, C, block_size, min_size, start_c + n*n, start_c + n*n + sq_size, start_c + n*n + 2 * sq_size);
            
            for (i = 0; i < sq_size ; ++i) {
                C[start_c + i+block_size*((2*sq_size + block_size + i)/block_size)] -=  C[start_c + n*n + 2 * sq_size + i];
            }
            
            for (i = 0; i < sq_size; ++i) {     
                C[start_c + n*n + i] = A[start_a + i+block_size*(i/block_size)] + A[start_a + i+block_size*((2*sq_size + i)/block_size)];
                C[start_c + n*n + sq_size + i] = B[start_b + i+block_size*((2*sq_size + i + block_size)/block_size)]; //B_22 
            }
            strassen(C, C, C, block_size, min_size, start_c + n*n, start_c + n*n + sq_size, start_c + n*n + 2 * sq_size);

            for (i = 0; i < sq_size ; ++i) {
                C[start_c + i+block_size*(i/block_size)] -= C[start_c + n*n + 2 * sq_size + i];
                C[start_c + i+block_size*((2*sq_size + i)/block_size)] = C[start_c + n*n + 2 * sq_size + i];
            }

            for (i = 0; i < sq_size; ++i) {     
                C[start_c + n*n + i] = A[start_a + i+block_size*(i/block_size)]; //A_11         
                C[start_c + n*n + sq_size + i] = B[start_b + i+block_size*((2*sq_size + i)/block_size)] - B[start_b + i+block_size*((2*sq_size + i + block_size)/block_size)];
            }
            strassen(C, C, C, block_size, min_size, start_c + n*n, start_c + n*n + sq_size, start_c + n*n + 2 * sq_size);

            for (i = 0; i < sq_size ; ++i) {
                C[start_c + i+block_size*((2*sq_size + i)/block_size)] += C[start_c + n*n + 2 * sq_size + i];
                C[start_c + i+block_size*((2*sq_size + block_size + i)/block_size)] += C[start_c + n*n + 2 * sq_size + i];
            }

            for (i = 0; i < sq_size; ++i) {     
                C[start_c + n*n + i] = A[start_a + i+block_size*((2*sq_size + i + block_size)/block_size)]; //A_22 
                C[start_c + n*n + sq_size + i] = B[start_b + i+block_size*((block_size + i)/block_size)] - B[start_b + i+block_size*(i/block_size)];
            }
            strassen(C, C, C, block_size, min_size, start_c + n*n, start_c + n*n + sq_size, start_c + n*n + 2 * sq_size);

            for (i = 0; i < sq_size ; ++i) {
                C[start_c + i+block_size*(i/block_size)] += C[start_c + n*n + 2 * sq_size + i];
                C[start_c + i+block_size*((block_size + i)/block_size)] = C[start_c + n*n + 2 * sq_size + i];
            }

            for (i = 0; i < sq_size; ++i) {     
                C[start_c + n*n + i] = A[start_a + i+block_size*((block_size + i)/block_size)] + A[start_a + i+block_size*((2*sq_size + i + block_size)/block_size)];
                C[start_c + n*n + sq_size + i] = B[start_b + i+block_size*(i/block_size)]; //B_11 
            }
            strassen(C, C, C, block_size, min_size, start_c + n*n, start_c + n*n + sq_size, start_c + n*n + 2 * sq_size);

            for (i = 0; i < sq_size ; ++i) {
                C[start_c + i+block_size*((block_size + i)/block_size)] += C[start_c + n*n + 2 * sq_size + i];
                C[start_c + i+block_size*((2*sq_size + block_size + i)/block_size)] -= C[start_c + n*n + 2 * sq_size + i];
            }

            if (old_n != n) {
                for (i = 0; i < old_n*old_n; ++i) {
                    C[start_c + i] = C[start_c + i+(n - old_n)*(i/old_n)];    
                }
            }
            n = old_n;


            return;
        }
    }

int main(int argc, char **argv)

{

    int n = atoi(argv[1]);//(1 << atoi(argv[1])); //Bit shift trick to find powers of 2
    int min_size = atoi(argv[2]);
    int sq_size = n*n;
    int alloc_size = 4*n*n;//max_size*max_size;
    //float A[sq_size];
    //float B[sq_size];
    //float C[sq_size];
    float alpha = 1.0;
    float beta = 0.0;
    clock_t start, end;
    double cpu_time_used_strassen, cpu_time_used_blas;
    float * A = (float *) malloc(alloc_size * sizeof(float));
    float * B = (float *) malloc(alloc_size * sizeof(float));
    float * C = (float *) malloc(alloc_size * sizeof(float));
    float * D = (float *) malloc(alloc_size * sizeof(float));
    int start_a = 0;
    int start_b = 0;
    int start_c = 0;

    srand48(1) ;

    for (int i=0; i<sq_size;i++){
        A[i] = drand48();
        B[i] = drand48();
    }
    start = clock();
    sgemm("N","N", &n, &n, &n, &alpha, A, &n, B, &n, &beta, D, &n);
    end = clock();

    cpu_time_used_blas = ((double)(end - start))/CLOCKS_PER_SEC;

    start = clock();
    strassen(A, B, C, n, min_size, start_a, start_b, start_c);
    end = clock();
    cpu_time_used_strassen = ((double)(end - start))/CLOCKS_PER_SEC;

    //If checking for residue
/*    float * R = (float *) malloc(alloc_size * sizeof(float));
    for (int i=0; i<sq_size;i++){
        R[i] = D[i] - C[start_c + i];
    }
*/

    printf(" %i, %g, %g \n", n, cpu_time_used_strassen, cpu_time_used_blas); //print for script
    //printf(" %i, %g \n", min_size, cpu_time_used_strassen); //print for block size test
    //printf("TIMES: \n Time_strassen: %g \n Time_blas: %g \n",cpu_time_used_strassen, cpu_time_used_blas);
/*    printf("A:  \n\n");
    imp_matriz(n,n,n,A);
    printf("B:  \n\n");
    imp_matriz(n,n,n,B);
    printf("D:  \n\n");
    imp_matriz(n,n,n,D);
    printf("C: \n\n");
    imp_matriz(n,n,n,C);
   printf("R:  \n\n");
    imp_matriz(n,n,n,R);
    printf("\n\n");
*/
    free(A);
    free(B);
    free(C);
}
