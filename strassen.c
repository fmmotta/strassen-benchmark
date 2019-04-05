#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <stdbool.h>
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

void strassen( float *A, float *B, float *C,int n, int min_size, int start, int start_orig, int r)
    {
        float alpha = 1.0;
        float beta = 0.0;
        int i;
        if (n < min_size) {
            sgemm("N","N", &n, &n, &n, &alpha, A, &n, B, &n, &beta, C, &n);
            return;
        }
        else if (n == 1){
            C[start_orig] = A[start_orig]*B[start_orig];
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
            if (old_n != n) {
                start = n*n;                    
            }
            //imp_matriz(n,n,n,A);

            //printf("sf 1 \n");
            int block_size = n/2;
            int sq_size = block_size*block_size;
            
            //printf("bs %i sq %i \n",block_size,sq_size);
            //printf("sf 3 \n");

            for (i = 0; i < sq_size; ++i) {
                A[i + start] = A[start_orig + i+block_size*((2*sq_size + i)/block_size)] - A[start_orig + i+block_size*((2*sq_size + i + block_size)/block_size)];
                B[i + start] = B[start_orig + i+block_size*((block_size + i)/block_size)] + B[start_orig + i+block_size*((2*sq_size + i + block_size)/block_size)]; 

                A[i + start + sq_size] = A[start_orig + i+block_size*(i/block_size)] + A[start_orig + i+block_size*((2*sq_size + i + block_size)/block_size)];  
                B[i + start + sq_size] = B[start_orig + i+block_size*(i/block_size)] + B[start_orig + i+block_size*((2*sq_size + i + block_size)/block_size)];

                A[i + start + 2*sq_size] = A[start_orig + i+block_size*(i/block_size)] - A[start_orig + i+block_size*((block_size + i)/block_size)];
                B[i + start + 2*sq_size] = B[start_orig + i+block_size*(i/block_size)] + B[start_orig + i+block_size*((2*sq_size + i)/block_size)];
           
                A[i + start + 3*sq_size] = A[start_orig + i+block_size*(i/block_size)] + A[start_orig + i+block_size*((2*sq_size + i)/block_size)];
                B[i + start + 3*sq_size] = B[start_orig + i+block_size*((2*sq_size + i + block_size)/block_size)]; //B_22 

                A[i + start + 4*sq_size] = A[start_orig + i+block_size*(i/block_size)]; //A_11         
                B[i + start + 4*sq_size] = B[start_orig + i+block_size*((2*sq_size + i)/block_size)] - B[start_orig + i+block_size*((2*sq_size + i + block_size)/block_size)];

                A[i + start + 5*sq_size] = A[start_orig + i+block_size*((2*sq_size + i + block_size)/block_size)]; //A_22 
                B[i + start + 5*sq_size] = B[start_orig + i+block_size*((block_size + i)/block_size)] - B[start_orig + i+block_size*(i/block_size)];

                A[i + start + 6*sq_size] = A[start_orig + i+block_size*((block_size + i)/block_size)] + A[start_orig + i+block_size*((2*sq_size + i + block_size)/block_size)];    
                B[i + start + 6*sq_size] = B[start_orig + i+block_size*(i/block_size)]; //B_11 
            }
            int seven_to_r = 1;
            for (i = 0; i < r; ++i) {
                seven_to_r *= 7;    
            }

            int orig_next_level = start_orig + (seven_to_r)*n*n;
            r += 1; // Recursion level
            //printf("sf 5 \n");
            for (i = 0; i < 7; ++i) {
                start = orig_next_level + (7 - i)*sq_size; //orig_next_level + (7^r - i)*sq_size;// + 7*i*(sq_size/4);
                strassen(A, B, C, block_size, min_size, start, orig_next_level,r);
            }

            printf("orig %i\n", start_orig);
            printf("orig nl %i\n", orig_next_level);
            printf("start %i\n", start);
            imp_matriz(10,10,10,C);
            //printf("sf 6 \n");
            for (i = 0; i < sq_size ; ++i) {
                C[start_orig + i+block_size*(i/block_size)] = C[orig_next_level + i] + C[orig_next_level + sq_size + i] - C[orig_next_level + 3*sq_size + i] + C[orig_next_level + 5*sq_size + i];//P_1[i] + P_2[i] - P_4[i] + P_6[i];
                C[start_orig + i+block_size*((2*sq_size + i)/block_size)] = C[orig_next_level + 3*sq_size + i] + C[orig_next_level + 4*sq_size + i]; //P_4[i] + P_5[i];
                C[start_orig + i+block_size*((block_size + i)/block_size)] = C[orig_next_level + 5*sq_size + i] + C[orig_next_level + 6*sq_size+i];//P_6[i] + P_7[i];
                C[start_orig + i+block_size*((2*sq_size + block_size + i)/block_size)] = C[orig_next_level + sq_size + i] - C[orig_next_level + 2*sq_size + i] + C[orig_next_level + 4*sq_size + i] - C[orig_next_level + 6*sq_size + i];//P_2[i] - P_3[i] + P_5[i] - P_7[i];
            }
            imp_matriz(10,10,10,C);
            //printf("sf 7 \n");
            return;
        }
    }

int main(int argc, char **argv)

{

    int n = atoi(argv[1]);//(1 << atoi(argv[1])); //Bit shift trick to find powers of 2
    int min_size = atoi(argv[2]);
    int sq_size = n*n;
    int size_c = n*n;
    int log_2_n = ceil(log2(n));
    int max_size = (1 << log_2_n);
    int alloc_size = max_size*max_size;

    for (int i = 0; i < log_2_n; ++i) {
        size_c *= 7;
    }
    //float A[sq_size];
    //float B[sq_size];
    //float C[sq_size];
    float alpha = 1.0;
    float beta = 0.0;
    clock_t time_1, end;
    double cpu_time_used_strassen, cpu_time_used_blas;
    float * A = (float *) malloc(alloc_size * sizeof(float));
    float * B = (float *) malloc(alloc_size * sizeof(float));
    float * C = (float *) malloc(4*n*n * sizeof(float));
    float * D = (float *) malloc(alloc_size * sizeof(float));
    
    int start = n*n;
    int start_orig = 0;
    int r = 0; // recursion level

    srand48(1) ;

    for (int i=0; i<sq_size;i++){
        A[i] = drand48();
        B[i] = drand48();
        C[i] = A[i];
        C[i+sq_size] = B[i];
    }
    //imp_matriz(n+1,n+1,n+1,A);
    //imp_matriz(n,n,n,B);
    //printf("n: %i \n", n);
    time_1 = clock();
    sgemm("N","N", &n, &n, &n, &alpha, A, &n, B, &n, &beta, D, &n);
    end = clock();

    cpu_time_used_blas = ((double)(end - time_1))/CLOCKS_PER_SEC;

    time_1 = clock();
    strassen(A, B, C, n, min_size, start, start_orig, r);
    end = clock();
    cpu_time_used_strassen = ((double)(end - time_1))/CLOCKS_PER_SEC;

    //If checking for residue
    float * R = (float *) malloc(alloc_size * sizeof(float));
    for (int i=0; i<sq_size;i++){
        R[i] = D[i] - C[i];
    }


//    printf(" %i, %g, %g \n", n, cpu_time_used_strassen, cpu_time_used_blas); //print for script
    //printf("TIMES: \n Time_strassen: %g \n Time_blas: %g \n",cpu_time_used_strassen, cpu_time_used_blas);
/*    printf("A:  \n\n");
    imp_matriz(n,n,n,A);
    printf("B:  \n\n");
    imp_matriz(n,n,n,B);
    printf("D:  \n\n");
    imp_matriz(n,n,n,D);
    printf("C: \n\n");
    imp_matriz(10,10,10,C);
    printf("R:  \n\n");
    imp_matriz(n,n,n,R);
    printf("\n\n");
*/
    /*free(A);
    free(B);
    free(C);
    free(D);
*/
}
