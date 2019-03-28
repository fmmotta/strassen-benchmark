#include <stdio.h>
#include <stdlib.h>
#include <time.h>
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

float strassen( float A, float B, int n, int min_size)
    {
        float C[n*n];
                
        if (n < min_size) {
            sgemm("N","N", &n, &n, &n, 1.0, A, &n, B, &n, 0, C, &n);
            return C;
        }
        else
        {
            int block_size = n/2;
            int sq_size = block_size*block_size;

            float A_11[sq_size];    
            float A_12[sq_size];    
            float A_21[sq_size];    
            float A_22[sq_size];    
            float B_11[sq_size];    
            float B_12[sq_size];    
            float B_21[sq_size];    
            float B_22[sq_size];    

            for (i = 0; i < sq_size ; ++i) {
                A_11[i] = A[i+block_size*(i/block_size)];        
                A_12[i] = A[i+block_size*((block_size + i)/block_size)];        
                A_21[i] = A[i+block_size*((4*block_size + i)/block_size)];        
                A_22[i] = A[i+block_size*((5*block_size + i)/block_size)];    

                B_11[i] = B[i+block_size*(i/block_size)];        
                B_12[i] = B[i+block_size*((block_size + i)/block_size)];        
                B_21[i] = B[i+block_size*((4*block_size + i)/block_size)];        
                B_22[i] = B[i+block_size*((5*block_size + i)/block_size)];        
            }

            float P_1[sq_size];
            float P_2[sq_size];
            float P_3[sq_size];
            float P_4[sq_size];
            float P_5[sq_size];
            float P_6[sq_size];
            float P_7[sq_size];

            float P_11[sq_size]; // A_12 - A22
            float P_12[sq_size]; // B_21 + B22

            float P_21[sq_size]; // A_11 + A22
            float P_22[sq_size]; // B_11 + B22

            float P_31[sq_size]; // A_11 - A21
            float P_32[sq_size]; // B_11 - B12

            float P_41[sq_size]; // A_11 + A12

            float P_52[sq_size]; // B_12 - B22

            float P_62[sq_size]; // B_21 - B11

            float P_71[sq_size]; // A_21 + A22

            for (i = 0; i < sq_size - 1; ++i) {
                P_11[i] = A_12[i] - A_22[i];        
                P_12[i] = B_21[i] + B_22[i];        

                P_21[i] = A_11[i] + A_22[i];        
                P_22[i] = B_11[i] + B_22[i];        

                P_31[i] = A_11[i] - A_21[i];        
                P_32[i] = B_11[i] + B_12[i];        

                P_41[i] = A_11[i] + A_12[i];        
                P_52[i] = B_12[i] - B_22[i];        
                P_62[i] = B_21[i] - B_11[i];        
                P_71[i] = A_21[i] + A_22[i];        
            }

            P_1 = strassen(P_11, P_12, block_size, min_size);
            P_2 = strassen(P_21, P_22, block_size, min_size);
            P_3 = strassen(P_31, P_32, block_size, min_size);
            P_4 = strassen(P_41, B_22, block_size, min_size);
            P_5 = strassen(A_11, P_52, block_size, min_size);
            P_6 = strassen(A_22, P_62, block_size, min_size);
            P_7 = strassen(P_71, B_11, block_size, min_size);


            for (i = 0; i < sq_size ; ++i) {
                C[i+block_size*(i/block_size)] = P_1[i] + P_2[i] - P_4[i] + P_6[i];
                C[i+block_size*((block_size + i)/block_size)] = P_4[i] + P_5[i];
                C[i+block_size*((4*block_size + i)/block_size)] = P_6[i] + P_7[i];
                C[i+block_size*((5*block_size + i)/block_size)] = P_2[i] - P_3[i] + P_5[i] - P_7[i];
            }

            return C;
        }

main(int argc, char **argv)

{

    int n = (1 << argv[1]); //Bit shift trick to find powers of 2
    int min_size[argv[2]];
    int sq_size = n*n;
    float A[sq_size];
    float B[sq_size];
    clock_t start, end;
    double cpu_time_used_strassen, cpu_time_used_blas;

    srand48(1) ;

    for (i=0; i<sq_size;i++)
        A[i] = drand48();
        B[i] = drand48();

    start = clock();
    strassen(&A, &B, n, min_size);
    end = clock();
    cpu_time_used_strassen = ((double)(end - start))/CLOCKS_PER_SEC;

    start = clock();
    sgemm("N","N", &n, &n, &n, 1.0, A, &n, B, &n, 0, C, &n);
    end = clock();
    cpu_time_used_blas = ((double)(end - start))/CLOCKS_PER_SEC;

}
