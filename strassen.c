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

void strassen( float *A, float *B, float *C,int n, int min_size)
    {
        float alpha = 1.0;
        float beta = 0.0;
        int i;
        int odd = 0;        
        if (n < min_size) {
            sgemm("N","N", &n, &n, &n, &alpha, A, &n, B, &n, &beta, C, &n);
            return;
        }
        else if (n == 1){
            C[0] = A[0]*B[0];
            return;
        }
        else
        {

            //imp_matriz(n,n,n,A);
            if (n%2){
                n += 1;
                odd = 1;
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
            //imp_matriz(n,n,n,A);

            //printf("sf 1 \n");
            int block_size = n/2;
            int sq_size = block_size*block_size;
            
            //printf("bs %i sq %i \n",block_size,sq_size);
            //printf("sf 3 \n");

            float * P_11 = (float *) malloc(sq_size * sizeof(float));
            float * P_12 = (float *) malloc(sq_size * sizeof(float));
            float * P_21 = (float *) malloc(sq_size * sizeof(float));
            float * P_22 = (float *) malloc(sq_size * sizeof(float));
            float * P_31 = (float *) malloc(sq_size * sizeof(float));
            float * P_32 = (float *) malloc(sq_size * sizeof(float));
            float * P_41 = (float *) malloc(sq_size * sizeof(float));
            float * P_42 = (float *) malloc(sq_size * sizeof(float));
            float * P_51 = (float *) malloc(sq_size * sizeof(float));
            float * P_52 = (float *) malloc(sq_size * sizeof(float));
            float * P_61 = (float *) malloc(sq_size * sizeof(float));
            float * P_62 = (float *) malloc(sq_size * sizeof(float));
            float * P_71 = (float *) malloc(sq_size * sizeof(float));
            float * P_72 = (float *) malloc(sq_size * sizeof(float));

            for (i = 0; i < sq_size; ++i) {
                P_11[i] = A[i+block_size*((2*sq_size + i)/block_size)] - A[i+block_size*((2*sq_size + i + block_size)/block_size)];
                P_12[i] = B[i+block_size*((block_size + i)/block_size)] + B[i+block_size*((2*sq_size + i + block_size)/block_size)]; 

                P_21[i] = A[i+block_size*(i/block_size)] + A[i+block_size*((2*sq_size + i + block_size)/block_size)];  
                P_22[i] = B[i+block_size*(i/block_size)] + B[i+block_size*((2*sq_size + i + block_size)/block_size)];

                P_31[i] = A[i+block_size*(i/block_size)] - A[i+block_size*((block_size + i)/block_size)];
                P_32[i] = B[i+block_size*(i/block_size)] + B[i+block_size*((2*sq_size + i)/block_size)];
            
                P_41[i] = A[i+block_size*(i/block_size)] + A[i+block_size*((2*sq_size + i)/block_size)];
                P_42[i] = B[i+block_size*((2*sq_size + i + block_size)/block_size)]; //B_22 

                P_51[i] = A[i+block_size*(i/block_size)]; //A_11         
                P_52[i] = B[i+block_size*((2*sq_size + i)/block_size)] - B[i+block_size*((2*sq_size + i + block_size)/block_size)];

                P_61[i] = A[i+block_size*((2*sq_size + i + block_size)/block_size)]; //A_22 
                P_62[i] = B[i+block_size*((block_size + i)/block_size)] - B[i+block_size*(i/block_size)];

                P_71[i] = A[i+block_size*((block_size + i)/block_size)] + A[i+block_size*((2*sq_size + i + block_size)/block_size)];    
                P_72[i] = B[i+block_size*(i/block_size)]; //B_11 
            }

            float * P_1 = (float *) malloc(sq_size * sizeof(float));
            float * P_2 = (float *) malloc(sq_size * sizeof(float));
            float * P_3 = (float *) malloc(sq_size * sizeof(float));
            float * P_4 = (float *) malloc(sq_size * sizeof(float));
            float * P_5 = (float *) malloc(sq_size * sizeof(float));
            float * P_6 = (float *) malloc(sq_size * sizeof(float));
            float * P_7 = (float *) malloc(sq_size * sizeof(float));

            //printf("sf 5 \n");
            strassen(P_11, P_12, P_1, block_size, min_size);
            strassen(P_21, P_22, P_2, block_size, min_size);
            strassen(P_31, P_32, P_3, block_size, min_size);
            strassen(P_41, P_42, P_4, block_size, min_size);
            strassen(P_51, P_52, P_5, block_size, min_size);
            strassen(P_61, P_62, P_6, block_size, min_size);
            strassen(P_71, P_72, P_7, block_size, min_size);

           

            //printf("sf 6 \n");
            for (i = 0; i < sq_size ; ++i) {
                C[i+block_size*(i/block_size)] = P_1[i] + P_2[i] - P_4[i] + P_6[i];
                C[i+block_size*((2*sq_size + i)/block_size)] = P_4[i] + P_5[i];
                C[i+block_size*((block_size + i)/block_size)] = P_6[i] + P_7[i];
                C[i+block_size*((2*sq_size + block_size + i)/block_size)] = P_2[i] - P_3[i] + P_5[i] - P_7[i];
            }
            //imp_matriz(n,n,n,C);
            if (odd){
                n -= 1;     
                for (i = 0; i < n*n; ++i) {
                    C[i] = C[i+(i/n)];    
                }
            }
            free(P_1);
/*            free(P_1);
            free(P_2);
            free(P_3);
            free(P_4);
            free(P_5);
            free(P_6);
            free(P_7);
            free(P_11);
            free(P_12);
            free(P_21);
            free(P_22);
            free(P_31);
            free(P_32);
            free(P_41);
            free(P_42);
            free(P_51);
            free(P_52);
            free(P_61);
            free(P_62);
            free(P_71);
            free(P_72);
            //printf("sf 7 \n");
*/            return;
        }
    }

int main(int argc, char **argv)

{

    int n = atoi(argv[1]);//(1 << atoi(argv[1])); //Bit shift trick to find powers of 2
    int min_size = atoi(argv[2]);
    int sq_size = n*n;
    int alloc_size = (n+1)*(n+1);
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

    srand48(1) ;

    for (int i=0; i<sq_size;i++){
        A[i] = drand48();
        B[i] = drand48();
    }
    //imp_matriz(n+1,n+1,n+1,A);
    //imp_matriz(n,n,n,B);
    //printf("n: %i \n", n);
    start = clock();
    sgemm("N","N", &n, &n, &n, &alpha, A, &n, B, &n, &beta, D, &n);
    end = clock();

    cpu_time_used_blas = ((double)(end - start))/CLOCKS_PER_SEC;

    start = clock();
    strassen(A, B, C, n, min_size);
    end = clock();
    cpu_time_used_strassen = ((double)(end - start))/CLOCKS_PER_SEC;

    //If checking for residue
/*    float * R = (float *) malloc(alloc_size * sizeof(float));
    for (int i=0; i<sq_size;i++){
        R[i] = D[i] - C[i];
    }
*/

    printf(" %i, %g, %g \n", n, cpu_time_used_strassen, cpu_time_used_blas); //print for script
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
}
