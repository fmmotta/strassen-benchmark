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

void strassen( float *A, float *B, float *C,int n, int min_size, int mem_tracker_a, int mem_tracker_b, int mem_tracker_c, int mem_tracker_aux1, int mem_tracker_aux2, int jump_size)
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
            int old_n = n;
            while (ceil(log2(n)) != floor(log2(n))){
                n += 1;
                for (i = (n*n - 1); i >= 0; --i){
                    if ((i > (n*n - n - 1)) || (i+1)%n == 0){
                         A[mem_tracker_a + i] = 0;
                         B[mem_tracker_b + i] = 0;
                    }
                    else {
                        A[mem_tracker_a + i] = A[mem_tracker_a + i - (i/n)];
                        B[mem_tracker_b + i] = B[mem_tracker_b + i - (i/n)];
                    }
                }            
            }
            //printf("sf 1 \n");
            int block_size = n/2;
            int sq_size = block_size*block_size;
            
            //float * P_i = (float *) malloc(sq_size * sizeof(float));
            //float * P_Aux1 = (float *) malloc(sq_size * sizeof(float));
            //float * P_Aux2 = (float *) malloc(sq_size * sizeof(float));
            
            printf("mem aux 1 0 %i\n", mem_tracker_aux1);
            printf("mem aux 2 0 %i\n", mem_tracker_aux2);
            mem_tracker_aux1 += 2*jump_size;
            mem_tracker_aux2 = mem_tracker_aux1 + sq_size;
            printf("mem aux 1 1 %i\n", mem_tracker_aux1);
            printf("mem aux 2 1 %i\n", mem_tracker_aux2);

            for (i = 0; i < sq_size; ++i) {
                C[mem_tracker_aux1 + i] = A[mem_tracker_a + i+block_size*((2*sq_size + i)/block_size)] - A[mem_tracker_a + i+block_size*((2*sq_size + i + block_size)/block_size)];
                C[mem_tracker_aux2 + i] = B[mem_tracker_b + i+block_size*((block_size + i)/block_size)] + B[mem_tracker_b + i+block_size*((2*sq_size + i + block_size)/block_size)]; 
            }
            int pos_c = mem_tracker_aux2 + sq_size;
            strassen(C, C, C, block_size, min_size, mem_tracker_aux1, mem_tracker_aux2, mem_tracker_aux2 + sq_size, mem_tracker_aux1, mem_tracker_aux2, sq_size);
            
            for (i = 0; i < sq_size ; ++i) {
                C[mem_tracker_c + i+block_size*(i/block_size)] = C[pos_c + i];
            }

            printf("mem aux 1 2 %i\n", mem_tracker_aux1);
            printf("mem aux 2 2 %i\n", mem_tracker_aux2);
            mem_tracker_aux1 += 2*sq_size;
            mem_tracker_aux2 = mem_tracker_aux1 + sq_size;
            printf("mem aux 1 3 %i\n", mem_tracker_aux1);
            printf("mem aux 2 3 %i\n", mem_tracker_aux2);

            for (i = 0; i < sq_size; ++i) {     
                C[mem_tracker_aux1 + i] = A[mem_tracker_a + i+block_size*(i/block_size)] + A[mem_tracker_a + i+block_size*((2*sq_size + i + block_size)/block_size)];  
                C[mem_tracker_aux2 + i] = B[mem_tracker_b + i+block_size*(i/block_size)] + B[mem_tracker_b + i+block_size*((2*sq_size + i + block_size)/block_size)];
            }
            pos_c = mem_tracker_aux2 + sq_size;
            //strassen(P_Aux1, P_Aux2, P_i, block_size, min_size, mem_tracker_a, mem_tracker_b, mem_tracker_c, mem_tracker_aux1, mem_tracker_aux2);
            strassen(C, C, C, block_size, min_size, mem_tracker_aux1, mem_tracker_aux2, mem_tracker_aux2 + sq_size, mem_tracker_aux1, mem_tracker_aux2, sq_size);
            
            for (i = 0; i < sq_size ; ++i) {
                C[mem_tracker_c + i+block_size*(i/block_size)] += C[pos_c + i];
                C[mem_tracker_c + i+block_size*((2*sq_size + block_size + i)/block_size)] = C[pos_c + i]; 
            }

            mem_tracker_aux1 += 2*sq_size;
            mem_tracker_aux2 = mem_tracker_aux1 + sq_size;

            for (i = 0; i < sq_size; ++i) {     
                C[mem_tracker_aux1 + i] = A[mem_tracker_a + i+block_size*(i/block_size)] - A[mem_tracker_a + i+block_size*((block_size + i)/block_size)];
                C[mem_tracker_aux2 + i] = B[mem_tracker_b + i+block_size*(i/block_size)] + B[mem_tracker_b + i+block_size*((2*sq_size + i)/block_size)];
            }
            pos_c = mem_tracker_aux2 + sq_size;
            //strassen(P_Aux1, P_Aux2, P_i, block_size, min_size, mem_tracker_a, mem_tracker_b, mem_tracker_c, mem_tracker_aux1, mem_tracker_aux2);
            strassen(C, C, C, block_size, min_size, mem_tracker_aux1, mem_tracker_aux2, mem_tracker_aux2 + sq_size, mem_tracker_aux1, mem_tracker_aux2, sq_size);
            
            for (i = 0; i < sq_size ; ++i) {
                C[mem_tracker_c + i+block_size*((2*sq_size + block_size + i)/block_size)] -=  C[pos_c + i];
            }

            mem_tracker_aux1 += 2*sq_size;
            mem_tracker_aux2 = mem_tracker_aux1 + sq_size;
            
            for (i = 0; i < sq_size; ++i) {     
                C[mem_tracker_aux1 + i] = A[mem_tracker_a + i+block_size*(i/block_size)] + A[mem_tracker_a + i+block_size*((2*sq_size + i)/block_size)];
                C[mem_tracker_aux2 + i] = B[mem_tracker_b + i+block_size*((2*sq_size + i + block_size)/block_size)]; //B_22 
            }
            //strassen(P_Aux1, P_Aux2, P_i, block_size, min_size, mem_tracker_a, mem_tracker_b, mem_tracker_c, mem_tracker_aux1, mem_tracker_aux2);
            pos_c = mem_tracker_aux2 + sq_size;
            strassen(C, C, C, block_size, min_size, mem_tracker_aux1, mem_tracker_aux2, mem_tracker_aux2 + sq_size, mem_tracker_aux1, mem_tracker_aux2, sq_size);

            for (i = 0; i < sq_size ; ++i) {
                C[mem_tracker_c + i+block_size*(i/block_size)] -= C[pos_c + i];
                C[mem_tracker_c + i+block_size*((2*sq_size + i)/block_size)] = C[pos_c + i];
            }

            mem_tracker_aux1 += 2*sq_size;
            mem_tracker_aux2 = mem_tracker_aux1 + sq_size;

            for (i = 0; i < sq_size; ++i) {     
                C[mem_tracker_aux1 + i] = A[mem_tracker_a + i+block_size*(i/block_size)]; //A_11         
                C[mem_tracker_aux2 + i] = B[mem_tracker_b + i+block_size*((2*sq_size + i)/block_size)] - B[mem_tracker_b + i+block_size*((2*sq_size + i + block_size)/block_size)];
            }
            //strassen(P_Aux1, P_Aux2, P_i, block_size, min_size, mem_tracker_a, mem_tracker_b, mem_tracker_c, mem_tracker_aux1, mem_tracker_aux2);
            pos_c = mem_tracker_aux2 + sq_size;
            strassen(C, C, C, block_size, min_size, mem_tracker_aux1, mem_tracker_aux2, mem_tracker_aux2 + sq_size, mem_tracker_aux1, mem_tracker_aux2, sq_size);

            for (i = 0; i < sq_size ; ++i) {
                C[mem_tracker_c + i+block_size*((2*sq_size + i)/block_size)] += C[pos_c + i];
                C[mem_tracker_c + i+block_size*((2*sq_size + block_size + i)/block_size)] += C[pos_c + i];
            }

            mem_tracker_aux1 += 2*sq_size;
            mem_tracker_aux2 = mem_tracker_aux1 + sq_size;

            for (i = 0; i < sq_size; ++i) {     
                C[mem_tracker_aux1 + i] = A[mem_tracker_a + i+block_size*((2*sq_size + i + block_size)/block_size)]; //A_22 
                C[mem_tracker_aux2 + i] = B[mem_tracker_b + i+block_size*((block_size + i)/block_size)] - B[mem_tracker_b + i+block_size*(i/block_size)];
            }
            //strassen(P_Aux1, P_Aux2, P_i, block_size, min_size, mem_tracker_a, mem_tracker_b, mem_tracker_c, mem_tracker_aux1, mem_tracker_aux2);
            pos_c = mem_tracker_aux2 + sq_size;
            strassen(C, C, C, block_size, min_size, mem_tracker_aux1, mem_tracker_aux2, mem_tracker_aux2 + sq_size, mem_tracker_aux1, mem_tracker_aux2, sq_size);

            for (i = 0; i < sq_size ; ++i) {
                C[mem_tracker_c + i+block_size*(i/block_size)] += C[pos_c + i];
                C[mem_tracker_c + i+block_size*((block_size + i)/block_size)] = C[pos_c + i];
            }

            mem_tracker_aux1 += 2*sq_size;
            mem_tracker_aux2 = mem_tracker_aux1 + sq_size;

            for (i = 0; i < sq_size; ++i) {     
                C[mem_tracker_aux1 + i] = A[mem_tracker_a + i+block_size*((block_size + i)/block_size)] + A[mem_tracker_a + i+block_size*((2*sq_size + i + block_size)/block_size)];
                C[mem_tracker_aux2 + i] = B[mem_tracker_b + i+block_size*(i/block_size)]; //B_11 
            }
            //strassen(P_Aux1, P_Aux2, P_i, block_size, min_size, mem_tracker_a, mem_tracker_b, mem_tracker_c, mem_tracker_aux1, mem_tracker_aux2);
            pos_c = mem_tracker_aux2 + sq_size;
            strassen(C, C, C, block_size, min_size, mem_tracker_aux1, mem_tracker_aux2, mem_tracker_aux2 + sq_size, mem_tracker_aux1, mem_tracker_aux2, sq_size);

            for (i = 0; i < sq_size ; ++i) {
                C[mem_tracker_c + i+block_size*((block_size + i)/block_size)] += C[pos_c + i];
                C[mem_tracker_c + i+block_size*((2*sq_size + block_size + i)/block_size)] -= C[pos_c + i];
            }

            if (old_n != n) {
                for (i = 0; i < old_n*old_n; ++i) {
                    C[mem_tracker_c + i] = C[mem_tracker_c + i+(n - old_n)*(i/old_n)];    
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
    int log_2_n = (int) ceil(log2(n));
    int max_size = (1 << log_2_n);//atoi(argv[3]);
    int size_c = n*n;
    int alloc_size = max_size*max_size;

    for (int i = 0; i < log_2_n; ++i) {
        size_c *= 7;
    }
    //float A[sq_size];
    //float B[sq_size];
    //float C[sq_size];
    float alpha = 1.0;
    float beta = 0.0;
    clock_t start, end;
    double cpu_time_used_strassen, cpu_time_used_blas;
    float * A = (float *) malloc(alloc_size * sizeof(float));
    float * B = (float *) malloc(alloc_size * sizeof(float));
    float * C = (float *) malloc(size_c * sizeof(float));
    float * D = (float *) malloc(alloc_size * sizeof(float));
    int mem_tracker_a = 0;
    int mem_tracker_b = 0;
    int mem_tracker_c = 0;
    int mem_tracker_aux1 = 0;
    int mem_tracker_aux2 = 0;

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
    strassen(A, B, C, n, min_size, mem_tracker_a, mem_tracker_b, mem_tracker_c, mem_tracker_aux1, mem_tracker_aux2, sq_size);
    end = clock();
    cpu_time_used_strassen = ((double)(end - start))/CLOCKS_PER_SEC;

    printf("mt1:  %i\n\n", mem_tracker_aux1);
    printf("mt2:  %i\n\n", mem_tracker_aux2);
    printf("ma:  %i\n\n", mem_tracker_a);
    printf("mb:  %i\n\n", mem_tracker_b);
    printf("mc:  %i\n\n", mem_tracker_c);
    //If checking for residue
    float * R = (float *) malloc(alloc_size * sizeof(float));
    for (int i=0; i<sq_size;i++){
        R[i] = D[i] - C[i];
    }


    //printf(" %i, %g, %g \n", n, cpu_time_used_strassen, cpu_time_used_blas); //print for script
    //printf(" %i, %g \n", min_size, cpu_time_used_strassen); //print for block size test
    //printf("TIMES: \n Time_strassen: %g \n Time_blas: %g \n",cpu_time_used_strassen, cpu_time_used_blas);
    printf("A:  \n\n");
    imp_matriz(n,n,n,A);
    printf("B:  \n\n");
    imp_matriz(n,n,n,B);
    printf("D:  \n\n");
    imp_matriz(n,n,n,D);
    printf("C: \n\n");
    imp_matriz(11,11,11,C);
   printf("R:  \n\n");
    imp_matriz(n,n,n,R);
    printf("\n\n");

    free(A);
    free(B);
    free(C);
    free(D);
}
