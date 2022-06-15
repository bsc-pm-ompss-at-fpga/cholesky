/*
* Copyright (c) 2020-2022, Barcelona Supercomputing Center
*                          Centro Nacional de Supercomputacion
*
* This program is free software: you can redistribute it and/or modify  
* it under the terms of the GNU General Public License as published by  
* the Free Software Foundation, version 3.
*
* This program is distributed in the hope that it will be useful, but 
* WITHOUT ANY WARRANTY; without even the implied warranty of 
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
* General Public License for more details.
*
* You should have received a copy of the GNU General Public License 
* along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <errno.h>
#include <assert.h>

#include "cholesky.h"
#include "cholesky.fpga.h"

#if defined(OPENBLAS_IMPL) || defined(POTRF_SMP)
#pragma omp target device(smp) copy_deps
#else
#pragma omp target device(fpga) copy_deps
#endif
#pragma omp task inout([ts*ts]A)
void omp_potrf(type_t *A)
{
#if defined(OPENBLAS_IMPL) || defined(POTRF_SMP)
   static const char L = 'L';
   int info;
   potrf(&L, &ts, A, &ts, &info);
#else
   #pragma HLS inline
   #pragma HLS array_partition variable=A cyclic factor=FPGA_PWIDTH/64
   for (int j = 0; j < ts; ++j) {
      type_t tmp = A[j*ts + j];
      for (int k = 0; k < j; ++k) {
         #pragma HLS pipeline II=1
         type_t Akj = A[k*ts + j];
         tmp -= Akj*Akj;
      }

      A[j*ts + j] = sqrtf(tmp);

      for (int i = j + 1; i < ts; ++i) {
         type_t tmp = A[j*ts + i];
         for (int k = 0; k < j; ++k) {
            #pragma HLS pipeline II=1
            tmp -= A[k*ts + i]*A[k*ts + j];
         }
         A[j*ts + i] = tmp/A[j*ts + j];
      }
   }
#endif
}

#ifdef OPENBLAS_IMPL
#pragma omp target device(smp) copy_deps
#else
#pragma omp target device(fpga) num_instances(TRSM_NUMACCS) copy_deps
#endif
#pragma omp task in([ts*ts]A) inout([ts*ts]B)
void omp_trsm(const type_t *A, type_t *B)
{
#ifdef OPENBLAS_IMPL
   trsm(CBLAS_MAT_ORDER, CBLAS_RI, CBLAS_LO, CBLAS_T, CBLAS_NU,
      ts, ts, 1.0, A, ts, B, ts);
#else
   #pragma HLS inline
   #pragma HLS array_partition variable=A cyclic factor=FPGA_PWIDTH/64
   #pragma HLS array_partition variable=B cyclic factor=ts/FPGA_OTHER_II
   #pragma HLS array_partition variable=tmp_row cyclic factor=ts/(2*FPGA_OTHER_II)
   type_t tmp_row[ts];

   for (int k = 0; k < ts; ++k) {
      type_t temp = 1. / A[k*ts + k];
      for (int i = 0; i < ts; ++i) {
         #pragma HLS unroll factor=ts/FPGA_OTHER_II
         #pragma HLS pipeline II=1
         //Sometimes Vivado HLS doesn't achieve II=1 because it detects
         //some false dependence on B, this fixes the issue. Same for the other loop
         #pragma HLS DEPENDENCE variable=B inter false
         B[k*ts + i] = tmp_row[i] = temp * B[k*ts + i];
      }

      for (int j = k + 1; j < ts; ++j) {
         #pragma HLS pipeline II=FPGA_OTHER_II
         #pragma HLS DEPENDENCE variable=B inter false
         for (int i = 0; i < ts; ++i) {
            B[j*ts + i] -= A[k*ts + j] * tmp_row[i];
         }
      }
   }
#endif
}

#ifdef OPENBLAS_IMPL
#pragma omp target device(smp) copy_deps
#else
#pragma omp target device(fpga) num_instances(SYRK_NUMACCS) copy_deps
#endif
#pragma omp task in([ts*ts]A) inout([ts*ts]B)
void omp_syrk(const type_t *A, type_t *B)
{
#ifdef OPENBLAS_IMPL
   syrk(CBLAS_MAT_ORDER, CBLAS_LO, CBLAS_NT,
      ts, ts, -1.0, A, ts, 1.0, B, ts);
#else
   #pragma HLS inline
   #pragma HLS array_partition variable=A cyclic factor=ts/FPGA_OTHER_II
   #pragma HLS array_partition variable=B cyclic factor=ts/FPGA_OTHER_II

   for (int k = 0; k < ts; ++k) {
      for (int i = 0; i < ts; ++i) {
         #pragma HLS pipeline II=FPGA_OTHER_II
         for (int j = 0; j < ts; ++j) {
            //NOTE: Instead of reduce the 'i' iterations, multiply by 0
            B[i*ts + j] += -A[k*ts + i] * (j < i ? 0 : A[k*ts + j]);
         }
      }
   }
#endif
}

#ifdef OPENBLAS_IMPL
#pragma omp target device(smp) copy_deps
#else
#pragma omp target device(fpga) num_instances(GEMM_NUMACCS) copy_deps
#endif
#pragma omp task in([ts*ts]A, [ts*ts]B) inout([ts*ts]C)
void omp_gemm(const type_t *A, const type_t *B, type_t *C)
{
#ifdef OPENBLAS_IMPL
   gemm(CBLAS_MAT_ORDER, CBLAS_NT, CBLAS_T,
      ts, ts, ts, -1.0, A, ts, B, ts, 1.0, C, ts);
#else
   #pragma HLS inline
   #pragma HLS array_partition variable=A cyclic factor=ts/(2*FPGA_GEMM_II)
   #pragma HLS array_partition variable=B cyclic factor=FPGA_PWIDTH/64
   #pragma HLS array_partition variable=C cyclic factor=ts/FPGA_GEMM_II
   #ifdef USE_URAM
   #pragma HLS resource variable=A core=XPM_MEMORY uram
   #pragma HLS resource variable=B core=XPM_MEMORY uram
   #endif

   for (int k = 0; k < ts; ++k) {
      for (int i = 0; i < ts; ++i) {
         #pragma HLS pipeline II=FPGA_GEMM_II
         for (int j = 0; j < ts; ++j) {
            C[i*ts + j] += A[k*ts + j] * -B[k*ts + i];
         }
      }
   }
#endif
}

#ifdef OPENBLAS_IMPL
#pragma omp target device(smp) copy_deps
#else
#pragma omp target device(fpga) copy_inout([nt*nt*ts*ts]A)
#endif
#pragma omp task
void cholesky_blocked(const int nt, type_t* A)
{
   for (int k = 0; k < nt; k++) {

      // Diagonal Block factorization
      omp_potrf( A + (k*nt + k)*ts*ts );

      // Triangular systems
      #ifdef OPENBLAS_IMPL
      for (int i = k+1; i < nt; i++) {
      #else
      // Create in inverse order because Picos wakes up ready
      // chain tasks in that order
      for (int i = nt-1; i >= k+1; i--) {
      #endif
         omp_trsm( A + (k*nt + k)*ts*ts,
                   A + (k*nt + i)*ts*ts );
      }

      // Update trailing matrix
      for (int i = k + 1; i < nt; i++) {
         for (int j = k + 1; j < i; j++) {
            omp_gemm( A + (k*nt + i)*ts*ts,
                      A + (k*nt + j)*ts*ts,
                      A + (j*nt + i)*ts*ts );
         }
         omp_syrk( A + (k*nt + i)*ts*ts,
                   A + (i*nt + i)*ts*ts );
      }
   }
   #pragma omp taskwait
}

// Robust Check the factorization of the matrix A2
static int check_factorization(int N, type_t *A1, type_t *A2, int LDA, char uplo)
{
#ifdef VERBOSE
   printf ("Checking result ...\n");
#endif

   char NORM = 'I', ALL = 'A', UP = 'U', LO = 'L', TR = 'T', NU = 'N', RI = 'R';
   type_t alpha = 1.0;
   type_t const b = 2.0;
#ifdef USE_DOUBLE
   const int t = 53;
#else
   const int t = 24;
#endif
   type_t const eps = pow_di( b, -t );

   type_t *Residual = (type_t *)malloc(N*N*sizeof(type_t));
   type_t *L1       = (type_t *)malloc(N*N*sizeof(type_t));
   type_t *L2       = (type_t *)malloc(N*N*sizeof(type_t));
   type_t *work     = (type_t *)malloc(N*sizeof(type_t));

   memset((void*)L1, 0, N*N*sizeof(type_t));
   memset((void*)L2, 0, N*N*sizeof(type_t));

   lacpy(&ALL, &N, &N, A1, &LDA, Residual, &N);

   /* Dealing with L'L or U'U  */
   if (uplo == 'U'){
      lacpy(&UP, &N, &N, A2, &LDA, L1, &N);
      lacpy(&UP, &N, &N, A2, &LDA, L2, &N);
      trmm(CBLAS_MAT_ORDER, CBLAS_LF, CBLAS_UP, CBLAS_T, CBLAS_NU,
         N, N, alpha, L1, N, L2, N);
   } else {
      lacpy(&LO, &N, &N, A2, &LDA, L1, &N);
      lacpy(&LO, &N, &N, A2, &LDA, L2, &N);
      trmm(CBLAS_MAT_ORDER, CBLAS_RI, CBLAS_LO, CBLAS_T, CBLAS_NU,
         N, N, alpha, L1, N, L2, N);
   }

   /* Compute the Residual || A -L'L|| */
   for (int i = 0; i < N; i++) {
      for (int j = 0; j < N; j++) {
         Residual[j*N+i] = L2[j*N+i] - Residual[j*N+i];
      }
   }

   type_t Rnorm = lange(&NORM, &N, &N, Residual, &N, work);
   type_t Anorm = lange(&NORM, &N, &N, A1, &N, work);

   printf("==================================================\n");
   printf("Checking the Cholesky Factorization \n");
#ifdef VERBOSE
   printf("-- Rnorm = %e \n", Rnorm);
   printf("-- Anorm = %e \n", Anorm);
   printf("-- Anorm*N*eps = %e \n", Anorm*N*eps);
   printf("-- ||L'L-A||_oo/(||A||_oo.N.eps) = %e \n",Rnorm/(Anorm*N*eps));
#endif

   const int info_factorization = isnan(Rnorm/(Anorm*N*eps)) ||
      isinf(Rnorm/(Anorm*N*eps)) || (Rnorm/(Anorm*N*eps) > 60.0);

   if ( info_factorization ){
      fprintf(stderr, "\n-- Factorization is suspicious ! \n\n");
   } else {
      printf("\n-- Factorization is CORRECT ! \n\n");
   }

   free(work);
   free(L2);
   free(L1);
   free(Residual);

   return info_factorization;
}

void initialize_matrix(const int n, type_t *matrix)
{
   int ISEED[4] = {0,0,0,1};
   int intONE=1;

#ifdef VERBOSE
   printf("Initializing matrix with random values ...\n");
#endif

   for (int i = 0; i < n*n; i+=n) {
      larnv(&intONE, &ISEED[0], &n, &matrix[i]);
   }

   type_t a = (type_t)n;
   for (int i=0; i<n; i++) {
      for (int j=0; j<n; j++) {
         matrix[j*n + i] = matrix[j*n + i] + matrix[i*n + j];
         matrix[i*n + j] = matrix[j*n + i];
      }
      //add_to_diag
      matrix[i*n + i] += a;
   }
}

int main(int argc, char* argv[])
{
   char *result[3] = {"n/a","sucessful","UNSUCCESSFUL"};

   if ( argc < 3 ) {
      fprintf( stderr, "USAGE:\t%s <matrix size> [<check>]\n", argv[0] );
      exit( -1 );
   }
   const int  n = atoi(argv[1]); // matrix size
   int check    = argc > 2 ? atoi(argv[2]) : 1; // check result?
   const int nt = n / ts; // number of tiles
   if ( n % ts != 0 ) {
      fprintf( stderr, "ERROR:\t<matrix size> is not multiple of <block size>\n" );
      exit( -1 );
   }

   // Allocate matrix
   type_t * const matrix = (type_t *) malloc(n * n * sizeof(type_t));
   assert(matrix != NULL);

   double tIniStart = wall_time();

   // Init matrix
   initialize_matrix(n, matrix);

   type_t * original_matrix = NULL;
   if ( check == 1 ) {
      // Allocate matrix
      original_matrix = (type_t *) malloc(n * n * sizeof(type_t));
      assert(original_matrix != NULL);
      memcpy(original_matrix, matrix, n * n * sizeof(type_t));
   }

   // Allocate blocked matrix
   type_t *Ah[nt][nt];
   type_t *Ab;
   const size_t s = ts * ts * sizeof(type_t);
   #ifdef USE_DMA_MEM
      Ab = (type_t *)nanos_fpga_malloc(s*nt*nt);
   #else
      Ab = malloc(s*nt*nt);
   #endif
   assert(Ab != NULL);

   for (int i = 0; i < nt; i++) {
      for (int j = 0; j < nt; j++) {
         Ah[i][j] = Ab + (i*nt + j)*ts*ts;
      }
   }

   const double tEndStart = wall_time();

#ifdef VERBOSE
   printf ("Executing ...\n");
#endif

   convert_to_blocks(nt, n, (type_t(*)[n]) matrix, Ah);

   const double tIniWarm = wall_time();

   //Warm up execution
   if (check == 2) {
       cholesky_blocked(nt, Ab);
       #pragma omp taskwait noflush
   }

   const double tEndWarm = wall_time();
   const double tIniExec = tEndWarm;

   //Performance execution
   cholesky_blocked(nt, Ab);

   #pragma omp taskwait noflush
   const double tEndExec = wall_time();
   const double tIniFlush = tEndExec;

   //The following TW will copy out the data moved to FPGA devices
   #pragma omp taskwait

   const double tEndFlush = wall_time();
   const double tIniToLinear = tEndFlush;

   convert_to_linear(nt, n, Ah, (type_t (*)[n]) matrix);

   const double tEndToLinear = wall_time();
   const double tIniCheck = tEndToLinear;

   if ( check == 1 ) {
      const char uplo = 'L';
      if ( check_factorization(n, original_matrix, matrix, n, uplo) ) check = 10;
      free(original_matrix);
   }

   const double tEndCheck = wall_time();

   // Print results
   float gflops = (float)n/1e3;
   gflops = (gflops*gflops*gflops/3.f)/(tEndExec - tIniExec);
   printf( "==================== RESULTS ===================== \n" );
   printf( "  Benchmark: %s (%s)\n", "Cholesky", "OmpSs" );
   printf( "  Elements type: %s\n", ELEM_T_STR );
#ifdef VERBOSE
   printf( "  Matrix size:           %dx%d\n", n, n);
   printf( "  Block size:            %dx%d\n", ts, ts);
#endif
   printf( "  Init. time (secs):     %f\n", tEndStart    - tIniStart );
   printf( "  Warm up time (secs):   %f\n", tEndWarm     - tIniWarm );
   printf( "  Execution time (secs): %f\n", tEndExec     - tIniExec );
   printf( "  Flush time (secs):     %f\n", tEndFlush    - tIniFlush );
   printf( "  Convert linear (secs): %f\n", tEndToLinear - tIniToLinear );
   printf( "  Checking time (secs):  %f\n", tEndCheck    - tIniCheck );
   printf( "  Performance (GFLOPS):  %f\n", gflops );
   printf( "================================================== \n" );

   // Free blocked matrix
   #ifdef USE_DMA_MEM
      nanos_fpga_free(Ab);
   #else
      free(Ab);
   #endif

   // Free matrix
   free(matrix);

   return check == 10 ? 1 : 0;
}
