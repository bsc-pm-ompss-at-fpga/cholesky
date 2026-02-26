/*
* Copyright (c) 2020, BSC (Barcelona Supercomputing Center)
* All rights reserved.
*
* Redistribution and use in source and binary forms, with or without
* modification, are permitted provided that the following conditions are met:
*     * Redistributions of source code must retain the above copyright
*       notice, this list of conditions and the following disclaimer.
*     * Redistributions in binary form must reproduce the above copyright
*       notice, this list of conditions and the following disclaimer in the
*       documentation and/or other materials provided with the distribution.
*     * Neither the name of the <organization> nor the
*       names of its contributors may be used to endorse or promote products
*       derived from this software without specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY BSC ''AS IS'' AND ANY
* EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
* WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
* DISCLAIMED. IN NO EVENT SHALL <copyright holder> BE LIABLE FOR ANY
* DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
* (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
* LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
* ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
* (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
* SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <errno.h>
#include <assert.h>

//#include "nanos6.h"
#include "nanos6/distributed.h"
//Needed for nanos6_get_num_cpus
#include "nanos6/debug.h"

#include "cholesky.h"
#include "cholesky.fpga.h"

const unsigned int FPGA_GEMM_II = FPGA_GEMM_LOOP_II;
const unsigned int FPGA_OTHER_II = FPGA_OTHER_LOOP_II;
const int ts = BLOCK_SIZE; // tile size
const unsigned int FPGA_PWIDTH = FPGA_MEMORY_PORT_WIDTH;
const unsigned int SYRK_NUMACCS = SYRK_NUM_ACCS;
const unsigned int GEMM_NUMACCS = GEMM_NUM_ACCS;
const unsigned int TRSM_NUMACCS = TRSM_NUM_ACCS;

static void gather_block(const int N, type_t *Alin, type_t *A)
{
   for (int i = 0; i < ts; i++) {
      for (int j = 0; j < ts; j++) {
         A[i*ts + j] = Alin[i*N + j];
      }
   }
}

static void scatter_block(const int N, const type_t *A, type_t *Alin)
{
   for (int i = 0; i < ts; i++) {
      for (int j = 0; j < ts; j++) {
         Alin[i*N + j] = A[i*ts + j];
      }
   }
}

static void scatter_transposed_block(const int N, const type_t *A, type_t *Alin)
{
   for (int i = 0; i < ts; i++) {
      for (int j = 0; j < ts; j++) {
         Alin[i*N + j] = A[j*ts + i];
      }
   }
}


static void convert_to_blocks(const int DIM, const int N, type_t (*Alin)[N], type_t *A[DIM][DIM])
{
   for (int i = 0; i < DIM; i++) {
      for (int j = 0; j < DIM; j++) {
         gather_block(N, &Alin[i*ts][j*ts], A[i][j]);
      }
   }
}

static void convert_to_linear(const int DIM, const int N, type_t *A[DIM][DIM], type_t (*Alin)[N])
{
   for (int i = 0; i < DIM; i++) {
      for (int j = 0; j < DIM; j++) {
         scatter_block(N, A[i][j], (type_t *) &Alin[i*ts][j*ts]);
      }
   }
}

static inline int get_block_idx(uint i, uint j, uint nt, unsigned char r, unsigned char s, const char cfidx, const uint ndiags, const unsigned char remstart)
{
#pragma HLS inline
    uint cf;
    uint globaln = i+j;
    uint n = globaln/s;
    if (cfidx == 0)
        cf = 0;
    else if (cfidx == 1)
        cf = (n+1)/2;
    else if (cfidx == 2)
        cf = (n+1)/4;
    else
        cf = (n+2)/4;
    uint rem = n - ndiags;
    uint remblocks = (rem+1)*remstart + s*((rem*(rem+1))/2) + (rem+1);
    uint nblocks = (s*n*(n+1) + 2*(n+1)*r)/4 + (n+1) - cf - (globaln >= nt ? remblocks : 0);
    uint block_idx = nblocks - ((r + n*s)/2 + 1) + (globaln >= nt ? remstart + rem*s + 1 : 0) + (globaln >= nt ? nt-1-i : j);
    return block_idx;
}

unsigned char calc_owner(int i, int j, unsigned char size) {
#pragma HLS inline
   return (i + j)%size;
}

static inline size_t get_total_blocks(int nt, int r, int s)
{
   int cf;
   int n = (2*nt-1)/s + (r < (2*nt-1)%s ? 1 : 0) - 1;
   int globaln = r + n*s;
   if (s%2 == 0 && r%2 == 0) {
      cf = 0;
   }
   else if (s%2 == 0 && r%2 == 1) {
      cf = (n+1)/2;
   }
   else if (s%2 == 1 && r%2 == 0) {
      cf = (n+1)/4;
   } else { //s%2 == 1 && r%2 == 1
      cf = (n+2)/4;
   }
   int rem = n - (nt/s + (r < nt%s ? 1 : 0));
   int remstart = r >= nt%s ? r-nt%s : r+s - nt%s;
   int remblocks = (rem+1)*remstart + s*((rem*(rem+1))/2) + (rem+1);
   int nblocks = (s*n*(n+1) + 2*(n+1)*r)/4 + (n+1) - cf - (globaln >= nt ? remblocks : 0);
   return nblocks;
}


#pragma oss task in([len]data)
void flushData(const type_t *data, int len) {
    //dummy task to pull data from fpga
}

#pragma oss task device(fpga) inout([ts*ts]A) copy_deps owner(taskOwner)
void omp_potrf(type_t *A, int taskOwner)
{
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
}

#pragma oss task device(fpga) num_instances(TRSM_NUMACCS) copy_deps in([ts*ts]A;ADowner) inout([ts*ts]B) owner(taskOwner)
void omp_trsm(const type_t *A, type_t *B, int ADowner, int taskOwner)
{
#ifdef TRSM_SMP
   trsm(CBLAS_MAT_ORDER, CBLAS_RI, CBLAS_LO, CBLAS_T, CBLAS_NU,
      ts, ts, 1.0, A, ts, B, ts);
#else
   #pragma HLS inline
   #pragma HLS array_partition variable=A cyclic factor=FPGA_PWIDTH/64
   #pragma HLS array_partition variable=B cyclic factor=ts/FPGA_OTHER_II
   type_t tmp_row[ts];
   #pragma HLS array_partition variable=tmp_row cyclic factor=ts/(2*FPGA_OTHER_II)

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

#pragma oss task device(fpga) num_instances(SYRK_NUMACCS) copy_deps in([ts*ts]A;Aowner) inout([ts*ts]B) owner(taskOwner)
void omp_syrk(const type_t *A, type_t *B, int Aowner, int taskOwner)
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

#pragma oss task device(fpga) num_instances(GEMM_NUMACCS) copy_deps in([ts*ts]A;Aowner) in([ts*ts]B;Bowner) inout([ts*ts]C) owner(taskOwner)
void omp_gemm(const type_t A[ts*ts], const type_t B[ts*ts], type_t C[ts*ts], int Aowner, int Bowner, int taskOwner)
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
   #if defined(__VITIS_HLS__)
      #pragma HLS bind_storage variable=A type=RAM_T2P impl=URAM
      #pragma HLS bind_storage variable=B type=RAM_T2P impl=URAM
   #else
      #pragma HLS resource variable=A core=XPM_MEMORY uram
      #pragma HLS resource variable=B core=XPM_MEMORY uram
   #endif
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

unsigned const int recv_buffer_size = 2;

//#ifdef GENERATE_BITSTREAM
//#pragma oss task device(fpga) inout([(nt*(nt+1)/2)*ts*ts]A) in([2*nt*ts*ts]recv_buffer)
//#else
//#pragma oss task device(broadcaster) inout([(nt*(nt+1)/2)*ts*ts]A) in([2*nt*ts*ts]recv_buffer)
//#endif
#pragma oss task device(fpga) inout([(nt*(nt+1)/2)*ts*ts]A) in([2*nt*ts*ts]recv_buffer) owner("all")
void cholesky_blocked(const int nt, type_t* A, type_t* recv_buffer)
{

    int r = OMPIF_Comm_rank();
    int s = OMPIF_Comm_size();
    //ap_uint<2> cfidx;
    //Select special cases based on rank and cluster sizes
    unsigned char cfidx;
    uint ndiags = nt/s + (r < nt%s ? 1 : 0);
    unsigned char remstart = r >= nt%s ? r-nt%s : r+s - nt%s;
    if (s%2 == 0 && r%2 == 0)
        cfidx = 0;
    else if (s%2 == 0 && r%2 == 1)
        cfidx = 1;
    else if (s%2 == 1 && r%2 == 0)
        cfidx = 2;
    else
        cfidx = 3;


   main_loop:
   for (int k = 0; k < nt; k++) {

      // Diagonal Block factorization
      const unsigned char task_owner = calc_owner(k, k, s);
      int A_idx;
      //Looks like this condition could be removed
      if (task_owner == r) {
         A_idx = get_block_idx(k, k, nt, r, s, cfidx, ndiags, remstart);
      }
      omp_potrf( A + A_idx*ts*ts, task_owner );

      // Triangular systems
      inner_loop:
      for (int i = k+1; i < nt; i++) {

         const unsigned char task_owner = calc_owner(i, k, s);
         const unsigned char data_owner0 = calc_owner(k, k, s);
         type_t *Ablock, *Bblock;

         char is_recv_buffer_A = (task_owner == r && data_owner0 != r);
         Ablock = is_recv_buffer_A ?
             recv_buffer + ((k%recv_buffer_size)*nt)*ts*ts:
             A + get_block_idx(k, k, nt, r, s, cfidx, ndiags, remstart)*ts*ts;

         Bblock = A + get_block_idx(i, k, nt, r, s, cfidx, ndiags, remstart)*ts*ts;

         //if (task_owner == r && data_owner0 != r)
         //   Ablock = recv_buffer + ((k%recv_buffer_size)*nt)*ts*ts;
         //else if (data_owner0 == r) {
         //   Ablock = A + get_block_idx(k, k, nt, r, s, cfidx, ndiags, remstart)*ts*ts;
         //}
         //if (task_owner == r) {
         //   Bblock = A + get_block_idx(i, k, nt, r, s, cfidx, ndiags, remstart)*ts*ts;
         //}

         omp_trsm( Ablock, Bblock, data_owner0, task_owner);
      }

      // Update trailing matrix
      inner_loop2:
      for (int i = k + 1; i < nt; i++) {
         gemm_loop:
         for (int j = k + 1; j < i; j++) {
            const unsigned char task_owner = calc_owner(i, j, s);
            //Data owner for first argument
            const unsigned char data_owner0 = calc_owner(i, k, s);
            //Data owner of second argument
            const unsigned char data_owner1 = calc_owner(j, k, s);
            type_t *Ablock, *Bblock, *Cblock;
            //Not owning data -> copy to scratchpad buffer

            char is_recv_buffer_A = (task_owner == r && data_owner0 != r);
            Ablock = is_recv_buffer_A ?
               recv_buffer + ((k%recv_buffer_size)*nt + i)*ts*ts:
               A + get_block_idx(i, k, nt, r, s, cfidx, ndiags, remstart)*ts*ts;
            //if (task_owner == r && data_owner0 != r)
            //   Ablock = recv_buffer + ((k%recv_buffer_size)*nt + i)*ts*ts;
            ////Owning the data -> use working matrix
            //else if (data_owner0 == r) {
            //   int A_idx = get_block_idx(i, k, nt, r, s, cfidx, ndiags, remstart);
            //   Ablock = A + A_idx*ts*ts;
            //}
            char is_recv_buffer_B = (task_owner == r && data_owner1 != r);
            Bblock = is_recv_buffer_B ?
               recv_buffer + ((k%recv_buffer_size)*nt + j)*ts*ts:
               A + get_block_idx(j, k, nt, r, s, cfidx, ndiags, remstart)*ts*ts;


            //if (task_owner == r && data_owner1 != r)
            //   Bblock = recv_buffer + ((k%recv_buffer_size)*nt + j)*ts*ts;
            //else if (data_owner1 == r) {
            //   uint B_idx = get_block_idx(j, k, nt, r, s, cfidx, ndiags, remstart);
            //   Bblock = A + B_idx*ts*ts;
            //}
            //Does not matter if not owning the task
            //if (task_owner == r) {
               uint C_idx = get_block_idx(i, j, nt, r, s, cfidx, ndiags, remstart);
               Cblock = A + C_idx*ts*ts;
            //}

            omp_gemm( Ablock, Bblock, Cblock, data_owner0, data_owner1, task_owner);
         }
         const unsigned char task_owner = calc_owner(i, i, s);
         const unsigned char data_owner0 = calc_owner(i, k, s);
         type_t *Ablock, *Bblock;
         char is_recv_buffer_A = (task_owner == r && data_owner0 != r);
         //This needs to be done as a conditiona assignment to work around a compiler issue
         Ablock = is_recv_buffer_A ?
            recv_buffer + ((k%recv_buffer_size)*nt + i)*ts*ts:
            A + get_block_idx(i, k, nt, r, s, cfidx, ndiags, remstart)*ts*ts;

         //if (task_owner == r && data_owner0 != r)
         //   Ablock = recv_buffer + ((k%recv_buffer_size)*nt + i)*ts*ts;
         //else if (data_owner0 == r) {
         //   int A_idx = get_block_idx(i, k, nt, r, s, cfidx, ndiags, remstart);
         //   Ablock = A + A_idx*ts*ts;
         //}
         //This will be always executed y the task owner
         //if (task_owner == r) {
            int B_idx = get_block_idx(i, i, nt, r, s, cfidx, ndiags, remstart);
            Bblock = A + B_idx*ts*ts;
         //}
         omp_syrk( Ablock, Bblock, data_owner0, task_owner );
      }
   }
   #pragma oss taskwait
}


static void triangular_to_linear(type_t* triangular, type_t* linear, int n, int nt)
{
   int diag_offset = 0;
   for (int d = 0; d < 2*nt-1; ++d) {
      int blocks = d/2 + 1 - (d >= nt ? d-nt + 1 : 0);
      int i = d >= nt ? nt-1 : d;
      int j = d >= nt ? d-nt+1 : 0;
      type_t* block_array;
      block_array = triangular + diag_offset*ts*ts;
      diag_offset += blocks;
      for (int b = 0; b < blocks; ++b) {
         const type_t* triangular_b = block_array + b*ts*ts;
         scatter_block(n, triangular_b, linear + j*ts*n + i*ts); //linear matrix is column major
         if (i != j) scatter_transposed_block(n, triangular_b, linear + i*ts*n + j*ts);
         --i, ++j;
      }
   }
}

void initialize_matrix_blocked_lower(const int n, type_t *matrix)
{
   // This code was designed for an MPI application, every rank would allocate and initialize its own
   // part of the matrix. In this version, the CPU allocates and initializes the whole matrix.
   // If the CPU doesn't have enough memory for the full matrix, a simple solution is to allocate just the space
   // for the FPGA with more blocks, and initialize the part corresponding to each FPGA by changing this
   // rank and size parameters.
   int rank = 0;
   int size = 1;

   //ISEED is INTEGER array, dimension (4)
   //On entry, the seed of the random number generator; the array
   //elements must be between 0 and 4095, and ISEED(4) must be odd.
   //On exit, the seed is updated.
   //int ISEED[4] = {0,0,0,1};
   const int intONE=1;

   const int nt = n/ts;

   const int tb = get_total_blocks(nt, rank, size);

   const long long msize = tb*ts*ts;
   int cpus = nanos6_get_num_cpus();
   const long long bsize = msize/cpus < 1024 ? 1024 : (msize/cpus > 2147483648ull ? 2147483648ull : msize/cpus);

   for (long long i = 0; i < msize; i += bsize) {
      #pragma oss task firstprivate(i)
      {
         int ISEED[4] = {0, i >> 32, i & 0xFFFFFFFF, 1};
         int final_size = msize-i > bsize ? bsize : msize-i; 
         larnv(&intONE, &ISEED[0], &final_size, matrix + i);
      }
   }
   #pragma oss taskwait

   type_t a = (type_t)n;
   int diag_offset = 0;
   for (int d = rank; d < 2*nt-1; d += size) {
      int blocks = d/2 + 1 - (d >= nt ? d-nt + 1 : 0);
      for (int b = 0; b < blocks; b++) {
         //#pragma oss task firstprivate(diag_offset, b, blocks, d)
         {
         type_t* bmat = matrix + (diag_offset + b)*ts*ts;
         for (int j = 0; j < ts; ++j) {
            for (int i = (d%2 == 0 && b == blocks-1 ? j : 0); i < ts; ++i) {
               bmat[j*ts + i] *= 2;
               if (d%2 == 0 && b == blocks-1) {//diagonal block
                  if (i == j)
                     bmat[j*ts + i] += a;
                  else
                     bmat[i*ts + j] = bmat[j*ts + i];
               }
            }
         }
         }
      }
      diag_offset += blocks;
   }
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

int main(int argc, char* argv[])
{
   char *result[3] = {"n/a","sucessful","UNSUCCESSFUL"};

   if ( argc < 3 ) {
      fprintf( stderr, "USAGE:\t%s <matrix size> [<check>]\n", argv[0] );
      fprintf( stderr, "Check valid values:\n" );
      fprintf( stderr, "\t0: Do not validate results\n" );
      fprintf( stderr, "\t1: Validate results\n" );
      fprintf( stderr, "\t2: Run an additional warm-up execution and validate results\n" );
      return 1;
   }
   const int  n = atoi(argv[1]); // matrix size
   int check    = argc > 2 ? atoi(argv[2]) : 1; // check result?
   const int nt = n / ts; // number of tiles
   if ( n % ts != 0 ) {
      fprintf( stderr, "ERROR:\t<matrix size> is not multiple of <block size>\n" );
      exit( -1 );
   }
   int nranks = nanos6_dist_num_devices();
   if (nranks <= 0) {
      fprintf(stderr, "No devices found");
   }

   // Allocate matrix
   type_t * matrix;
   // Allocate blocked matrix
   type_t *Ab;
   type_t *recv_buffer;
   const size_t s = ts * ts * sizeof(type_t);
   //const size_t tb = get_total_blocks(nt, 0, 1);
   const size_t tb = nt*(nt+1)/2;
   Ab = malloc(tb*s);
   recv_buffer = malloc(2*nt*s);
   if (Ab == NULL || recv_buffer == NULL) {
      fprintf(stderr, "Could not allocate matrix\n");
      return 1;
   }

   double tIniStart = wall_time();

   // Init matrix
   initialize_matrix_blocked_lower(n, Ab);

   type_t * original_matrix = NULL;
   if ( check == 1 ) {
      // Allocate matrix to check
      original_matrix = (type_t *) malloc(n * n * sizeof(type_t));
      if (original_matrix == NULL) {
         fprintf(stderr, "Could not allocate original matrix\n");
         return 1;
      }
      triangular_to_linear(Ab, original_matrix, n, nt);
      //print_matrix(original_matrix, n);
   }

   //compute max needed memory in the fpga
   unsigned int max_tb = get_total_blocks(nt, 0, nranks);
   for (int i = 1; i < nranks; ++i) {
	  unsigned int blocks = get_total_blocks(nt, i, nranks);
	  if (blocks > max_tb)
		 max_tb = blocks;
   }

//   //map buffers to remote nodes
   nanos6_dist_map_address(Ab, max_tb*s);
   nanos6_dist_map_address(recv_buffer, 2*nt*s);
   //copy data to the devices
   //I think nt copy descriptors should be enough
   nanos6_dist_memcpy_info_t *memcpy_infos =
      (nanos6_dist_memcpy_info_t*)malloc((2*nt-1)*sizeof(nanos6_dist_memcpy_info_t));

   const double tBeginCopy = wall_time();

   unsigned int diag_offset = 0;
   unsigned int *r_diag_offset = (unsigned int*)malloc(nranks*sizeof(unsigned int));
   memset(r_diag_offset, 0, nranks*sizeof(unsigned int));
   for (int d = 0; d < 2*nt-1; ++d) {
      int blocks = d/2 + 1 - (d >= nt ? d-nt + 1 : 0);
      int r = d%nranks;
      memcpy_infos[d].size = blocks*s;
      memcpy_infos[d].sendOffset = diag_offset*s;
      memcpy_infos[d].recvOffset = r_diag_offset[r]*s;
      memcpy_infos[d].devId = r;
      r_diag_offset[r] += blocks;
      diag_offset += blocks;
   }
   nanos6_dist_memcpy_vector(Ab, 2*nt-1, memcpy_infos, NANOS6_DIST_COPY_TO);
   const double tEndCopy = wall_time();
   const double tElapsedCopy = tEndCopy-tBeginCopy;
   free(r_diag_offset);
   free(memcpy_infos);


   const double tEndStart = wall_time();

#ifdef VERBOSE
   printf ("Executing ...\n");
#endif

   //we already have a matrix of blocks copy stuff to remote nodes

   const double tIniWarm = wall_time();

   //Warmup removed for distributed versions
   //TODO: remove warmup printing

   const double tEndWarm = wall_time();
   const double tIniExec = tEndWarm;

   //Performance execution
   cholesky_blocked(nt, Ab, recv_buffer);

   //Noflush is not currently supported, use nanos api calls instead
   //#pragma oss taskwait noflush([nt*nt*ts*ts]Ab)
   //Noflush is not needed as distributed copies are done manually
   //nanos6_set_noflush(Ab, (nt*(nt+1)/2)*ts*ts*sizeof(*Ab));
   #pragma oss taskwait
   const double tEndExec = wall_time();
   const double tIniFlush = tEndExec;

   //flushData(Ab, (nt*(nt+1)/2)*ts*ts);

   #pragma oss taskwait

   const double tEndFlush = wall_time();
   const double tIniToLinear = tEndFlush;


   const double tEndToLinear = wall_time();
   double tIniCheck;

   if ( check == 1 ) {
      //tIniToLinear = wall_time();
      unsigned int diag_offset = 0;
      unsigned int *r_diag_offset = (unsigned int*)malloc(nranks*sizeof(unsigned int));
      memset(r_diag_offset, 0, nranks*sizeof(unsigned int));
      for (unsigned int d = 0; d < 2*nt-1; ++d) {
         unsigned int blocks = d/2 + 1 - (d >= nt ? d-nt + 1 : 0);
         unsigned int r = d%nranks;
         nanos6_dist_memcpy_from_device(r, Ab, blocks*s, r_diag_offset[r]*s, diag_offset*s);
         r_diag_offset[r] += blocks;
         diag_offset += blocks;
      }
      free(r_diag_offset);
      matrix = (type_t*)malloc(n*n*sizeof(type_t));
      if (matrix == NULL) {
          fprintf(stderr, "Could not allocate auxiliar matrix\n");
          return 1;
      }
      triangular_to_linear(Ab, matrix, n, nt);
      //tEndToLinear = wall_time();
      tIniCheck = tEndToLinear;
      const char uplo = 'L';
      if ( check_factorization(n, original_matrix, matrix, n, uplo) ) check = 10;
      free(original_matrix);
      free(matrix);
      //tEndCheck = wall_time();
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
   free(Ab);

   return check == 10 ? 1 : 0;
}
