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

#include <sys/time.h>
#include <time.h>
#include <math.h>
#if USE_MKL
# include <mkl_lapacke.h>
# include <mkl.h>
#elif USE_OPENBLAS
# include <lapacke.h>
# include <cblas.h>
#else
//NOTE: Cannot throw the error as Vivado HLS will not compile
//# error No backend library found. See README for more information
#endif

#ifndef BLOCK_SIZE
#  error BLOCK_SIZE not defined
#endif
#ifndef FPGA_MEMORY_PORT_WIDTH
#  error FPGA_MEMORY_PORT_WIDTH variable not defined
#endif
#ifndef SYRK_NUM_ACCS
#  error SYRK_NUM_ACCS variable not defined
#endif
#ifndef GEMM_NUM_ACCS
#  error GEMM_NUM_ACCS variable not defined
#endif
#ifndef TRSM_NUM_ACCS
#  error TRSM_NUM_ACCS variable not defined
#endif
#ifndef FPGA_OTHER_LOOP_II
#  error FPGA_OTHER_LOOP_II variable not defined
#endif
#ifndef FPGA_GEMM_LOOP_II
#  error FPGA_GEMM_LOOP_II variable not defined
#endif

#pragma omp target device(fpga)
const unsigned int FPGA_GEMM_II = FPGA_GEMM_LOOP_II;
#pragma omp target device(fpga)
const unsigned int FPGA_OTHER_II = FPGA_OTHER_LOOP_II;
const int ts = BLOCK_SIZE; // tile size
#pragma omp target device(fpga)
const unsigned int FPGA_PWIDTH = FPGA_MEMORY_PORT_WIDTH;
const unsigned int SYRK_NUMACCS = SYRK_NUM_ACCS;
const unsigned int GEMM_NUMACCS = GEMM_NUM_ACCS;
const unsigned int TRSM_NUMACCS = TRSM_NUM_ACCS;

#if defined(USE_DOUBLE)
#  define type_t     double
#  define ELEM_T_STR "double"
#  define gemm       cblas_dgemm
#  define trsm       cblas_dtrsm
#  define trmm       cblas_dtrmm
#  define syrk       cblas_dsyrk
#  if USE_MKL
#    define potrf    dpotrf
#    define lacpy    dlacpy
#    define lange    dlange
#    define larnv    dlarnv
#  else
#    define potrf    LAPACK_dpotrf
#    define lacpy    LAPACK_dlacpy
#    define lange    LAPACK_dlange
#    define larnv    LAPACK_dlarnv
#  endif
#else
#  define type_t     float
#  define ELEM_T_STR "float"
#  define gemm       cblas_sgemm
#  define trsm       cblas_strsm
#  define trmm       cblas_strmm
#  define syrk       cblas_ssyrk
#  if USE_MKL
#    define potrf    spotrf
#    define lacpy    slacpy
#    define lange    slange
#    define larnv    slarnv
#  else
#    define potrf    LAPACK_spotrf
#    define lacpy    LAPACK_slacpy
#    define lange    LAPACK_slange
#    define larnv    LAPACK_slarnv
#  endif
#endif
#define CBLAS_MAT_ORDER   CblasColMajor
#define CBLAS_T           CblasTrans
#define CBLAS_NT          CblasNoTrans
#define CBLAS_LO          CblasLower
#define CBLAS_UP          CblasUpper
#define CBLAS_LF          CblasLeft
#define CBLAS_RI          CblasRight
#define CBLAS_NU          CblasNonUnit

double wall_time () {
   struct timespec ts;
   clock_gettime(CLOCK_MONOTONIC,&ts);
   return (double) (ts.tv_sec) + (double) ts.tv_nsec * 1.0e-9;
}

static type_t pow_di(type_t x, int n)
{
   type_t rv = 1.0;

   if (n < 0) {
      n = -n;
      x = 1.0 / x;
   }

   for (; n; n >>= 1, x *= x) {
      if (n & 1) rv *= x;
   }

   return rv;
}

static void gather_block(const int N, type_t *Alin, type_t *A)
{
   for (int i = 0; i < ts; i++) {
      for (int j = 0; j < ts; j++) {
         A[i*ts + j] = Alin[i*N + j];
      }
   }
}

static void scatter_block(const int N, type_t *A, type_t *Alin)
{
   for (int i = 0; i < ts; i++) {
      for (int j = 0; j < ts; j++) {
         Alin[i*N + j] = A[i*ts + j];
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
