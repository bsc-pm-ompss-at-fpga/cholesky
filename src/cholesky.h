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

