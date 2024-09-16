# Cholesky

**Name**: Cholesky Factorization Kernel  
**Contact Person**: OmpSs@FPGA Team, ompss-fpga-support@bsc.es  
**License Agreement**: GPL  
**Platform**: OmpSs@FPGA


### Description
This application performs a cholesky decomposition/factorization over a square matrix.
The matrix is distributed by blocks of contiguous memory.

The task implementation requires support for one external library that implements the mathematical operations. Supported ones are:
 - [Intel MKL](https://software.intel.com/en-us/mkl)
 - [Open BLAS](http://www.openblas.net/)

### Build instructions
Clone the repository:
```
git clone https://pm.bsc.es/gitlab/ompss-at-fpga/benchmarks/cholesky.git
cd cholesky
```

Build the application binaries:
```
make BOARD=zedboard CROSS_COMPILE=arm-linux-gnueabihf-
```
##### Build variables
You can change the build process defining or modifying some environment variables.
The supported ones are:
  - `CFLAGS`
    - `-DUSE_DMA_MEM`. Defining the `USE_DMA_MEM` variable the blocked matrix is allocated in kernel memory instead of user-space memory.
    - `-DUSE_DOUBLE`. Defining the `USE_DOUBLE` variable the matix elements are of type `double` instead of `float`.
    - `-DVERBOSE`. Defining the `VERBOSE` variable the application steps are shown meanwhile executed.
  - `LDFLAGS`
  - `GCC`. If not defined, the default value is: `gcc`.
  - `MCC`. If not defined, the default value is: `mcc`. However, for SMP machines we recommend the use of `smpcc`.
  - `CROSS_COMPILE`
  - `MKL_DIR`. Installation directory of MKL library. The default value is: `$MKLROOT`.
    - `MKL_INC_DIR`. Installation directory of includes for MKL library. The default value is: `$MKL_DIR/include`.
    - `MKL_LIB_DIR`. Installation directory of OS libraries for MKL library. The default value is: `$MKL_DIR/lib`.
  - `OPENBLAS_DIR`. Installation directory of OpenBLAS library. The default value is: `$OPENBLAS_HOME`.
    - `OPENBLAS_INC_DIR`. Installation directory of includes for OpenBLAS library. The default value is: `$OPENBLAS_DIR/include`.
    - `OPENBLAS_LIB_DIR`. Installation directory of OS libraries for OpenBLAS library. The default value is: `$OPENBLAS_DIR/lib`.
  - `BOARD`. Board option used when generating the bitstreams.
  - `FPGA_CLOCK`. Target frequency of FPGA accelerators in the bitstreams. The default value is: `200`.
  - `FPGA_MEMORY_PORT_WIDTH`. Bit-width of accelerators memory port to access main memory. The default value is: `128`.
  - `BLOCK_SIZE`. Dimension of matrix blocks that FPGA accelerators deal with. The default value is: `32`.
  - `SYRK_NUM_ACCS`. Number of FPGA accelerators for syrk task. The default value is: `1`.
  - `GEMM_NUM_ACCS`. Number of FPGA accelerators for gemm task. The default value is: `1`.
  - `TRSM_NUM_ACCS`. Number of FPGA accelerators for trsm task. The default value is: `1`.
  - `FPGA_GEMM_II`. Initiation interval, in cycles, for gemm middle loop. The default value is: `1`.
  - `FPGA_OTHER_II`. Initiation interval, in cycles, for syrk and trsm middle loop. The default value is: `1`.
  - `POTRF_SMP`. Use SMP arch for potrf tasks. The default value is: `1`.

Note that in order to compile the application either `MKL_DIR` or `OPENBLAS_DIR` (or the derivate variables) must point to a valid installation.

For example, the build step to cross-compile the application for ARM using the `smpcc` profile and OpenBLAS library may be:
```
export MCC=smpcc
export CROSS_COMPILE=arm-linux-gnueabihf-
export OPENBLAS_DIR=/opt/install-arm/openblas
make
```


### Run instructions
The name of each binary file created by build step ends with a suffix which determines the version:
 - program-p: performance version
 - program-i: instrumented version
 - program-d: debug version

All versions use the same arguments structure:
```
./cholesky <matrix size> [<check>]
```
where:
 - `matrix size` is the dimension of the matrices. (Mandatory)
 - `check` defines if the result must be checked. Default is: TRUE. (Optional)
