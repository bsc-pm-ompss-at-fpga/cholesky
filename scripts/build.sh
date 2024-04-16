#!/bin/bash -e

BUILD_TARGET=$1

if [ "$BOARD" == "" ]; then
  echo "BOARD environment variable not defined"
  exit 1
elif [ "$FPGA_HWRUNTIME" == "" ]; then
  echo "FPGA_HWRUNTIME environment variable not defined"
  exit 1
elif [ "$FPGA_CLOCK" == "" ]; then
  echo "FPGA_CLOCK environment variable not defined. Using default: 100"
  export FPGA_CLOCK=${FPGA_CLOCK:-100}
fi

PROG_NAME=cholesky
OUT_DIR=$(pwd -P)/build
RES_FILE=$(pwd -P)/resources.json

if ! module load openblas ; then
  echo "CI_NODE ('$CI_NODE') not supported. Cannot load OpenBLAS."
  exit 1
fi

# Cleanup
make clean
mkdir -p $OUT_DIR

if [ "$BUILD_TARGET" == "binary" ]; then
  #Only build the binaries
  make ${PROG_NAME}-p ${PROG_NAME}-i ${PROG_NAME}-d
  mv ${PROG_NAME}-p ${PROG_NAME}-i ${PROG_NAME}-d $OUT_DIR
elif [ "$BUILD_TARGET" == "design" ]; then
  #Only generate the design
  make design-p design-i design-d

  #Remove OUT_DIR directory since we are not generating output products
  rm -rf $OUT_DIR
else
  make bitstream-p

  if [ -e ${PROG_NAME}_ait/${PROG_NAME}.bin ] ; then
    mv ${PROG_NAME}_ait/${PROG_NAME}.bin $OUT_DIR/bitstream.bin
  fi
  mv ${PROG_NAME}_ait/${PROG_NAME}.bit $OUT_DIR/bitstream.bit
  mv ${PROG_NAME}_ait/${PROG_NAME}.xtasks.config $OUT_DIR/xtasks.config

  printf "{\"benchmark\": \"${PROG_NAME}\", " >>$RES_FILE
  printf "\"toolchain\": \"ompss-2\", " >>$RES_FILE
  printf "\"hwruntime\": \"${FPGA_HWRUNTIME}\", " >>$RES_FILE
  printf "\"board\": \"${BOARD}\", " >>$RES_FILE
  printf "\"builder\": \"${CI_NODE}\", " >>$RES_FILE
  printf "\"version\": \"${SYRK_NUM_ACCS}syrk ${GEMM_NUM_ACCS}gemm ${TRSM_NUM_ACCS}trsm ${BLOCK_SIZE}BS memport_128 noflush\", " >>$RES_FILE
  printf "\"accels_freq\": \"${FPGA_CLOCK}\", " >>$RES_FILE
  printf "\"memory_port_width\": \"${FPGA_MEMORY_PORT_WIDTH}" >>$RES_FILE
  for PARAM in BRAM DSP FF LUT; do
    printf "\", \"${PARAM}_HLS\": \"" >>$RES_FILE
    grep "$PARAM" ${PROG_NAME}_ait/${PROG_NAME}.resources-hls.txt | awk '{printf $2}' >>$RES_FILE
  done
  if grep -q URAM ${PROG_NAME}_ait/${PROG_NAME}.resources-hls.txt; then
    printf "\", \"URAM_HLS\": \"" >>$RES_FILE
    grep "URAM" ${PROG_NAME}_ait/${PROG_NAME}.resources-hls.txt | awk '{printf $2}' >>$RES_FILE
  fi
  for PARAM in BRAM DSP FF LUT; do
    printf "\", \"${PARAM}_IMPL\": \"" >>$RES_FILE
    grep "${PARAM}" ${PROG_NAME}_ait/${PROG_NAME}.resources-impl.txt | awk '{printf $2}' >>$RES_FILE
  done
  if grep -q URAM ${PROG_NAME}_ait/${PROG_NAME}.resources-impl.txt; then
    printf "\", \"URAM_IMPL\": \"" >>$RES_FILE
    grep "URAM" ${PROG_NAME}_ait/${PROG_NAME}.resources-impl.txt | awk '{printf $2}' >>$RES_FILE
  fi
  for PARAM in WNS TNS NUM_ENDPOINTS NUM_FAIL_ENDPOINTS; do
    printf "\", \"${PARAM}\": \"" >>$RES_FILE
    grep "${PARAM}" ${PROG_NAME}_ait/${PROG_NAME}.timing-impl.txt | awk '{printf $2}' >>$RES_FILE
  done
  printf "\"}\n" >>$RES_FILE
fi
