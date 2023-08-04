#!/bin/bash -el

RES_FILE=$(pwd -P)/test_results.json

if [ "$BOARD" == "alveo_u200" ]; then
	for EXEC_MODE in d p; do
	  for MATRIX_SIZE in 5120 10240; do
	    echo "=== Check mode: ${EXEC_MODE}, msize: ${MATRIX_SIZE} ==="
	    ##NOTE: Check == 2 -> enables the warm-up mode
	    CHECK=$([ "$MATRIX_SIZE" == "5120" ] && echo 1 || echo 2)
		NANOS6_CONFIG_OVERRIDE="devices.fpga.reverse_offload=true" \
            timeout --preserve-status 150s \
            taskset --cpu-list "0-4" \
            ./build/cholesky-${EXEC_MODE} ${MATRIX_SIZE} ${CHECK}
	    cat test_result.json >>$RES_FILE
	    echo "," >>$RES_FILE
	  done
	done
elif [ "$BOARD" == "zcu102" ]; then
	for EXEC_MODE in d p; do
	  for MATRIX_SIZE in 2048 4096; do
	    echo "=== Check mode: ${EXEC_MODE}, msize: ${MATRIX_SIZE} ==="
	    ##NOTE: Check == 2 -> enables the warm-up mode
	    CHECK=$([ "$MATRIX_SIZE" == "2048" ] && echo 1 || echo 2)
		NANOS6_CONFIG_OVERRIDE="devices.fpga.reverse_offload=true" \
            timeout --preserve-status 150s \
            taskset --cpu-list "0-4" \
            ./build/cholesky-${EXEC_MODE} ${MATRIX_SIZE} ${CHECK}
	    cat test_result.json >>$RES_FILE
	    echo "," >>$RES_FILE
	  done
	done
fi
