#!/bin/bash -el

RES_FILE=$(pwd -P)/test_results.json

MATRIX_SIZES=(${MATRIX_SIZES-1024 2048})

declare -A NANOS6_CONFIG_EXEC_MODE
NANOS6_CONFIG_EXEC_MODE['d']="devices.fpga.reverse_offload=true, version.debug=true"
NANOS6_CONFIG_EXEC_MODE['p']="devices.fpga.reverse_offload=true"

declare -A RUNTIME_MODE_EXEC_MODE
RUNTIME_MODE_EXEC_MODE['d']="debug"
RUNTIME_MODE_EXEC_MODE['p']="perf"

# Do not override NANOS6_CONFIG_OVERRIDE, as it might contain node-specific default values

for EXEC_MODE in d p; do
  for IDX in ${!MATRIX_SIZES[@]}; do
    MATRIX_SIZE=${MATRIX_SIZES[$IDX]}
    echo "=== Check mode: ${EXEC_MODE}, msize: ${MATRIX_SIZE} ==="
    #NOTE: Check == 1 -> Enable output result check
    #      Check == 2 -> Run an additional warm-up execution before the performance one
    CHECK=$([ "$IDX" == "0" ] && echo 1 || echo 2)
    NANOS6_CONFIG_OVERRIDE="$NANOS6_CONFIG_OVERRIDE,${NANOS6_CONFIG_EXEC_MODE[$EXEC_MODE]}" \
    RUNTIME_MODE=${RUNTIME_MODE_EXEC_MODE[$EXEC_MODE]} \
      taskset --cpu-list "0-4" \
      ./build/cholesky-${EXEC_MODE} ${MATRIX_SIZE} ${CHECK}
    cat test_result.json >>$RES_FILE
    echo "," >>$RES_FILE
  done
done
