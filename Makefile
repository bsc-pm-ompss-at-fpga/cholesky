.PHONY: clean info
all: help

PROGRAM_     = cholesky

MCC         ?= fpgacc
MCC_         = $(CROSS_COMPILE)$(MCC)
GCC_         = $(CROSS_COMPILE)gcc
CFLAGS_      = $(CFLAGS) -O3 -std=gnu99 -Wall -Wno-unused -Wno-unknown-pragmas
MCC_FLAGS_   = $(MCC_FLAGS) --ompss -DRUNTIME_MODE=\"perf\"
MCC_FLAGS_I_ = $(MCC_FLAGS_) --instrument -DRUNTIME_MODE=\"instr\"
MCC_FLAGS_D_ = $(MCC_FLAGS_) --debug -g -k -DRUNTIME_MODE=\"debug\"
LDFLAGS_     = $(LDFLAGS) -lm

# FPGA bitstream Variables
FPGA_HWRUNTIME         ?= pom
FPGA_CLOCK             ?= 200
FPGA_MEMORY_PORT_WIDTH ?= 128
INTERCONNECT_OPT       ?= performance
SYRK_NUM_ACCS          ?= 1
GEMM_NUM_ACCS          ?= 1
TRSM_NUM_ACCS          ?= 1
BLOCK_SIZE             ?= 32
POTRF_SMP              ?= 1
FPGA_GEMM_II           ?= 1
FPGA_OTHER_II          ?= 1

ifeq ($(POTRF_SMP),1)
	CFLAGS_ += -DPOTRF_SMP
endif

## Check if copies are needed
ifneq (,$(findstring USE_DMA_MEM,$(CFLAGS_)))
	# When using Kernel memory no copies are needed
	MCC_FLAGS_ += --no-copy-deps
else
	# Copy deps is the default
endif

## MKL Variables
MKL_DIR      ?= $(MKLROOT)
MKL_INC_DIR  ?= $(MKL_DIR)/include
MKL_LIB_DIR  ?= $(MKL_DIR)/lib
MKL_SUPPORT_ = $(if $(and $(wildcard $(MKL_INC_DIR)/mkl.h ), \
               $(wildcard $(MKL_LIB_DIR)/libmkl_sequential.so )),YES,NO)

## Open Blas Variables
OPENBLAS_DIR      ?= $(OPENBLAS_HOME)
OPENBLAS_INC_DIR  ?= $(OPENBLAS_DIR)/include
OPENBLAS_LIB_DIR  ?= $(OPENBLAS_DIR)/lib
OPENBLAS_SUPPORT_ = $(if $(and $(wildcard $(OPENBLAS_INC_DIR)/lapacke.h ), \
                    $(wildcard $(OPENBLAS_LIB_DIR)/libopenblas.so )),YES,NO)

ifeq ($(MKL_SUPPORT_),YES)
	CFLAGS_  += -I$(MKL_INC_DIR) -DUSE_MKL
#	LDFLAGS_ += -L$(MKL_LIB_DIR) -lmkl_sequential -lmkl_core -lmkl_intel_lp64
	LDFLAGS_ += -L$(MKL_LIB_DIR) -lmkl_sequential -lmkl_core -lmkl_rt
else ifeq ($(OPENBLAS_SUPPORT_),YES)
	CFLAGS_  += -I$(OPENBLAS_INC_DIR) -DUSE_OPENBLAS
	LDFLAGS_ += -L$(OPENBLAS_LIB_DIR) -lopenblas -Wl,-rpath=$(OPENBLAS_LIB_DIR)
endif

CFLAGS_ += -DFPGA_OTHER_LOOP_II=$(FPGA_OTHER_II) -DFPGA_GEMM_LOOP_II=$(FPGA_GEMM_II) -DBLOCK_SIZE=$(BLOCK_SIZE) -DFPGA_MEMORY_PORT_WIDTH=$(FPGA_MEMORY_PORT_WIDTH) -DSYRK_NUM_ACCS=$(SYRK_NUM_ACCS) -DGEMM_NUM_ACCS=$(GEMM_NUM_ACCS) -DTRSM_NUM_ACCS=$(TRSM_NUM_ACCS) -DFPGA_HWRUNTIME=\"$(FPGA_HWRUNTIME)\"
FPGA_LINKER_FLAGS_ =--Wf,--name=$(PROGRAM_),--board=$(BOARD),-c=$(FPGA_CLOCK),--hwruntime=$(FPGA_HWRUNTIME)
ifdef FPGA_MEMORY_PORT_WIDTH
	MCC_FLAGS_ += --variable=fpga_memory_port_width:$(FPGA_MEMORY_PORT_WIDTH)
endif
ifdef MEMORY_INTERLEAVING_STRIDE
	FPGA_LINKER_FLAGS_ += --Wf,--memory_interleaving_stride=$(MEMORY_INTERLEAVING_STRIDE)
endif
ifdef SIMPLIFY_INTERCONNECTION
	FPGA_LINKER_FLAGS_ += --Wf,--simplify_interconnection
endif
ifdef INTERCONNECT_OPT
	FPGA_LINKER_FLAGS_ += --Wf,--interconnect_opt=$(INTERCONNECT_OPT)
endif
ifdef INTERCONNECT_REGSLICE
	FPGA_LINKER_FLAGS_ += --Wf,--interconnect_regslice,$(INTERCONNECT_REGSLICE)
endif
ifdef FLOORPLANNING_CONSTR
	FPGA_LINKER_FLAGS_ += --Wf,--floorplanning_constr=$(FLOORPLANNING_CONSTR)
endif
ifdef SLR_SLICES
	FPGA_LINKER_FLAGS_ += --Wf,--slr_slices=$(SLR_SLICES)
endif
ifdef PLACEMENT_FILE
	FPGA_LINKER_FLAGS_ += --Wf,--placement_file=$(PLACEMENT_FILE)
endif
ifeq ($(FPGA_HWRUNTIME),som)
	## Ignore the deps when spawning tasks inside the FPGA (only with SOM)
	FPGA_LINKER_FLAGS_ += --variable=fpga_ignore_deps_task_spawn:1
else ifeq ($(FPGA_HWRUNTIME),pom)
	FPGA_LINKER_FLAGS_ += --Wf,--picos_tm_size=256,--picos_dm_size=645,--picos_vm_size=775,--picos_max_args_per_task=3,--picos_max_deps_per_task=3,--picos_max_copies_per_task=3
endif
ifdef USE_URAM
	CFLAGS += -DUSE_URAM
endif

help:
	@echo 'Supported targets:       $(PROGRAM_)-p, $(PROGRAM_)-i, $(PROGRAM_)-d, $(PROGRAM_)-seq, design-p, design-i, design-d, bitstream-p, bitstream-i, bitstream-d, clean, help'
	@echo 'Environment variables:   CFLAGS, CROSS_COMPILE, LDFLAGS, MCC, MCC_FLAGS'
	@echo 'FPGA env. variables:     BOARD, FPGA_HWRUNTIME, FPGA_CLOCK, FPGA_MEMORY_PORT_WIDTH, SYRK_NUM_ACCS, GEMM_NUM_ACCS, TRSM_NUM_ACCS, BLOCK_SIZE, POTRF_SMP'
	@echo 'MKL env. variables:      MKLROOT, MKL_DIR, MKL_INC_DIR, MKL_LIB_DIR'
	@echo 'OpenBLAS env. variables: OPENBLAS_HOME, OPENBLAS_DIR, OPENBLAS_INC_DIR, OPENBLAS_LIB_DIR, OPENBLAS_IMPL'

$(PROGRAM_)-p: ./src/$(PROGRAM_)_$(FPGA_HWRUNTIME).c
	$(MCC_) $(CFLAGS_) $(MCC_FLAGS_) $^ -o $@ $(LDFLAGS_)

$(PROGRAM_)-i: ./src/$(PROGRAM_)_$(FPGA_HWRUNTIME).c
	$(MCC_) $(CFLAGS_) $(MCC_FLAGS_I_) $^ -o $@ $(LDFLAGS_)

$(PROGRAM_)-d: ./src/$(PROGRAM_)_$(FPGA_HWRUNTIME).c
	$(MCC_) $(CFLAGS_) $(MCC_FLAGS_D_) $^ -o $@ $(LDFLAGS_)

$(PROGRAM_)-seq: ./src/$(PROGRAM_)_$(FPGA_HWRUNTIME).c
	$(GCC_) $(CFLAGS_) -DRUNTIME_MODE=\"seq\" $^ -o $@ $(LDFLAGS_)

design-p: ./src/$(PROGRAM_)_$(FPGA_HWRUNTIME).c
	$(eval TMPFILE := $(shell mktemp))
	$(MCC_) $(CFLAGS_) $(MCC_FLAGS_) --bitstream-generation $(FPGA_LINKER_FLAGS_) \
		--Wf,--to_step=design \
		$^ -o $(TMPFILE) $(LDFLAGS_)
	rm $(TMPFILE)

design-i: ./src/$(PROGRAM_)_$(FPGA_HWRUNTIME).c
	$(eval TMPFILE := $(shell mktemp))
	$(MCC_) $(CFLAGS_) $(MCC_FLAGS_I_) --bitstream-generation $(FPGA_LINKER_FLAGS_) \
		--Wf,--to_step=design \
		$^ -o $(TMPFILE) $(LDFLAGS_)
	rm $(TMPFILE)

design-d: ./src/$(PROGRAM_)_$(FPGA_HWRUNTIME).c
	$(eval TMPFILE := $(shell mktemp))
	$(MCC_) $(CFLAGS_) $(MCC_FLAGS_D_) --bitstream-generation $(FPGA_LINKER_FLAGS_) \
		--Wf,--to_step=design,--debug_intfs=both \
		$^ -o $(TMPFILE) $(LDFLAGS_)
	rm $(TMPFILE)

bitstream-p: ./src/$(PROGRAM_)_$(FPGA_HWRUNTIME).c
	$(eval TMPFILE := $(shell mktemp))
	$(MCC_) $(CFLAGS_) $(MCC_FLAGS_) --bitstream-generation $(FPGA_LINKER_FLAGS_) \
		$^ -o $(TMPFILE) $(LDFLAGS_)
	rm $(TMPFILE)

bitstream-i: ./src/$(PROGRAM_)_$(FPGA_HWRUNTIME).c
	$(eval TMPFILE := $(shell mktemp))
	$(MCC_) $(CFLAGS_) $(MCC_FLAGS_I_) --bitstream-generation $(FPGA_LINKER_FLAGS_) \
		$^ -o $(TMPFILE) $(LDFLAGS_)
	rm $(TMPFILE)

bitstream-d: ./src/$(PROGRAM_)_$(FPGA_HWRUNTIME).c
	$(eval TMPFILE := $(shell mktemp))
	$(MCC_) $(CFLAGS_) $(MCC_FLAGS_D_) --bitstream-generation $(FPGA_LINKER_FLAGS_) \
		--Wf,--debug_intfs=both \
		$^ -o $(TMPFILE) $(LDFLAGS_)
	rm $(TMPFILE)

info:
	@echo "========== OPENBLAS =========="
	@echo "  SUPPORT enabled:  $(OPENBLAS_SUPPORT_)"
	@echo "  OPENBLAS_DIR:     $(OPENBLAS_DIR)"
	@echo "  OPENBLAS_INC_DIR: $(OPENBLAS_INC_DIR)"
	@echo "  OPENBLAS_LIB_DIR: $(OPENBLAS_LIB_DIR)"
	@echo "  Headers:          $(if $(wildcard $(OPENBLAS_INC_DIR)/lapacke.h ),YES,NO)"
	@echo "  Lib files (.so):  $(if $(wildcard $(OPENBLAS_LIB_DIR)/libopenblas.so ),YES,NO)"
	@echo "=============================="
	@echo "============= MKL ============"
	@echo "  SUPPORT enabled:  $(MKL_SUPPORT_)"
	@echo "  MKL_DIR:          $(MKL_DIR)"
	@echo "  MKL_INC_DIR:      $(MKL_INC_DIR)"
	@echo "  MKL_LIB_DIR:      $(MKL_LIB_DIR)"
	@echo "  Headers:          $(if $(wildcard $(MKL_INC_DIR)/mkl.h ),YES,NO)"
	@echo "  Lib files (.so):  $(if $(wildcard $(MKL_LIB_DIR)/libmkl_sequential.so ),YES,NO)"
	@echo "=============================="

clean:
	rm -fv *.o $(PROGRAM_)-? $(MCC_)_$(PROGRAM_)*.c *_ompss.cpp ait_$(PROGRAM_)*.json
	rm -fr $(PROGRAM_)_ait
