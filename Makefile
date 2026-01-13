.PHONY: all help clean
.PHONY: design-p design-i design-d
.PHONY: bitstream-p bitstream-i bitstream-d
.PHONY: $(PROGRAM)-seq $(PROGRAM)-p $(PROGRAM)-i $(PROGRAM)-d
all: help

PROGRAM_ = cholesky

PROGRAM_SRC = \
	src/cholesky.c

# Parameters
FPGA_CLOCK             ?= 200
FPGA_MEMORY_PORT_WIDTH ?= 128
SYRK_NUM_ACCS          ?= 1
GEMM_NUM_ACCS          ?= 1
TRSM_NUM_ACCS          ?= 1
BLOCK_SIZE             ?= 32
POTRF_SMP              ?= 1
FPGA_GEMM_II           ?= 1
FPGA_OTHER_II          ?= 1

help:
	@echo 'Supported targets:           $(PROGRAM_)-p, $(PROGRAM_)-i, $(PROGRAM_)-d, $(PROGRAM_)-seq, design-p, design-i, design-d, bitstream-p, bitstream-i, bitstream-d, clean, help'
	@echo 'FPGA env. variables:         BOARD, FPGA_CLOCK, FPGA_MEMORY_PORT_WIDTH, MEMORY_INTERLEAVING_STRIDE, SIMPLIFY_INTERCONNECTION, INTERCONNECT_OPT, INTERCONNECT_REGSLICE, FLOORPLANNING_CONSTR, SLR_SLICES, PLACEMENT_FILE'
	@echo 'Benchmark env. variables:    SYRK_NUM_ACCS, GEMM_NUM_ACCS, TRSM_NUM_ACCS, BLOCK_SIZE, POTRF_SMP, FPGA_GEMM_II, FPGA_OTHER_II'
	@echo 'MKL env. variables:          MKLROOT, MKL_DIR, MKL_INC_DIR, MKL_LIB_DIR'
	@echo 'OpenBLAS env. variables:     OPENBLAS_HOME, OPENBLAS_DIR, OPENBLAS_INC_DIR, OPENBLAS_LIB_DIR, OPENBLAS_IMPL'
	@echo 'Compiler env. variables:     CFLAGS, CROSS_COMPILE, LDFLAGS'

CLANG_TARGET =
ifdef CROSS_COMPILE
	CLANG_TARGET += -target $(CROSS_COMPILE)
endif

COMPILER_         = clang
COMPILER_FLAGS_   = $(CFLAGS) $(CLANG_TARGET) -fompss-2 -fompss-fpga-wrapper-code
COMPILER_FLAGS_I_ = -fompss-fpga-instrumentation
COMPILER_FLAGS_D_ = -g -fompss-fpga-hls-tasks-dir $(PWD)
LINKER_FLAGS_     = $(LDFLAGS)

AIT_FLAGS__        = --name=$(PROGRAM_) --board=$(BOARD) -c=$(FPGA_CLOCK)
AIT_FLAGS_DESIGN__ = --to_step=design
AIT_FLAGS_D__      = --debug_intfs=both -k -i -v

# Optional FPGA optimization variables
ifdef FPGA_MEMORY_PORT_WIDTH
	COMPILER_FLAGS_ += -fompss-fpga-memory-port-width $(FPGA_MEMORY_PORT_WIDTH)
endif
ifdef USER_CONFIG
	AIT_FLAGS__ += --user_config=config_files/$(USER_CONFIG).json
endif
ifdef MEMORY_INTERLEAVING_STRIDE
	AIT_FLAGS__ += --memory_interleaving_stride=$(MEMORY_INTERLEAVING_STRIDE)
endif
ifdef INTERCONNECT_PRIORITIES
	AIT_FLAGS__ += --interconnect_priorities
endif
ifdef INTERCONNECT_REGSLICES
	AIT_FLAGS__ += --interconnect_regslices
endif
ifdef DISABLE_STATIC_CONSTRAINTS
	AIT_FLAGS__ += --disable_static_constraints
endif
ifdef AIT_EXTRA_FLAGS
	AIT_FLAGS__ += $(AIT_EXTRA_FLAGS)
endif

## Cholesky specifics
# Picos configuration
AIT_FLAGS__ += --max_deps_per_task=3 --max_args_per_task=3 --max_copies_per_task=3 --picos_tm_size=256 --picos_dm_size=645 --picos_vm_size=775

# Preprocessor flags
COMPILER_FLAGS_ += -DBOARD=\"$(BOARD)\" -DFPGA_MEMORY_PORT_WIDTH=$(FPGA_MEMORY_PORT_WIDTH) -DFPGA_CLOCK=$(FPGA_CLOCK)
COMPILER_FLAGS_ += -DFPGA_OTHER_LOOP_II=$(FPGA_OTHER_II) -DFPGA_GEMM_LOOP_II=$(FPGA_GEMM_II) -DBLOCK_SIZE=$(BLOCK_SIZE) -DFPGA_MEMORY_PORT_WIDTH=$(FPGA_MEMORY_PORT_WIDTH) -DSYRK_NUM_ACCS=$(SYRK_NUM_ACCS) -DGEMM_NUM_ACCS=$(GEMM_NUM_ACCS) -DTRSM_NUM_ACCS=$(TRSM_NUM_ACCS)

ifdef USE_URAM
	COMPILER_FLAGS_ += -DUSE_URAM
endif

ifeq ($(POTRF_SMP),1)
	COMPILER_FLAGS_ += -DPOTRF_SMP
endif

# Linker flags
LINKER_FLAGS_ += -lm

## MKL variables
MKL_DIR      ?= $(MKLROOT)
MKL_INC_DIR  ?= $(MKL_DIR)/include
MKL_LIB_DIR  ?= $(MKL_DIR)/lib
MKL_SUPPORT_ = $(if $(and $(wildcard $(MKL_INC_DIR)/mkl.h ), \
               $(wildcard $(MKL_LIB_DIR)/libmkl_sequential.so )),YES,NO)

## Open Blas variables
OPENBLAS_DIR      ?= $(OPENBLAS_HOME)
OPENBLAS_INC_DIR  ?= $(OPENBLAS_DIR)/include
OPENBLAS_LIB_DIR  ?= $(OPENBLAS_DIR)/lib
OPENBLAS_SUPPORT_ = $(if $(and $(wildcard $(OPENBLAS_INC_DIR)/lapacke.h ), \
                    $(wildcard $(OPENBLAS_LIB_DIR)/libopenblas.so )),YES,NO)

ifeq ($(MKL_SUPPORT_),YES)
	COMPILER_FLAGS_  += -I$(MKL_INC_DIR) -DUSE_MKL
	LINKER_FLAGS_    += -L$(MKL_LIB_DIR) -lmkl_sequential -lmkl_core -lmkl_rt
else ifeq ($(OPENBLAS_SUPPORT_),YES)
	COMPILER_FLAGS_  += -I$(OPENBLAS_INC_DIR) -DUSE_OPENBLAS
	LINKER_FLAGS_    += -L$(OPENBLAS_LIB_DIR) -lopenblas -Wl,-rpath=$(OPENBLAS_LIB_DIR)
endif

AIT_FLAGS_        = -fompss-fpga-ait-flags "$(AIT_FLAGS__)"
AIT_FLAGS_DESIGN_ = -fompss-fpga-ait-flags "$(AIT_FLAGS_DESIGN__)"
AIT_FLAGS_D_      = -fompss-fpga-ait-flags "$(AIT_FLAGS_D__)"

.PHONY: info
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

$(PROGRAM_)-p: $(PROGRAM_SRC)
	$(COMPILER_) $(COMPILER_FLAGS_) $^ -o $@ $(LINKER_FLAGS_)

$(PROGRAM_)-i: $(PROGRAM_SRC)
	$(COMPILER_) $(COMPILER_FLAGS_) $(COMPILER_FLAGS_I_) $^ -o $@ $(LINKER_FLAGS_)

$(PROGRAM_)-d: $(PROGRAM_SRC)
	$(COMPILER_) $(COMPILER_FLAGS_) $(COMPILER_FLAGS_D_) $^ -o $@ $(LINKER_FLAGS_)

$(PROGRAM_)-seq: $(PROGRAM_SRC)
	$(COMPILER_) $(COMPILER_FLAGS_) $(COMPILER_FLAGS_SEQ_) $^ -o $@ $(LINKER_FLAGS_)

design-p: $(PROGRAM_SRC)
	$(eval TMPFILE := $(shell mktemp))
	$(COMPILER_) $(COMPILER_FLAGS_) \
		$(AIT_FLAGS_) $(AIT_FLAGS_DESIGN_) \
		$^ -o $(TMPFILE) $(LINKER_FLAGS_)
	rm $(TMPFILE)

design-i: $(PROGRAM_SRC)
	$(eval TMPFILE := $(shell mktemp))
	$(COMPILER_) $(COMPILER_FLAGS_) $(COMPILER_FLAGS_I_) \
		$(AIT_FLAGS_) $(AIT_FLAGS_DESIGN_) \
		$^ -o $(TMPFILE) $(LINKER_FLAGS_)
	rm $(TMPFILE)

design-d: $(PROGRAM_SRC)
	$(eval TMPFILE := $(shell mktemp))
	$(COMPILER_) $(COMPILER_FLAGS_) $(COMPILER_FLAGS_D_) \
		$(AIT_FLAGS_) $(AIT_FLAGS_DESIGN_) $(AIT_FLAGS_D_) \
		$^ -o $(TMPFILE) $(LINKER_FLAGS_)
	rm $(TMPFILE)

bitstream-p: $(PROGRAM_SRC)
	$(eval TMPFILE := $(shell mktemp))
	$(COMPILER_) $(COMPILER_FLAGS_) \
		$(AIT_FLAGS_) \
		$^ -o $(TMPFILE) $(LINKER_FLAGS_)
	rm $(TMPFILE)

bitstream-i: $(PROGRAM_SRC)
	$(eval TMPFILE := $(shell mktemp))
	$(COMPILER_) $(COMPILER_FLAGS_) $(COMPILER_FLAGS_I_) \
		$(AIT_FLAGS_) \
		$^ -o $(TMPFILE) $(LINKER_FLAGS_)
	rm $(TMPFILE)

bitstream-d: $(PROGRAM_SRC)
	$(eval TMPFILE := $(shell mktemp))
	$(COMPILER_) $(COMPILER_FLAGS_) $(COMPILER_FLAGS_D_) \
		$(AIT_FLAGS_) $(AIT_FLAGS_D_) \
		$^ -o $(TMPFILE) $(LINKER_FLAGS_)
	rm $(TMPFILE)

clean:
	rm -fv *.o $(PROGRAM_)-? $(PROGRAM_)_hls_automatic_clang.cpp ait_extracted.json
	rm -frv $(PROGRAM_)_ait
