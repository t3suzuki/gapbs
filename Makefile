# See LICENSE.txt for license details.


MEMKIND_PATH = /home/tomoya-s/memkind/install
ABT_PATH = /home/tomoya-s/argobots/install
#CXX = clang++
CXX_FLAGS += -std=c++11 -O3 -Wall -g
CXX_FLAGS += -DDAX=1
CXX_FLAGS += -DABT=1
#PAR_FLAG = -fopenmp
CXX_FLAGS += -I$(ABT_PATH)/include -I$(MEMKIND_PATH)/include
LD_FLAGS += -L$(ABT_PATH)/lib -labt -L$(MEMKIND_PATH)/lib -lmemkind

ifneq (,$(findstring icpc,$(CXX)))
	PAR_FLAG = -openmp
endif

ifneq (,$(findstring sunCC,$(CXX)))
	CXX_FLAGS = -std=c++11 -xO3 -m64 -xtarget=native
	PAR_FLAG = -xopenmp
endif

ifneq ($(SERIAL), 1)
	CXX_FLAGS += $(PAR_FLAG)
endif

KERNELS = bc bfs cc cc_sv pr pr_spmv sssp tc
SUITE = $(KERNELS) converter

.PHONY: all
all: $(SUITE)

% : src/%.cc src/*.h
	$(CXX) $(CXX_FLAGS) $< -o $@ $(LD_FLAGS)

# Testing
include test/test.mk

# Benchmark Automation
include benchmark/bench.mk


.PHONY: clean
clean:
	rm -f $(SUITE) test/out/*
