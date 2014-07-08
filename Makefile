CC=gcc
#CXX=g++
CXX=clang++

COMMON=-g -O3 -Wall # -g for VTune Analysis
LDFLAGS=-ltbb

CFLAGS = $(COMMON)
CXXFLAGS = $(COMMON)

# Modify for your platform
TBB21_INSTALL_DIR=/opt/intel/tbb/2.1
TBB_ARCH_PLATFORM=em64t/cc4.1.0_libc2.4_kernel2.6.16.21

all: BarnesHut

install: BarnesHut
	cp src/BarnesHut driver

BarnesHut: src/BarnesHut.cpp
	# To fix the load path into the executable (ELF header) so that we don't need to set the LD_LIBRARY_PATH at runtime
	# LD_RUN_PATH only works on Linux. On OS X we need something else .... 
	# To check if the path is embedded properly, run 'readelf -d' on the executable and look for the rpath setting 
	# LD_RUN_PATH=$(TBB21_INSTALL_DIR)/$(TBB_ARCH_PLATFORM)/lib $(CXX) $(CXXFLAGS) -I$(TBB21_INSTALL_DIR)/include -L$(TBB21_INSTALL_DIR)/$(TBB_ARCH_PLATFORM)/lib $(LDFLAGS) src/BarnesHut.cpp -o BarnesHut 

	$(CXX) $(CXXFLAGS) $(LDFLAGS) src/BarnesHut.cpp -o BarnesHut 
	mv BarnesHut src

clean:
	rm -f driver/BarnesHut src/BarnesHut
