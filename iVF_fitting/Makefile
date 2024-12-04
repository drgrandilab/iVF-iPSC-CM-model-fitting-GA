
CXX=icpc

# CXX=mpiicpc

INC =  -Ilib -#I/opt/OpenBLAS/include/ -L/opt/OpenBLAS/lib/
CFLAGS =  -Ilib  -std=c++11 -g
CXXFLAGS = -Ilib -std=c++11 -O3 -lsundials_cvode  -xcore-avx2  -axcore-avx2 -g# -D USE_V1 #-lopenblas # -axcommon-avx512 -g #-lm


INCLUDE=-I/usr/local/include -Ilib
LIB    =-L/usr/local/lib

# LDFLAGS= ${INCLUDE} ${LIB}  #-vec-report -profile-functions -profile-loops=all -profile-loops-report=2
# LDLIBS= -lsundials_cvode#-lopenblas -lpthread


# cvode_main_ipsc_baseline.cpp
# cvode_main_ipsc_baseline:	cvode_main_ipsc_baseline.cpp
# 	${CC} ${CFLAGS} ${INCLUDE} ${LDFLAGS}  cvode_main_ipsc_baseline.cpp -o main_ipsc_baseline -lsundials_cvode #  -lm

BUILD_DIR=lib/build
INC_PARAMS=-Ilib -I/usr/local/include

# CXX = g++
srcs = $(wildcard lib/*.cpp)
objs := $(srcs:%=$(BUILD_DIR)/%.o)
deps = $(srcs:.cpp=.d)


cvode_main_ipsc_baseline_05Hz: cvode_main_ipsc_baseline_05Hz.cpp $(objs)  
	$(CXX) $(CXXFLAGS) $(INC_PARAMS) $^ -o $@


cvode_main_ipsc_baseline_10Hz: cvode_main_ipsc_baseline_10Hz.cpp $(objs)  
	$(CXX) $(CXXFLAGS) $(INC_PARAMS) $^ -o $@

cvode_main_ipsc_ramp: cvode_main_ipsc_ramp.cpp $(objs)  
	$(CXX) $(CXXFLAGS) $(INC_PARAMS) $^ -o $@



clean:
	$(RM) $(objs) $(deps) cvode_main_ipsc_baseline


$(BUILD_DIR)/%.cpp.o: %.cpp
	$(MKDIR_P) $(dir $@)
	$(CXX) $(INC_PARAMS)   $(CXXFLAGS)  -MMD -MP -c $< -o $@ 
.PHONY: clean

MKDIR_P ?= mkdir -p
