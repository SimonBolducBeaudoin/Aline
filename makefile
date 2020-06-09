tempo1 = $(wildcard ../Multi_array/*.o)
tempo2 = $(wildcard ../FFTW_extra/obj/*.o)
tempo3 = $(wildcard ../Interpolation/obj/*.o)
tempo4 = $(wildcard ../Omp_extra/obj/*.o)
EXTERNAL_OBJ = $(tempo1) $(tempo2) $(tempo3) $(tempo4)

TARGET_NAME = Aline
TARGET_STATIC = $(TARGET_NAME).a
PYLIB_EXT = $(if $(filter $(OS),Windows_NT),.pyd,.so)
TARGET_PYLIB = ../Python/$(TARGET_NAME)$(PYLIB_EXT)

IDIR = includes
ODIR = obj
LDIR = lib
SDIR = src

SRC  = $(wildcard $(SDIR)/*.cpp)
OBJ  = $(patsubst $(SDIR)/%.cpp,$(ODIR)/%.o,$(SRC))
DEPS = $(OBJ:.o=.d)

CC = $(OS:Windows_NT=x86_64-w64-mingw32-)g++
CFLAGS = -Ofast -march=native -std=c++14 -MMD -MP -Wall $(OS:Windows_NT=-DMS_WIN64 -D_hypot=hypot)
OMPFLAGS = -fopenmp -fopenmp-simd
FFTWFLAGS = -lfftw3
MATHFLAGS = -lm
SHRFLAGS = -fPIC -shared

PY = $(OS:Windows_NT=/c/Anaconda2/)python
ifeq ($(USERNAME),Sous-sol)
    PY = $(OS:Windows_NT=/c/ProgramData/Anaconda2/)python
endif

PYINCL := $(shell $(PY) -m pybind11 --includes)
ifneq ($(OS),Windows_NT)
    PYINCL += -I /usr/include/python2.7/
endif

PYLIBS = $(OS:Windows_NT=-L /c/Anaconda2/ -l python27)
ifeq ($(USERNAME),Sous-sol)
    PYLIBS = $(OS:Windows_NT=-L /c/ProgramData/Anaconda2/ -l python27)
endif

$(TARGET_PYLIB): $(OBJ)
	@ echo " "
	@ echo "---------Compile library $(TARGET_PYLIB)---------"
	$(CC) $(SHRFLAGS) -o $(TARGET_PYLIB) $(OBJ) $(EXTERNAL_OBJ) $(CFLAGS) $(OMPFLAGS) $(FFTWFLAGS) $(MATHFLAGS) $(PYLIBS)
	
$(ODIR)/%.o : $(SDIR)/%.cpp
	@ echo " "
	@ echo "---------Compile object $@ from $<--------"
	$(CC) $(SHRFLAGS) -c -Wall -o $@ $< $(CFLAGS) $(OMPFLAGS) $(PYINCL) 

-include $(DEPS)

clean:
	@rm -f $(TARGET_PYLIB) $(TARGET_STATIC) $(OBJ) $(DEPS)
	 	 
.PHONY: clean, dummy