# Basics
CC		= icc
CPP     = icc
LD		= icc
CXX     = icc

SHOC_COMMON=./common

MICROOT=./bin

BINDIR=./bin

SUBDIRS = stencil2d

VPATH=$(SHOC_COMMON) level0 md reduction scan triad spmv sort fft gemm s3d 

COMMON_OBJS = main.o \
			  Option.o \
			  OptionParser.o \
			  Timer.o \
			  ResultDatabase.o \
			  ProgressBar.o

BENCH_OBJS = BusSpeedReadback.o \
			 BusSpeedDownload.o \
			 MD.o \
			 Reduction.o \
			 Scan.o \
			 Triad.o \
			 Spmv.o\
			 Sort.o\
			 FFT.o\
			 GEMM.o \
			 S3D.o\
			 DeviceMemory.o\
			 MaxFlops.o


BENCHMARKPROG = $(patsubst %.o,$(BINDIR)/%,$(BENCH_OBJS))

# Flags to include MIC MKL library
MIC_MKL_LIBS=-lmkl_core

CFLAGS          = -offload-option,mic,compiler,"-mP2OPT_hpo_vec_check_dp_trip=F -fimf-precision=low -fimf-domain-exclusion=15  -opt-report 1" -mP2OPT_hlo_pref_issue_second_level_prefetch=F -mP2OPT_hlo_pref_issue_first_level_prefetch=F -vec-report0  -O3  -openmp -intel-extensions -I$(SHOC_COMMON) -opt-report-phase:offload -openmp-report


LDFLAGS = -vec-report0 -parallel -mkl

%.o: %.cpp
	$(CC) -c $< $(CFLAGS) $(CXXFLAGS)

$(BINDIR)/%: %.o 
	$(CC) -o $@ $(CFLAGS) $(CXXFLAGS) $(LDFLAGS) $(COMMON_OBJS) $< $(LIBS) 

all : $(BENCHMARKPROG)
	for d in $(SUBDIRS); do (cd $$d; $(MAKE) ); done

$(BENCHMARKPROG) : $(COMMON_OBJS)

clean :
	rm *.o *.out *.xnlog
	rm ./bin/*.o
	rm ./bin/*.out

# individual targets
fft : $(BINDIR)/FFT

$(BINDIR)/FFT : $(COMMON_OBJS)

gemm : $(BINDIR)/GEMM

$(BINDIR)/GEMM : $(COMMON_OBJS)

md : $(BINDIR)/MD

$(BINDIR)/MD : $(COMMON_OBJS)

reduction : $(BINDIR)/Reduction

$(BINDIR)/Reduction : $(COMMON_OBJS)

scan: $(BINDIR)/Scan

$(BINDIR)/Scan : $(COMMON_OBJS)

spmv: $(BINDIR)/Spmv

$(BINDIR)/Sort : $(COMMON_OBJS)

sort: $(BINDIR)/Sort

$(BINDIR)/S3D : $(COMMON_OBJS)

s3d : $(BINDIR)/S3D
