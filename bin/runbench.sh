#!/bin/bash
echo "Running  SHOC Benchmarks";
export MIC_ENV_PREFIX=MIC
export MIC_USE_2MB_BUFFERS=32K
export MIC_KMP_AFFINITY=granularity=fine,balanced
export MIC_BUFFERSIZE=128M
export MIC_OMP_NUM_THREADS=240
export MIC_MKL_DYNAMIC=false
echo "Running BusSpeedDownload";
./BusSpeedDownload -s 4&>busspeeddownload.log
echo "Running BusSpeedReadback";
./BusSpeedReadback -s 4&>busspeedreadback.log
export MIC_OMP_NUM_THREADS=240
echo "Running DeviceMemory";
./DeviceMemory &>devicememory.log
echo "Running FFT";
export MIC_OMP_NUM_THREADS=240
./FFT --MB 256 &>fft.log
echo "Running GEMM";
export MIC_OMP_NUM_THREADS=240
./GEMM --N 4096 &>gemm.log
echo "Running MaxFlops";
export MIC_OMP_NUM_THREADS=240
./MaxFlops -s 4&>maxflops.log
echo "Running Reduction";
export MIC_OMP_NUM_THREADS=120
./Reduction -s 4&>reduction.log
echo "Running S3D";
export MIC_OMP_NUM_THREADS=240
./S3D -s 4&>s3d.log
echo "Running Sort";
./Sort -s 4&>sort.log
echo "Running Spmv";
export MIC_OMP_NUM_THREADS=236
./Spmv -s 4 &>spmv.log
echo "Running Triad";
export MIC_OMP_NUM_THREADS=240
./Triad -s 4 &>triad.log
echo "Running MonteCarlo";
export MIC_OMP_NUM_THREADS=240
OFFLOAD_INIT=on_start
./MC -s 4 &>MC.log

echo "Running MD";
export MIC_OMP_NUM_THREADS=236
export MIC_KMP_AFFINITY=compact
export MIC_KMP_BLOCKTIME=2s
export MIC_USE_2MB_BUFFERS=128K
./MD -s 4 &> md.log
echo "Running Scan";
export MIC_OMP_NUM_THREADS=240
export KMP_DYNAMIC_MODE=load_balance
./Scan -s 4&>scan.log
echo "Running Stencil2D";
export MIC_OMP_NUM_THREADS=236
export KMP_DYNAMIC_MODE=load_balance
export MIC_USE_2MB_BUFFERS=128K
export MIC_OMP_NUM_THREADS=240
./Stencil2Dmain -s 4 &>stencil2d.log
