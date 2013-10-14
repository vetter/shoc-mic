// This example from an alpha release of the Scalable HeterOgeneous Computing
// (SHOC) Benchmark Suite Alpha v1.1.4a-mic for Intel MIC architecture
// Contact: Kyle Spafford <kys@ornl.gov>
//          Rezaur Rahman <rezaur.rahman@intel.com>
//
// Copyright (c) 2011, UT-Battelle, LLC
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of Oak Ridge National Laboratory, nor UT-Battelle, LLC,
//    nor the names of its contributors may be used to endorse or promote
//    products derived from this software without specific prior written
//    permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY,
// OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
// THE POSSIBILITY OF SUCH DAMAGE.

#if defined(__APPLE__)
#if _GLIBCXX_ATOMIC_BUILTINS == 1
#undef _GLIBCXX_ATOMIC_BUILTINS
#endif // _GLIBCXX_ATOMIC_BUILTINS
#endif // __APPLE__

#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <stdio.h>

#include "omp.h"
#include "math.h"
#include "offload.h"
#include "Timer.h"
#include "MICStencil.cpp"

#define LINESIZE    64
#define ALLOC       alloc_if(1)
#define FREE        free_if(1)
#define RETAIN      free_if(0)
#define REUSE       alloc_if(0)

////////////////////////////////////////////////////////////////
// TODO: Tune Threads, Partitions according to card's parameters
////////////////////////////////////////////////////////////////

template <class T> void
MICStencil<T>::operator()( Matrix2D<T>& mtx, unsigned int nIters )
{
    unsigned int uDimWithHalo    = mtx.GetNumRows();
    unsigned int uHaloWidth      = LINESIZE / sizeof(T);
    unsigned int uImgElements    = uDimWithHalo * uDimWithHalo;

    __declspec(target(mic), align(LINESIZE)) T* pIn = mtx.GetFlatData();

    __declspec(target(mic), align(sizeof(T)))    T wcenter      = this->wCenter;
    __declspec(target(mic), align(sizeof(T)))    T wdiag        = this->wDiagonal;
    __declspec(target(mic), align(sizeof(T)))    T wcardinal    = this->wCardinal;

    #pragma offload target(mic) in(pIn:length(uImgElements) ALLOC RETAIN)
    {
        // Just copy pIn to compute the copy transfer time
    }

    #pragma offload target(mic) in(pIn:length(uImgElements) REUSE RETAIN)    \
                                in(uImgElements) in(uDimWithHalo)            \
                                in(wcenter) in(wdiag) in(wcardinal)
    {
        unsigned int uRowPartitions = sysconf(_SC_NPROCESSORS_ONLN) / 4 - 1;
        unsigned int uColPartitions = 4;    // Threads per core for KNC

        unsigned int uRowTileSize    = (uDimWithHalo - 2 * uHaloWidth) / uRowPartitions;
        unsigned int uColTileSize    = (uDimWithHalo - 2 * uHaloWidth) / uColPartitions;

        uRowTileSize = ((uDimWithHalo - 2 * uHaloWidth) % uRowPartitions > 0) ? (uRowTileSize + 1) : (uRowTileSize);

        // Should use the "Halo Val" when filling the memory space
        T *pTmp     = (T*)pIn;
        T *pCrnt = (T*)memset((T*)_mm_malloc(uImgElements * sizeof(T), LINESIZE), 0, uImgElements * sizeof(T));

        #pragma omp parallel firstprivate(pTmp, pCrnt, uRowTileSize, uColTileSize, uHaloWidth, uDimWithHalo)
        {
            unsigned int uThreadId = omp_get_thread_num();

            unsigned int uRowTileId = uThreadId / uColPartitions;
            unsigned int uColTileId = uThreadId % uColPartitions;

            unsigned int uStartLine = uRowTileId * uRowTileSize + uHaloWidth;
            unsigned int uStartCol  = uColTileId * uColTileSize + uHaloWidth;

            unsigned int uEndLine = uStartLine + uRowTileSize;
            uEndLine = (uEndLine > (uDimWithHalo - uHaloWidth)) ? uDimWithHalo - uHaloWidth : uEndLine;

            unsigned int uEndCol    = uStartCol  + uColTileSize;
            uEndCol  = (uEndCol  > (uDimWithHalo - uHaloWidth)) ? uDimWithHalo - uHaloWidth : uEndCol;

            T    cardinal0 = 0.0;
            T    diagonal0 = 0.0;
            T    center0   = 0.0;

            unsigned int cntIterations, i, j;

            for (cntIterations = 0; cntIterations < nIters; cntIterations ++)
            {
                // Do Stencil Operation
                for (i = uStartLine; i < uEndLine; i++)
                {
                    T * pCenter      = &pTmp [ i * uDimWithHalo];
                    T * pTop         = pCenter - uDimWithHalo;
                    T * pBottom      = pCenter + uDimWithHalo;
                    T * pOut         = &pCrnt[ i * uDimWithHalo];

                    __assume_aligned(pCenter, 64);
                    __assume_aligned(pTop,    64);
                    __assume_aligned(pBottom, 64);
                    __assume_aligned(pOut,    64);

                    #pragma simd vectorlengthfor(float)
                    for (j = uStartCol; j < uEndCol; j++)
                    {
                        cardinal0   = pCenter[j - 1] + pCenter[j + 1] + pTop[j] + pBottom[j];
                        diagonal0   = pTop[j - 1] + pTop[j + 1] + pBottom[j - 1] + pBottom[j + 1];
                        center0     = pCenter[j];

                        pOut[j]     = wcardinal * cardinal0 + wdiag * diagonal0 + wcenter * center0;
                    }
                }

                #pragma omp barrier
                ;

                // Switch pointers
                T* pAux    = pTmp;
                pTmp     = pCrnt;
                pCrnt    = pAux;
            } // End For

        } // End Parallel

        _mm_free(pCrnt);
    } // End Offload

    #pragma offload target(mic) out(pIn:length(uImgElements) REUSE FREE)
    {
        // Just copy back pIn
    }
}

void
EnsureStencilInstantiation( void )
{
    MICStencil<float> csf( 0, 0, 0, 0 );
    Matrix2D<float> mf( 2, 2 );
    csf( mf, 0);

    MICStencil<double> csd( 0, 0, 0, 0 );
    Matrix2D<double> md( 2, 2 );
    csd( md, 0);
}
