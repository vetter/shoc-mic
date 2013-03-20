// This example from an alpha release of the Scalable HeterOgeneous Computing
// (SHOC) Benchmark Suite Alpha v1.1.4a-mic for Intel MIC architecture
// Contact: Kyle Spafford <kys@ornl.gov>
//          Rezaur Rahman <rezaur.rahman@intel.com>
//
// Copyright (c) 2011, UT-Battelle, LLC
// Copyright (c) 2013, Intel Corporation
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

#include <stdlib.h>
#include <string.h>
#include <cassert>
#include "InitializeMatrix2D.h"



template<class T>
void
Initialize<T>::operator()( Matrix2D<T>& mtx )
{
    srand48( seed );

    int nTileRows = mtx.GetNumRows() - 2 * haloWidth;
    if( (rowPeriod != -1) && (rowPeriod < nTileRows) )
    {
        nTileRows = rowPeriod;
    }

    int nTileCols = mtx.GetNumColumns() - 2 * haloWidth;
    if( (colPeriod != -1) && (colPeriod < nTileCols) )
    {
        nTileCols = colPeriod;
    }


    // initialize first tile
    for( unsigned int i = 0; i < nTileRows; i++ )
    {
        for( unsigned int j = 0; j < nTileCols; j++ )
        {
#ifndef READY
            mtx.GetData()[i+haloWidth][j+haloWidth] = i * j;
#else
            mtx.GetData()[i+haloWidth][j+haloWidth] = (T)drand48();
#endif // READY
        }
    }

    // initialize any remaining tiles
    // first we fill along rows a tile at a time,
    // then fill out along columns a row at a time
    if( colPeriod != -1 )
    {
        int nTiles = (mtx.GetNumColumns() - 2*haloWidth) / colPeriod;
        if( (mtx.GetNumColumns() - 2*haloWidth) % colPeriod != 0 )
        {
            nTiles += 1;
        }

        for( unsigned int t = 1; t < nTiles; t++ )
        {
            for( unsigned int i = 0; i < nTileRows; i++ )
            {
                memcpy( &(mtx.GetData()[haloWidth + i][haloWidth + t*nTileCols]),
                        &(mtx.GetData()[haloWidth + i][haloWidth]),
                        nTileCols * sizeof(T) );
            }
        }
    }
    if( rowPeriod != -1 )
    {
        int nTiles = (mtx.GetNumRows() - 2*haloWidth) / rowPeriod;
        if( (mtx.GetNumRows() - 2*haloWidth) % rowPeriod != 0 )
        {
            nTiles += 1;
        }

        for( unsigned int t = 1; t < nTiles; t++ )
        {
            for( unsigned int i = 0; i < nTileRows; i++ )
            {
                memcpy( &(mtx.GetData()[haloWidth + t*nTileRows + i][haloWidth]),
                        &(mtx.GetData()[haloWidth + i][haloWidth]),
                        (mtx.GetNumColumns() - 2*haloWidth) * sizeof(T) );
            }
        }
    }

    // initialize halo
    for( unsigned int i = 0; i < mtx.GetNumRows(); i++ )
    {
        for( unsigned int j = 0; j < mtx.GetNumColumns(); j++ )
        {
            bool inHalo = false;

            if( (i < haloWidth) || (i > mtx.GetNumRows() - 1 - haloWidth) )
            {
                inHalo = true;
            }
            else if( (j < haloWidth) || (j > mtx.GetNumColumns() - 1 - haloWidth) )
            {
                inHalo = true;
            }

            if( inHalo )
            {
                mtx.GetData()[i][j] = haloVal;
            }
        }
    }
}

