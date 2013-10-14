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

#include <iostream>
#include <sstream>
#include <assert.h>
#include <math.h>
#include "omp.h"

#include "OptionParser.h"
#include "ResultDatabase.h"
#include "Timer.h"
#include "BadCommandLine.h"
#include "InvalidArgValue.h"
#include "Matrix2D.h"
#include "HostStencilFactory.h"
#include "HostStencil.h"
#include "MICStencilFactory.h"
#include "MICStencil.h"
#include "InitializeMatrix2D.h"
#include "ValidateMatrix2D.h"
#include "StencilUtil.h"

#include "InitializeMatrix2D.cpp"
#include "ValidateMatrix2D.cpp"
#include "StencilUtil.cpp"
#include "StencilFactory.cpp"
#include "CommonMICStencilFactory.cpp"
#include "HostStencil.cpp"
#include "MICStencil.cpp"
#include "HostStencilFactory.cpp"
#include "MICStencilFactory.cpp"

// prototypes of auxiliary functions defined in this file or elsewhere
void CheckOptions( const OptionParser& opts );
void EnsureStencilInstantiation( void );

#define LINESIZE 64

template<class T> void
MICValidate(const Matrix2D<T>& s,
            const Matrix2D<T>& t,
            double valErrThreshold,
            unsigned int nValErrsToPrint)
{
    assert( (s.GetNumRows() == t.GetNumRows()) &&
            (s.GetNumColumns() == t.GetNumColumns()) );
    unsigned int uHaloWidth = LINESIZE / sizeof(T);

    for( unsigned int i = uHaloWidth; i < s.GetNumRows() - uHaloWidth; i++ )
    {
        for( unsigned int j = uHaloWidth; j < s.GetNumColumns() - uHaloWidth; j++ )
        {
            T expVal    = s.GetConstData()[i][j];
            T actualVal = t.GetConstData()[i][j];
            T delta     = fabsf( actualVal - expVal );
            T relError  = (expVal != 0.0f) ? delta / expVal : 0.0f;

            if( relError > valErrThreshold )
            {
                std::cout<<"Failed\n";
                return;
            }
        }
    }

    std::cout<<"Passed\n";
}

template<class T> void
DoTest( const char* timerDesc, ResultDatabase& resultDB, OptionParser& opts )
{
    StencilFactory<T>*  stdStencilFactory  = NULL;
    Stencil<T>*         stdStencil         = NULL;
    StencilFactory<T>*  testStencilFactory = NULL;
    Stencil<T>*         testStencil        = NULL;

    stdStencilFactory   = new HostStencilFactory<T>;
    testStencilFactory  = new MICStencilFactory<T>;
    assert( (stdStencilFactory != NULL) && (testStencilFactory != NULL) );

    // Do a sanity check on option values
    CheckOptions( opts );
    stdStencilFactory->CheckOptions( opts );
    testStencilFactory->CheckOptions( opts );

    // Extract and validate options
    std::vector<long long> arrayDims = opts.getOptionVecInt( "customSize" );
    if( arrayDims.size() != 2 )
        cerr << "Dim size: " << arrayDims.size() << "\n";

    if (arrayDims[0] == 0) // User has not specified a custom size
    {
        const int probSizes[4] = { 768, 1408, 2048, 4096 };
        int sizeClass = opts.getOptionInt("size");

        if (!(sizeClass >= 0 && sizeClass < 5))
        {
            //throw InvalidArgValue( "Size class must be between 1-4" );
        }

        arrayDims[0] = arrayDims[1] = probSizes[sizeClass - 1];
    }

    long int seed = (long)opts.getOptionInt( "seed" );
    bool beVerbose = opts.getOptionBool( "verbose" );
    unsigned int nIters = (unsigned int)opts.getOptionInt( "num-iters" );
    double valErrThreshold = (double)opts.getOptionFloat( "val-threshold" );
    unsigned int nValErrsToPrint = (unsigned int)opts.getOptionInt( "val-print-limit" );

    // Define Halo
    unsigned int haloWidth = LINESIZE / sizeof(T);
    float haloVal          = (float)opts.getOptionFloat( "haloVal" );

    // Build a description of this experiment
    std::ostringstream experimentDescriptionStr;
    experimentDescriptionStr<< nIters << ':'<< arrayDims[0] << 'x' << arrayDims[1];

    unsigned int  nPasses = (unsigned int)opts.getOptionInt( "passes" );
    unsigned long npts    = (arrayDims[0] + 2 * haloWidth - 2) * (arrayDims[1] + 2*haloWidth - 2);
    unsigned long nflops  = npts * 11 * nIters;
    cout<<"FLOP are = "<< nflops <<endl;

    Matrix2D<T> exp(arrayDims[0] + 2 * haloWidth, arrayDims[1] + 2 * haloWidth);
    Initialize<T> init(seed, haloWidth, haloVal);

    init(exp);
    if(beVerbose)
        std::cout << "initial state:\n" << exp << std::endl;

    stdStencil = stdStencilFactory->BuildStencil(opts);

    (*stdStencil)(exp, nIters);

    if( beVerbose )
        std::cout << "expected result:\n" << exp << std::endl;

    // Compute the result on the Xeon Phi device
    Matrix2D<T> data(arrayDims[0] + 2 * haloWidth, arrayDims[1] + 2 * haloWidth);
    testStencil = testStencilFactory->BuildStencil( opts );

    std::cout<<"Passes:"<<nPasses<<endl;
    for( unsigned int pass = 0; pass < nPasses; pass++ )
    {
        init(data);

        double start         = curr_second();
        (*testStencil)(data, nIters);
        double elapsedTime     = curr_second() - start;

        double gflopsPCIe     = (nflops / elapsedTime) / 1e9;

        resultDB.AddResult(timerDesc, experimentDescriptionStr.str(), "GFLOPS_PCIe", gflopsPCIe);

        if( beVerbose )
            std::cout << "observed result, pass " << pass << ":\n"<< data<< std::endl;

        MICValidate(exp, data, valErrThreshold, nValErrsToPrint);
    }

    // clean up - normal termination
    delete stdStencil;
    delete stdStencilFactory;
    delete testStencil;
    delete testStencilFactory;
}

void RunBenchmark(OptionParser& opts, ResultDatabase& resultDB )
{
    std::cout << "Running Single Precision test :" << std::endl;
    DoTest<float>( "SP_Sten2D", resultDB, opts);

    std::cout << "Running Double Precision test :" << std::endl;
    DoTest<double>( "DP_Sten2D", resultDB, opts );
}

// Adds command line options to given OptionParser
void addBenchmarkSpecOptions( OptionParser& opts )
{
    opts.addOption( "customSize",      OPT_VECINT, "0,0",   "specify custom problem size");
    opts.addOption( "num-iters",       OPT_INT,    "1000",  "number of stencil iterations" );
    opts.addOption( "weight-center",   OPT_FLOAT,  "0.25",  "center value weight" );
    opts.addOption( "weight-cardinal", OPT_FLOAT,  "0.15",  "cardinal values weight" );
    opts.addOption( "weight-diagonal", OPT_FLOAT,  "0.05",  "diagonal values weight" );
    opts.addOption( "seed",            OPT_INT,    "71594", "random number generator seed" );
    opts.addOption( "val-threshold",   OPT_FLOAT,  "0.1",   "validation error threshold" );
    opts.addOption( "val-print-limit", OPT_INT,    "15",    "number of validation errors to print" );
    opts.addOption( "haloVal",         OPT_FLOAT,  "0.0",   "value to use for halo data" );
}

// validate stencil-independent values
void CheckOptions( const OptionParser& opts )
{
    // check matrix dimensions - must be 2d, must be positive
    std::vector<long long> arrayDims = opts.getOptionVecInt( "customSize" );
    if( arrayDims.size() != 2 )
    {
        throw InvalidArgValue( "overall size must have two dimensions" );
    }

    if( (arrayDims[0] < 0) || (arrayDims[1] < 0) )
    {
        throw InvalidArgValue( "each size dimension must be positive" );
    }

    // validation error threshold must be positive
    float valThreshold = opts.getOptionFloat( "val-threshold" );
    if( valThreshold <= 0.0f )
    {
        throw InvalidArgValue( "validation threshold must be positive" );
    }

    // number of validation errors to print must be non-negative
    int nErrsToPrint = opts.getOptionInt( "val-print-limit" );
    if( nErrsToPrint < 0 )
    {
        throw InvalidArgValue( "number of validation errors to print must be non-negative" );
    }
}
