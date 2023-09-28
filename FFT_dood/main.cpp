// FFT_dood.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include "FFT_lib.h"

int main()
{
    uint32_t N = 8;

    // test sample values
    double x_real[8] = {
        // sinusoid with f = 1000 Hz, values taken from MATLAB
        0,
        -6.428332918551266e-05,
        -1.285666583710253e-04,
        -1.019005173792452e-04,
        -2.571333167420507e-04,
        -4.123661161048562e-04,
        -2.038010347584904e-04,
        -3.590338341212959e-04
    };

    FFT_lib::FFTAlgorithm_Types algorithm = FFT_lib::FFTAlgorithm_Types::DECIMATION_IN_FREQUENCY;

    FFT_lib testFFT(N, algorithm); // TODO need to remove x_real as an input since this will cause an object to be created each time we want to do an FFT.

    std::complex<double>* fft_results;
    double* magnitude;

    testFFT.init();
    fft_results = testFFT.calculateFFT(x_real);

    magnitude = testFFT.calculateMagnitude(fft_results);

    int endbreakpoint = 0;
}
