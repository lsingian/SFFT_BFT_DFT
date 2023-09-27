// FFT_dood.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include "FFT_lib.h"

int main()
{
    uint32_t N = 8;

    // test sample values
    double x_real[8] = {
        // sinusoid with f = 1000 Hz, amplitude multiplied by 10^8. values taken from MATLAB
        0,
        -6.428332918551266e-05,
        -1.285666583710253e-04,
        -1.019005173792452e-04,
        -2.571333167420507e-04,
        -4.123661161048562e-04,
        -2.038010347584904e-04,
        -3.590338341212959e-04
    };

    double x_real2[8] = {
        // sinusoid with f = 5000 Hz, amplitude multiplied by 10^8. values taken from MATLAB
        0,
        -4.856823539568490e-05,
        -9.713647079136980e-05,
        -1.457047061870547e-04,
        -1.942729415827396e-04,
        0.001212350345858,
        -2.914094123741094e-04,
        -0.001795169170606
    };

    FFT_lib::FFTAlgorithm_Types algorithm = FFT_lib::FFTAlgorithm_Types::DECIMATION_IN_FREQUENCY;

    FFT_lib testFFT(N, algorithm); // TODO need to remove x_real as an input since this will cause an object to be created each time we want to do an FFT.

    std::complex<double>* fft_results;

    testFFT.init();
    fft_results = testFFT.calculateFFT(x_real);

    double real = fft_results[3].real();
    double imaginary = fft_results[3].imag();


    int endbreakpoint = 0;
}
