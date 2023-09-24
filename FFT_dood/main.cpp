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

    FFT_lib testFFT(x_real, N);

    std::complex<double> a = 1.0;
    std::complex<double> b = 2.0;
    std::complex<double> A;
    std::complex<double> B;
    std::complex<double>* A_ptr = &A;
    std::complex<double>* B_ptr = &B;

    testFFT.init();
    testFFT.testFunc();
    testFFT.calculateFFT();
    testFFT.two_pointDFT(a, b, A_ptr, B_ptr);

    std::cout << "Hello World!\n";
}

// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
