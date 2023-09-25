#include "FFT_lib.h"
#include <iostream>
#include <cmath>
#include "math.h"

#define PI 3.14159265358979323846


void FFT_lib::init() 
{
	calculateButterflyPairsPerStage();
	calculateTwiddlePerStage();
}
void FFT_lib::testFunc() 
{
	std::cout << "This is a test func dood! " << std::endl;
}

void FFT_lib::resetButterflyIndicesState() 
{
	for (int i = 0; i < N; i++)
	{
		butterFlyIndicesState[i] = false;

		//butt.indexArray[i].index = (uint32_t*) malloc(sizeof(uint32_t)); // TODO address the mallocs later
		//butt.indexArray[i].index = 0;
		//butt.indexArray[i].paired = false;
	}
}

void FFT_lib::calculateButterflyPairsPerStage()
{
	// row iteration of the 2D array 
	for (int stageIdx = 0; stageIdx < numStages; stageIdx++) 
	{
		// The butterfly pair is the index located 2^(numStages - stageIdx - 1) away.
		// Ex. 8-point FFT: N = 8 = 2^4 
		// For the first stage, 
		// the pair for 0 is 4, which is 2^(4 - 1) away from 0.
		// the pair for 1 is 5,  which is 2^(4 - 1) away from 1.
		uint32_t pairIndexDelta = pow(2, numStages - stageIdx - 1);
		uint32_t pairIdx = 0; // used for indexing the 2D array that contains the butterfly pairs for each stage
		
		resetButterflyIndicesState();

		// column iteration of the 2D array
		for (int i = 0; i < N; i++) 
		{   
			// index has already been paired, no need to pair
			if (butterFlyIndicesState[i]) 
			{
				continue;
			}

			// Find the index of the butterfly pair
			uint32_t pairIndex = i + butterFlyIndices[pairIndexDelta];

			// Store the butterfly pair in the 2D array of stage rows and pairs columns.
			stageButterflyPairs[stageIdx][pairIdx++] = butterFlyIndices[i];
			butterFlyIndicesState[i] = true;

			stageButterflyPairs[stageIdx][pairIdx++] = pairIndex;
			butterFlyIndicesState[i + pairIndexDelta] = true; 
		}
	}

	int test = 0;

	//// Construct final stage pairs - combination of odd and even sequence DFTs
	//for (uint32_t n = 0; n < N / 2; n++)
	//{
	//	butt.A[n] = 2 * n;
	//	butt.B[n] = butt.A[n] + 1;
	//}

	//// Construct preceding stage pairs, starting from the final stage pairs.
	//for (uint32_t iStage = 1; iStage < numStages + 1; iStage++)
	//{

	//}
};

void FFT_lib::calculateTwiddlePerStage() 
{
	// For radix-2 Decimation in time (DIT) N-point FFT algorithm,
	// Each stage has N/2 butterfly pairs.
	// Therefore, each stage has N/2 twiddle factors
	// The k-th twiddle factor is [ (k * 2^P) / N ], where...
	// P is stage, P = 1, 2, ... , log2(N) 
	// k = 0, 1, 2, ... (N/2) - 1
	// N = points in FFT
	// [] is the bit-reverse operation. Length of a bit is log2(N) - 1.
	// Reference: https://www.dsprelated.com/showcode/232.php#:~:text=The%20N%2Dpoint%20DIT%20FFT,of%20Figure%201(a).
	const uint32_t t_numStages = 4;
	const uint32_t t_N = 16;
	// const uint32_t t_kMax = t_N / 2; //  kmax for decimation in time FFT.

	std::complex<double> W_table[t_numStages][t_N / 2]; // N/2 twiddle factors per stage
	uint32_t k_test_table[t_numStages][t_N / 2]; // N/2 twiddle factors per stage

	// row iteration of twiddle factor 2D array
	for (int stageIdx = 0; stageIdx < t_numStages; stageIdx++)
	{
		uint32_t t_kMax = (t_N / pow(2,stageIdx + 1)); // kmax is number of unique twiddle factors. It is calculated every stage for decimation in frequency FFT.
		
		// decimation in frequency k-th Twiddle factor is simply, k * (2^P / 2)
		uint32_t kFactor = pow(2, stageIdx + 1) / 2;

		// column iteration of twiddle factor 2D array
		for (int kIndex = 0; kIndex < t_kMax; kIndex++) 
		{
			// TODO debug when switching to decimation in time FFT algo, for now, using decimation in frequency
			//uint32_t tempNum = (kp * pow(2, stageIdx + 1)) / t_N;
			//uint32_t k = bitReverse(tempNum); // TODO - bit reverse
			//std::complex<double> W_k = std::polar(1.0, -2 *  PI * k / N);
			//W_table[stageIdx][kIndex] = W_k;

			uint32_t k = kIndex * kFactor;
			std::complex<double> W_k = std::polar(1.0, -2 * PI * k / N);
			//k_test_table[stageIdx][kIndex] = k; // this is correct, but its a duplicated step

			// populate N/2 columns of 2D array with twiddle factors
			for (int i = 0; i < t_N / 2; i++)
			{
				if (kIndex + (i * t_kMax) < t_N / 2) // ensure that we don't access out of array indices
				{
					//k_test_table[stageIdx][kIndex + (i * t_kMax)] = k; // twiddle factor at kIndex repeats every kMax = N/(2^P) indices
					W_table[stageIdx][kIndex + (i * t_kMax)] = W_k; // twiddle factor at kIndex repeats every kMax = N/(2^P) indices
				}
				else
				{
					break; // succeeding i * t_kMax will cause array index access to go out of bounds
				}
			}
		}
	}

	int test = 1;

	//std::complex<float> _Wnk_Nc(float n, float k)
	//{
	//	std::complex<float> _Wnk_Ncomp;
	//	_Wnk_Ncomp.real(cosf(-2.0f * _Pi * (float)n * k / sampleRate));
	//	_Wnk_Ncomp.imag(sinf(-2.0f * _Pi * (float)n * k / sampleRate));
	//	return _Wnk_Ncomp;
	//}

	//std::complex<float> _Wk_Nc(float k)
	//{
	//	std::complex<float> _Wk_Ncomp;
	//	_Wk_Ncomp.real(cosf(-2.0f * _Pi * k / sampleRate));
	//	_Wk_Ncomp.imag(sinf(-2.0f * _Pi * k / sampleRate));
	//	return _Wk_Ncomp;
	//}
}

uint32_t FFT_lib::bitReverse(uint32_t inputNum) 
{
	unsigned int reversedNum = 0;
	unsigned int numOfBits = sizeof(inputNum) * 8; // Number of bits in 'inputNum'

	for (int i = 0; i < numOfBits; i++) {
		reversedNum = (reversedNum << 1) | (inputNum & 1);
		inputNum >>= 1;
	}
	return reversedNum;

	// Code taken from: https://stackoverflow.com/questions/932079/in-place-bit-reversed-shuffle-on-an-array

	//// Code taken from: https://saturncloud.io/blog/efficient-algorithm-for-bit-reversal-in-c-a-data-scientists-perspective/
	//// Will do my own implementation later
	//unsigned int numOfBits = sizeof(inputNum) * 8; // Number of bits in 'inputNum'
	//unsigned int reversedNum = 0;

	//for (int i = 0; i < numOfBits; i++) {
	//	if ((inputNum & (1 << i))) {
	//		reversedNum |= 1 << ((numOfBits - 1) - i);
	//	}
	//}
	//return reversedNum;
}


void FFT_lib::calculateFFT() // TODO add samples as input here
{
	// TODO add assert that N-point must be a power of 2.	

	// calculate bit-reversed indices pairs for each stage - add this in the init function
	// 
	const uint32_t maxStageIter = numStages;
	const uint32_t maxEvenIndex = N - 2;
	const uint32_t maxOddIndex = N - 1;

	//for (uint32_t iStage = numStages; iStage > 0; iStage--)
	//{
	//	if (1 == iStage) 
	//	{
	//		// final stage combines both even and odd DFTs
	//	}
	//	else 
	//	{
	//		// compute the even pairs
	//		while (pow(2, iStage) - 2) 
	//		{
	//		
	//		}
	//		// odd pairs are just the even pairs + 1
	//	}


	//}
}

void FFT_lib::two_pointDFT(std::complex<double> a, std::complex<double> b, std::complex<double>* A, std::complex<double>* B)
{
	const int N = 8;

	//std::complex<double> j;
	//j = -1;
	//j =  sqrt(j);
	
	//std::complex<double> W_dood;
	//W_dood = exp(j * (2 * PI / N));

	////std::complex<double> L_dood = std::polar(1.0, -2 * PI / N); // need to add k in the exponent later

	//*A = a + (W_dood * b);
	//*B = a - (W_dood * b);
}
