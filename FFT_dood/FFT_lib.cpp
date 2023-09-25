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
	}
}

void FFT_lib::calculateButterflyPairsPerStage()
{
	// row iteration of the 2D array 
	for (int stageIdx = 0; stageIdx < numStages; stageIdx++) 
	{
		// The butterfly pair is the index located 2^(numStages - stageIdx - 1) away.
		// Ex. 8-point FFT: N = 8 = 2^4 
		// For the first stage (note Stages use 1-based indexing), 
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
	// const uint32_t t_kMax = t_N / 2; //  kmax for decimation in time FFT.

	//std::complex<double> W_table[t_numStages][t_N / 2]; // N/2 twiddle factors per stage

	#define DEBUG
	#ifdef DEBUG
		const uint32_t t_N = 8;
		const uint32_t t_numStages = t_N / 2;
		uint32_t k_test_table[t_numStages][t_N / 2]; // DEBUG only N/2 twiddle factors per stage
	#endif

	// row iteration of twiddle factor 2D array
	for (int stageIdx = 0; stageIdx < numStages; stageIdx++)
	{
		uint32_t kMax = (N / pow(2,stageIdx + 1)); // kmax is the number of UNIQUE twiddle factors. It is calculated every stage for decimation in frequency FFT.
		
		// decimation in frequency k-th Twiddle factor is simply, k * (2^P / 2)
		uint32_t kFactor = pow(2, stageIdx + 1) / 2;
		
		// TODO clean this uip by switching the two for loops (kIndex and i loops)

		// column iteration of twiddle factor 2D array
		for (int kIndex = 0; kIndex < kMax; kIndex++) 
		{
			// TODO debug when switching to decimation in time FFT algo, for now, using decimation in frequency
			//uint32_t tempNum = (kp * pow(2, stageIdx + 1)) / t_N;
			//uint32_t k = bitReverse(tempNum); // TODO - bit reverse
			//std::complex<double> W_k = std::polar(1.0, -2 *  PI * k / N);
			//W_table[stageIdx][kIndex] = W_k;

			uint32_t k = kIndex * kFactor;
			std::complex<double> W_k = std::polar(1.0, -2 * PI * k / N);
			//std::complex<double> W_k(cos(2 * PI * k / N), -sin( 2 * PI * k / N));

			// populate N/2 columns of 2D array with twiddle factors
			for (int i = 0; i < N / 2; i++)
			{
				if (kIndex + (i * kMax) < N / 2) // ensure that we don't access out of array indices
				{
					#ifdef DEBUG
						k_test_table[stageIdx][kIndex + (i * kMax)] = k; // twiddle factor at kIndex repeats every kMax = N/(2^P) indices
					#endif
					W_table[stageIdx][kIndex + (i * kMax)] = W_k; // twiddle factor at kIndex repeats every kMax = N/(2^P) indices
				}
				else
				{
					break; // succeeding i * t_kMax will cause array index access to go out of bounds
				}
			}
		}
	}

	//int test1 = 0;

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


std::complex<double>* FFT_lib::calculateFFT() // TODO add samples as input here
{

	// calculate bit-reversed indices pairs for each stage - add this in the init function
	// 
	const uint32_t maxStageIter = numStages;
	const uint32_t maxEvenIndex = N - 2;
	const uint32_t maxOddIndex = N - 1;
	const uint32_t butterfliesPerStage = N / 2;

	for (uint32_t stageIdx = 0; stageIdx < numStages; stageIdx++)
	{
		// sizes:
		// stageButterflyPairs[numStages][N]
		// W_table[numStages][N/2]
		// inputSamples[N];

		uint32_t twiddleIdx = 0;

		// iterate through the butterfly pairs for the current stage
		for (int i = 0; i < butterfliesPerStage * 2; i++)
		{
			// determine indices of the butterfly pairs
			uint32_t A_index = stageButterflyPairs[stageIdx][i++];
			uint32_t B_index = stageButterflyPairs[stageIdx][i];

			two_pointDFT(inputSamples, A_index, B_index);

			// multiply the "B" output of the butterfly by the twiddle factor
			inputSamples[B_index] = inputSamples[B_index] * W_table[stageIdx][twiddleIdx++];
		}
	}

	return inputSamples;
}

void FFT_lib::two_pointDFT(std::complex<double>* samples, uint32_t A_index, uint32_t B_index)
{
	// grab the values of the indices
	std::complex<double> x_A = samples[A_index];
	std::complex<double> x_B = samples[B_index];

	// TODO see if we can do this in place
	samples[A_index] = x_A + x_B;
	samples[B_index] = x_A - x_B;



	//std::complex<double> j;
	//j = -1;
	//j =  sqrt(j);
	
	//std::complex<double> W_dood;
	//W_dood = exp(j * (2 * PI / N));

	////std::complex<double> L_dood = std::polar(1.0, -2 * PI / N); // need to add k in the exponent later

	//*A = a + (W_dood * b);
	//*B = a - (W_dood * b);
}
