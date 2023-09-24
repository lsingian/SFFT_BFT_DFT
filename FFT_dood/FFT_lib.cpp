#include "FFT_lib.h"
#include <iostream>
#include <cmath>
#include "math.h"

#define PI 3.14159265358979323846

void FFT_lib::init() 
{
	calculateIndexPairsPerStage();
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

void FFT_lib::calculateIndexPairsPerStage()
{
	// stageButterflyPairs[numStages][N]

	for (int stageIdx = 0; stageIdx < numStages; stageIdx++) 
	{
		uint32_t powpow = pow(2, numStages - stageIdx - 1);
		uint32_t pairIdx = 0; // used for the 2D array that contains the pairs for each stage
		resetButterflyIndicesState();

		for (int i = 0; i < N; i++) 
		{   
			// index has already been paired, no need to pair
			if (butterFlyIndicesState[i]) 
			{
				continue;
			}

			// Find the pair
			// Pair is the index located 2^(numStages - stageIdx - 1) away.
			uint32_t pairIndex = i + butterFlyIndices[powpow];

			stageButterflyPairs[stageIdx][pairIdx++] = butterFlyIndices[i];
			butterFlyIndicesState[i] = true;
			stageButterflyPairs[stageIdx][pairIdx++] = pairIndex;
			butterFlyIndicesState[i + powpow] = true; 

			//stageButterflyPairs[stageIdx][i] = butterFlyIndices[i];
			//butterFlyIndices[i] = 0;
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
	// calculate the twiddle factors - add this in the init function

	for (int iStage = 0; iStage < numStages; iStage++)
	{
		// get values of pairs for bit-reverse indices

		// calculate two point FFT for each pair

	}
}


void FFT_lib::calculateFFT() 
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

	twiddle j;
	j = -1;
	j =  sqrt(j);

	twiddle W_dood;
	W_dood = exp(j * (2 * PI / N));

	//twiddle L_dood = std::polar(1.0, -2 * PI / N); // need to add k in the exponent later

	*A = a + (W_dood * b);
	*B = a - (W_dood * b);
}
