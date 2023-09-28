/*
	Library that implements Decimation-In_Frequency (DIF) radix-2 FFT algorithm
	Author: Larry Singian
	Date: 9/23/2023
*/
#include <complex>
#include "math.h"

class FFT_lib
{
	public:
		enum class FFTAlgorithm_Types
		{
			INVALID_ALGO = -1,
			DECIMATION_IN_TIME,
			DECIMATION_IN_FREQUENCY,
			MAX_ALGO
		};

	private:
		const uint32_t MAX_N = 2048; 

		uint32_t N;	// N-point FFT
		uint32_t numStages; // stages in the FFTs
		FFTAlgorithm_Types algorithm;
		std::complex<double>* inputSamples; // stores the samples that the FFT will be performed on 
		double* magnitude; // stores the magnitude result of calculateMagnitude;

		// TODO - using these two arrays for now, switch to struct later
		uint32_t* butterFlyIndices;
		bool* butterflyIndexHasPair; // false if index is not paired yet, true if its already been paired

		uint32_t** stageButterflyPairs; // contains the indices of the Butterfly pairs for each stage. Rows are the stage, Columns are the indices.
		std::complex<double>** W_table; // contains the twiddle factors each stage

		void calculateButterflyPairsPerStage();
		void calculateTwiddlePerStage();
		void computeButterfly_DIF(std::complex<double>* samples, uint32_t A_index, uint32_t B_index, std::complex<double> W);
		void computeButterfly_DIT(std::complex<double>* samples, uint32_t A_index, uint32_t B_index, std::complex<double> W);
		void resetbutterflyIndexHasPair();
		template <typename T> void bitReverse(T* inputArray);

	public:
		void init();
		std::complex<double>* calculateFFT(double* realSamples);
		double* calculateMagnitude(std::complex<double>* inputSample);

		FFT_lib(uint32_t N, FFTAlgorithm_Types algorithm)
		{
			// limit the FFTs to size MAX_N for now. Future implementation should remove this restriction
			if (N > MAX_N) 
			{
				this->N = MAX_N;
			}
			else 
			{
				this->N = N;
			}

			if (FFTAlgorithm_Types::INVALID_ALGO < algorithm && algorithm < FFTAlgorithm_Types::MAX_ALGO) 
			{
				this->algorithm = algorithm;
			}
			else 
			{
				this->algorithm = FFTAlgorithm_Types::DECIMATION_IN_TIME;
			}

			// TODO if realSamples is not a power of 2, then pad samples with zeros

			numStages = round(log2(this->N));

			inputSamples = new std::complex<double>[this->N];
			magnitude = new double[this->N];
			butterFlyIndices = new uint32_t[this->N];
			butterflyIndexHasPair = new bool[this->N];
			stageButterflyPairs = new uint32_t*[numStages]; // allocate memory for each row
			W_table = new std::complex<double>*[numStages]; // allocate memory for each row

			// allocate memory for each column of stageButterFlyPairs. N columns
			for (int i = 0; i < numStages; i++)
			{
				stageButterflyPairs[i] = new uint32_t[this->N];
			}

			// allocate memory for each column of W_table. N/2 columns
			for (int i = 0; i < numStages; i++)
			{
				W_table[i] = new std::complex<double>[this->N/2];
			}

			// butterFlyIndices are simply the array indices - TODO this array will not be needed in the future.
			for (int i = 0; i < this->N; i++) 
			{
				butterFlyIndices[i] = i;
				butterflyIndexHasPair[i] = false;
			}
		}

		~FFT_lib()
		{
			//  Free up memory allocated!! not needed for embedded software if malloc occurs during initialization, but just good practice.
			delete[] inputSamples;
			delete[] magnitude;
			delete[] butterFlyIndices;
			delete[] butterflyIndexHasPair;

			// delete each row
			for (int i = 0; i < numStages; i++) 
			{
				delete[] stageButterflyPairs[i];
				delete[] W_table[i];
			}

			delete[] stageButterflyPairs;
			delete[] W_table;
		}
};