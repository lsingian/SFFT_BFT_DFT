#pragma once
/*
	Library that implements Decimation-In_Frequency (DIF) radix-2 FFT algorithm
	Author: Larry Singian
	Date: 9/23/2023
*/
#include <complex>

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
		FFTAlgorithm_Types algorithm;
		void calculateButterflyPairsPerStage();
		void calculateTwiddlePerStage();
		void resetButterflyIndicesState();
		void bitReverse(std::complex<double>* inputArray);


	public:
		void init();
		std::complex<double>* calculateFFT();
		void computeButterfly_DIF(std::complex<double>* samples, uint32_t A_index, uint32_t B_index, std::complex<double> W);
		void computeButterfly_DIT(std::complex<double>* samples, uint32_t A_index, uint32_t B_index, std::complex<double> W);

		struct butterflyIndex
		{
			uint32_t index;
			bool paired;
		};

		butterflyIndex* indexArray;

		// TODO - using these two arrays for now, switch to struct later
		uint32_t* butterFlyIndices;
		bool* butterFlyIndicesState; // false if index is not paired yet, true if its already been paired
		std::complex<double>* inputSamples;

		uint32_t** stageButterflyPairs; // contains the indices of the Butterfly pairs for each stage. Rows are the stage, Columns are the indices.
		std::complex<double>** W_table; 

		struct butterflyValues
		{
			double* A;
			double* B;
			size_t numButterflyElements;
			// butterflyIndex* indexArray;
		};
	
		butterflyValues butt;
		uint32_t N;
		const uint32_t MAX_N = 2048; 
		uint32_t numStages;
		uint32_t currentStage;

		FFT_lib(double* realSamples, uint32_t N, FFTAlgorithm_Types algorithm)
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

			currentStage = 0;
			numStages = log2(this->N);

			butt.numButterflyElements = this->N;
			butt.A = realSamples;
			butt.B = realSamples;
			//butt.indexArray = (butterflyIndex*) malloc(this->N * sizeof(butterflyIndex));

			inputSamples = new std::complex<double>[this->N];
			butterFlyIndices = new uint32_t[this->N]; // TODO address this allocation later
			butterFlyIndicesState = new bool[this->N];
			stageButterflyPairs = new uint32_t*[numStages]; // allocate memory for each row
			W_table = new std::complex<double>*[numStages]; // allocate memory for each row

			// TODO change to std::copy later
			for (int i = 0; i < this->N; i++)
			{
				inputSamples[i] = realSamples[i];
				// TODO eventually add imaginary samples, real samples are sufficient for now (for audio applications)
			}

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

			for (int i = 0; i < this->N; i++) 
			{
				butterFlyIndices[i] = i;
				butterFlyIndicesState[i] = false;
			}
		}

		~FFT_lib()
		{
			//  Free up memory allocated!! not needed for embedded software, but just good practice.
			delete[] butterFlyIndices;
			delete[] butterFlyIndicesState;

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