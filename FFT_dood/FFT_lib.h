#pragma once
/*
	Library that implements Decimation-In_Frequency (DIF) radix-2 FFT algorithm
	Author: Larry Singian
	Date: 9/23/2023
*/
#include <complex>

class FFT_lib
{
	private:
		// uint32_t N; // number of points in DFT
		void calculateButterflyPairsPerStage();
		void calculateTwiddlePerStage();
		void resetButterflyIndicesState();
		uint32_t bitReverse(uint32_t inputNum);


	public:
		void two_pointDFT(std::complex<double>* samples, uint32_t A_index, uint32_t B_index);

		//typedef std::complex<double> twiddle;

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

		FFT_lib(double* realSamples, uint32_t N)
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
				// TODO eventually add imaginary samples
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

			// C way of allocating memory...TODO delete later
			//butterFlyIndices = (uint32_t*) malloc(this->N * sizeof(uint32_t)); // TODO address the malloc later
			//butterFlyIndicesState = (bool*) malloc(this->N * sizeof(bool));
			//stageButterflyPairs = (uint32_t**) malloc(numStages * sizeof(uint32_t) * this->N * sizeof(uint32_t));

			// indexArray = (butterflyIndex*)malloc(this->N * sizeof(butterflyIndex));

			for (int i = 0; i < this->N; i++) 
			{
				butterFlyIndices[i] = i;
				butterFlyIndicesState[i] = false;

				//butt.indexArray[i].index = (uint32_t*) malloc(sizeof(uint32_t)); // TODO address the mallocs later
				//butt.indexArray[i].index = 0;
				//butt.indexArray[i].paired = false;
			}


			//butt.A = (double*) malloc(N); // TODO address malloc later
			//butt.B = (double*) malloc(N);
		}

		~FFT_lib()
		{
			// TODO free up memory allocated!! not needed for embedded software, but just good practice.
			delete[] butterFlyIndices;
			delete[] butterFlyIndicesState;

			// delete each row
			for (int i = 0; i < numStages; i++) 
			{
				delete[] stageButterflyPairs[i];
			}

			delete[] stageButterflyPairs;
		}

		void init();
		void testFunc();
		std::complex<double>* calculateFFT();
};