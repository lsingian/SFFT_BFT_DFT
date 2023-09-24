#pragma once
/*
	Library that implements radix-2 FFT algorithm
*/
#include <complex>

class FFT_lib
{
	private:
		// uint32_t N; // number of points in DFT
		void calculateIndexPairsPerStage();
		void calculateTwiddlePerStage();
		void resetButterflyIndicesState();

	public:
		void two_pointDFT(std::complex<double> a, std::complex<double> b, std::complex<double>* A, std::complex<double>* B);

		typedef std::complex<double> twiddle;

		struct butterflyIndex
		{
			uint32_t index;
			bool paired;
		};

		butterflyIndex* indexArray;

		// TODO - using these two arrays for now, switch to struct later
		uint32_t* butterFlyIndices;
		bool* butterFlyIndicesState; // false if index is not paired yet, true if its already been paired

		uint32_t** stageButterflyPairs; // contains the Butterfly pairs for each stage. rows are the stage, columns are the indices.

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

			currentStage = 0;
			numStages = log2(this->N);

			butt.numButterflyElements = this->N;
			butt.A = realSamples;
			butt.B = realSamples;
			//butt.indexArray = (butterflyIndex*) malloc(this->N * sizeof(butterflyIndex));

			butterFlyIndices = new uint32_t[this->N]; // TODO address this allocation later
			butterFlyIndicesState = new bool[this->N];
			stageButterflyPairs = new uint32_t*[numStages]; // allocate memory for each row

			// allocate memory for each column of stageButterFlyPairs
			for (int i = 0; i < numStages; i++)
			{
				stageButterflyPairs[i] = new uint32_t[this->N];
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
			// TODO free up memory allocated!! probably not needed for embedded software, but just good practice.
			delete[] butterFlyIndices;
			delete[] butterFlyIndicesState;

		}

		void init();
		void testFunc();
		void calculateFFT();








};