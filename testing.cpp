// digitRecognitionV3General.cpp : Defines the entry point for the console application.
//


/*
Model folder naming according to the words
1. digit(0-9) ---> ./resource/final model/model-digit/
2. by(10) --->./resource/final model/model-10/
3. into(11) --->./resource/final model/model-11/
4. minus(12) --->./resource/final model/model-12/
5. plus(13) --->./resource/final model/model-13/


Total words : 
	0-9		----    10
	+,-,*,/ ---		4

-------------------------
	total    ---	14 (words)
*/
//#include "stdafx.h"
#include <stdio.h>
#include <iostream>
#include <limits.h>
#include <conio.h>
#include <stdlib.h>
#include <string.h>
#include <windows.h>
#include <direct.h>

# define READ "r"
# define WRITE "w"
# define APPEND "a"

# define ARRAY_SIZE_L 70000
# define ARRAY_SIZE_S 5000
# define FRAME_SIZE 320
# define WINDOW_SIZE 320
# define TRAIN_SIZE 10

#define P 12
#define Q 12

# define T 85
# define CODEBOOK_SIZE 32

# define DIGITS 14
# define ROLLNO "194101022"

using namespace std;



// Create 2D-Array
long double * * create2DArray(int numOfRow, int numOfCol) {
	int i;
	long double * *p = (long double * *) malloc(sizeof(long double*)*numOfRow); // memory allocation for row
	for(i = 0; i < numOfRow; i++) {
		p[i] = (long double *) malloc(sizeof(long double) * numOfCol); // memory allocation for column
	}
	return p;
}
// Display the content of 1D-array
void display1DArray(long double * arr, int numOfCol) {
	for(int i = 0; i < numOfCol; i++)
		printf("%Lf ",arr[i]);
	printf("\n");
}

// Display the content of 1D-integer array 
void display1DIntArray(int * arr, int numOfCol) {
	for(int i = 0; i < numOfCol; i++)
		printf("%d ",arr[i]);
	printf("\n");
}
// Display content of the 2D-Array
void display2DArray(long double * *p, int numOfRow, int numOfCol) {
	int i, j;
	for(i = 0; i < numOfRow; i++) {
		printf("row-%d: ",i);
		for(j = 0; j < numOfCol; j++) {
			printf("%Lf ",p[i][j]);
		}
		printf("\n");
	}
}
// Display content of the 2D-Array
void scientific_display2DArray(long double * *p, int numOfRow, int numOfCol) {
	int i, j;
	
	for(i = 0; i < numOfRow; i++) {
		printf("row-%d: \n",i);
		for(j = 0; j < numOfCol; j++) {
			//printf(%0.4G,p[i][j]);
			//printf("%0.30Lf ",p[i][j]);
			cout << std::scientific;
			cout << p[i][j]<<" ";
		}
		printf("\n");
	}
}
// Parameter : file name and pointer
// returns the number of samples in the file
int cpFromFileToArray(char * fname, long double * arr) {
	int index = 0;
	FILE * f = fopen(fname, READ);
	while(fscanf(f,"%Lf",&arr[index++]) != EOF);
	fclose(f);
	return index - 1;
}
// copy the content of an array to the file
// Here seperator is ('\n') newline
void cpFromArrayToFile(char * fname, long double * arr, int arr_len) {
	FILE * f = fopen(fname, WRITE);
	for(int i = 0; i < arr_len; i++)
		fprintf(f,"%Lf\n",arr[i]);
	fclose(f);
}
// copy the content of an array to the file
// Here seperator is ('\t') 4 whitespace (tab)
void cpFromArrayToFileV2(char * fname, long double * arr, int arr_len) {
	FILE * f = fopen(fname, WRITE);
	for(int i = 0; i < arr_len; i++)
		fprintf(f,"%Lf\t",arr[i]);
	fclose(f);
}
// Apend the content of an array to the file
void appendArrayToFile(char * fname, long double * arr, int low, int high) {
	FILE * f = fopen(fname, APPEND);
	for(int i = low; i <= high; i++)
		fprintf(f,"%Lf\n",arr[i]);
	fclose(f);
}
// Apend the content(seperated by comman) of an array to the file
void appendArrayToFileV2(char * fname, long double * arr, int low, int high) {
	FILE * f = fopen(fname, APPEND);
	for(int i = low; i <= high - 1; i++)
		fprintf(f,"%Lf,",arr[i]);
	fprintf(f,"%Lf\n",arr[high]);
	fclose(f);
}

// Calculation of DCShift
long double getDCShift(long double * arr, int len) {
	long double shift = 0;
	for(int i = 0; i < len; i++)
		shift += arr[i];
	shift = shift / len;
	return shift;
}
// Perform DC Shift on the data
void shiftData(long double *arr, int len, long double shift) {
	for(int i = 0; i < len; i++)
		arr[i] = arr[i] - shift;	
}
// Get absolute max from the array
long double getAbsMax(long double * arr, int len) {
	long double max = INT_MIN,temp;	
	for(int i = 0; i < len;i++) {
		temp = arr[i];
		if(temp < 0)
			temp = temp * -1;
		if(temp > max)
			max = temp;
	}
	return max;
}
// Normalisation of the data
void normaliseData(long double * arr, int len) {
	long double max = getAbsMax(arr,len);
	for(int i = 0; i < len; i++)
		arr[i] = (arr[i] / max) * 5000;
}
// Find hightest amplitude index
int getHighestAmplitudeIndex(long double * arr, int len) {
	long double max = INT_MIN,temp;
	int index;
	for(int i = 0; i < len;i++) {
		temp = arr[i];
		if(temp < 0)
			temp = temp * -1;
		if(temp > max) {
			max = temp;
			index = i;
		}
	}
	return index;
}
// calculate front and end frame numbers of the data
void getFrontAndEndMarkerOfNFrames(int n, int hightestAmplitudeIndex, int * front, int * end) {
	int n1; // Number of frames to the left of hightestAmplitudeIndex
	int n2; // Number of frames to the right of hightestAmplitudeIndex
	if(n % 2 == 0) { // Number of frames are even
		n1 = n2 = n / 2;
	}
	else { // Number of frames are odd
		n1 = n / 2;
		n2 = n1 + 1;
	}
	printf("n1: %d and n2: %d\n", n1, n2);
	*front = hightestAmplitudeIndex - n1 * (FRAME_SIZE - 240) + 1;
	*end = hightestAmplitudeIndex + n2 * (FRAME_SIZE - 240);
}

//--------------------------------------------------------------Modules taken from yes-no-------------------------------------------------------
// Calculation of energy per frame
long double getEnergyPerFrame(int frame_size, long double * arr, int low, int high) {
	long double e = 0,temp;
	for(int i = low; i <= high; i++){
		temp = arr[i];
		e += temp * temp;
	}
	e = e / frame_size;
	return e;
}
// Calculation of energy
// Returns the number of frames
int storeEnergy(int frame_size, long double * arr, int len_arr, long double * e) {
	int frame_count = 0;
	for(int i = 0; i + frame_size - 1 < len_arr; i = i + (frame_size - 240)) {
		e[frame_count] = getEnergyPerFrame(frame_size,arr,i,i + frame_size - 1);
		frame_count = frame_count + 1;
	}
	return frame_count;
}
// Here Frame Number starts from zero
int getFrontMarker(long double * energy, int arr_len) {
	long double pre_energy, curr_energy;
	pre_energy = energy[0]; 
	int i;
	for(i = 1; i < arr_len; i++) {
		curr_energy = energy[i];
		if(curr_energy > 1.1 * pre_energy)
			break;
		pre_energy = curr_energy;
	}
	if(i == arr_len) // Error checking
		return -1;
	else
		return i;
}
// Calculate the end marker
int getEndMarker(long double * energy, int arr_len) {
	long double pre_energy, curr_energy;
	pre_energy = energy[arr_len - 1];
	int i;
	for(i = arr_len - 2; i>= 0; i--) {
		curr_energy = energy[i];
		if(curr_energy > 1.1 * pre_energy)
			break;
		pre_energy = curr_energy;
	}
	if(i == 0)
		return -1;
	else
		return i;
}

//------------------------------------------Modules Taken From Vowel Recognition-----------------------------------------------------------------

// Calculation of Hamming values
long double * getHamVal() {
	long double * ham = (long double *) malloc (sizeof(long double)* WINDOW_SIZE);
	for(int i = 0; i < WINDOW_SIZE; i++) {
		long double angle =(long double) 2 * 22 * i / (7 * (WINDOW_SIZE - 1));
		ham[i] = 0.54 - 0.46 * (long double)cos(angle);
	}
	return ham;
}
// Applying Hamming window on a frame(data)
long double * applyHamWinOnData(long double * ham,long double * data, int p) {
	long double * dataAfterHam = (long double *) malloc(sizeof(long double) * WINDOW_SIZE);
	for(int i = 0; i < WINDOW_SIZE;i++) {
		dataAfterHam[i] = ham[i] * data[i];
	}
	return dataAfterHam;
}
// Calculation of r_i
long double * getRi(long double * arr, int p) {
	long double * r = (long double *) malloc(sizeof(long double)*(p+1)); // r(0),r(1),....,r(p)
	for(int i = 0; i <= p; i++){
		r[i] = 0;
		for(int j = 0; j <= (WINDOW_SIZE - 1 - i);j++) {
			r[i] += arr[j] * arr[j + i];
		}
	}
	return r;
}
// Calculation of a_i using Durbin's algorithm
long double * getAi(long double * r , int p) {
	long double * a_i = (long double *) malloc (sizeof(long double)* (p + 1));
    long double * energy = (long double *) malloc (sizeof(long double)* (p + 1));
	long double * k = (long double *) malloc (sizeof(long double)* (p + 1));
	long double x;
	long double a[P + 1][P + 1];
    energy[0] = r[0]; // Step 1
    for(int i = 1; i <= p; i++) {
        // Step 2
        x = 0;
        for(int j = 1; j <= i - 1; j++)
           x +=  a[i - 1][j] * r[i - j];
        k[i] = (r[i] - x) / energy[i - 1];
        a[i][i] = k[i]; // Step 3
        for(int j = 1; j <= i - 1; j++) // Step 4
            a[i][j] = a[i - 1][j] - (k[i] * a[i - 1][i - j]);
        energy[i] = (1 - k[i] * k[i]) * energy[i - 1]; // Step 5
    }
    for(int i = 1; i <= p; i++) // Storing the values of ai
        a_i[i] = a[p][i];
    return a_i;
}
// Calculation of c_i
// Note: Indexing of c_i is started from 1
long double * getCi(long double * a, int p) {
	long double * c = (long double *) malloc (sizeof(long double)* (p + 1));
	long double temp;
	c[0] = 0;
	for(int m = 1; m <=p ; m++) {
		temp = 0;
		for(int k = 1; k <= m - 1;k++) {
			temp += ((long double) k / m) * c[k] * a[m - k];
		}
		c[m] = a[m] + temp;
	}
	return c;
}
// Calculation of raised sine values
long double * getRaisedSinVal() {
	// When q = (3/2) * p then sin values varies from 0 to q - 1 
	// So memory required is q instead WINDOW_SIZE times data type
	long double * sine = (long double *) malloc(sizeof(long double) * Q); 
	for(int i = 1; i <= Q; i++) {
		long double angle =(long double) 22 * i / (7 * Q);
		sine[i-1] = 1 + (Q / 2) * (long double) sin(angle);
	}
	return sine;
}
// Applying raised sine window on Ci
long double * applyRaisedSinWinOnCi(long double * sine, long double *c_i) {
	long double * modified_c_i = (long double *) malloc(sizeof(long double)*(Q+1));
	modified_c_i[0] = 0;
	for(int i = 0; i < Q; i++) {
		modified_c_i[i+1] = sine[i] * c_i[i+1]; // Indexing of c_i is started from 0
	} 
	return modified_c_i;
}

//------------------------------------------------------------------------------------------------------------------ 
// Calculate ci for a frame and reture the raised sine ci's
// arr: starting address of the frame 
// p : dimension of the ceptral coefficient
long double * getCiForAFrame(long double * arr, int p) {

		long double * r = getRi(arr, p); // Calculation of Ri
		long double * a_i = getAi(r,p);// Calculation of Ai
		long double * c_i = getCi(a_i,p); // Calculation of c_i
		long double * raised_sin_val = getRaisedSinVal(); // calculation of raised sine value
		long double * raised_sine_ci = applyRaisedSinWinOnCi(raised_sin_val,c_i); // Appyling raised sine values on Ci's

		free(r);
		free(a_i);
		free(c_i);
		free(raised_sin_val);
		return raised_sine_ci;
}
// Calculation of ci's for single utterance of a digit which contains n frames
// Indexing of ci's start from 1
// Here Frames are overlapping ones
long double * *getCiOfTFrames(long double * arr, int t, int start, int p) {
	long double * initAddr = arr;
	long double * ham = getHamVal();
	long double * *tCiVector = create2DArray(t, p + 1); // indexing of ci's starts from 1
	for(int i = 0; i < t ; i++) {
		initAddr = arr + start + i * (FRAME_SIZE - 240);
		initAddr = applyHamWinOnData(ham,initAddr, p); // applying hamming window
		tCiVector[i] = getCiForAFrame(initAddr, p);
		//display1DArray(tCiVector[i], p + 1);
	}
	return tCiVector;
}

//----------------------------------------------------------Modules taken from BinarySplitAlgorithm-----------------------------
// Convert CSV file into 2D-Array
long double ** convertCSVTo2DArray(char * fcsv, int numOfRow, int numOfCol, int * numOfActualRow) {
	int count = 0, numOfVector = 0;
	long double ** arr = create2DArray(numOfRow, numOfCol); // Memory allocation for 2D-Array
	FILE * f = fopen(fcsv, READ);
	while(fscanf(f,"%Lf,",&arr[numOfVector][count]) != EOF) {
		count = count + 1;
		if(count == P) {
			count = 0;
			numOfVector = numOfVector + 1;
		}
	}
	fclose(f);
	*numOfActualRow = numOfVector;
	return arr;
}
// convert codebook in .txt to 2-D Array
long double ** getCodebook(char * f_txt, int numOfRow, int numOfCol, int * numOfActualRow) {
	int count = 0, numOfVector = 0;
	long double ** arr = create2DArray(numOfRow, numOfCol); // Memory allocation for 2D-Array
	FILE * f = fopen(f_txt, READ);
	while(fscanf(f,"%Lf",&arr[numOfVector][count]) != EOF) {
		count = count + 1;
		if(count == P) {
			count = 0;
			numOfVector = numOfVector + 1;
		}
	}
	fclose(f);
	*numOfActualRow = numOfVector;
	return arr;
}
//------------------------------------------------------------Modules taken from HMM---------------------------------------------
// Read initial model from the files
long double ** *getInitialModelFromFile(char * fileNameOfPIMatrix, char * fileNameOfAMatrix, char * fileNameOfBMatrix, int n, int m) {

	int i, j;

	FILE * f_pi = fopen(fileNameOfPIMatrix, READ);
	FILE * f_a = fopen(fileNameOfAMatrix, READ);
	FILE * f_b = fopen(fileNameOfBMatrix, READ);

	long double ** *model = (long double ***) malloc(sizeof(long double **) * 3);
	model[0] = create2DArray(1, n); // Memory allocation for PI matrix
	model[1] = create2DArray(n, n); // Memory allocation for A matrix
	model[2] = create2DArray(n, m); // Memory allcoation for B matrix

	
	i = 0;
	while(fscanf(f_pi, "%Lf", &model[0][0][i++]) != EOF); // store PI values from file to array

	i = 0;
	j = 0;
	while(fscanf(f_a, "%Lf", &model[1][i][j]) != EOF) { // store A matrix values from file to array
		j = j + 1;
		if(j == n) {
			i = i + 1;
			j = 0;
		}
	}

	i = 0;
	j = 0;
	while(fscanf(f_b, "%Lf", &model[2][i][j]) != EOF) { // store B matrix values from file to array
		j = j + 1;
		if(j == m) {
			i = i + 1;
			j = 0;
		}
	}

	fclose(f_pi);
	fclose(f_a);
	fclose(f_b);

	return model;
}
// Calculation of forward variable
// It returns probability of obseravation sequence given a model
long double calculateForwardVariable(long double * *alpha, long double ** *model, int *obsSeq, int n, int t) {
	long double probOfObsSeqGivenModel = 0, temp;
	int i, j, time_t;

	// Initialisation
	for(i = 0; i < n; i++) { // At time t = 1
		alpha[i][0] = model[0][0][i] * model[2][i][obsSeq[0] - 1];
	}

	// Induction
	for(time_t = 0; time_t <= t - 2; time_t++) { // For time t = 2 to t = t 
		for(j = 0; j < n; j++) { // For each state 
			temp = 0;
			for(i = 0; i < n; i++) {
				temp += alpha[i][time_t] * model[1][i][j];
			}
	//...........................................................ERROR PRONE.........................................................
			alpha[j][time_t + 1] = temp * model[2][j][obsSeq[time_t + 1] - 1];
		}
	}

	// Termination
	probOfObsSeqGivenModel = 0;
	for(i = 0; i < n; i++)
		probOfObsSeqGivenModel += alpha[i][t - 1];
	return probOfObsSeqGivenModel;
}
// Calculation of backward variable
long double calculateBackwardVariable(long double * *beta, long double ** *model, int *obsSeq, int n, int t) {
	int i, j, time_t;
	long double temp = 0;

	// Initialization
	for(i = 0; i < n; i++) // At time t = t
		beta[i][t - 1] = 1; 
	
	// Induction
	for(time_t = t - 2; time_t >= 0; time_t--) { // For time t = t - 1 to t = 1
		for(i = 0; i < n; i++) { // For each state
			temp = 0;
			for(j = 0; j < n; j++) {
				temp += model[1][i][j] * model[2][j][obsSeq[time_t + 1] - 1] * beta[j][time_t + 1];
			}
			beta[i][time_t] = temp;
		}
	}

	return 0;
}
// Calculation of Gamma variable
// It returns the array of state sequence
int * calculateGamma(long double * *alpha, long double * *beta, long double * *gamma, long double probOfObsSeqGivenModel, int n, int t) {
	int * stateSeq = (int *) malloc(sizeof(int ) * t);
	int max_state = -1, time_t, i, j;
	long double temp_y = 0, temp_prob;

	for(time_t = 0; time_t < t; time_t++) { // At each time instance
		max_state = -1;
		temp_y = 0;
		for(i = 0; i < n; i++) { // For each state
			//gamma[i][time_t] = alpha[i][time_t] * beta[i][time_t] / probOfObsSeqGivenModel; // calculate gamma
			//------------------------------------------------------------------------------------------------
			temp_prob = 0;
			for(j = 0; j < n; j++) {
				temp_prob += alpha[j][time_t] * beta[j][time_t];
			}
			gamma[i][time_t] = alpha[i][time_t] * beta[i][time_t] / temp_prob; // calculate gamma
			//-------------------------------------------------------------------------------------------------
			if(gamma[i][time_t] > temp_y) {
				temp_y = gamma[i][time_t];
				max_state = i + 1;
			}
			stateSeq[time_t] = max_state; // State which has highest probability for particular instance of time
		}
	}


	return stateSeq;
}
// Calculation of delta variable using Viterbi Algorithm 
int * calculateDelta(long double ** *model, int * obsSeq, long double * *delta, long double * *psi, int n, int t, long double * pstar) {
	int * stateSeq = (int *) malloc(sizeof(int ) * t);
	int i, j, time_t, max_state = -1;
	long double temp_d = 0, temp = 0;

	// Initialization
	for(i = 0; i < n; i++) {
		delta[i][0] = model[0][0][i] * model[2][i][obsSeq[0] - 1];
		psi[i][0] = 0;
	}


	// Recursion
	for(time_t = 1; time_t < t; time_t++) {
		for(j = 0; j < n; j++) {
			temp_d = 0; 
			max_state = -1;
			for(i = 0; i < n; i++) {
				temp = delta[i][time_t - 1] * model[1][i][j];
				if(temp >= temp_d) {
					temp_d = temp;
					max_state = i + 1;
				}
			}
			delta[j][time_t] = temp_d * model[2][j][obsSeq[time_t] - 1];
			psi[j][time_t] = max_state;
		}
	}

	// Termination
	temp_d = 0;
	max_state = -1;

	for(i = 0; i < n; i++) {
		if(delta[i][t - 1] >= temp_d) {
			temp_d = delta[i][t - 1];
			max_state = i + 1;
		}
	}
	stateSeq[t - 1] = max_state; // And temp_d is the max state probability

	// State sequence (or path) backtracking
	for(time_t = t - 2; time_t >= 0; time_t--) {
		stateSeq[time_t] = psi[stateSeq[time_t + 1] - 1][time_t + 1];
	}

	printf("p star : %0.300Lf\n q star T: %d\n", temp_d, max_state);
	*pstar = temp_d;
	printf("State Sequence: \n");
	for(i = 0; i < t; i++)
		printf("%d ", stateSeq[i]);
	printf("\n\n");

	return stateSeq;
}
// Calculation of xi variable
long double ** *calculateXi(long double ** *model, int *obsSeq,long double * *alpha, long double * *beta, long double probOfObsSeqGivenModel, int n, int t) {
	long double ** *xi = (long double ** *) malloc(sizeof(long double * *) * t - 1);
	int i, j, time_t;
	long double temp_prob = 0;

	for(time_t = 0; time_t < t - 1; time_t++) {
		xi[time_t] = create2DArray(n, n);
	}

	for(time_t = 0; time_t < t - 1; time_t++) {
		for(i = 0; i < n; i++) {
			for(j = 0; j < n; j++) {

				xi[time_t][i][j] = (alpha[i][time_t] * model[1][i][j] * model[2][j][obsSeq[time_t + 1] - 1] * beta[j][time_t + 1]) / probOfObsSeqGivenModel;
			}
		}

	}

	


	return xi;
}
// Re-estimation of the model
long double ** * reEstimateModel(long double ** *model, int * obsSeq, long double * *gamma, long double ** *xi, int n, int m, int t) {
	int i, j, k, time_t;
	long double expNoOfTransitionFromSiToSj, expNoOfTransitionFromSi;
	long double expTimeSjAndObsSymK, expTimeSj;
	long double ** *reEstimatedModel = (long double ***) malloc(sizeof(long double **) * 3);
	reEstimatedModel[0] = create2DArray(1, n); // Memory allocation for PI matrix
	reEstimatedModel[1] = create2DArray(n, n); // Memory allocation for A matrix
	reEstimatedModel[2] = create2DArray(n, m); // Memory allcoation for B matrix




	for(i = 0; i < n; i++) { // Re-estimation of PI matrix
		reEstimatedModel[0][0][i] = gamma[i][0]; 
	}

	for(i = 0; i < n; i++) { // Re-estimation of A matrix
		for(j = 0; j < n; j++) {
			expNoOfTransitionFromSiToSj = 0; 
			expNoOfTransitionFromSi = 0;
			for(time_t = 0; time_t < t - 1; time_t++) {
				expNoOfTransitionFromSiToSj += xi[time_t][i][j];
				expNoOfTransitionFromSi += gamma[i][time_t];
			}
			reEstimatedModel[1][i][j] = expNoOfTransitionFromSiToSj / expNoOfTransitionFromSi; 
		}
	}

	for(k = 0; k < m; k++){ // Re-estimation of B matrix
		for(j = 0; j < n; j++) {
			expTimeSjAndObsSymK = 0; 
			expTimeSj = 0;
			for(time_t = 0; time_t < t; time_t++) {
				if((obsSeq[time_t] - 1) == k)
					expTimeSjAndObsSymK += gamma[j][time_t];
				expTimeSj += gamma[j][time_t];
			}
			reEstimatedModel[2][j][k] = expTimeSjAndObsSymK / expTimeSj;
		}
	}
	return reEstimatedModel;
}
// Modification of B matrix
void modifyBMatrix(long double ** *model, int n, int m, long double threshold) {
	int i, j, max_row, max_col, zeroCount;
	long double max = 0, rowSum, temp;
	//printf("Threshold: %0.300Lf\n",threshold);


	for(i = 0; i < n; i++) {
		zeroCount = 0;
		max = 0.0;
		rowSum = 0.0;
		temp = 0;
		for(j = 0; j < m; j++) {
			if(model[2][i][j] >= max) {
				max = model[2][i][j];
				max_row = i;
				max_col = j;
			}
			if(model[2][i][j] <= threshold) {
				temp += model[2][i][j];
				model[2][i][j] = threshold;
				zeroCount += 1;
				//printf("row:%d column:%d modified Value: %0.30Lf\n", i , j, model[2][i][j]);
				//rowSum += model[2][i][j];
			}
			rowSum += model[2][i][j];
		}

		model[2][max_row][max_col] -= (threshold * zeroCount - temp); 

		
	}


}






// Display the Model
void displayModel(long double ** *model, int n, int m) {
	printf("PI matrix:\n");
	//getch();
	display1DArray(model[0][0], n);
	//getch();
	printf(" A matrix:\n");
	//getch();
	scientific_display2DArray(model[1], n, n);
	//getch();
	printf(" B matrix:\n");
	//getch();
	scientific_display2DArray(model[2], n, m);
	//getch();
}
// Calculation of Tokhura's Distance
// Indexing of ci_ref starts from 0 where as indexing of ci_test starts from 1
// Note : precision is not matching with the data in excel sheet
long double getTokhuraDistanceV2(long double * ci_ref,long double * ci_test, long double * tokhuraWeights) {
	long double dist = 0, diff;
	for(int i = 1; i <= Q; i++) {
		diff = ci_test[i] - ci_ref[i - 1]; // indexing of ci_ref starts from 0 where as indexing of ci_test starts from 1
		dist += diff * diff * tokhuraWeights[i];
	}
	return dist;
}
// Generate Observation Sequence
int * generateObsSeq(long double * *codeBook,int numOfRowsInCodeBook,long double * *tCiVector, int t, int p) {
	long double min, dist;
	long double tokhuraWeights[13] = {0, 1, 3, 7, 13, 19, 22, 25, 33, 42, 50, 56, 61}; // indexing starts from 1
	int * obsSeq = (int *) malloc(sizeof(int) * t);
	int min_index;
	for(int i = 0; i < t; i++) { // For each ci vector
		min = INT_MAX;
		for(int j = 0; j < numOfRowsInCodeBook; j++) {
			// index of codeBook[j] starts from 0 where as index of tCiVector starts from 1
			// that's why getTokhuraDistance() is modified into getTokhuraDistanceV2()
			dist = getTokhuraDistanceV2(codeBook[j], tCiVector[i], tokhuraWeights); 
			if(dist < min) {
				min = dist;
				min_index = j;
			}
		}
		obsSeq[i] = min_index + 1; // obsSeq contains values lies between 1 to numOfRowsInCodeBook;
	}
	return obsSeq;
}
long double ** *getModelUsingHMM(long double ** *model, int * obsSeq, int n, int m, int t) {
	long double * *alpha = create2DArray(n, t); // Forward Variable
	long double * *beta = create2DArray(n, t); // Backward Variable
	long double * *gamma= create2DArray(n, t); // Gamma variable
	long double * *delta = create2DArray(n, t); // Delta variable
	long double * *psi = create2DArray(n, t); // Psi variable

	long double ** *xi; // xi variable 
	long double ** *reEstimatedModel; // re-estimated model will be obtained from solution of problem-3

	int *stateSeq;

	long double pstar = 0; // probability that shows quality of the model
	long double probOfObsSeqGivenModel = 0; // Probability of observation sequence given a model
	long double Next_probOfObsSeqGivenModel = 0; // Next Probability of observation sequence given a model 
	int count = 0;

	long double threshold; // Threshold to adjust B matrix
	long double base = 10;
	long double exponent = -30;

	threshold = pow(base,exponent); 
	count = 0;
	while(1) {
		count = count + 1;

		//-----------------------------------------------------Solution of Problem 1---------------------------------------------//
		probOfObsSeqGivenModel = calculateForwardVariable(alpha, model, obsSeq, n, t); // Calculation of forward (alpha) variable
		printf("Iteration : %d \nP(O/model):%0.200Lf\n\n", count, probOfObsSeqGivenModel);

		calculateBackwardVariable(beta, model, obsSeq, n, t); // Calculation of backward (beta) variable
		//printf("Backward Variable: \n");
		//display2DArray(beta, n, t);
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		//-----------------------------------------------------Solution of Problem 2---------------------------------------------//
		stateSeq = calculateGamma(alpha, beta, gamma, probOfObsSeqGivenModel, n, t); // Calculation of gamma variable
		//printf("Gamma Variable: \n");
		//display2DArray(gamma, n, t);
		stateSeq = calculateDelta(model, obsSeq, delta, psi, n, t, &pstar); // Calculation of delta variable

		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		//-----------------------------------------------------Solution of Problem 3---------------------------------------------//
		xi = calculateXi(model, obsSeq, alpha, beta, probOfObsSeqGivenModel, n, t); // Calculation of xi variable
		reEstimatedModel = reEstimateModel(model, obsSeq, gamma, xi, n, m, t); // Calculation of re-estimate model
		//printf("Modified B matrix:\n");
		//display2DArray(reEstimatedModel[2], n, m);
		modifyBMatrix(reEstimatedModel, n, m, threshold); // Modifiy the B matrix i.e. give probs for zero probs
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


		if(count == 20)
			break;
		model = reEstimatedModel;
	}

	return model;
}
// Calculate avgModel from models of all utterances of a digit
long double ** *getAvgModel(long double *** *modelsForADigit,int numOfUtterances, int n, int m) {
	long double ** *avgModel = (long double ***) malloc(sizeof(long double **) * 3);
	avgModel[0] = create2DArray(1, n); // Memory allocation for PI matrix
	avgModel[1] = create2DArray(n, n); // Memory allocation for A matrix
	avgModel[2] = create2DArray(n, m); // Memory allcoation for B matrix

	// Initialization of the avg model
	for(int i = 0; i < n; i++) // Initialisation of pi matrix
		avgModel[0][0][i] = 0;

	for(int i = 0; i < n; i++) // Initialisation of A matrix
		for(int j = 0; j < n; j++)
			avgModel[1][i][j] = 0;

	for(int i = 0; i < n; i++) // Initialisation of B matrix
		for(int j = 0; j < m; j++)
			avgModel[2][i][j] = 0;

	// Above initialisation can be done using only two loops => Do this in optimisation phase(It will be helpful in the final project)
	// Same for the avg 

	int temp_m_n = (n > m?n:m); // max(n,m)
	for(int utterance = 1; utterance <= numOfUtterances; utterance++) {
		// Calculation of pi matrix of avgModel (No need to calculate)
		for(int i = 0; i < n; i++) {
			avgModel[0][0][i] += modelsForADigit[utterance - 1][0][0][i]; // Calculation of pi matrix of avgModel (No need to calculate)
			for(int j = 0 ; j < temp_m_n; j++) {
				if(j < n)
					avgModel[1][i][j] += modelsForADigit[utterance - 1][1][i][j]; // Calculation of A matrix of avgModel
				if(j < m)
					avgModel[2][i][j] += modelsForADigit[utterance - 1][2][i][j]; // Calculation of B matrix of avgModel
			}
		}
	}

	// divide all the sum values by number of values
	for(int i = 0; i < n; i++) {
		avgModel[0][0][i] /= numOfUtterances;
		for(int j = 0 ; j < temp_m_n; j++) {
			if(j < n)
				avgModel[1][i][j] /= numOfUtterances; // Calculation of A matrix of avgModel
			if(j < m)
				avgModel[2][i][j] /= numOfUtterances; // Calculation of B matrix of avgModel
		}
	}
	return avgModel;
}
// Find the model with highest probability of P(O/model)
int getModelWithHighProb(long double *** *modelsForAllDigit, int startDigit, int endDigit ,int numOfDigits,int * obsSeq, int n,int t) {
	int recognizedDigit = -1;
	long double * *alpha = create2DArray(n, t); // Forward Variable
	long double probOfObsSeqGivenModel = 0; // Probability of observation sequence given a model
	long double maxProb = 0;
	printf("Number of Digit: %d\n", numOfDigits);
	for(int i = startDigit; i <= endDigit; i++) {
		//printf("Digit: %d\n", i);
		probOfObsSeqGivenModel = calculateForwardVariable(alpha, modelsForAllDigit[i], obsSeq, n, t); // Calculation of forward (alpha) variable
		//cout << std::scientific;
		//cout << "Prob: " << i <<" "<<probOfObsSeqGivenModel << endl;
		//printf("%0.3Lf\n", probOfObsSeqGivenModel);
		if(probOfObsSeqGivenModel > maxProb) {
			maxProb = probOfObsSeqGivenModel;
			recognizedDigit = i;
		}
	}

	return recognizedDigit;
}
// copy a model to the file
// Parameters: model, path to the folder in which model is going to be copied, number of states, number of observations(or codebook size)
// Model is copied in three seperated files namely PI.txt, A.txt, B.txt
void cpModelToFile(long double ** *model, char * path, int n, int m) {
	char fileName[100];
	strcpy(fileName, path);
	strcat(fileName, "PI.txt");
	printf("%s\n", fileName);
	//cpFromArrayToFileV2(fileName, model[0][0], n);
	// copy content of the PI matrix to the PI.txt file
	FILE * f = fopen(fileName, WRITE);
	for(int i = 0; i < n; i++) 
		fprintf(f,"%Lf ",model[0][0][i]);
	fclose(f);

	strcpy(fileName, path);
	strcat(fileName, "A.txt");
	printf("%s\n", fileName);
	// copy content of the A matrix to the A.txt file
	f = fopen(fileName, WRITE);
	for(int i = 0; i < n; i++) {
		for(int j = 0; j < n; j++) {
			fprintf(f,"%0.300Lf ",model[1][i][j]);
		}
		fprintf(f,"\n",1);
	}
	fclose(f);

	strcpy(fileName, path);
	strcat(fileName, "B.txt");
	printf("%s\n", fileName);

	// copy content of the B matrix to the B.txt file
	f = fopen(fileName, WRITE);
	for(int i = 0; i < n; i++) {
		for(int j = 0; j < m; j++) {
			fprintf(f,"%0.300Lf ",model[2][i][j]);
		}
		fprintf(f,"\n",1);
	}
	fclose(f);
}
long double ** *cpFromFileToModel(char * path, int n, int m) {
	char fileName[100];
	long double ** *model = (long double ***) malloc(sizeof(long double **) * 3);
	int i, j;
	FILE * f;
	model[0] = create2DArray(1, n); // Memory allocation for PI matrix
	model[1] = create2DArray(n, n); // Memory allocation for A matrix
	model[2] = create2DArray(n, m); // Memory allcoation for B matrix

	
	strcpy(fileName, path);
	strcat(fileName, "PI.txt");
	printf("%s\n", fileName);
	
	i = 0;
	f = fopen(fileName, READ);
	while(fscanf(f, "%Lf ", &model[0][0][i++]) != EOF); // store PI values from file to array
	fclose(f);

	strcpy(fileName, path);
	strcat(fileName, "A.txt");
	printf("%s\n", fileName);

	i = 0;
	j = 0;
	f = fopen(fileName, READ);
	while(fscanf(f, "%Lf ", &model[1][i][j]) != EOF) { // store A matrix values from file to array
		j = j + 1;
		if(j == n) {
			i = i + 1;
			j = 0;
		}
	}
	fclose(f);

	strcpy(fileName, path);
	strcat(fileName, "B.txt");
	printf("%s\n", fileName);

	i = 0;
	j = 0;
	f = fopen(fileName, READ);
	while(fscanf(f, "%Lf ", &model[2][i][j]) != EOF) { // store B matrix values from file to array
		j = j + 1;
		if(j == m) {
			i = i + 1;
			j = 0;
		}
	}
	fclose(f);
	return model;
}
// Calculate average energy of the data
long double getAvgEnergy(long double * energy, int arr_len) {
	long double avgEnergy = 0;
	for(int i = 0; i < arr_len; i++)
		avgEnergy += energy[i];

	avgEnergy /= arr_len;
	return avgEnergy;
}
// Calculate front marker based on average energy
int getFrontMarker_avg_energy_based(long double * energy, int arr_len, long double avgEnergy) {
	int i;
	long double precentIncrement = 1.1;
	for(i = 0; i < arr_len; i++) 
		if(energy[i] > precentIncrement * avgEnergy) 
			break;
	if(i == arr_len)
		return -1;
	else
		return i;
}

// Calculate end marker based on average energy
int getEndMarker_avg_energy_based(long double * energy, int arr_len, long double avgEnergy) {
	int i;
	long double precentIncrement = 1.1;
	for(i = arr_len - 1; i>= 0; i--) 
		if(energy[i] > precentIncrement * avgEnergy) 
			break;
	if(i == 0)
		return -1;
	else
		return i;
}

void training(char * rollNum) {

	
	//char rollno[50] = ROLLNO;
	char rollno[100];
	strcpy(rollno, rollNum);
	//printf("Model training is ahead\n");
	//getch();
	char fileNameOfData[100] = ".\\resource\\data.txt"; // it stores the recorded data
	char * fileNameOfSilence = ".\\resource\\silence.txt"; // it stores the recorded silence
	char * fileNameOfEnergy = ".\\resource\\energy.txt";
	char * fileNameOfTrimmedData = ".\\resource\\trimmedData.txt";
	char * fileNameOfCodeBook = ".\\resource\\codebook.txt";
	char * fileNameOfPIMatrix = ".\\resource\\initial model\\PI.txt";
	char * fileNameOfAMatrix = ".\\resource\\initial model\\A.txt";
	char * fileNameOfBMatrix = ".\\resource\\initial model\\B.txt";
	char pathForFinalModel[100];

	
	long double arr[ARRAY_SIZE_L]; // arr contains the data of the original signal (which is loaded from file)
	long double energy[ARRAY_SIZE_S]; // energy contains the energy of all the frames in the data


	int numOfDigits; // Number of digits in the training dataset
	int numOfUtterances; // Number of utterances in the training dataset 
	
	int * t_array; // array of number of observations of all the utterances of the digit
	long double * *tCiVector; // Contains t Ci vectors for a utterance of the digit


	int * obsSeq; // Observation sequence
	long double * *codeBook; 
	int numOfRowsInCodeBook; // Number of rows in the code book


	int * *obsSeqOfAllUtterances;
	long double *** *modelsForAllDigit; // it contains 10 models as (0,1,2,...,9) 10 digits are used.

	int m; // Number of possible observations or size of the code book
	int n; // Number of possible states
	int t; // Number of observations or Number of frames for particular utterance 

	long double ** *initialModel; // Initial Inertia Model
	long double *** *modelsForADigit; // it contains number of models equal to number of utterances
	long double ** *avgModel; // average model obtained by taking average of the models of each utterance of a digit


	char path[100];

	sprintf(path, "./resource/%s_digitRecording", rollNum);
	//printf("Path: %s\n", path);
	int err = mkdir(path);
	if(err == -1) {

		//return ;
	}
	else {
		printf("New user is created.\n");
	}

	sprintf(pathForFinalModel, "./resource/final model_%s", rollno);
	printf("Path: %s\n", pathForFinalModel);
	err = mkdir(pathForFinalModel);

	if(err == -1) {
		printf("Models are already formed.\nPlease proceed to use the calcy\n");
		return ;
	}
	else {
		// Initialization and memeory allocation

	m = CODEBOOK_SIZE;
	n = 5;

	numOfUtterances = 1; // number of utterances 20
	numOfDigits = 14; // number of digits 10

	numOfRowsInCodeBook = CODEBOOK_SIZE;
	printf("Model training is ahead\n");
	
	t_array = (int *) malloc(sizeof(int) * numOfUtterances);

	numOfRowsInCodeBook = CODEBOOK_SIZE;


	codeBook = getCodebook(fileNameOfCodeBook, CODEBOOK_SIZE, P, &numOfRowsInCodeBook); // Read codebook from file
	obsSeqOfAllUtterances = (int * *) malloc(sizeof(int *) * numOfUtterances);
	initialModel = getInitialModelFromFile(fileNameOfPIMatrix, fileNameOfAMatrix, fileNameOfBMatrix, n, m); // Read initial model from file
	modelsForADigit = (long double *** *) malloc(sizeof(long double ***) * numOfUtterances);
	modelsForAllDigit = (long double *** *) malloc(sizeof(long double ***) * numOfDigits);
	

	for(int digit = 0; digit < numOfDigits; digit++) { // digits 0,9 then 10,11,12,13,14 (by, into, minus, plus, calcy)
		for(int utterance = 1; utterance <= numOfUtterances; utterance++) {
			sprintf(fileNameOfData, ".\\resource\\%s_digitRecording\\%s_%d_%d.txt", rollno, rollno, digit, utterance);
			
			//------------------------------------------------Getting DC shift-----------------------------------
			int len1 = cpFromFileToArray(fileNameOfSilence, arr);

	
			long double shift = getDCShift(arr,len1);

			//---------------------------------------------------------------------------------------------------

			//------------------------------------------------Applying DC Shift----------------------------------
			printf("%s\n", fileNameOfData);
			int len2 = cpFromFileToArray(fileNameOfData,arr);

			shiftData(arr,len2,shift);  // Shift the data


			//-----------------------------------------------------------------------------------------------------
			//-------------------------------------------------Normalisation of data------------------------------
			normaliseData(arr,len2);  // Normalize the data


			//----------------------------------------------------------------------------------------------------
			int len3 = storeEnergy(FRAME_SIZE, arr, len2, energy); // Get the energy 	
			cpFromArrayToFile(fileNameOfEnergy,energy,len3);

			//-----------------------------------------------------------------------------------------------------

			//-------------------------------------------------Get Markers----------------------------------------
			int front = 0; // frame number
			int end = len3 - 1;
			int t = end - front + 1; // number of frames
			


			// Based on 10 % of avg energy
			long double avgEnergy = getAvgEnergy(energy, len3);
			int front_2 = getFrontMarker_avg_energy_based(energy, len3, avgEnergy); // frame number
			int end_2 = getEndMarker_avg_energy_based(energy, len3, avgEnergy);

			int numberOfFrameToBeAdded = 10; // here n = 5

			front_2 = front_2 - numberOfFrameToBeAdded;
			end_2 = end_2 + numberOfFrameToBeAdded;

			if(t >= 200) {
				front = front_2;
				end = end_2;
			}
			t = end - front + 1; // number of frames in the digit



			printf("Number of frames: %d\n", t);

			// store the trimmed signal to analyse the the signal using cool edit software
			int startOfData = front * 80; // due to sliding of 80
			int endOfData = end * 80;
			printf("Front Sample: %d End Sample: %d\n", front * 80, end * 80);
			printf("Front time: %Lf End time: %Lf\n", (long double)front * 80 / 16000, (long double) end * 80 / 16000);


			//----------------------------------------------------------------------------------------------------------

			//-------------------------------------------------Get Ci's of t frames-------------------------------------

			t_array[utterance - 1] = t;

			tCiVector = getCiOfTFrames(arr, t, front, P);

	
			//----------------------------------------------------------------------------------------------------------
			//--------------------------------------------------Generate Observation Sequence----------------------------

			obsSeqOfAllUtterances[utterance - 1] = generateObsSeq(codeBook, numOfRowsInCodeBook, tCiVector, t, P);

			//-----------------------------------------------------------------------------------------------------------
			//---------------------------------------------Obtain a model (using obsSeq and intial inertia model)--------

			modelsForADigit[utterance - 1] = getModelUsingHMM(initialModel, obsSeqOfAllUtterances[utterance - 1], n, m, t);
			//------------------------------------------------------------------------------------------------------------
		}
		//printf("Digit: %d is completed\n\n", digit);
		//getch();
		avgModel = getAvgModel(modelsForADigit, numOfUtterances, n, m);


		//-----------------------------------------------------Convergence of avgModel--------------------------------------------
		for(int converge_index = 1; converge_index <= 0; converge_index++) {
			for(int utter_i = 1; utter_i <= numOfUtterances; utter_i++) {
				modelsForADigit[utter_i - 1] = getModelUsingHMM(avgModel, obsSeqOfAllUtterances[utter_i - 1], n, m, t_array[utter_i - 1]);
			}
			// free(avgModel); 
			avgModel = getAvgModel(modelsForADigit, numOfUtterances, n, m);
		}

		//--------------------------------------------------------------------------------------------------------------------------

		modelsForAllDigit[digit] = avgModel;
		
		
		sprintf(pathForFinalModel, "./resource/final model_%s/model-%d", rollno,digit);
		mkdir(pathForFinalModel);
		printf("Directory is created: %s", pathForFinalModel);

		
		sprintf(pathForFinalModel, ".\\resource\\final model_%s\\model-%d\\", rollno,digit);
		cpModelToFile(avgModel, pathForFinalModel, n, m);
		printf("Digit: %d is completed and its model is formed\n\n", digit);	
	}

	}
	
	
}

void prerecordedTesting() {
	char rollno[50] = ROLLNO;

	int numOfDigits;
	int numOfUtterances;

	char * fileNameOfSilence = ".\\resource\\silence.txt"; // it stores the recorded silence
	char * fileNameOfEnergy = ".\\resource\\energy.txt";
	char * fileNameOfCodeBook = ".\\resource\\codebook.txt";
	char pathForFinalModel[100]; 
	char fileNameOfData[100];

	long double arr[ARRAY_SIZE_L]; // arr contains the data of the original signal (which is loaded from file)
	long double energy[ARRAY_SIZE_S]; // energy contains the energy of all the frames in the data

	long double *** *modelsForAllDigit; // it contains 10 models as (0,1,2,...,9) 10 digits are used.
	int * t_array; // array of number of observations of all the utterances of the digit
	long double * *tCiVector; // Contains t Ci vectors for a utterance of the digit

	int m; // Number of possible observations or size of the code book
	int n; // Number of possible states
	//int t; // Number of observations or Number of frames for particular utterance 

	int * obsSeq; // Observation sequence
	int correctlyIdentifiedDigit;
	int totalNumberOfTestData;

	long double * accuracy;
	long double * *codeBook; 
	int numOfRowsInCodeBook; // Number of rows in the code book
	// Initialization 
	


	numOfDigits = 15;
	numOfUtterances = 10;

	m = CODEBOOK_SIZE;
	n = 5;

	numOfRowsInCodeBook = CODEBOOK_SIZE;

	codeBook = getCodebook(fileNameOfCodeBook, CODEBOOK_SIZE, P, &numOfRowsInCodeBook); // Read codebook from file
	modelsForAllDigit = (long double *** *) malloc(sizeof(long double ***) * numOfDigits);
	t_array = (int *) malloc(sizeof(int) * numOfUtterances);
	accuracy = (long double *) malloc(sizeof(long double) * numOfDigits);
	correctlyIdentifiedDigit = 0;
	totalNumberOfTestData = 0;
	for(int i = 0; i < numOfDigits; i++)
		accuracy[i] = 0;


	// Load all 10 models obtained by training
	for(int digit = 0; digit < numOfDigits; digit++) {
		sprintf(pathForFinalModel, ".\\resource\\final model\\model-%d\\", digit);
		modelsForAllDigit[digit] = cpFromFileToModel(pathForFinalModel, n, m);
	}






	for(int digit = 0; digit < numOfDigits; digit++) {
		for(int utterance = 1; utterance <= numOfUtterances; utterance++) { 
			//--------------------------------------------------File Name formation (For Testing Part)------------------
			//e.g. ".\\resource\\194101017_digitRecording\\194101017_0_1.txt"
			sprintf(fileNameOfData, ".\\resource\\%s_digitRecording\\%s_%d_%d.txt", rollno,rollno,digit, utterance + 20);

			//-----------------------------------------------------------------------------------------------------------

			//------------------------------------------------Getting DC shift-----------------------------------
			int len1 = cpFromFileToArray(fileNameOfSilence, arr);
			printf("Samples in the silence: %d\n",len1);

	
			long double shift = getDCShift(arr,len1);
			printf("Shift: %Lf\n", shift);
			//getch();
			//---------------------------------------------------------------------------------------------------
			//------------------------------------------------Applying DC Shift----------------------------------
			int len2 = cpFromFileToArray(fileNameOfData,arr);
			printf("Samples in the data: %d\n",len2);

			shiftData(arr,len2,shift);  // Shift the data


			//-----------------------------------------------------------------------------------------------------
			//-------------------------------------------------Normalisation of data------------------------------
			normaliseData(arr,len2);  // Normalize the data

			//----------------------------------------------------------------------------------------------------
			int len3 = storeEnergy(FRAME_SIZE, arr, len2, energy); // Get the energy 	
			cpFromArrayToFile(fileNameOfEnergy,energy,len3);

			//-----------------------------------------------------------------------------------------------------
			//-------------------------------------------------Get Markers----------------------------------------
			int front = 0; // frame number
			int end = len3 - 1;
			int t = end - front + 1; // number of frames
			
			//printf("Markers without trimming:\nfront:%d\nend:%d\n", front, end);

			// Based on 10 % of avg energy
			long double avgEnergy = getAvgEnergy(energy, len3);
			int front_2 = getFrontMarker_avg_energy_based(energy, len3, avgEnergy); // frame number
			int end_2 = getEndMarker_avg_energy_based(energy, len3, avgEnergy);

			int numberOfFrameToBeAdded = 10; // here n = 5

			front_2 = front_2 - numberOfFrameToBeAdded;
			end_2 = end_2 + numberOfFrameToBeAdded;

			if(t >= 200) {
				front = front_2;
				end = end_2;
			}
			t = end - front + 1; // number of frames in the digit
			//printf("Markers with trimming:\nfront:%d\nend:%d\n", front_2, end_2);


			printf("Number of frames: %d\n", t);

			// store the trimmed signal to analyse the the signal using cool edit software
			int startOfData = front * 80; // due to sliding of 80
			int endOfData = end * 80;
			printf("Front Sample: %d End Sample: %d\n", front * 80, end * 80);
			printf("Front time: %Lf End time: %Lf\n", (long double)front * 80 / 16000, (long double) end * 80 / 16000);


			//----------------------------------------------------------------------------------------------------------
			//-------------------------------------------------Get Ci's of t frames-------------------------------------

			t_array[utterance - 1] = t;
			//printf("Number of frames in the observation: %d\n", t);
			//getch();
			tCiVector = getCiOfTFrames(arr, t, front, P);

	
			//----------------------------------------------------------------------------------------------------------
			//--------------------------------------------------Generate Observation Sequence----------------------------
			obsSeq = generateObsSeq(codeBook, numOfRowsInCodeBook, tCiVector, t, P);

			//-----------------------------------------------------------------------------------------------------------
			//---------------------------------------------------Choose model with high probability of P(O/model)-------------------------------
			int recognizedDigit = getModelWithHighProb(modelsForAllDigit, 0, DIGITS - 1, DIGITS, obsSeq, n, t);
			if(recognizedDigit == digit) {
				accuracy[digit] += 1;
				correctlyIdentifiedDigit += 1;
			}
			totalNumberOfTestData += 1;
			printf("File under Testing: \n%s\n", fileNameOfData);
			printf("Recognized Digit: %d\n\n\n",recognizedDigit);
			//------------------------------------------------------------------------------------------------------------
		}
	}

	
	// Print Individual Accuracy
	for(int i = 0 ; i < numOfDigits; i++) {
		accuracy[i] /= numOfDigits;
		printf("Accuracy of digit - %d : %Lf\n", i, accuracy[i] * 100);
	}

	// Print overall accuracy (=correctly recognized digit / total number of test data)
	long double avg_acc = (long double) correctlyIdentifiedDigit / totalNumberOfTestData;
	printf("correctlyIdentifiedDigit: %d\ntotalNumberOfTestData: %d\n", correctlyIdentifiedDigit, totalNumberOfTestData);
	printf("Average Accuracy: %Lf\n", avg_acc * 100);
	
}

void prerecordedSingleUtterance() {
	char rollno[50] = ROLLNO;

	int numOfDigits;
	int numOfUtterances;

	//char * fileNameOfSilence = ".\\resource\\silence.txt"; // it stores the recorded silence
	char * fileNameOfSilence = "silence.txt"; // it stores the recorded silence
	char * fileNameOfEnergy = ".\\resource\\energy.txt";
	char * fileNameOfCodeBook = ".\\resource\\codebook.txt";
	char pathForFinalModel[100]; 
	char fileNameOfData[100] = "trimmedData.txt";

	long double arr[ARRAY_SIZE_L]; // arr contains the data of the original signal (which is loaded from file)
	long double energy[ARRAY_SIZE_S]; // energy contains the energy of all the frames in the data

	long double *** *modelsForAllDigit; // it contains 10 models as (0,1,2,...,9) 10 digits are used.
	int * t_array; // array of number of observations of all the utterances of the digit
	long double * *tCiVector; // Contains t Ci vectors for a utterance of the digit

	int m; // Number of possible observations or size of the code book
	int n; // Number of possible states
	//int t; // Number of observations or Number of frames for particular utterance 

	int * obsSeq; // Observation sequence
	int correctlyIdentifiedDigit;
	int totalNumberOfTestData;

	long double * accuracy;
	long double * *codeBook; 
	int numOfRowsInCodeBook; // Number of rows in the code book
	// Initialization 
	


	numOfDigits = 15;
	numOfUtterances = 10;

	m = CODEBOOK_SIZE;
	n = 5;

	numOfRowsInCodeBook = CODEBOOK_SIZE;


	codeBook = getCodebook(fileNameOfCodeBook, CODEBOOK_SIZE, P, &numOfRowsInCodeBook); // Read codebook from file
	modelsForAllDigit = (long double *** *) malloc(sizeof(long double ***) * numOfDigits);
	t_array = (int *) malloc(sizeof(int) * numOfUtterances);
	accuracy = (long double *) malloc(sizeof(long double) * numOfDigits);
	correctlyIdentifiedDigit = 0;
	totalNumberOfTestData = 0;
	for(int i = 0; i < numOfDigits; i++)
		accuracy[i] = 0;


	// Load all 10 models obtained by training
	for(int digit = 0; digit < DIGITS; digit++) {
		sprintf(pathForFinalModel, ".\\resource\\final model\\model-%d\\", digit);
		modelsForAllDigit[digit] = cpFromFileToModel(pathForFinalModel, n, m);
	}





	numOfDigits = 1;
	numOfUtterances = 1;
	for(int digit = 0; digit < numOfDigits; digit++) {
		for(int utterance = 1; utterance <= numOfUtterances; utterance++) { 

			//-----------------------------------------------------------------------------------------------------------

			//------------------------------------------------Getting DC shift-----------------------------------
			int len1 = cpFromFileToArray(fileNameOfSilence, arr);
			printf("Samples in the silence: %d\n",len1);
			//getch();
			//display1DArray(arr,len1);
	
			long double shift = getDCShift(arr,len1);
			printf("Shift: %Lf\n", shift);
			//getch();
			//---------------------------------------------------------------------------------------------------
			//------------------------------------------------Applying DC Shift----------------------------------
			int len2 = cpFromFileToArray(fileNameOfData,arr);
			printf("Samples in the data: %d\n",len2);

			shiftData(arr,len2,shift);  // Shift the data


			//-----------------------------------------------------------------------------------------------------
			//-------------------------------------------------Normalisation of data------------------------------
			normaliseData(arr,len2);  // Normalize the data


			//----------------------------------------------------------------------------------------------------
			int len3 = storeEnergy(FRAME_SIZE, arr, len2, energy); // Get the energy 	
			cpFromArrayToFile(fileNameOfEnergy,energy,len3);

			//-----------------------------------------------------------------------------------------------------
			//-------------------------------------------------Get Markers----------------------------------------
			int front = 0; // frame number
			int end = len3 - 1;
			int t = end - front + 1; // number of frames
			
			//printf("Markers without trimming:\nfront:%d\nend:%d\n", front, end);

			// Based on 10 % of avg energy
			long double avgEnergy = getAvgEnergy(energy, len3);
			int front_2 = getFrontMarker_avg_energy_based(energy, len3, avgEnergy); // frame number
			int end_2 = getEndMarker_avg_energy_based(energy, len3, avgEnergy);

			int numberOfFrameToBeAdded = 10; // here n = 5

			front_2 = front_2 - numberOfFrameToBeAdded;
			end_2 = end_2 + numberOfFrameToBeAdded;

			if(t >= 200) {
				front = front_2;
				end = end_2;
			}
			t = end - front + 1; // number of frames in the digit
			//printf("Markers with trimming:\nfront:%d\nend:%d\n", front_2, end_2);

			printf("Markers with trimming:\nfront:%d\nend:%d\nDiff: %d\n", front, end, end - front);
			printf("Number of frames: %d\n", t);
			getch();
			// store the trimmed signal to analyse the the signal using cool edit software
			int startOfData = front * 80; // due to sliding of 80
			int endOfData = end * 80;
			printf("Front Sample: %d End Sample: %d\n", front * 80, end * 80);
			printf("Front time: %Lf End time: %Lf\n", (long double)front * 80 / 16000, (long double) end * 80 / 16000);


			//----------------------------------------------------------------------------------------------------------
			//-------------------------------------------------Get Ci's of t frames-------------------------------------

			t_array[utterance - 1] = t;

			tCiVector = getCiOfTFrames(arr, t, front, P);

	
			//----------------------------------------------------------------------------------------------------------
			//--------------------------------------------------Generate Observation Sequence----------------------------
			obsSeq = generateObsSeq(codeBook, numOfRowsInCodeBook, tCiVector, t, P);

			//-----------------------------------------------------------------------------------------------------------
			//---------------------------------------------------Choose model with high probability of P(O/model)-------------------------------
			int recognizedDigit = getModelWithHighProb(modelsForAllDigit, 0, DIGITS - 1, DIGITS, obsSeq, n, t);

			printf("Recognized Digit: %d\n\n\n",recognizedDigit);
			//------------------------------------------------------------------------------------------------------------
		}
	}
}


void prerecordedThreeSignal(int * op1, int * op, int * op2) {
	char rollno[50] = ROLLNO;

	int recDigit[3];
	char * threeSignal[3] = {"operand1.txt", "operator.txt", "operand2.txt"};
	int numOfDigits;
	int numOfUtterances;

	//char * fileNameOfSilence = ".\\resource\\silence.txt"; // it stores the recorded silence
	char * fileNameOfSilence = "silence.txt"; // it stores the recorded silence
	char * fileNameOfEnergy = ".\\resource\\energy.txt";
	char * fileNameOfCodeBook = ".\\resource\\codebook.txt";
	char pathForFinalModel[100]; 
	char fileNameOfData[100] = "trimmedData.txt";

	long double arr[ARRAY_SIZE_L]; // arr contains the data of the original signal (which is loaded from file)
	long double energy[ARRAY_SIZE_S]; // energy contains the energy of all the frames in the data

	long double *** *modelsForAllDigit; // it contains 10 models as (0,1,2,...,9) 10 digits are used.
	int * t_array; // array of number of observations of all the utterances of the digit
	long double * *tCiVector; // Contains t Ci vectors for a utterance of the digit

	int m; // Number of possible observations or size of the code book
	int n; // Number of possible states
	//int t; // Number of observations or Number of frames for particular utterance 

	int * obsSeq; // Observation sequence
	int correctlyIdentifiedDigit;
	int totalNumberOfTestData;

	long double * accuracy;
	long double * *codeBook; 
	int numOfRowsInCodeBook; // Number of rows in the code book
	// Initialization 
	


	numOfDigits = 3;
	numOfUtterances = 1;

	m = CODEBOOK_SIZE;
	n = 5;

	numOfRowsInCodeBook = CODEBOOK_SIZE;


	codeBook = getCodebook(fileNameOfCodeBook, CODEBOOK_SIZE, P, &numOfRowsInCodeBook); // Read codebook from file
	modelsForAllDigit = (long double *** *) malloc(sizeof(long double ***) * numOfDigits);
	t_array = (int *) malloc(sizeof(int) * numOfUtterances);
	accuracy = (long double *) malloc(sizeof(long double) * numOfDigits);
	correctlyIdentifiedDigit = 0;
	totalNumberOfTestData = 0;
	for(int i = 0; i < numOfDigits; i++)
		accuracy[i] = 0;


	// Load all 10 models obtained by training
	for(int digit = 0; digit < DIGITS; digit++) {
		sprintf(pathForFinalModel, ".\\resource\\final model\\model-%d\\", digit);
		modelsForAllDigit[digit] = cpFromFileToModel(pathForFinalModel, n, m);
	}





	numOfDigits = 3;
	numOfUtterances = 1;
	for(int digit = 0; digit < numOfDigits; digit++) {
		for(int utterance = 1; utterance <= numOfUtterances; utterance++) { 


			//------------------------------------------------Getting DC shift-----------------------------------
			int len1 = cpFromFileToArray(fileNameOfSilence, arr);
			printf("Samples in the silence: %d\n",len1);
			//getch();
			//display1DArray(arr,len1);
	
			long double shift = getDCShift(arr,len1);
			printf("Shift: %Lf\n", shift);
			//getch();
			//---------------------------------------------------------------------------------------------------
			//------------------------------------------------Applying DC Shift----------------------------------
			int len2 = cpFromFileToArray(threeSignal[digit],arr);
			printf("Samples in the data: %d\n",len2);

			shiftData(arr,len2,shift);  // Shift the data


			//-----------------------------------------------------------------------------------------------------
			//-------------------------------------------------Normalisation of data------------------------------
			normaliseData(arr,len2);  // Normalize the data


			//----------------------------------------------------------------------------------------------------
			int len3 = storeEnergy(FRAME_SIZE, arr, len2, energy); // Get the energy 	
			cpFromArrayToFile(fileNameOfEnergy,energy,len3);

			//-----------------------------------------------------------------------------------------------------
			//-------------------------------------------------Get Markers----------------------------------------
			int front = 0; // frame number
			int end = len3 - 1;
			int t = end - front + 1; // number of frames
			


			// Based on 10 % of avg energy
			long double avgEnergy = getAvgEnergy(energy, len3);
			int front_2 = getFrontMarker_avg_energy_based(energy, len3, avgEnergy); // frame number
			int end_2 = getEndMarker_avg_energy_based(energy, len3, avgEnergy);

			int numberOfFrameToBeAdded = 10; // here n = 5

			front_2 = front_2 - numberOfFrameToBeAdded;
			end_2 = end_2 + numberOfFrameToBeAdded;

			if(t >= 200) {
				front = front_2;
				end = end_2;
			}
			t = end - front + 1; // number of frames in the digit
			//printf("Markers with trimming:\nfront:%d\nend:%d\n", front_2, end_2);

			printf("Markers with trimming:\nfront:%d\nend:%d\nDiff: %d\n", front, end, end - front);
			printf("Number of frames: %d\n", t);
			//getch();
			// store the trimmed signal to analyse the the signal using cool edit software
			int startOfData = front * 80; // due to sliding of 80
			int endOfData = end * 80;
			printf("Front Sample: %d End Sample: %d\n", front * 80, end * 80);
			printf("Front time: %Lf End time: %Lf\n", (long double)front * 80 / 16000, (long double) end * 80 / 16000);


			//----------------------------------------------------------------------------------------------------------
			//-------------------------------------------------Get Ci's of t frames-------------------------------------

			t_array[utterance - 1] = t;
			//printf("Number of frames in the observation: %d\n", t);
			//getch();
			tCiVector = getCiOfTFrames(arr, t, front, P);

	
			//----------------------------------------------------------------------------------------------------------
			//--------------------------------------------------Generate Observation Sequence----------------------------
			obsSeq = generateObsSeq(codeBook, numOfRowsInCodeBook, tCiVector, t, P);

			//-----------------------------------------------------------------------------------------------------------
			//---------------------------------------------------Choose model with high probability of P(O/model)-------------------------------
			int recognizedDigit;
			if(digit == 1) {  // compare operator with operator models
				recognizedDigit = getModelWithHighProb(modelsForAllDigit, 10, 13, DIGITS, obsSeq, n, t);
				if(recognizedDigit == -1)
					recognizedDigit = 13;

			}
			else { // compare operand with operand models
				recognizedDigit = getModelWithHighProb(modelsForAllDigit, 0, 9, DIGITS, obsSeq, n, t);
				if(recognizedDigit == -1)
					recognizedDigit = 7;
			}	
			
			

			printf("Recognized Digit: %d\n\n\n",recognizedDigit);
			recDigit[digit] = recognizedDigit;
			//getch();
			//------------------------------------------------------------------------------------------------------------
		}
	}
	*op1 = recDigit[0];
	*op = recDigit[1];
	*op2 = recDigit[2];
}


int liveTesting() {

	char rollno[100];
	printf("Enter username: ");
	scanf("%s", rollno);

	char * fileNameOfSilence = "silence.txt"; // For live record testing
	char * fileNameOfData = "data.txt"; // For live record testing
	char * fileNameOfEnergy = ".\\resource\\energy.txt";
	char * fileNameOfCodeBook = ".\\resource\\codebook.txt";
	char * trimmedData = "trimmedData.txt";

	char pathForFinalModel[100]; 

	long double arr[ARRAY_SIZE_L]; // arr contains the data of the original signal (which is loaded from file)
	long double energy[ARRAY_SIZE_S]; // energy contains the energy of all the frames in the data
	long double * arr_temp;

	long double *** *modelsForAllDigit; // it contains 10 models as (0,1,2,...,9) 10 digits are used.
	long double * *tCiVector; // Contains t Ci vectors for a utterance of the digit

	int numOfDigits;
	int m; // Number of possible observations or size of the code book
	int n; // Number of possible states

	int * obsSeq; // Observation sequence
	long double * *codeBook;
	int numOfRowsInCodeBook; // Number of rows in the code book

	// Initialization
	numOfDigits = 15;
	numOfRowsInCodeBook = CODEBOOK_SIZE;

	m = CODEBOOK_SIZE;
	n = 5;

	codeBook = getCodebook(fileNameOfCodeBook, CODEBOOK_SIZE, P, &numOfRowsInCodeBook); // Read codebook from file
	modelsForAllDigit = (long double *** *) malloc(sizeof(long double ***) * numOfDigits);


	sprintf(pathForFinalModel, "./resource/final model_%s", rollno);
	int err = mkdir(pathForFinalModel);

	if(err == -1) {
		printf("Valid user.\n");
	}
	else {
		printf("Invalid user\t");
		rmdir(pathForFinalModel);
		return -1;
	}

	/*
	// Load all 15(numOfDigits) models obtained by training
	for(int digit = 0; digit < DIGITS; digit++) {
		//sprintf(pathForFinalModel, ".\\resource\\final model\\model-%d\\", digit);
		sprintf(pathForFinalModel, "./resource/final model_%s/model-%d/", rollno,digit);
		modelsForAllDigit[digit] = cpFromFileToModel(pathForFinalModel, n, m);
	}

	*/
	
		////--------------------------------------------------Integrated Module for silence----------------------------
		printf("Record Silence:\n");
		//getch();
		system("Recording_Module.exe 3 input_file.wav silence.txt");
		//-----------------------------------------------------------------------------------------------------------
		//--------------------------------------------------Integrated Module for digit------------------------------
		printf("Record Digit:\n");
		//getch();
		system("Recording_Module.exe 3 input_file.wav data.txt");
		//------------------------------------------------------------------------------------------------------------

	/*
		//------------------------------------------------Getting DC shift-----------------------------------
		int len1 = cpFromFileToArray(fileNameOfSilence, arr);
		printf("Samples in the silence: %d\n",len1);
		//getch();
		//display1DArray(arr,len1);
	
		long double shift = getDCShift(arr,len1);
		printf("Shift: %Lf\n", shift);
		//getch();
		//---------------------------------------------------------------------------------------------------
		//------------------------------------------------Applying DC Shift----------------------------------
		int len2 = cpFromFileToArray(fileNameOfData,arr);
		arr_temp = (long double *) malloc(sizeof(long double) * ARRAY_SIZE_L);
		cpFromFileToArray(fileNameOfData,arr_temp);
		printf("Samples in the data: %d\n",len2);
		//getch();
		//display1DArray(arr,len2);
		//printf("\n");
		shiftData(arr,len2,shift);  // Shift the data


		//-----------------------------------------------------------------------------------------------------
		//-------------------------------------------------Normalisation of data------------------------------
		normaliseData(arr,len2);  // Normalize the data

		//printf("After norm: \n");
		//getch();
		//display1DArray(arr,len2);
		//----------------------------------------------------------------------------------------------------
		int len3 = storeEnergy(FRAME_SIZE, arr, len2, energy); // Get the energy 	
		cpFromArrayToFile(fileNameOfEnergy,energy,len3);

		//-----------------------------------------------------------------------------------------------------
		//-------------------------------------------------Get Markers----------------------------------------
		int front = 0; // frame number
		int end = len3 - 1;
		int t = end - front + 1; // number of frames
			
		//printf("Markers without trimming:\nfront:%d\nend:%d\n", front, end);

		// Based on 10 % of avg energy
		long double avgEnergy = getAvgEnergy(energy, len3);
		int front_2 = getFrontMarker_avg_energy_based(energy, len3, avgEnergy); // frame number
		int end_2 = getEndMarker_avg_energy_based(energy, len3, avgEnergy);

		int numberOfFrameToBeAdded = 10; // here n = 5

		front_2 = front_2 - numberOfFrameToBeAdded;
		end_2 = end_2 + numberOfFrameToBeAdded;

		if(t >= 200) {
			front = front_2;
			end = end_2;
		}
		t = end - front + 1; // number of frames in the digit
		printf("Markers with trimming:\nfront:%d\nend:%d\nDiff: %d\n", front, end, end - front);
		//getch();


		printf("Number of frames: %d\n", t);

		// store the trimmed signal to analyse the the signal using cool edit software
		int startOfData = front * 80; // due to sliding of 80
		int endOfData = end * 80;
		printf("Front Sample: %d End Sample: %d\n", front * 80, end * 80);
		printf("Front time: %Lf End time: %Lf\n", (long double)front * 80 / 16000, (long double) end * 80 / 16000);


		//----------------------------------------------------------------------------------------------------------


		
		//-------------------------------------------------Get Ci's of t frames-------------------------------------

		tCiVector = getCiOfTFrames(arr, t, front, P);

		//----------------------------------------------------------------------------------------------------------

		//--------------------------------------------------Generate Observation Sequence----------------------------
		obsSeq = generateObsSeq(codeBook, numOfRowsInCodeBook, tCiVector, t, P);
		printf("\nObservation Sequence: \n");
		display1DIntArray(obsSeq, t);

		//-----------------------------------------------------------------------------------------------------------
		//---------------------------------------------------Choose model with high probability of P(O/model)-------------------------------
		int recognizedDigit = getModelWithHighProb(modelsForAllDigit, 0, DIGITS - 1, DIGITS, obsSeq, n, t);

		//------------------------------------------------------------------------------------------------------------
		*/

		return 0;
}
void joys() {
	cout << "Reached" << endl;
}


