# include "record.h"
//# include "recordUserData.h"
#include"testing.h"
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

# define ARRAY_SIZE_L 60000
# define ARRAY_SIZE_S 1000
# define FRAME_SIZE 320
# define WINDOW_SIZE 320

using namespace std;


void record(char * rollNo)
{
	system("Recording_Module.exe 3 input_file.wav input_file.txt");

}

void trainingModule(char * rollNo) {

	char * fileNameOfData = ".\\data.txt"; // This file contains the data of the recorded signal (word/digit)
	char * fileNameOfSilence = ".\\silence.txt"; // This file contains the data of the recorded silence
	char * fileNameOfEnergy = ".\\energy.txt";

	char path[100];

	

	long double arr[ARRAY_SIZE_L]; // arr contains the data of the original signal (which is loaded from file)
	long double energy[ARRAY_SIZE_S]; // energy contains the energy of all the frames in the data

	FILE * f;

	char rollNum[50] = "194101017";
	char word[50];
	char fileNameOfDataEntry[100];
	int utterance;
	printf("Enter Roll number: ");
	scanf("%s", rollNum);

	strcpy(rollNo, rollNum);
	sprintf(path, "./resource/%s_digitRecording", rollNum);
	printf("Path: %s\n", path);
	int err = mkdir(path);
	if(err == -1) {
		printf("User is already exists\n");
		printf("Please close the terminal and retrain using new user name.\n");

		return ;
	}
	else {
		printf("New user is created.\n");
	}

	getch();

	printf("\n\nNote : Recording will be done for all the words\n\n");
	////--------------------------------------------------Integrated Module for silence---------------------------
	printf("Record Silence:\n");
	system("Recording_Module.exe 3 input_file.wav silence.txt");
	//------------------------------------------------------------------------------------------------------------
	for(int i_l = 0; i_l < 14; i_l++) {
			sprintf(path, ".\\resource\\%s_digitRecording\\", rollNum); //.\\resource\\%s_digitRecording\\

			printf("Enter word to be spoken: ");
			scanf("%s", word);
			printf("Enter the number of utterance: ");
			scanf("%d", &utterance);


	


		//*******************************************************RECORDING WORD ALONG WITH THE SILENCE***********************************//



			for(int j = 0; j < utterance; j++) {

				sprintf(fileNameOfDataEntry, "%s%s_%s_%d.txt", path, rollNum, word, j + 1);
				printf("File Name: %s\n", fileNameOfDataEntry);
				////--------------------------------------------------Integrated Module for digit-----------------------------
				printf("Record Digit:\n");
				system("Recording_Module.exe 3 input_file.wav data.txt");
				//------------------------------------------------------------------------------------------------------------


				//*******************************************************PREPROCESSING THE WORD SIGNAL*******************************************//

				int len1 = cpFromFileToArray(fileNameOfSilence, arr);
				long double shift = getDCShift(arr,len1); // calculate DC-shift from the silence


				int len2 = cpFromFileToArray(fileNameOfData,arr);

				shiftData(arr,len2,shift);  // Shift the data

				normaliseData(arr,len2);  // Normalize the data

				//*******************************************************ENERGY CALCULATION**************************************************//
				int len3 = storeEnergy(FRAME_SIZE, arr, len2, energy); // Get the energy 	
				cpFromArrayToFile(fileNameOfEnergy,energy,len3);

				//*******************************************************VOICE DETECTION**************************************************//


				// Approach 2: get the frame which has 10% more energy than the average energy
					long double avgEnergy = getAvgEnergy(energy, len3);
					//printf("\n\nAvg Energy: %Lf\n", avgEnergy);

					int front_2 = getFrontMarker_avg_energy_based(energy, len3, avgEnergy);
					int end_2 = getEndMarker_avg_energy_based(energy, len3, avgEnergy);

					int numberOfFrameToBeAdded = 10; // here n = 5

					front_2 = front_2 - numberOfFrameToBeAdded;
					end_2 = end_2 + numberOfFrameToBeAdded;


					printf("Number of frames: %d\n", end_2 - front_2 + 1);


					f = fopen(fileNameOfDataEntry, WRITE);
					for(int i = front_2 * 80; i <= end_2 * 80; i++)
						fprintf(f,"%Lf\n",arr[i]);
					fclose(f);
					getch();
			}


			printf("Completed recording for the %s\n\n\n", word);
			getch();
	}
 
}
