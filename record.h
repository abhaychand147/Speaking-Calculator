
#include <iostream>
#include<conio.h>
#include<fstream>
#include<cmath>
#include<math.h>
#include<string>

using namespace std;


void record(char * name);
int cpFromFileToArray(char * fname, long double * arr);
void cpFromArrayToFile(char * fname, long double * arr, int arr_len);
long double getDCShift(long double * arr, int len);
void shiftData(long double *arr, int len, long double shift);
long double getAbsMax(long double * arr, int len);
void normaliseData(long double * arr, int len);
long double getEnergyPerFrame(int frame_size, long double * arr, int low, int high);
int storeEnergy(int frame_size, long double * arr, int len_arr, long double * e);
int getFrontMarker(long double * energy, int arr_len);
int getEndMarker(long double * energy, int arr_len);
long double getAvgEnergy(long double * energy, int arr_len);
int getFrontMarker_avg_energy_based(long double * energy, int arr_len, long double avgEnergy);
int getEndMarker_avg_energy_based(long double * energy, int arr_len, long double avgEnergy);
void trainingModule(char * rollNo);
void trainingModule1();

