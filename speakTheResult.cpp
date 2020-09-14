# include "speakTheResult.h"
# include <iostream>
# include <Windows.h>
# include <cmath>

#pragma comment(lib, "winmm.lib")

using namespace std;


void speakIntegerResult(int result) {
	LPCWSTR singleDigit[10] = {L"0.wav", L"1.wav", L"2.wav", L"3.wav", L"4.wav", L"5.wav", L"6.wav", L"7.wav", L"8.wav", L"9.wav"}; // L"one.wav";
	LPCWSTR doubleDigit[9] = {L"10.wav", L"20.wav", L"30.wav", L"40.wav", L"50.wav", L"60.wav", L"70.wav", L"80.wav", L"90.wav"}; 
	LPCWSTR elevenTo19[9] = {L"11.wav", L"12.wav", L"13.wav", L"14.wav", L"15.wav", L"16.wav", L"17.wav", L"18.wav", L"19.wav"};
	LPCWSTR minus = L"minus.wav";
	int unitPlace = result % 10;

	if(result < 0) {
		PlaySound(minus, NULL, SND_SYNC);
		result *= -1;
		unitPlace *= -1;
	}
	if(result < 10) {
		PlaySound(singleDigit[unitPlace], NULL, SND_SYNC);
	}
	else if(result > 10 && result < 20) {
		int tenthPlace = result / 10;
		PlaySound(elevenTo19[unitPlace - 1], NULL, SND_SYNC);
	}
	else if(result < 100){
		int tenthPlace = result / 10;
		PlaySound(doubleDigit[tenthPlace - 1], NULL, SND_SYNC);
		if(unitPlace != 0) { // number can be from (10,20,...,90)
			PlaySound(singleDigit[unitPlace], NULL, SND_SYNC);
		}
	}
}
void speakResult(long double num) {
	LPCWSTR point = L"point.wav";
	int x;
	int precision = 3;
	int integral_part = floor(num);
	long double fractional_part = num - integral_part;
	int count = 0;


	//cout << integral_part << "\n"<<fractional_part << endl;
	fractional_part += 0.0000000001; // to maintain the precision of long double
	speakIntegerResult(integral_part);
	if(fractional_part != 0) {
		PlaySound(point, NULL, SND_SYNC);

		while(count != precision) {
			fractional_part *= 10;
			//cout << "Floor of " << fractional_part;
			//cout << "is " << floor(fractional_part) << endl;
			integral_part = floor(fractional_part);
			fractional_part = fractional_part - integral_part;
			//cout << "Fraction: "<<fractional_part << endl;
			//cout << "Integral: " << integral_part << endl;
			speakIntegerResult(integral_part);
			
			//cout << fractional_part << endl;
			count++;
		}
	}
	
}