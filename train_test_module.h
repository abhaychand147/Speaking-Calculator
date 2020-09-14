#ifndef TRAIN_TEST_MODULE_H
#define TRAIN_TEST_MODULE_H

struct markers
{
	int start1, end1, start2, end2, start3, end3;
};


#endif



void normalize(long double maxVal, long double minVal);
void removeDcShift(long double avg);
void fetchData(int digit_no);
int setStartMarkerUtterence();
int setEndMarkerUtterence(int startMarker);
void trimUtterence(int startFrameNo,int endFrameNo);
markers setMarkers();
void trimSignal(markers m);
void findEnergy();
void overlapping();
void overlapping_test(int dc);
void hammingWindow();
void durbin();
void findRi();
void findCi();
void applyRSW();
void step1();
void calculateObsSeq(int choseCodebook);
void step2();
void setMatrices(int digit_no, int utterence_no, int iteration_no);
void forwardProcedure();
void backwardProcedure();
void viterbiAlgorithm();
void reestimation();
void initialization();
void step3();
void setMatrices_test(int digit_no_file, int digit_no_model,int utterence_no);
void step4();
void setMatrices_test_live(int digit_no_model, int e);
void live_train();
void step5(int* op1, char* op, int * op2);