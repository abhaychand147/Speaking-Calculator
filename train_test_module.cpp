#include "train_test_module.h"
#include <iostream>
#include <sstream>
#include <fstream>
#include <conio.h>
#include <vector>
#include <stdio.h>
#include "math.h"

using namespace std;

fstream f1,f2,f3;

int digits = 14, utterences = 20, dataCount = 0, frameSize = 320, maxRange = 100000, iteration1 = 3, iteration2 = 20, testUtterences = 10;
long double energyPerFrame[100000], operand1[100000], operand2[100000], operator1[100000];
long double data[100000], Ri[13], Ai[13], frame[320], Ci[13] ;
int dcOperand1= 0, dcOperand2 = 0, dcOperator=0;
char result2;
int result1, result3;
static int p = 12;



fstream fi,fo;

//input variables
int N=5;
int M=32;
int T=0;
long double a[6][6];
long double a_avg[6][6];
long double b[6][33];
long double b_avg[6][33];
long double pi[6] = { 0.0, 1.0, 0.0, 0.0, 0.0, 0.0 };
int obs[500];

//solution to problem 1 variables
long double alpha[500][6];
long double beta[500][6];
long double prob1 = 0.0;

//solution to problem 2 variables
long double gamma[500][6];
long double delta[500][6];
long double psy[500][6];
long double pStar;
int qStar[500];

//solution to problem 3 variables
long double et[500][6][6];
long double newpi[6];
long double newa[6][6];
long double newb[6][33];
long double XiSum[6][6];
long double gammaSum[6];

//tokhura weights
long double tw[12] = {1.0, 3.0, 7.0, 13.0, 19.0, 22.0, 25.0, 33.0, 42.0, 50.0, 56.0, 61.0};

//to normalize the data within the range +- 5000;
void normalize(long double maxVal, long double minVal)
{
	if( abs(minVal) > maxVal )
		maxVal = abs(minVal);
	long double range = 5000;
	for(int i = 0; i < dataCount; ++i)
		data[i] = ((data[i])/maxVal)*range;
}
//To remove the dc shift that exists due to ambient noise
void removeDcShift(long double avg)
{
	int nc=0;
	int noise[200000];
	long double value = 0.0;
	long double max = 0.0;
	FILE* fp = fopen("noise.txt","r");
	if(fp) {
		while(fscanf(fp,"%d",&noise[nc++]) != EOF);
		nc--;
		for(int i=0 ; i<nc ; i++)
			value += noise[i];
		value /= nc;
		for(int i=0 ; i<dataCount ; i++){
			data[i] -= value;
			if(max < abs(data[i])) max = abs(data[i]);
		}
		for(int i=0 ; i<dataCount ; i++){
			data[i] = (data[i]*5000)/max;
		}
	}
	else printf("cant open file");
}
// fetching signal from txt file intially
void fetchData(int digit_no)
{
	string s;
	if(digit_no <=9)
	{
		for(int i=0; i<10; i++)
			f1>>s;
	}
	for(int i=0; i<maxRange; i++)
		data[i]=0;
	long double val =0, maxVal =0, minVal =0, sum =0, avg =0;
	int flag = 0;
	dataCount = 0;
	while ( f1>>val )
	{
		dataCount++;
		data[dataCount-1] = val;
		sum += val;
		if(flag = 0)
		{
			maxVal = val;
			minVal = val;
			flag = 1;
		}
		if( val > maxVal)
			maxVal = val;
		else if( val < minVal)
			minVal = val;
	}

	//padding zeros if number of frames is not integer
	float i = floorf(dataCount/frameSize);
	if( (dataCount/frameSize) - i != 0)
	{
		int noOfZerosPadded = frameSize - (dataCount - i*frameSize);
		for(int j = 0; j<noOfZerosPadded; j++)
			data[dataCount +j] = 0;
	}
	avg = sum/dataCount;
	removeDcShift(avg);
	normalize(maxVal, minVal);
}
int setStartMarkerUtterence()	//PLACING THE MARKER AT THE END OF THE VOICE SIGNAL
{
	int markend=0;
	int start=0;
	int k = ceilf(dataCount/frameSize);
	for(int i=0;i<=k;i++)
	{
		if(((2*(energyPerFrame[i] ))<=(energyPerFrame[i+1])) && energyPerFrame[i+1] > 5000 )
		{
			start=i+1;
			return start;
		}
	}
	return 0;
}
int setEndMarkerUtterence(int startMarker)
{
	int k = 0;
	int end=ceilf(dataCount/frameSize)-1;
	for(int i=end;i>=k;i--)
	{
		if(((1.5*(energyPerFrame[i]))<=(energyPerFrame[i-1])) && energyPerFrame[i-1] > 5000)
		{
			end=i;
			if(end <= startMarker)
				return ceilf(dataCount/frameSize)-1;
			return end;
		}
	}
	return end; 
}
void trimUtterence(int startFrameNo,int endFrameNo)
{
	dataCount = 0;
	long double temp[120000]; 
	for(int i=startFrameNo*frameSize; i<(endFrameNo + 1)*frameSize; i++)
	{
		temp[dataCount]=data[i];
		dataCount++;
	}
	for(int i=0; i<dataCount; i++)
		data[i] = temp[i];
	for(int i = dataCount ; i< 20000; i++)
		data[i]=0;
}
inline markers setMarkers()
{
	markers m;
	m.end1 = 0;m.end2=0;m.end3=0;m.start1=0;m.start2=0;m.start3=0;

	//start1
	int k = ceilf(dataCount/frameSize);
	for(int i=0;i<k;i++)
	{
		//if((4*(energyPerFrame[i] ))<=(energyPerFrame[i+1]))
		if((4*(energyPerFrame[i] ))<=(energyPerFrame[i+1]) && ((i*frameSize) > 1500))
		{
			m.start1=i;
			if(energyPerFrame[i+1] > 5000)
			{
				m.start1=i+1;
				break;
			}
		}
		else
			m.start1 = 30;	//hardcoded
	}

	//end1
	m.end1 = m.start1+int(8500/320);
	for(int i = m.end1; i>m.start1; i--)
	{
		if((1.5*(energyPerFrame[i]))<=(energyPerFrame[i-1]))
		{
			m.end1 =i;
			if(energyPerFrame[i-1] > 5000)
				break;
		}
		else
			m.end1 = m.start1 + int(8500/320);
	}

	//end3
	k = 0;
	int end=ceilf(dataCount/frameSize)-1;
	for(int i=end;i>=k;i--)
	{
		if((1.5*(energyPerFrame[i]))<=(energyPerFrame[i-1]))
		{
			m.end3 =i;
			if(energyPerFrame[i-1] > 5000)
				break;
		}
		else
			m.end3 = ceilf(dataCount/frameSize);
	}

	//start3
	m.start3=m.end3-int(8500/320);
	for(int i = m.start3; i<m.end3; i++)
	{
		if((4*(energyPerFrame[i]))<=(energyPerFrame[i+1]))
		{
			m.start3 =i;
			if(energyPerFrame[i+1] > 5000)

				break;
		}
		else
			m.start3 = m.end3 - int(8500/320);
	}

	//middle part
	//start2
	m.start2 = m.end1;
	for(int i = m.start2; i<m.start3; i++)
	{
		if((4*(energyPerFrame[i]))<=(energyPerFrame[i+1]))
		{
			m.start2 =i;
			if(energyPerFrame[i+1] > 5000)
				break;
		}
		else
			m.start2 = m.end1;
	}

	//end2
	m.end2 =m.start3;
	for(int i = m.end2; i>m.start2; i--)
	{
		if((1.5*(energyPerFrame[i]))<=(energyPerFrame[i-1]))
		{
			m.end2 =i;
			if(energyPerFrame[i-1] > 5000)
			{
				m.end2 =i;
				break;
			}
		}
		else
			m.end2 = m.start3;
	}
	return m;
} 
void trimSignal(markers m)
{
	dcOperand1= 0, dcOperand2 = 0, dcOperator=0;
	long double temp[100000];

	//for operand1
	for(int i =0; i<100000; i++)
		temp[i]=0;
	for(int i=m.start1*frameSize; i<(m.end1 + 1)*frameSize; i++)
	{
		temp[dcOperand1]=data[i];
		dcOperand1++;
	}
	for(int i=0; i<dcOperand1; i++)
		operand1[i] = temp[i];
	for(int i = dcOperand1 ; i< 100000; i++)
		operand1[i]=0;

	//for operator
	for(int i =0; i<100000; i++)
		temp[i]=0;
	for(int i=m.start2*frameSize; i<(m.end2 + 1)*frameSize; i++)
	{
		temp[dcOperator]=data[i];
		dcOperator++;
	}
	for(int i=0; i<dcOperator; i++)
		operator1[i] = temp[i];
	for(int i = dcOperator ; i< 100000; i++)
		operator1[i]=0;

	//for operand2
	for(int i =0; i<100000; i++)
		temp[i]=0;
	for(int i=m.start3*frameSize; i<(m.end3 + 1)*frameSize; i++)
	{
		temp[dcOperand2]=data[i];
		dcOperand2++;
	}
	for(int i=0; i<dcOperand2; i++)
		operand2[i] = temp[i];
	for(int i = dcOperand2 ; i< 100000; i++)
		operand2[i]=0;
}
//finds energy for each frame, stores it and returns frame number whose energy was maximum
void findEnergy()
{
	for(int i=0; i<maxRange; i++)
		energyPerFrame[i] = 0;
	int remainder = dataCount%frameSize;
	for( int i=0; i<(frameSize - remainder) ; i++)
		data[dataCount + i] = 0;
	int noOfFrames = ceilf(dataCount/frameSize);
	int maxEnergyFrame = 0, flag = 0;

	for (int i = 0; i < noOfFrames; i++)
	{
		long double energy = 0;
		for(int j = ((i == 0) ? 0 : (frameSize*i)) ;
			j < ((i== (dataCount/frameSize - 1))? dataCount : (frameSize*(i+1)));
			j++)
		{
			energy = (energy + (data[j])*(data[j]));
		}
		energy = energy/frameSize;
		energyPerFrame[i] = energy;
	}
}
//overlapping signal
void overlapping()
{
	int temp[120000], count = 0, k = 0, frameValCount = 0;
	while(k < dataCount)
	{
		if(count >=120000)
			break;
		if(frameValCount==320)
		{
			frameValCount = 0;
			k = k-240;
		}
		temp[count] = data[k];
		count++;
		frameValCount++;
		k++;
	}
	dataCount = count;
	for( k =0; k< dataCount ; k++)
		data[k] = temp[k];
	for( k =dataCount; k< maxRange ; k++)
		data[k] = 0;
}
void overlapping_test(int dc)
{
	dataCount = dc;
	int temp[120000], count = 0, k = 0, frameValCount = 0;
	while(k < dataCount)
	{
		if(count >=120000)
			break;
		if(frameValCount==320)
		{
			frameValCount = 0;
			k = k-240;
		}
		temp[count] = data[k];
		count++;
		frameValCount++;
		k++;
	}
	dataCount = count;
	for( k =0; k< dataCount ; k++)
		data[k] = temp[k];
	for( k =dataCount; k< maxRange ; k++)
		data[k] = 0;
}
void hammingWindow()
{
	long double hw[320];
	for(int i=0; i<320; i++)
		hw[i] = 0.54 - (0.46*(cos((2*3.14159265*i)/319)));
	int count = 0;
	for(int i= 0; i<ceilf(dataCount/frameSize); i++)
	{
		for(int j=0; j<320; j++)
		{
			data[count] = data[count]*hw[j];
			count++;
		}
	}
}
void durbin()
{
	long double E[13];
	long double k[13];
	long double alpha[13][13];
	for(int i = 0; i<=p;i++)
	{	
		E[i]=0;
		k[i]=0;
		for(int j = 0; j<=p; j++)
			alpha[i][j] = 0;
	}
	E[0]=Ri[0];
	int flag = 0;
	for(int i = 1; i<=p ; i++)
	{
		long double sum = 0;
		int j = 1;
		while( j<=(i-1))
		{
			sum += alpha[j][i-1]*Ri[i-j];
			j++;
		}
		k[i] = (Ri[i] - sum )/E[i-1];
		alpha[i][i] = k[i];

		j=1;
		while( j<=(i-1))
		{
			alpha[j][i]= alpha[j][i-1] - (k[i]*alpha[i-j][i-1]);
			j++;
		}
		E[i]= (1- (k[i]*k[i]))*E[i-1];
		flag = 1;
	}
	for(int i=1; i<=p;i++)
	{
		for(int j=1; j<=p;j++)
		{
			if(j==p)
				Ai[i] = alpha[i][p];
		}
	}
}
void findRi()
{
	int lag = 0;
	for(int i=0;i<=p;i++)
		Ri[i]=0;
	for(int i=0;i<=p;i++)
	{
		for(int j=0; j<(320-i);j++)
			Ri[i] += frame[j]*frame[j+i];
	}
}
void findCi()
{
	Ci[0]=Ri[0];
	for(int m = 1; m<=p;m++)
	{
		long double sum=0;
		for(int k=1; k<=(m-1);k++)
			sum += ((float)k/(float)m)*Ci[k]*Ai[m-k];
		Ci[m] = Ai[m] + sum;
	}
}
void applyRSW()
{
	for( int l = 1 ; l<=p; l++)
	{
		Ci[l] = Ci[l]*(1 + (p/2)*(sin((3.14159265*l)/p)));
		f1<<Ci[l]<<" ";
	}
	f1<<endl;
}
//opens each digit and furthers calls functions for removal of dc shift, normalization, applying hamming window, triming signal, finding ri, ai, ci, and applying raised sine window.  
void step1()
{
	int i = 0, j = 0;
	for( i = 0; i< digits ; i++)
	{
		cout<<i<<" ";
		stringstream k1; //converting integer digit to char using k.str()
		k1 << i;
		string folderName = "194101022_" + k1.str();

		//if it is an operator
		if(i>9)
		{
			string operatorName;
			if(i==10)
				operatorName = "plus";
			if(i==11)
				operatorName = "minus";
			if(i==12)
				operatorName = "by";
			if(i==13)
				operatorName = "into";
			folderName = "194101022_" + operatorName;
		}

		for( j=1; j<= utterences; j++)
		{
			stringstream k2;
			k2 << j ;
			string m;
			//string fileName = "194101007_digitsRecordings/" + folderName + "/" + folderName + "_" + k2.str() + ".txt";	
			//string fileName = "Uchiha Itachi/" + folderName + "_" + k2.str() + ".txt";
			string fileName = "194101022_digit/" + folderName + "_" + k2.str() + ".txt";
			//fetching the data signal initially
			f1.open(fileName,ios::in);
			fetchData(i);
			f1.close();

			//finding energy and triming the signal
			int startFrameNo = 0, endFrameNo = 0;
			findEnergy();

			startFrameNo = setStartMarkerUtterence();
			endFrameNo = setEndMarkerUtterence(startFrameNo);
			trimUtterence(startFrameNo ,endFrameNo );

			//apply overlapping
			overlapping();
			hammingWindow();
			string b = "Ci/Ci_" + k1.str() + "_" + k2.str() + ".txt";
			if(i>9)
			{
				string operatorName;
				if(i==10)
					operatorName = "plus";
				if(i==11)
					operatorName = "minus";
				if(i==12)
					operatorName = "by";
				if(i==13)
					operatorName = "into";
				b = "Ci/Ci_" + operatorName + "_" + k2.str() + ".txt";
			}
			f1.open(b, ios:: out);
			long double val = 0;
			for(int l = 0; l<dataCount/frameSize; l++)
			{
				int m = 0;
				while(m<320)
				{
					frame[m] = data[l*frameSize + m];
					m++;
				}
				findRi();
				durbin();
				findCi();
				applyRSW();
			}
			f1.close();
		}
	}
}
void calculateObsSeq(int choseCodebook)
{
	int obs = 0;
	string codebookName;
	if(choseCodebook == 0)	//digitCodebook
		codebookName = "codebook.txt";
	else if(choseCodebook == 1)	//OperatorCodebook
		codebookName = "codebook.txt";
	//codebookName = "joys_operator_codebook.txt";
	f3.open(codebookName, ios::in);
	long double val1, disNew = 0, min = 0;
	int flag = 0;
	for(int i =1 ; i<= 32; i++)
	{
		disNew = 0;
		for(int j=0; j<12; j++)
		{
			f3 >> val1;
			disNew += tw[j]*((Ci[j] - val1)*(Ci[j] - val1)) ;
		}
		if(flag == 0)
		{
			flag = 1;
			obs = i;
			min = disNew;
		}
		else
		{
			if(disNew < min)
			{
				min = disNew;
				obs = i;
			}
		}
	}
	f2<<obs<<" ";
	f3.close();
}
//for making observation sequences
void step2()
{
	for(int i=0; i<digits; i++)
	{
		stringstream k1; //converting integer digit to char using k.str()
		k1 << i;
		for(int j=1; j<=utterences; j++)
		{
			for(int l=1; l<=p; l++)
				Ci[l]=0;
			stringstream k2;
			k2 << j ;
			string fileName1 = "Ci/Ci_" + k1.str() + "_" + k2.str() + ".txt";
			string fileName2 = "obs_seq/obs_" + k1.str() + "_" + k2.str() + ".txt";
			if(i>9)
			{
				string operatorName;
				if(i==10)
					operatorName = "plus";
				if(i==11)
					operatorName = "minus";
				if(i==12)
					operatorName = "by";
				if(i==13)
					operatorName = "into";
				fileName1 = "Ci/Ci_" + operatorName + "_" + k2.str() + ".txt";
				fileName2 = "obs_seq/obs_" + operatorName + "_" + k2.str() + ".txt";
			}
			f1.open(fileName1, ios::in);
			f2.open(fileName2, ios::out);
			long double val;
			int m=0;
			while(f1 >> val)
			{
				if(m == 12)
				{
					if(i>9)
						calculateObsSeq(1);
					else
						calculateObsSeq(0);

					f2<<endl;
					m=0;
				}
				Ci[m] = val;
				m++;
			}
			f1.close();
			f2.close();
		}
	}
}
//fetching initial model
void setMatrices(int digit_no, int utterence_no, int iteration_no)
{
	stringstream k1;
	k1 <<digit_no;
	string fileName;
	if(iteration_no ==0)
		fileName  = "hmm_training/a.txt";
	else
	{
		fileName  = "hmm_training/a_" + k1.str() + ".txt";
		if(digit_no>9)
		{
			string operatorName;
			if(digit_no==10)
				operatorName = "plus";
			if(digit_no==11)
				operatorName = "minus";
			if(digit_no==12)
				operatorName = "by";
			if(digit_no==13)
				operatorName = "into";
			fileName  = "hmm_training/a_" + operatorName + ".txt";
		}
	}
	long double val=0;
	int count = 0;
	fi.open(fileName,ios::in);
	int i=1;
	int j=1;

	//a matrix
	while(fi>>val)
	{
		if(j==6)
		{
			j=1;
			i++;
		}
		if(count==(N*N))
			break;
		count++;
		a[i][j]=val;
		j++;
	}
	fi.close();

	//b matrix
	if(iteration_no ==0)
		fileName  = "hmm_training/b.txt";
	else
	{
		fileName  = "hmm_training/b_" + k1.str() + ".txt";
		if(digit_no>9)
		{
			string operatorName;
			if(digit_no==10)
				operatorName = "plus";
			if(digit_no==11)
				operatorName = "minus";
			if(digit_no==12)
				operatorName = "by";
			if(digit_no==13)
				operatorName = "into";
			fileName  = "hmm_training/b_" + operatorName + ".txt";
		}
	}

	fi.open(fileName,ios::in);
	i=1;
	j=1;
	count=0;
	while(fi>>val)
	{
		if(j==M+1)
		{
			j=1;
			i++;
		}
		if(count==(N*M))
			break;
		count++;
		b[i][j]=val;
		j++;
	}
	fi.close();

	//obs values
	stringstream k2;
	k2<<utterence_no;
	fileName  = "obs_seq/obs_" + k1.str() + "_" + k2.str() + ".txt";
	if(digit_no>9)
	{
		string operatorName;
		if(digit_no==10)
			operatorName = "plus";
		if(digit_no==11)
			operatorName = "minus";
		if(digit_no==12)
			operatorName = "by";
		if(digit_no==13)
			operatorName = "into";
		fileName  = "obs_seq/obs_" + operatorName + "_" + k2.str() + ".txt";
	}
	fi.open(fileName,ios::in);
	count = 1;
	int val1=0;
	while(fi>>val1)
	{
		if(count >=160)
			break;
		obs[count]=val1;
		count++;
	}
	T = count-1;
	fi.close();
}
//forward procedure for calculation of alpha
void forwardProcedure()
{
	int i,j,t;
	long double sum =0.0;
	//initialization
	for(i=1; i<=N; i++)
	{
		alpha[1][i]= pi[i]*b[i][obs[1]];
	}
	//recursion
	for(t=1;t<T;t++)
	{
		for(j = 1; j<=N; j++)
		{
			sum = 0;
			for(i =1; i<=N; i++)
			{
				sum += alpha[t][i]*a[i][j];
			}
			alpha[t+1][j]=sum*b[j][obs[t+1]];
		}
	}
	for(i=1; i<=N; i++)
	{
		prob1 += alpha[T][i];
	}
	cout<<"prob : "<<prob1<<endl;
}
//backward procedure for calculation of beta
void backwardProcedure()
{
	int i,j,t;
	for(i = 1; i<=N; i++)
		beta[T][i]=1;
	long double sum = 0.0;
	for(t = T-1; t>0;t--)
	{
		for(i=1;i<=N;i++)
		{
			sum=0;
			for(j=1; j<=N; j++)
			{
				sum +=a[i][j]*b[j][obs[t+1]]*beta[t+1][j];
			}
			beta[t][i]=sum;
		}
	}
}
//calculating gamma and viterbi algorithm for p* and q*
void viterbiAlgorithm()
{
	int i,j,t;
	long double m = 1.0;
	//gamma calculation
	for(t =1; t<=T; t++)
	{
		for(j=1;j<=N;j++)
		{
			m = alpha[t][j]*beta[t][j];
			gamma[t][j] = m / prob1;
		}
	}

	//initialization
	for(i=1;i<=N;i++)
	{
		delta[1][i]=pi[i]*b[i][obs[1]];
		psy[1][i] = 0;
	}
	long double maxIndex = 0;
	long double maxVal = 0;

	//recursion
	for(t=2;t<=T;t++)
	{
		for(j =1;j<=N;j++)
		{
			maxVal=-1;
			for(i=1;i<=N;i++)
			{
				if((delta[t-1][i]*a[i][j]) > maxVal)
				{
					maxVal = delta[t-1][i]*a[i][j];
					maxIndex = i;
				}
			}
			delta[t][j]=maxVal*b[j][obs[t]];
			psy[t][j] = maxIndex;
		}
	}

	maxIndex = 0;
	maxVal = 0;
	for(i=1;i<=N;i++)
	{
		if(delta[T][i] > maxVal)
		{
			maxVal = delta[T][i];
			maxIndex=i;
		}
	}
	pStar=maxVal;
	qStar[T] = maxIndex;
	fo<<endl<<"p* : "<<pStar<<endl;
	printf("p* : %E\n",pStar);
	for(t=T-1; t>0;t--)
	{
		qStar[t] = psy[t+1][qStar[t+1]];
	}
	cout<<"q* = ";
	fo<<"q* : ";
	for(t=1;t<=T;t++)
	{
		fo<<qStar[t]<<" ";
		cout<<qStar[t]<<" ";
	}
	cout<<endl<<endl;
}
//re-estimation for calculation of a and b matrices
void reestimation()
{
	int i,j,t;
	for(t=1;t<T;t++)
	{
		for(i = 1; i<=N; i++)
		{
			for(j=1; j<=N; j++)
			{
				et[t][i][j]=(alpha[t][i]*a[i][j]*b[j][obs[t+1]]*beta[t+1][j])/prob1;
			}
		}
	}
	for(i =1; i<N; i++)
	{
		newpi[i] = gamma[1][i];
	}
	long double sum = 0;

	//Computing expected number of transitions from state Si
	for(i=1;i<=N;i++)
	{
		sum = 0;
		for(t=1; t<T;t++)
		{
			sum +=gamma[t][i];
		}
		gammaSum[i] = sum;
	}

	//Computing expected number of transitions from state Si and state Sj
	for(i=1;i<=N;i++)
	{
		for(j=1;j<=N;j++)
		{
			sum=0;
			for(t=1; t<T;t++)
			{
				sum +=et[t][i][j];
			}
			XiSum[i][j] = sum;
		}
	}
	for(i =1; i<=N; i++)
	{
		for(j=1; j<=N; j++)
		{
			newa[i][j] = XiSum[i][j]/gammaSum[i];
		}
	}
	long double threshold = pow(10.0,-30);
	long double rowsum = 0;
	long double maxvalue=0 ;
	int maxvalindex = 1;
	int count =0;
	for(j=1; j<=N; j++)
	{
		maxvalindex = 1;
		maxvalue = -1.0;
		count = 0;
		for(i = 1; i<=M; i++)
		{
			long double gSum = 0;
			long double sum=0.0;
			for(t=1; t<=T; t++)
			{
				if((i)==obs[t])
					sum+= gamma[t][j];
				gSum += gamma[t][j];
			}
			if( sum/gSum < threshold)
			{
				sum = threshold;
				gSum=1;
				count++;
			}
			if( (sum/gSum) > maxvalue )
			{
				maxvalue =  (sum/gSum) ;
				maxvalindex = i;
			}
			newb[j][i] = sum/gSum;
		}
		newb[j][maxvalindex] = newb[j][maxvalindex] - count*threshold;
	}

	//updation of a matrix
	for(i =1; i<=N; i++)
	{
		for(j=1; j<=N; j++)
		{
			a[i][j]=newa[i][j];
			//printf("%E ", a[j][i]);	//uncomment to print b matrix after every iteration
		}
		//cout<<endl;
	}

	//updation of b matrix
	//cout<<endl;
	for(j=1; j<=N; j++)
	{
		for(i = 1; i<=M; i++)
		{
			b[j][i] = newb[j][i];
			//printf("%E ", b[j][i]);	//uncomment to print b matrix after every iteration
		}
		//cout<<endl;
	}
}
void initialization()
{
	for(int i =0; i< 500; i++)
	{
		for(int j = 0; j<6; j++)
		{
			alpha[i][j]=0;
			beta[i][j]=0;
			delta[i][j]=0;
			psy[i][j]=0;
			gamma[i][j]=0;
		}
		qStar[i]=0;
	}
	for(int k = 0; k<500; k++)
	{
		for(int i =0; i< 6; i++)
		{
			for(int j = 0; j<6; j++)
			{
				et[k][i][j]=0;
				newa[i][j] = 0;
				XiSum[i][j]= 0;
			}
			newpi[i]=0;
			for(int j = 0; j<33; j++)
			{
				newb[i][j]=0;
			}
			gammaSum[i] = 0;
		}
	}
	prob1 = 0.0;
	pStar = 0;
}
//training hmm
void step3()
{
	int loopCount = 0;
	for(int i=0; i<digits; i++)
	{
		stringstream ki;
		ki<< i;
		for(int j = 0; j<iteration1 ; j++)
		{	
			stringstream kj;
			kj<< j;
			for(int m =0; m<=N; m++)
				for(int n=0; n<=N; n++)
					a_avg[m][n]=0;
			for(int m=0; m<=N; m++)
				for(int n = 0; n<=M; n++)
					b_avg[m][n] = 0;
			for(int k = 1; k<=utterences; k++)
			{
				cout<<endl<<"********************"<<i<<"_"<<j<<"_"<<k<<"********************"<<endl;
				stringstream kk;
				kk<< k;
				for(int l=0; l<20; l++)
				{
					initialization();
					loopCount++;

					setMatrices(i,k, l );
					forwardProcedure();
					backwardProcedure();
					viterbiAlgorithm();
					reestimation();

					string fileName1 = "hmm_training/a_" + ki.str() + ".txt";
					if(i>9)
					{
						string operatorName;
						if(i==10)
							operatorName = "plus";
						if(i==11)
							operatorName = "minus";
						if(i==12)
							operatorName = "by";
						if(i==13)
							operatorName = "into";
						fileName1 = "hmm_training/a_" + operatorName + ".txt";
					}
					fo.open(fileName1,ios::out);
					fo.precision(6);
					for(int m =1; m<=N; m++)
					{
						for(int n=1; n<=N; n++)
							fo << a[m][n]<<"	";
						fo.scientific;
						fo<<endl;
					}
					fo.close();

					fileName1 = "hmm_training/b_" + ki.str() + ".txt";
					if(i>9)
					{
						string operatorName;
						if(i==10)
							operatorName = "plus";
						if(i==11)
							operatorName = "minus";
						if(i==12)
							operatorName = "by";
						if(i==13)
							operatorName = "into";
						fileName1 = "hmm_training/b_" + operatorName + ".txt";
					}
					fo.open(fileName1,ios::out);
					fo.precision(6);
					for(int m=1; m<=N; m++)
					{
						for(int n = 1; n<=M; n++)
							fo << b[m][n]<< "	";
						fo.scientific;
						fo<<endl;
					}
					fo.close();
				}
				for(int m =1; m<=N; m++)
					for(int n=1; n<=N; n++)
						a_avg[m][n] += a[m][n];
				for(int m=1; m<=N; m++)
					for(int n = 1; n<=M; n++)
						b_avg[m][n] += b[m][n];
			}
			string fileName1 = "hmm_training/a_" + ki.str() + ".txt";
			if(i>9)
			{
				string operatorName;
				if(i==10)
					operatorName = "plus";
				if(i==11)
					operatorName = "minus";
				if(i==12)
					operatorName = "by";
				if(i==13)
					operatorName = "into";
				fileName1 = "hmm_training/a_" + operatorName + ".txt";
			}
			fo.open(fileName1,ios::out);
			fo.precision(6);
			for(int m =1; m<=N; m++)
			{
				for(int n=1; n<=N; n++)
				{
					a_avg[m][n] /= utterences;
					a[m][n] = a_avg[m][n];
					fo.scientific;
					fo << a[m][n]<<"	";
				}
				fo<<endl;
			}
			fo.close();

			fileName1 = "hmm_training/b_" + ki.str() + ".txt";
			if(i>9)
			{
				string operatorName;
				if(i==10)
					operatorName = "plus";
				if(i==11)
					operatorName = "minus";
				if(i==12)
					operatorName = "by";
				if(i==13)
					operatorName = "into";
				fileName1 = "hmm_training/b_" + operatorName + ".txt";
			}
			fo.open(fileName1,ios::out);
			fo.precision(6);
			for(int m=1; m<=N; m++)
			{
				for(int n = 1; n<=M; n++)
				{
					b_avg[m][n] /= utterences;
					b[m][n] = b_avg[m][n];
					fo.scientific;
					fo << b[m][n]<<"	";
				}
				fo<<endl;
			}
			fo.close();
		}
		string fileName = "final/final_model_a_" + ki.str() + ".txt";
		if(i>9)
		{
			string operatorName;
			if(i==10)
				operatorName = "plus";
			if(i==11)
				operatorName = "minus";
			if(i==12)
				operatorName = "by";
			if(i==13)
				operatorName = "into";
			fileName = "final/final_model_a_" + operatorName + ".txt";
		}
		f1.open(fileName, ios::out);
		for(int i =1; i<=N; i++)
		{
			for(int j=1; j<=N; j++)
				f1 << a_avg[i][j]<<" ";
			f1.scientific;
			f1<<endl;
		}
		f1.close();
		fileName = "final/final_model_b_" + ki.str() + ".txt";
		if(i>9)
		{
			string operatorName;
			if(i==10)
				operatorName = "plus";
			if(i==11)
				operatorName = "minus";
			if(i==12)
				operatorName = "by";
			if(i==13)
				operatorName = "into";
			fileName = "final/final_model_b_" + operatorName + ".txt";
		}
		f1.open(fileName, ios::out);
		for(int i =1; i<=N; i++)
		{
			for(int j=1; j<=M; j++)
				f1 << b_avg[i][j]<<"	";
			f1.scientific;
			f1 <<endl;
		}
		f1.close();
	}
}
void setMatrices_test(int digit_no_file, int digit_no_model,int utterence_no)
{
	stringstream k1;
	k1 <<digit_no_model;
	string fileName;

	for(int  i=0; i<6; i++)
	{
		for (int j = 0; j<6; j++)
			a[i][j]=0;
		for (int j = 0; j<33; j++)
			b[i][j]=0;
	}
	fileName  = "final/final_model_a_" + k1.str() + ".txt";
	if(digit_no_model>9)
	{
		string operatorName;
		if(digit_no_model==10)
			operatorName = "plus";
		if(digit_no_model==11)
			operatorName = "minus";
		if(digit_no_model==12)
			operatorName = "by";
		if(digit_no_model==13)
			operatorName = "into";
		fileName  = "final/final_model_a_" + operatorName + ".txt";
	}
	long double val=0;
	int count = 0;
	fi.open(fileName,ios::in);
	int i=1;
	int j=1;

	//a matrix
	while(fi>>val)
	{
		if(j==6)
		{
			j=1;
			i++;
		}
		if(count==(N*N))
			break;
		count++;
		a[i][j]=val;
		j++;
	}
	fi.close();

	//b matrix
	fileName  = "final/final_model_b_" + k1.str() + ".txt";
	if(digit_no_model>9)
	{
		string operatorName;
		if(digit_no_model==10)
			operatorName = "plus";
		if(digit_no_model==11)
			operatorName = "minus";
		if(digit_no_model==12)
			operatorName = "by";
		if(digit_no_model==13)
			operatorName = "into";
		fileName  = "final/final_model_b_" + operatorName + ".txt";
	}
	fi.open(fileName,ios::in);
	i=1;
	j=1;
	count=0;
	while(fi>>val)
	{
		if(j==M+1)
		{
			j=1;
			i++;
		}
		if(count==(N*M))
			break;
		count++;
		b[i][j]=val;
		j++;
	}
	fi.close();

	for(int  i=0; i<500; i++)
		obs[i]=0;
	//obs values
	stringstream k2, k3;

	k3 << digit_no_file;
	k2<<utterence_no+20;
	fileName  = "test_obs_seq/obs_" + k3.str() + "_" + k2.str() + ".txt";
	if(digit_no_file>9)
	{
		string operatorName;
		if(digit_no_file==10)
			operatorName = "plus";
		if(digit_no_file==11)
			operatorName = "minus";
		if(digit_no_file==12)
			operatorName = "by";
		if(digit_no_file==13)
			operatorName = "into";
		fileName  = "test_obs_seq/obs_" + operatorName + "_" + k2.str() + ".txt";
	}
	fi.open(fileName,ios::in);
	count = 1;
	int val1=0;
	while(fi>>val1)
	{
		if(count >=160) 
			break;
		obs[count]=val1;
		count++;
	}
	T = count-1;
	fi.close();
}
//testing recorded
void step4()
{
	int i = 0, j = 0;
	for( i = 0; i< digits ; i++)
	{
		stringstream k1; //converting integer digit to char using k.str()
		k1 << i;
		string folderName = "194101022_" + k1.str();
		if(i>9)
		{
			string operatorName;
			if(i==10)
				operatorName = "plus";
			if(i==11)
				operatorName = "minus";
			if(i==12)
				operatorName = "by";
			if(i==13)
				operatorName = "into";
			folderName = "194101022_" + operatorName;
		}
		for( j=1; j<= testUtterences; j++)
		{
			stringstream k2;
			k2 << j+20 ;
			//string fileName = "194101007_digitsRecordings/" + folderName + "/" + folderName + "_" + k2.str() + ".txt";
			//string fileName = "Uchiha Itachi/" + folderName + "_" + k2.str() + ".txt";
			string fileName = "194101022_digit/" + folderName + "_" + k2.str() + ".txt";

			//fetching the data signal initially
			f1.open(fileName,ios::in);
			fetchData(i);
			f1.close();

			//finding energy and triming the signal
			int startFrameNo = 0, endFrameNo = 0;
			findEnergy();
			startFrameNo = setStartMarkerUtterence();
			endFrameNo = setEndMarkerUtterence(startFrameNo);
			trimUtterence(startFrameNo ,endFrameNo );

			//apply overlapping
			overlapping();

			hammingWindow();
			string b = "test_Ci/Ci_" + k1.str() + "_" + k2.str() + ".txt";
			if(i>9)
			{
				string operatorName;
				if(i==10)
					operatorName = "plus";
				if(i==11)
					operatorName = "minus";
				if(i==12)
					operatorName = "by";
				if(i==13)
					operatorName = "into";
				b = "test_Ci/Ci_" + operatorName + "_" + k2.str() + ".txt";
			}

			f1.open(b, ios:: out);
			long double val = 0;
			for(int l = 0; l<dataCount/frameSize; l++)
			{
				int m = 0;
				while(m<320)
				{
					frame[m] = data[l*frameSize + m];
					m++;
				}
				findRi();
				durbin();
				findCi();
				applyRSW();
			}
			f1.close();
		}
	}

	for( i=0; i<digits; i++)
	{
		stringstream k1; //converting integer digit to char using k.str()
		k1 << i;
		for(j=1; j<=testUtterences; j++)
		{
			for(int l=1; l<=p; l++)
				Ci[l]=0;
			stringstream k2;
			k2 << j+20 ;
			string fileName1 = "test_Ci/Ci_" + k1.str() + "_" + k2.str() + ".txt";
			string fileName2 = "test_obs_seq/obs_" + k1.str() + "_" + k2.str() + ".txt";
			if(i>9)
			{
				string operatorName;
				if(i==10)
					operatorName = "plus";
				if(i==11)
					operatorName = "minus";
				if(i==12)
					operatorName = "by";
				if(i==13)
					operatorName = "into";
				fileName1 = "test_Ci/Ci_" + operatorName + "_" + k2.str() + ".txt";
				fileName2 = "test_obs_seq/obs_" + operatorName + "_" + k2.str() + ".txt";
			}
			f1.open(fileName1, ios::in);
			f2.open(fileName2, ios::out);
			long double val;
			int n=0;
			while(f1 >> val)
			{
				if(n == 12)
				{
					calculateObsSeq(0);
					f2<<endl;
					n=0;
				}
				Ci[n] = val;
				n++;
			}
			f1.close();
			f2.close();
		}
	}
	cout<<endl;
	long double p[14];
	int digitPredicted = -1;
	int predictedCount = 0;
	for(i=0; i<digits ; i++)
	{
		stringstream ki;
		ki<< i;
		for(int k = 1; k<=testUtterences; k++)
		{
			for(j=0; j<10; j++)
				p[j]=0;
			digitPredicted =-1;
			if(i>9)
			{
				string operatorName;
				if(i==10)
					operatorName = "plus";
				if(i==11)
					operatorName = "minus";
				if(i==12)
					operatorName = "by";
				if(i==13)
					operatorName = "into";
				cout<<endl<<"file: digit_"<<operatorName<<"_utterence_"<<k<<" - "<<endl;
			}
			else
				cout<<endl<<"file: digit_"<<i<<"_utterence_"<<k<<" - "<<endl;
			stringstream kk;
			kk<< k;
			for(j = 0; j<digits; j++)
			{
				initialization();
				setMatrices_test(i,j,k );
				cout<<j<<" ";
				forwardProcedure();
				p[j]= prob1;
				cout<<"";
			}
			long double maxP = -1.0;
			for(j=0; j<digits;j++)
			{
				if(p[j]>maxP)
				{
					maxP = p[j];
					digitPredicted = j;
				}
			}
			if(digitPredicted>9)
			{
				string operatorName;
				if(digitPredicted==10)
					operatorName = "plus";
				if(digitPredicted==11)
					operatorName = "minus";
				if(digitPredicted==12)
					operatorName = "by";
				if(digitPredicted==13)
					operatorName = "into";
				cout<<"detected: "<<operatorName<<endl;
			}
			else
				cout<<"detected: "<<digitPredicted<<endl;
			if(digitPredicted == i)
				predictedCount++;
		}
	}
	float accuracy = (predictedCount*100)/140;
	cout<<endl<<"Accuracy: "<<accuracy<<"%"<<endl;
}
void setMatrices_test_live(int digit_no_model, int e)
{
	stringstream k1;
	k1 <<digit_no_model;
	string fileName;

	for(int  i=0; i<6; i++)
	{
		for (int j = 0; j<6; j++)
			a[i][j]=0;
		for (int j = 0; j<33; j++)
			b[i][j]=0;
	}
	if(e==1)
	{
		string operatorName;
		if(digit_no_model==0)
			operatorName = "plus";
		if(digit_no_model==1)
			operatorName = "minus";
		if(digit_no_model==2)
			operatorName = "by";
		if(digit_no_model==3)
			operatorName = "into";
		fileName  = "final/final_model_a_" + operatorName + ".txt";
	}
	else 
	{
		fileName  = "final/final_model_a_" + k1.str() + ".txt";
	}
	long double val=0;
	int count = 0;
	fi.open(fileName,ios::in);
	int i=1;
	int j=1;

	//a matrix
	while(fi>>val)
	{
		if(j==6)
		{
			j=1;
			i++;
		}
		if(count==(N*N))
			break;
		count++;
		a[i][j]=val;
		j++;
	}
	fi.close();

	//b matrix
	fileName  = "final/final_model_b_" + k1.str() + ".txt";
	if(digit_no_model>9)
	{
		string operatorName;
		if(digit_no_model==10)
			operatorName = "plus";
		if(digit_no_model==11)
			operatorName = "minus";
		if(digit_no_model==12)
			operatorName = "by";
		if(digit_no_model==13)
			operatorName = "into";
		fileName  = "final/final_model_b_" + operatorName + ".txt";
	}
	fi.open(fileName,ios::in);
	i=1;
	j=1;
	count=0;
	while(fi>>val)
	{
		if(j==M+1)
		{
			j=1;
			i++;
		}
		if(count==(N*M))
			break;
		count++;
		b[i][j]=val;
		j++;
	}
	fi.close();

	for(int  i=0; i<500; i++)
		obs[i]=0;
	//obs values
	fileName  = "test_obs_seq.txt";
	fi.open(fileName,ios::in);
	count = 1;
	int val1=0;
	while(fi>>val1)
	{
		if(count >=160) 
			break;
		obs[count]=val1;
		count++;
	}
	T = count-1;
	fi.close();
}
void step5(int* op1, char* op, int * op2)
{

	getchar();
	for(int i=0; i<1; i++)
	{
		initialization();
		system("recording_module.exe 3 input_file.wav input_file.txt");
		f1.open("input_file.txt",ios::in);
		fetchData(12);
		f1.close();
		//finding energy and triming the signal
		int startFrameNo = 0, endFrameNo = 0;
		findEnergy();
		markers m = setMarkers();
		trimSignal(m);


		//apply overlapping
		for(int e=0; e<3; e++)
		{
			for(int i=0; i<maxRange; i++)
				data[i]=0;
			if(e==0)
			{
				for(int i=0; i<dcOperand1; i++)
					data[i] = operand1[i];
				overlapping_test(dcOperand1);
			}
			else if(e==1)
			{
				for(int i=0; i<dcOperator; i++)
					data[i] = operator1[i];
				overlapping_test(dcOperator);
			}
			else if(e==2)
			{
				for(int i=0; i<dcOperand2; i++)
					data[i] = operand2[i];
				overlapping_test(dcOperand2);
			}


			hammingWindow();
			f1.open("test.txt", ios:: out);
			long double val = 0;
			for(int l = 0; l<dataCount/frameSize; l++)
			{
				int m = 0;
				while(m<320)
				{
					frame[m] = data[l*frameSize + m];
					m++;
				}
				findRi();
				durbin();
				findCi();
				applyRSW();
			}
			f1.close();

			for(int l=1; l<=p; l++)
				Ci[l]=0;
			string fileName1 = "test.txt";
			string fileName2 = "test_obs_seq.txt";
			f1.open(fileName1, ios::in);
			f2.open(fileName2, ios::out);
			val=0;
			int n=0;
			while(f1 >> val)
			{
				if(n ==12)
				{
					calculateObsSeq(0);
					f2<<endl;
					n=0;
				}
				Ci[n] = val;
				n++;
			}
			f1.close();
			f2.close();

			cout<<endl;
			long double pOperand[10], pOperator[4];
			int digitPredicted;
			for(int j=0; j<10; j++)
				pOperand[j]=0;
			for(int j=0; j<4; j++)
				pOperator[j]=0;
			digitPredicted =-1;

			if(e==1)
			{
				for(int j = 0; j<digits-10; j++)
				{
					initialization();
					setMatrices_test_live(j,e);
					forwardProcedure();
					pOperator[j]= prob1;
				}
				long double maxP = -1.0;
				for(int j=0; j<digits-10;j++)
				{
					if(pOperator[j]>maxP)
					{
						maxP = pOperator[j];
						digitPredicted = j;
					}
				}

				char operatorName;
				if(digitPredicted==0)
					operatorName = '+';
				if(digitPredicted==1)
					operatorName = '-';
				if(digitPredicted==2)
					operatorName = '/';
				if(digitPredicted==3)
					operatorName = '*';
				result2 = operatorName;
			}
			else
			{
				for(int j = 0; j<digits-4; j++)
				{
					initialization();
					setMatrices_test_live(j,e);
					forwardProcedure();
					pOperand[j]= prob1;
				}
				long double maxP = -1.0;
				for(int j=0; j<digits-4;j++)
				{
					if(pOperand[j]>maxP)
					{
						maxP = pOperand[j];
						digitPredicted = j;
					}
				}
				if(e==0)
					result1 = digitPredicted;
				else if(e==2)
					result3 = digitPredicted;
			}
		}
		*op1=result1;
		*op=result2;
		*op2=result3;
		cout<<endl<<result1<<" "<<result2<<" "<<result3<<endl;

	}
}

void live_train()
{
	int choice;
	cout<<endl<<"Enter digit: 0-9/plus-10/minus-11/by-12/into-13: ";
	cin>>choice;
	stringstream kc;
	kc<<choice;
	string fileName;
	if(choice<=9)
	{
		string path = "194101022_digit/194101022_" + kc.str() + "_";
		for(int i=0; i<10; i++)
		{

			cout<<endl<<"utterence: "<<i;
			initialization();
			stringstream ki;
			ki<<i+1;
			fileName = path + ki.str() + ".txt";
			system(("recording_module.exe 3 input_file.wav "+fileName).c_str());
		}
	}
	else
	{
		string operatorName;
		if(choice==10)
			operatorName = "plus";
		if(choice==11)
			operatorName = "minus";
		if(choice==12)
			operatorName = "by";
		if(choice==13)
			operatorName = "into";
		string path = "new/194101022_" + operatorName + "_";
		for(int i=0; i<10; i++)
		{
			cout<<endl<<"utterence: "<<i;
			initialization();
			stringstream ki;
			ki<<i+1;
			fileName = path + ki.str() + ".txt";
			system(("recording_module.exe 3 input_file.wav "+fileName).c_str());
		}
	}
	int i = 0, j = 0;
	for( i = choice; i< choice+1 ; i++)
	{
		cout<<i<<" ";
		stringstream k1; //converting integer digit to char using k.str()
		k1 << i;
		string folderName = "194101022_" + k1.str();

		//if it is an operator
		if(i>9)
		{
			string operatorName;
			if(i==10)
				operatorName = "plus";
			if(i==11)
				operatorName = "minus";
			if(i==12)
				operatorName = "by";
			if(i==13)
				operatorName = "into";
			folderName = "194101022_" + operatorName;
		}

		for( j=1; j<= 10; j++)
		{
			stringstream k2;
			k2 << j ;
			string m;
			//string fileName = "194101007_digitsRecordings/" + folderName + "/" + folderName + "_" + k2.str() + ".txt";	
			//string fileName = "Uchiha Itachi/" + folderName + "_" + k2.str() + ".txt";
			string fileName = "194101022_digit/" + folderName + "_" + k2.str() + ".txt";
			//fetching the data signal initially
			f1.open(fileName,ios::in);
			fetchData(i);
			f1.close();

			//finding energy and triming the signal
			int startFrameNo = 0, endFrameNo = 0;
			findEnergy();

			startFrameNo = setStartMarkerUtterence();
			endFrameNo = setEndMarkerUtterence(startFrameNo);
			trimUtterence(startFrameNo ,endFrameNo );

			//apply overlapping
			overlapping();
			hammingWindow();
			string b = "Ci/Ci_" + k1.str() + "_" + k2.str() + ".txt";
			if(i>9)
			{
				string operatorName;
				if(i==10)
					operatorName = "plus";
				if(i==11)
					operatorName = "minus";
				if(i==12)
					operatorName = "by";
				if(i==13)
					operatorName = "into";
				b = "Ci/Ci_" + operatorName + "_" + k2.str() + ".txt";
			}
			f1.open(b, ios:: out);
			long double val = 0;
			for(int l = 0; l<dataCount/frameSize; l++)
			{
				int m = 0;
				while(m<320)
				{
					frame[m] = data[l*frameSize + m];
					m++;
				}
				findRi();
				durbin();
				findCi();
				applyRSW();
			}
			f1.close();
		}
	}

	for(int i=choice; i<choice+1; i++)
	{
		stringstream k1; //converting integer digit to char using k.str()
		k1 << i;
		for(int j=1; j<=10; j++)
		{
			for(int l=1; l<=p; l++)
				Ci[l]=0;
			stringstream k2;
			k2 << j ;
			string fileName1 = "Ci/Ci_" + k1.str() + "_" + k2.str() + ".txt";
			string fileName2 = "obs_seq/obs_" + k1.str() + "_" + k2.str() + ".txt";
			if(i>9)
			{
				string operatorName;
				if(i==10)
					operatorName = "plus";
				if(i==11)
					operatorName = "minus";
				if(i==12)
					operatorName = "by";
				if(i==13)
					operatorName = "into";
				fileName1 = "Ci/Ci_" + operatorName + "_" + k2.str() + ".txt";
				fileName2 = "obs_seq/obs_" + operatorName + "_" + k2.str() + ".txt";
			}
			f1.open(fileName1, ios::in);
			f2.open(fileName2, ios::out);
			long double val;
			int m=0;
			while(f1 >> val)
			{
				if(m == 12)
				{
					if(i>9)
						calculateObsSeq(1);
					else
						calculateObsSeq(0);

					f2<<endl;
					m=0;
				}
				Ci[m] = val;
				m++;
			}
			f1.close();
			f2.close();
		}
	}

	for(int i=choice; i<choice+1; i++)
	{
		stringstream ki;
		ki<< i;
		for(int j = 0; j<iteration1 ; j++)
		{	
			stringstream kj;
			kj<< j;
			for(int m =0; m<=N; m++)
				for(int n=0; n<=N; n++)
					a_avg[m][n]=0;
			for(int m=0; m<=N; m++)
				for(int n = 0; n<=M; n++)
					b_avg[m][n] = 0;
			for(int k = 1; k<=10; k++)
			{
				cout<<endl<<"********************"<<i<<"_"<<j<<"_"<<k<<"********************"<<endl;
				stringstream kk;
				kk<< k;
				for(int l=0; l<20; l++)
				{
					initialization();
					setMatrices(i,k, l );
					forwardProcedure();
					backwardProcedure();
					viterbiAlgorithm();
					reestimation();

					string fileName1 = "hmm_training/a_" + ki.str() + ".txt";
					if(i>9)
					{
						string operatorName;
						if(i==10)
							operatorName = "plus";
						if(i==11)
							operatorName = "minus";
						if(i==12)
							operatorName = "by";
						if(i==13)
							operatorName = "into";
						fileName1 = "hmm_training/a_" + operatorName + ".txt";
					}
					fo.open(fileName1,ios::out);
					fo.precision(6);
					for(int m =1; m<=N; m++)
					{
						for(int n=1; n<=N; n++)
							fo << a[m][n]<<"	";
						fo.scientific;
						fo<<endl;
					}
					fo.close();

					fileName1 = "hmm_training/b_" + ki.str() + ".txt";
					if(i>9)
					{
						string operatorName;
						if(i==10)
							operatorName = "plus";
						if(i==11)
							operatorName = "minus";
						if(i==12)
							operatorName = "by";
						if(i==13)
							operatorName = "into";
						fileName1 = "hmm_training/b_" + operatorName + ".txt";
					}
					fo.open(fileName1,ios::out);
					fo.precision(6);
					for(int m=1; m<=N; m++)
					{
						for(int n = 1; n<=M; n++)
							fo << b[m][n]<< "	";
						fo.scientific;
						fo<<endl;
					}
					fo.close();
				}
				for(int m =1; m<=N; m++)
					for(int n=1; n<=N; n++)
						a_avg[m][n] += a[m][n];
				for(int m=1; m<=N; m++)
					for(int n = 1; n<=M; n++)
						b_avg[m][n] += b[m][n];
			}
			string fileName1 = "hmm_training/a_" + ki.str() + ".txt";
			if(i>9)
			{
				string operatorName;
				if(i==10)
					operatorName = "plus";
				if(i==11)
					operatorName = "minus";
				if(i==12)
					operatorName = "by";
				if(i==13)
					operatorName = "into";
				fileName1 = "hmm_training/a_" + operatorName + ".txt";
			}
			fo.open(fileName1,ios::out);
			fo.precision(6);
			for(int m =1; m<=N; m++)
			{
				for(int n=1; n<=N; n++)
				{
					a_avg[m][n] /= utterences;
					a[m][n] = a_avg[m][n];
					fo.scientific;
					fo << a[m][n]<<"	";
				}
				fo<<endl;
			}
			fo.close();

			fileName1 = "hmm_training/b_" + ki.str() + ".txt";
			if(i>9)
			{
				string operatorName;
				if(i==10)
					operatorName = "plus";
				if(i==11)
					operatorName = "minus";
				if(i==12)
					operatorName = "by";
				if(i==13)
					operatorName = "into";
				fileName1 = "hmm_training/b_" + operatorName + ".txt";
			}
			fo.open(fileName1,ios::out);
			fo.precision(6);
			for(int m=1; m<=N; m++)
			{
				for(int n = 1; n<=M; n++)
				{
					b_avg[m][n] /= utterences;
					b[m][n] = b_avg[m][n];
					fo.scientific;
					fo << b[m][n]<<"	";
				}
				fo<<endl;
			}
			fo.close();
		}
		string fileName = "final/final_model_a_" + ki.str() + ".txt";
		if(i>9)
		{
			string operatorName;
			if(i==10)
				operatorName = "plus";
			if(i==11)
				operatorName = "minus";
			if(i==12)
				operatorName = "by";
			if(i==13)
				operatorName = "into";
			fileName = "final/final_model_a_" + operatorName + ".txt";
		}
		f1.open(fileName, ios::out);
		for(int i =1; i<=N; i++)
		{
			for(int j=1; j<=N; j++)
				f1 << a_avg[i][j]<<" ";
			f1.scientific;
			f1<<endl;
		}
		f1.close();
		fileName = "final/final_model_b_" + ki.str() + ".txt";
		if(i>9)
		{
			string operatorName;
			if(i==10)
				operatorName = "plus";
			if(i==11)
				operatorName = "minus";
			if(i==12)
				operatorName = "by";
			if(i==13)
				operatorName = "into";
			fileName = "final/final_model_b_" + operatorName + ".txt";
		}
		f1.open(fileName, ios::out);
		for(int i =1; i<=N; i++)
		{
			for(int j=1; j<=M; j++)
				f1 << b_avg[i][j]<<"	";
			f1.scientific;
			f1 <<endl;
		}
		f1.close();
	}
}

//int main()
//{
//	//Uncomment to train again pre-recorded
//	//step1();	//preprocessing
//	//step2();	//make obs seq
//	//step3();	//training
//	//cout<<"1. Test on pre-recorded data."<<endl<<"2. Live testing."<<endl;
//	//cout<<"Enter choice"<<endl;
//	//int choice;
//	//cin>>choice;
//	//if(choice == 1)
//	//	step4();	//recorded testing
//	//else if( choice == 2)
//	step5();	//live recording
//	//live_train();
//	getch();
//	getch();
//	return 0;
//}

