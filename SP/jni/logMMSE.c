/*
 * logMMSE.c
 *
 *  Created on: Apr 19, 2014
 *      Author: Balaji Prasanna
 */

#include<logMMSE.h>
#include<math.h>
#include<stdlib.h>
#include"Traindata.h"

#define MAXIT 100 					//Maximum allowed number of iterations.
#define EULER 0.5772156649 			//Euler's constant
#define FPMIN 1.0e-30 				//Close to smallest representable floating-point number.
#define EPS 1.0e-7 					//Desired relative error, not smaller than the machine precision.

int framecounter, nFFT,numWinPmin,SafetyNetInd;
float alphas;

inline float max(float a,float b);
inline float min(float a,float b);
inline void expint(float *in,float *out);

LogMMSE* initialLogMMSE(int windowsize, int fs, int nfft)
{
	LogMMSE* logmmse = (LogMMSE*)malloc(sizeof(LogMMSE));

	framecounter = 1;
	nFFT = nfft;
	numWinPmin = (int) (0.8*fs)/windowsize;
	SafetyNetInd = 0;

	logmmse->X = (float*)calloc(nfft,sizeof(float));
	logmmse->Yprev = (float*)calloc(nfft,sizeof(float));
	logmmse->noiseMu= (float*)calloc(nfft,sizeof(float));
	logmmse->noiseVar = (float*)calloc(nfft,sizeof(float));
	logmmse->SafetyNetP = (float**)calloc(nfft,sizeof(float*));

	int i=nfft-1;
	do{
		logmmse->SafetyNetP[i] = (float*)calloc(numWinPmin,sizeof(float*));
		i--;
	}while(i>=0);

	return logmmse;
}
void LogMmse(LogMMSE* logmmse,float* Y,float* Xest)
{

		int i=nFFT-1,j;
		int psiNTdB,gammadB;
		float temp,vk,eivk;
		float gamma[nFFT],psi[nFFT],psiNT[nFFT],SafetyNetPmin[nFFT],NoiseEst[nFFT],smoothpost[nFFT];
		double Speechpresenceprob[nFFT];

		//*******Initial Noise variance Estimation assuming 1st 6 frames are noise********//
		if(framecounter <= 6)
		{
			do
			{
				logmmse->noiseMu[i]+= sqrt(Y[i]);          // Noise mean
				i--;
			}while(i>=0);

			i=nFFT-1;
			do
			{
				logmmse->noiseVar[i] = (float)logmmse->noiseMu[i]*logmmse->noiseMu[i]*2.777777777777778e-02; 	//(abs(fft)/6)^2
				i--;
			}while(i>=0);
		}

		i=nFFT-1;
		do
		{
			gamma[i] = min((Y[i]/logmmse->noiseVar[i]),40);

			if(framecounter == 1)
			{
				psi[i] = 0.98 + 0.02*max((gamma[i]-1),0);
				psiNT[i] = max((0.98 + 0.02*gamma[i]),0.003162277660168);
				logmmse->SafetyNetP[i][SafetyNetInd] = 0.9*Y[i];
			}
			else
			{
				psi[i] = 0.98*(double)logmmse->X[i]/logmmse->noiseVar[i] + 0.02*max((gamma[i]-1),0);
				psi[i] = max(psi[i],0.003162277660168);  														// limiting psi to -25dB
				psiNT[i] = max((0.98*logmmse->Yprev[i]/logmmse->noiseVar[i]+0.02*gamma[i]),0.003162277660168);
				if(SafetyNetInd > 0) logmmse->SafetyNetP[i][SafetyNetInd] = 0.1*logmmse->SafetyNetP[i][SafetyNetInd-1]+0.9*Y[i];
				else logmmse->SafetyNetP[i][SafetyNetInd] = 0.1*logmmse->SafetyNetP[i][numWinPmin]+0.9*Y[i];
			}

			//*******Setting SafetynetPmin********//
			temp = logmmse->SafetyNetP[i][0];
			j=numWinPmin-1;
			do{
				if(temp>logmmse->SafetyNetP[i][j]) temp = logmmse->SafetyNetP[i][j];
				j--;
			}while(j>=0);
			SafetyNetPmin[i] = temp;

			//*******Fetching Noise Estimate********//
			logmmse->Yprev[i] = Y[i];
			psiNTdB = round(max(min(10.0*log10f(psiNT[i]),40),-19))+19;
			gammadB = round(max(min(10.0*log10f(gamma[i]),40),-30))+30;
			NoiseEst[i] = GainTable[classdecision][71*psiNTdB+gammadB]*Y[i];

			i--;
		}while(i >= 0);

		//*******Frequency Smoothing********//
		smoothpost[0] = (gamma[0]+gamma[1])/2;
		i = nFFT-2;
		do
		{
			smoothpost[i] = (gamma[i-1]+gamma[i]+gamma[i+1])/3;
			i--;
		}while(i>0);
		smoothpost[nFFT-1] = (gamma[nFFT-1]+gamma[nFFT-2])/2;

		i=nFFT-1;
		do
		{
			//*******Updating NoiseVariance for next frame********//
			if(framecounter == 1) Speechpresenceprob[i] = (double)0.9;
			else Speechpresenceprob[i] =  (double)0.1*Speechpresenceprob[i]+(double)0.9*(smoothpost[i]>4);

			alphas = 0.85 + 0.15 * Speechpresenceprob[i];
			logmmse->noiseVar[i] = alphas*logmmse->noiseVar[i] + (1-alphas)*NoiseEst[i];

			if((1.5*SafetyNetPmin[i]) > logmmse->noiseVar[i])
			{
				Speechpresenceprob[i] = 0;
				logmmse->noiseVar[i] = 1.5*SafetyNetPmin[i];
			}

			//*******log MMSE Estimator********//
			temp = psi[i]/(1+psi[i]);
			vk = gamma[i]*temp;
			expint(&vk,&eivk);                                  //Evaluating the exponential integral part
			Xest[i] = temp*exp(eivk);
			logmmse->X[i] = Y[i]*Xest[i]*Xest[i];

			i--;
		}while(i>=0);

		SafetyNetInd = SafetyNetInd+1;
		if(SafetyNetInd>=numWinPmin) SafetyNetInd = 0;
		framecounter++;
}

inline void expint(float *in,float *out)
//***********Evaluates the exponential integral En(x)************//
{
	float x=*in;
	int i,ii,nm1,n=1;
	float a,b,c,d,del,fact,h,psi,ans;

	nm1=n-1;
	if (n < 0 || x < 0.0 || (x==0.0 && (n==0 || n==1)))
		return;
	else {
		if (n == 0) ans=(exp(-x)/x);// Special case.
		else {
			if (x == 0.0) ans=(float)1.0/nm1;// Another special case.
			else {
				if (x > 1.0) {// Lentz's algorithm (x5.2).
					b=x+n;
					c=(float)1.0/FPMIN;
					d=(float)1.0/b;
					h=d;
					for (i=1;i<=MAXIT;i++) {
						a = -i*(nm1+i);
						b += 2.0;
						d=(float)1.0/(a*d+b); //Denominators cannot be zero.
						c=b+a/c;
						del=c*d;
						h *= del;
						if (fabs(del-1.0) < EPS) {
						*out=h*exp(-x)*0.5;
						return ;
						}
					}
					return;
				} else { //Evaluate series.
					ans = (nm1!=0 ? 1.0/nm1 : -log(x)-EULER);// Set first term.
					fact=1.0;
					for (i=1;i<=MAXIT;i++) {
						fact *= -x/i;
						if (i != nm1) del = -fact/(i-nm1);
						else {
							psi = -EULER; //Compute (n).
							for (ii=1;ii<=nm1;ii++) psi += 1.0/ii;
							del=fact*(-log(x)+psi);
						}
						ans += del;
						if (fabs(del) < fabs(ans)*EPS) {
							*out=ans*0.5;
							return;
						}
					}
					return;
				}
			}
		}
	}
return ;
}

inline float max(float a,float b)
{
	return a>b?a:b;
}
inline float min(float a,float b)
{
	return a<b?a:b;
}


void destroyLogMMSE(LogMMSE** logmmse)
{
	if(* logmmse!= NULL){
			if((*logmmse)->X != NULL){
				free((*logmmse)->X);
				(*logmmse)->X = NULL;
			}
			if((*logmmse)->Yprev != NULL){
				free((*logmmse)->Yprev);
				(*logmmse)->Yprev= NULL;
			}
			if((*logmmse)->noiseVar != NULL){
				free((*logmmse)->noiseVar);
				(*logmmse)->noiseVar = NULL;
			}
			if((*logmmse)->SafetyNetP != NULL){
				int i = nFFT-1;
				do{
					free((*logmmse)->SafetyNetP[i]);
					i--;
				}while(i>=0);
				free((*logmmse)->SafetyNetP);
				(*logmmse)->SafetyNetP= NULL;
			}
			free(*logmmse);
			*logmmse = NULL;
		}
}
