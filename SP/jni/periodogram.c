/* * periodogram.c
 *
 *  Created on: Apr 5, 2014
 *      Author: Balaji Prasanna
 */

#include "periodogram.h"
void hanning(float* hann, int windowsize);

Periodogram*
newPeriodogram(int nFFT, int windowsize)
{
	Periodogram* newperiodogram = (Periodogram*)malloc(sizeof(Periodogram));

	newperiodogram->nFFT = nFFT;
	newperiodogram->real = (float*)calloc(nFFT,sizeof(float));
	newperiodogram->imaginary = (float*)calloc(nFFT,sizeof(float));
	newperiodogram->sine = (float*)malloc((nFFT/2)*sizeof(float));
	newperiodogram->cosine = (float*)malloc((nFFT/2)*sizeof(float));
	newperiodogram->windowing = (float*)malloc(windowsize*sizeof(float));
	hanning(newperiodogram->windowing,windowsize);												// Hanning window precomputing

	//precompute twiddle factors for fft
		float arg;
		int i;
		for (i=0;i<nFFT/2;i++)
		{
			arg = -2*3.14159265358979*i/nFFT;
			newperiodogram->cosine[i] = cos(arg);
			newperiodogram->sine[i] = sin(arg);
		}

	return newperiodogram;
}

void hanning(float* hann, int windowsize)  // Calculates Hanning window for windowsize
{
	int i;
	for(i=0;i<windowsize;i++)
		hann[i] = 0.5*(1-cos(2*3.14159265358979*(float)(i+1)/((float)windowsize+1)));
}

void FFT(Periodogram* fft,float* input, int windowsize)
{
	int i,j,k,L,m,n,o,p,q,r;
	float tempReal,tempImaginary,cos,sin,xt,yt;
	k = fft->nFFT;
	i=k-1;
	do
	{
		if(i<windowsize) fft->real[i] = input[i]*fft->windowing[i];       // Taking input samples from 0 to Windowsize
		else fft->real[i] = 0;											  // Zero Padding remaining samples till nFFT points
		fft->imaginary[i] = 0;
		i--;
	}while(i>=0);

	j=0;
	m=k/2;

	//bit reversal
	for(i=1;i<(k-1);i++)
	{
		L=m;

		while(j>=L)
		{
			j=j-L;
			L=L/2;
		}

		j=j+L;

		if(i<j)
		{
			tempReal=fft->real[i];
			tempImaginary=fft->imaginary[i];
			fft->real[i]=fft->real[j];
			fft->imaginary[i]=fft->imaginary[j];
			fft->real[j]=tempReal;
			fft->imaginary[j]=tempImaginary;
		}
	}

	L=0;
	m=1;
	n=k/2;

	//computation
	for(i=k;i>1;i=(i>>1))
	{
		L=m;
		m=2*m;
		o=0;

		for(j=0;j<L;j++)
		{
			cos=fft->cosine[o];
			sin=fft->sine[o];
			o=o+n;

			for(p=j;p<k;p=p+m)
			{
				q=p+L;

				xt=cos*fft->real[q]-sin*fft->imaginary[q];
				yt=sin*fft->real[q]+cos*fft->imaginary[q];
				fft->real[q]=(fft->real[p]-xt);
				fft->imaginary[q]=(fft->imaginary[p]-yt);
				fft->real[p]=(fft->real[p]+xt);
				fft->imaginary[p]=(fft->imaginary[p]+yt);
			}
		}
		n=n>>1;
	}
}

void iFFT(Periodogram* fft, float* coefficients)								//IFFT based on (FFT(X*))*
{
	int i,j,k,L,m,n,o,p,q,r;
	float tempReal,tempImaginary,cos,sin,xt,yt;
	k = fft->nFFT;

	i=k-1;
	do{
		fft->real[i] = coefficients[i] * fft->real[i];							//Multiplying Coefficients from logMMSE
		fft->imaginary[i] = -1*coefficients[i] * fft->imaginary[i];				//Imaginary part negated for conjugating
		i--;
	}while(i>=0);

	j=0;
	m=k/2;

	//bit reversal
	for(i=1;i<(k-1);i++)
	{
		L=m;

		while(j>=L)
		{
			j=j-L;
			L=L/2;
		}

		j=j+L;

		if(i<j)
		{
			tempReal=fft->real[i];
			tempImaginary=fft->imaginary[i];
			fft->real[i]=fft->real[j];
			fft->imaginary[i]=fft->imaginary[j];
			fft->real[j]=tempReal;
			fft->imaginary[j]=tempImaginary;
		}
	}

	L=0;
	m=1;
	n=k/2;

	//computation
	for(i=k;i>1;i=(i>>1))
	{
		L=m;
		m=2*m;
		o=0;

		for(j=0;j<L;j++)
		{
			cos=fft->cosine[o];
			sin=fft->sine[o];
			o=o+n;

			for(p=j;p<k;p=p+m)
			{
				q=p+L;

				xt=cos*fft->real[q]-sin*fft->imaginary[q];
				yt=sin*fft->real[q]+cos*fft->imaginary[q];
				fft->real[q]=(fft->real[p]-xt);
				fft->imaginary[q]=(fft->imaginary[p]-yt);
				fft->real[p]=(fft->real[p]+xt);
				fft->imaginary[p]=(fft->imaginary[p]+yt);
			}
		}
		n=n>>1;
	}
}

void spectralDensity(Periodogram* periodogram, float* output)
{
	int n = (periodogram->nFFT) -1;
	do
	{
		output[n] = periodogram->real[n]*periodogram->real[n]+periodogram->imaginary[n]*periodogram->imaginary[n];
		n--;
	}while(n >= 0);
}

void
destroyPeriodogram(Periodogram** periodogram)
{
	if(*periodogram != NULL){
		if((*periodogram)->windowing!= NULL){
			free((*periodogram)->windowing);
			(*periodogram)->windowing = NULL;
		}
		if((*periodogram)->cosine != NULL){
			free((*periodogram)->cosine);
			(*periodogram)->cosine = NULL;
		}
		if((*periodogram)->sine != NULL){
			free((*periodogram)->sine);
			(*periodogram)->sine = NULL;
		}
		if((*periodogram)->real != NULL){
			free((*periodogram)->real);
			(*periodogram)->real = NULL;
		}
		if((*periodogram)->imaginary != NULL){
			free((*periodogram)->imaginary);
			(*periodogram)->imaginary = NULL;
		}
		free(*periodogram);
		*periodogram = NULL;
	}
}
