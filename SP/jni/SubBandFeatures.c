/*
 * SubBandFeatures.C
 *
 *  Created on: Apr 28, 2014
 *      Author: Balaji Prasanna
 */
#include "SubBandFeatures.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <android/log.h>

inline float entropy(float* in,int len);

int nbands = 8,nfft,sbfframecounter,bandPoints,OneSecondLongNumFrames;
float Normalize;

SubBandFeatures* initialSubBandFeatures(int nFFT,int fs,int stepsize)
{
	int i;
	nfft = nFFT;
	sbfframecounter = 1;
	bandPoints = nFFT/2/nbands;
	OneSecondLongNumFrames = fs/stepsize;
	Normalize = 1/((bandPoints+1)*OneSecondLongNumFrames);


	SubBandFeatures* SubBandFeature = (SubBandFeatures*)malloc(sizeof(SubBandFeatures));
	SubBandFeature->features= (float*)calloc(9,sizeof(float*));
	SubBandFeature->tempPeriodicity = (float*)calloc(nbands,sizeof(float*));
	SubBandFeature->tempEntropy= (float*)calloc(nbands,sizeof(float*));
	SubBandFeature->PrevSubBands = (float**)calloc(nbands,sizeof(float*));        // Stores previous sub-bands and its norm value in the last column of every band
	for(i=0;i<nbands;i++) SubBandFeature->PrevSubBands[i] = (float*)calloc((nfft/2/nbands)+2,sizeof(float*));

	return SubBandFeature;
}

void ComputeSubBandFeatures(SubBandFeatures* SubBandFeature, float* input)
{
	int i,j;
	float in[nfft],subBands[nbands][bandPoints+1+1],temp,max = 0;

	i=nbands-1;
	do
	{
	 	subBands[i][bandPoints+1] = 0;
		j=bandPoints-1;
		do
		{
			subBands[i][j] = sqrt(input[j+i*bandPoints]);								//Creating Sub-Band
			subBands[i][bandPoints+1] += input[j+i*bandPoints];							//For norm calculation
			j--;
		}
		while(j>=0);

		if(i != nbands-1)
		{
			subBands[i][bandPoints] = 0.0;
			subBands[i][bandPoints+1] = 1e-11 + sqrt(subBands[i][bandPoints+1]);
		}
		i--;
	}
	while(i>=0);

	subBands[nbands-1][bandPoints]= sqrt(input[nfft/2]);
	subBands[nbands-1][bandPoints+1] += input[nfft/2];
	subBands[nbands-1][bandPoints+1] = 1e-11 + sqrt(subBands[nbands-1][bandPoints+1]);

	i=nbands-1;
	do
	{
		j = bandPoints;
		do
		{
			temp = subBands[i][j]*SubBandFeature->PrevSubBands[i][j];
			if (abs(temp) > max) max = abs(temp);
			j--;
		}
		while(j>=0);
		if(subBands[i][bandPoints+1]*SubBandFeature->PrevSubBands[i][bandPoints+1] != 0)
		SubBandFeature->tempPeriodicity[i] += max/(subBands[i][bandPoints+1]*SubBandFeature->PrevSubBands[i][bandPoints+1]);		//Calculating Periodicity
		SubBandFeature->tempEntropy[i] += entropy(subBands[i],bandPoints+1);															// Caluclating Entropy;

		i--;
	}
	while(i>=0);
	i = nbands-1;
	do
	{
		j=bandPoints+1;
		do{ SubBandFeature->PrevSubBands[i][j] = subBands[i][j]; j--; }while(j>=0);
		i--;
	}
	while(i>=0);

	sbfframecounter++;
	if(sbfframecounter > OneSecondLongNumFrames)    											   //updating sub-band features only once in 1 second
	{
		i = nbands-1;
		do
		{
			if(i<5) SubBandFeature->features[i] = SubBandFeature->tempPeriodicity[i]*Normalize;    //Normalizing and taking only 1st 5 values
			if(i<4) SubBandFeature->features[i+5] = SubBandFeature->tempPeriodicity[i]*Normalize;  //Normalizing and taking only 1st 4 values
			SubBandFeature->tempPeriodicity[i] = 0.0;
			SubBandFeature->tempEntropy[i] = 0.0;

			i--;
		}
		while(i>=0);
		sbfframecounter = 1;
	}
}


inline float entropy(float* in,int len)
{
	int i;
	float H=0.0,interval;
	float max  = in[len-1];
	float min = in[len-1];
	int* p = (int*)calloc(len,sizeof(int));

	i=len-2;
	do
	{
		if(max > in[i]) max = in[i];
		if(min < in[i]) min = in[i];
		i--;
	}
	while(i>=0);
	interval = (float)(max-min)/(len-1);

	i=len-1;
	do
	{
	    p[(int)((in[i]-min)/interval)]++;              //Just counting their occurrence because we dont have PDF of data
	    i--;
	}
	while(i>=0);

	i=len-1;
	do
	{
		if (p[i] != 0)	H -= (float)p[i]/(float)len*logf((float)p[i]/(float)len)/0.6391;    //Estimates entropy
		i--;
	}
	while(i>=0);
	free(p);
	return H;
}

void destroySubBandFeatures(SubBandFeatures** SubBandFeature)
{
	if(* SubBandFeature!= NULL){
			if((*SubBandFeature)->Periodicity != NULL){
				free((*SubBandFeature)->Periodicity);
				(*SubBandFeature)->Periodicity = NULL;
			}
			if((*SubBandFeature)->Entropy != NULL){
				free((*SubBandFeature)->Entropy);
				(*SubBandFeature)->Entropy= NULL;
			}
			if((*SubBandFeature)->tempPeriodicity  != NULL){
				free((*SubBandFeature)->tempPeriodicity );
				(*SubBandFeature)->tempPeriodicity  = NULL;
			}
			if((*SubBandFeature)->tempEntropy != NULL){
				free((*SubBandFeature)->tempEntropy );
				(*SubBandFeature)->tempEntropy = NULL;
			}
			if((*SubBandFeature)->PrevSubBands != NULL){
				int i = nbands-1;
				do{
					free((*SubBandFeature)->PrevSubBands[i]);
					i--;
				}while(i>=0);
				free((*SubBandFeature)->PrevSubBands);
				(*SubBandFeature)->PrevSubBands= NULL;
			}
			free(*SubBandFeature);
			*SubBandFeature = NULL;
		}
}
