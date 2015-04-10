/*=====================================
 * vad.c
 *=====================================
 *
 * Includes VAD function to Detect voice activity and consider the voice-noise cropping gaurd band
 *
 *  Created on: Apr 5, 2014
 *      Author: Balaji Prasanna
 */

#include "vad.h"
#include <stdlib.h>
#include <math.h>
//Wavlet Symlet '6' filters
static float lpf[] = {0.0154041093270274,0.00349071208421747,-0.117990111148191,-0.0483117425856330,0.491055941926747,0.787641141030194,0.337929421727622,-0.0726375227864625,-0.0210602925123006,0.0447249017706658,0.00176771186424280,-0.00780070832503415};
static float hpf[] = {0.00780070832503415,0.00176771186424280,-0.0447249017706658,-0.0210602925123006,0.0726375227864625,0.337929421727622,-0.787641141030194,0.491055941926747,0.0483117425856330,-0.117990111148191,	-0.00349071208421747,0.0154041093270274};
int firsttime;
void pop(float* x,float k,int Nf);   //Stack Definition
void bpsort(float* x,int Nf);		 //Stack based Insertion Sort Definition

Vad* initialVAD(int freq, int stepsize)
{
Vad* vad = (Vad*)malloc(sizeof(Vad));
vad->inbuf = (float*) calloc(stepsize+12,sizeof(float));
vad->Dsbuf = (float*) calloc((freq/stepsize),sizeof(float));
vad->Dsold = (float*) calloc((freq/stepsize),sizeof(float));		// Buffer to store
vad->noVoiceCount=0;
vad->flag=0;
vad->tqb=0.0;
firsttime=1;
return vad;
}

int VAD(Vad* vad, short* in, int freq, int stepsize)
{
	//declarations for vad
	int i,j,idx,qb,Nf,decision;
	float temp1,temp2,D=0.0,Px=0.0,Dw=0.0,Dc=0.0,Ds=0.0;
	float* lsb = (float*) calloc(stepsize,sizeof(float));
	float* hsb = (float*) calloc(stepsize,sizeof(float));
	Nf=freq/stepsize;

	//Input buffering for vad
	for(i=0; i<12; i++) {
		vad->inbuf[i] = vad->inbuf[stepsize+ i];
		}
	for(i=0; i<stepsize; i++) {
		vad->inbuf[12+i] = in[i]*3.051757812500000e-05;  //Short to float Q(15) conversion  1/32768 = 3.051757812500000e-05;
		Px+=(vad->inbuf[12+i]*vad->inbuf[12+i]);		 //Calculating Total power in input frame
	}

	//Low frequency side-band and High frequency sub-band decomposition from input frame
	for(i=0;i<stepsize;i++){
		temp1 = 0.0;
		temp2 = 0.0;
		for(j=0;j<12;j++){
			idx = 12 + (i - j);
			temp1 += vad->inbuf[idx]*lpf[j];
			temp2 += vad->inbuf[idx]*hpf[j];
		}
		lsb[i] = temp1 ;
		hsb[i] = temp2 ;
	}

														//Finding Spectral Power difference of this frame
	temp1 = 0.0;
	temp2 = 0.0;
	for(i=0;i<stepsize;){
		temp1+=hsb[i]*hsb[i]; 							//Down-sampling & Total Power calculations
		temp2+=lsb[i]*lsb[i];
		i+=2;
	}

	if(temp1>temp2) D = (temp1-temp2)/(stepsize/2);
	else  D = (temp2-temp1)/(stepsize/2);				//Absolute D value
	Dw = -2.0*D*(0.5 + (16.0/logf(2.0))*logf(1.0+2.0*Px));
	Dc = (1-exp(Dw))/(1+exp(Dw));
	Ds = Dc + vad->Dsold[Nf-1]*0.65;


	pop(vad->Dsbuf,vad->Dsold[0],Nf);					//Popping out old Ds from main stack
	vad->Dsbuf[Nf-1]=Ds;								//Pushing in new Ds to main stack
	bpsort(vad->Dsbuf,Nf);								//Bp-sorting the main stack
	pop(vad->Dsold,0.0,Nf);								//Popping Ds memory stack
	vad->Dsold[Nf-1]=Ds;								//Updating Ds memory stack

		qb=4;
		while(((vad->Dsbuf[qb]- vad->Dsbuf[qb-4])<0.001) && (qb<Nf-1)){
				qb++;
			}

		vad->tqb = 0.975*vad->tqb + 0.025*vad->Dsbuf[qb];

		if(firsttime == 1)
		{
			vad->tqb = vad->Dsbuf[qb];
			firsttime = 0;
		}

	if(Ds > vad->tqb)
	{   //Voice detection
		vad->noVoiceCount = 0;
		vad->flag = 0;
		decision = 1;
	}
	else{
		//Noise detection
		if(vad->flag)
			decision = 0;
		else{
			decision = 1;
			vad->noVoiceCount++;
			if(vad->noVoiceCount > (Nf/5)) //guard time = 200ms; 1sec = Nf frames => 200ms = Nf/5 frames
				vad->flag = 1;
		}
	}
	return decision;
}

inline void pop(float* x,float k,int Nf)
{
	int i=0,j;
	while (x[i] < k) i=i+1;						//Locating the candidate
	for(j=i;j<Nf-1;j++) x[j]=x[j+1];			//Popping it out
}

inline void bpsort(float* x,int Nf)
{
	    if (x[Nf-1] > x[Nf-2])					//Already sorted
	        return;
	    else
	    {
	        int i=0,j;
	        while (x[Nf-1] > x[i]) i=i+1;		//Locate incoming candidate
	        float temp = x[Nf-1];
	        for(j=Nf-1;j>i;j--) x[j]=x[j-1];    //Make space for incoming candidate
	        x[i]=temp;							//Inserting Candidate
	    }
}


void
destroyVAD(Vad** vad)
{
	if(*vad != NULL){
		if((*vad)->inbuf != NULL){
			free((*vad)->inbuf);
			(*vad)->inbuf = NULL;
		}
		if((*vad)->Dsbuf != NULL){
			free((*vad)->Dsbuf);
			(*vad)->Dsbuf = NULL;
		}
		if((*vad)->Dsold != NULL){
			free((*vad)->Dsold);
			(*vad)->Dsold = NULL;
		}
		free(*vad);
		*vad = NULL;
	}
}
