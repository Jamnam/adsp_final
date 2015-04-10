/*
 * vad.h
 *
 *  Created on: Apr 5, 2014
 *      Author: Balaji Prasanna
 */
#ifndef VAD_H_
#define VAD_H_
typedef struct Vad{
		float* inbuf;
		float* Dsbuf;
		float* Dsold;
		int noVoiceCount,flag;
		float Ds,tqb;
}Vad;
Vad* initialVAD(int freq, int stepsize);
int VAD(Vad* vad, short* in, int freq, int stepsize);
void destroyVAD(Vad** vad);
#endif
