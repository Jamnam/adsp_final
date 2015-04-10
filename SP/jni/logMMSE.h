/*
 * logMMSE.h
 *
 *  Created on: Apr 19, 2014
 *      Author: Balaji Prasanna
 */

#ifndef LOGMMSE_H_
#define LOGMMSE_H_
#define _USE_MATH_DEFINES
typedef struct LogMMSE
{
float* X;
float* Yprev;
float* noiseMu;
float* noiseVar;
float** SafetyNetP;
}LogMMSE;
LogMMSE* initialLogMMSE(int windowsize, int samplingFreq, int nfft);
void LogMmse(LogMMSE* logmmse,float* Y,float* Xest);
void destroyLogMMSE(LogMMSE** logmmse);

#endif /* LOGMMSE_H_ */
