#include <jni.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "Timer.h"
#include "vad.h"
#include "periodogram.h"
#include "Traindata.h"
#include "logMMSE.h"
#include "SubBandFeatures.h"
#include "RandomForestClassifier.h"

typedef struct Variables {
	Timer* timer;
	Vad* vad;
	Periodogram* periodogram;
	SubBandFeatures* SubBandFeature;
	LogMMSE* logmmse;
	float* inputBuffer;
	float* outputBuffer;
	short* originalInput;
	float* coefficientBuffer;
	float* Xold;
	float* XveryOld;
	int frequency;
	int stepSize;
	int windowSize;
	int overlap;
	int nfft;
} Variables;
float divisor;
int mainframecounter = 1;
static void
compute(JNIEnv *env, jobject thiz,  jlong memoryPointer, jshortArray input)
{
	Variables* inParam = (Variables*) memoryPointer;
	startTimer(inParam->timer);
	
	short *_in = (*env)->GetShortArrayElements(env, input, NULL);
	
	int i, overlap, stepsize, Nf;
	overlap = inParam->overlap;
	stepsize = inParam->stepSize;
	float spd[inParam->nfft],features[9];
	Nf = inParam->frequency/stepsize;

	for(i=0; i<overlap; i++)
	{
		inParam->inputBuffer[i] = inParam->inputBuffer[stepsize + i];
	}

	for (i=0; i<stepsize; i++)
	{
		inParam->originalInput[i] = _in[i];
		inParam->inputBuffer[overlap + i] = _in[i]*3.051757812500000e-05;
	}

	FFT(inParam->periodogram,inParam->inputBuffer ,inParam->windowSize);

	spectralDensity(inParam->periodogram, spd);  //SPD = Spectral power density. *Not Special police delta

	if(!VAD(inParam->vad, _in,inParam->frequency,stepsize))
	{
		ComputeSubBandFeatures(inParam->SubBandFeature,spd);
		RandomForestClassifier(inParam->SubBandFeature->features);
	}

	(*env)->ReleaseShortArrayElements(env, input, _in, 0);

	LogMmse(inParam->logmmse,spd,inParam->coefficientBuffer);

	iFFT(inParam->periodogram, inParam->coefficientBuffer);

	//************Overlap Add and Synthesis*************//
	i = stepsize-1;
	do
	{
		inParam->outputBuffer[i] = inParam->Xold[i] + inParam->periodogram->real[i]*divisor*inParam->periodogram->windowing[i];
		inParam->Xold[i] = inParam->periodogram->real[i+stepsize]*divisor*inParam->periodogram->windowing[i+stepsize];
		i--;
	}while(i>=0);

	if(overlap > stepsize && 2*stepsize > overlap)
	{
		i = overlap - stepsize - 1;
		do
		{
			inParam->outputBuffer[i+2*stepsize-overlap] += inParam->XveryOld[i];
			inParam->XveryOld[i] = inParam->XveryOld[i+overlap-stepsize];
			inParam->XveryOld[i+overlap-stepsize] = inParam->periodogram->real[i+2*stepsize]*divisor*inParam->periodogram->windowing[i+2*stepsize];
			i--;
		}while(i>=0);
	}

	mainframecounter++;

	//**********The END***********//

	stopTimer(inParam->timer);
}

static jlong
initialize(JNIEnv* env, jobject thiz, jint frequency, jint stepsize, jint windowsize, jint decisionBufferLength)
{
	Variables* inParam = (Variables*) malloc(sizeof(Variables));
	inParam->timer = newTimer();
	__android_log_print(ANDROID_LOG_ERROR, "Speech Processing", "Initialization: Sampling Frequency %d, Step Size %d, Window Size %d, Decision Buffer %d", frequency, stepsize, windowsize, decisionBufferLength);
	int i;
	inParam->frequency = frequency;
	inParam->stepSize = stepsize;
	inParam->windowSize = windowsize;
	inParam->overlap = windowsize-stepsize;
	inParam->nfft = 1;
	for(i=0;i<ceil(log(windowsize)/log(2));i++) inParam->nfft *= 2;
	divisor = (float)1/inParam->nfft;
	inParam->vad = initialVAD(frequency,stepsize);
	inParam->periodogram = newPeriodogram(inParam->nfft,windowsize);
	inParam->SubBandFeature = initialSubBandFeatures(inParam->nfft,frequency,stepsize);
	inParam->logmmse = initialLogMMSE(windowsize,frequency,inParam->nfft);
	inParam->inputBuffer = (float*)calloc(windowsize,sizeof(float));
	inParam->outputBuffer = (float*)calloc(stepsize,sizeof(float));
	inParam->originalInput = (short*)malloc(stepsize*sizeof(short));
	inParam->coefficientBuffer = (float*)calloc(inParam->nfft,sizeof(float));
	inParam->Xold = (float*)calloc(stepsize,sizeof(float));
	if(inParam->overlap > stepsize) inParam->XveryOld = (float*)calloc(2*(inParam->overlap - stepsize),sizeof(float));
	i = inParam->nfft - 1;
	do{inParam->coefficientBuffer[i] = (float)1.0; i--;}while(i>=0);
	mainframecounter = 1;
	return (jlong)inParam;
}

static void
finish(JNIEnv* env, jobject thiz, jlong memoryPointer)
{
	Variables* inParam = (Variables*) memoryPointer;
	//cleanup memory
	if(inParam != NULL){
		tellTimerTime(inParam->timer);
		destroyTimer(&(inParam->timer));
		destroyVAD(&(inParam->vad));
		destroyPeriodogram(&(inParam->periodogram));
		destroyLogMMSE(&(inParam->logmmse));
		destroySubBandFeatures(&(inParam->SubBandFeature));
		if(inParam->inputBuffer != NULL){
			free(inParam->inputBuffer);
			inParam->inputBuffer = NULL;
		}
		if(inParam->outputBuffer != NULL){
			free(inParam->outputBuffer);
			inParam->outputBuffer = NULL;
		}
		if(inParam->originalInput != NULL){
			free(inParam->originalInput);
			inParam->originalInput = NULL;
		}
		if(inParam->coefficientBuffer != NULL){
			free(inParam->coefficientBuffer);
			inParam->coefficientBuffer = NULL;
		}
		if(inParam->Xold!= NULL){
			free(inParam->Xold);
			inParam->Xold = NULL;
		}
		if(inParam->XveryOld!= NULL){
			free(inParam->XveryOld);
			inParam->XveryOld = NULL;
		}
		free(inParam);
		inParam = NULL;
	}
}

static jfloat
getTime(JNIEnv* env, jobject thiz, jlong memoryPointer)
{
	Variables* inParam = (Variables*) memoryPointer;
	return getTimerMS(inParam->timer);
}

static jfloatArray
getOutput(JNIEnv* env, jobject thiz, jlong memoryPointer, jint outputSelect)
{
	Variables* inParam = (Variables*) memoryPointer;

	jshortArray output = (*env)->NewShortArray(env, inParam->stepSize);
	short *_output = (*env)->GetShortArrayElements(env, output, NULL);

	if(outputSelect == 0) { //Case 1 - Original input signal
		int i;
		for(i=0;i<inParam->stepSize;i++)
		{
			_output[i] = inParam->originalInput[i];
		}

	} else {				//Case 2 - Processed output signal
		int i;
		for(i=0;i<inParam->stepSize;i++)
		{
			_output[i] = (short)(inParam->outputBuffer[i]*32768.0f);
		}
	}

	(*env)->ReleaseShortArrayElements(env, output, _output, 0);
	return output;
}

static jfloatArray
getDebug(JNIEnv* env, jobject thiz, jlong memoryPointer, jint debugSelect)
{
	Variables* inParam = (Variables*) memoryPointer;

	jfloatArray debugOutput = NULL;

	if(debugSelect == 0) {

		//Test Case 1 - inputBuffer contents

		debugOutput = (*env)->NewFloatArray(env, inParam->windowSize);
		float *_debugOutput = (*env)->GetFloatArrayElements(env, debugOutput, NULL);

		int i;
		for (i=0; i<inParam->windowSize;i++)
		{
			_debugOutput[i] = inParam->inputBuffer[i];
		}

		(*env)->ReleaseFloatArrayElements(env, debugOutput, _debugOutput, 0);

	} else if (debugSelect == 1) {

		//Test Case 2 - outputBuffer contents

		debugOutput = (*env)->NewFloatArray(env, inParam->stepSize);
		float *_debugOutput = (*env)->GetFloatArrayElements(env, debugOutput, NULL);

		int i;
		for (i=0; i<inParam->stepSize;i++)
		{
			_debugOutput[i] = inParam->outputBuffer[i];
		}

		(*env)->ReleaseFloatArrayElements(env, debugOutput, _debugOutput, 0);

	}

	return debugOutput;
}

////////////////////////////////////////////////////////////////////////////////////////////
// JNI Setup - Functions and OnLoad
////////////////////////////////////////////////////////////////////////////////////////////

static JNINativeMethod nativeMethods[] =
	{//		Name							Signature											Pointer
			{"compute", 					"(J[S)V",											(void *)&compute				},
			{"initialize",					"(IIII)J",											(void *)&initialize				},
			{"finish",						"(J)V",												(void *)&finish					},
			{"getTime",						"(J)F",												(void *)&getTime				},
			{"getOutput",					"(JI)[S",											(void *)&getOutput				},
			{"getDebug",					"(JI)[F",											(void *)&getDebug				}
	};

jint
JNI_OnLoad(JavaVM* vm, void* reserved)
{
	JNIEnv* env;
	jint result;
	//get a hook to the environment
	result = (*vm)->GetEnv(vm, (void**) &env, JNI_VERSION_1_6);
	if (result == JNI_OK) {
		//find the java class to hook the native methods to
		jclass filters = (*env)->FindClass(env, "com/dsp/speechpipeline/SpeechProcessing");
		if (filters != NULL) {
			result = (*env)->RegisterNatives(env, filters, nativeMethods, sizeof(nativeMethods)/sizeof(nativeMethods[0]));
			(*env)->DeleteLocalRef(env, filters);
			if(result == JNI_OK){
				return JNI_VERSION_1_6;
			} else {
				//something went wrong with the method registration
				return JNI_ERR;
			}
		} else {
			//class wasn't found
			return JNI_ERR;
		}
	} else {
		//could not get environment
		return JNI_ERR;
	}
}
