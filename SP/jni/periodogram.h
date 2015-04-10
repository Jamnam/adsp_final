/*
 * periodogram.h
 *
 *  Created on: Apr 5, 2014
 *      Author: Balaji Prasanna
 */

#ifndef PERIODOGRAM_H
#define PERIODOGRAM_H
#define _USE_MATH_DEFINES

#include <stdlib.h>
#include <math.h>

typedef struct Periodogram {
		int nFFT;
		float* windowing;
		float* sine;
		float* cosine;
		float* real;
		float* imaginary;
} Periodogram;

Periodogram* newPeriodogram(int points,int windowsize);
void FFT(Periodogram* periodogram,float* input, int windowsize);
void spectralDensity(Periodogram* periodogram,float* output);
void iFFT(Periodogram* fft, float* coefficients);
void destroyTransform(Periodogram** periodogram);
#endif
