/*
 * SubBandFeatures.h
 *
 *  Created on: Apr 28, 2014
 *      Author: Balaji Prasanna
 */

#ifndef SUBBANDFEATURES_H_
#define SUBBANDFEATURES_H_
#define _USE_MATH_DEFINES

typedef struct SubBandFeatures{

		float** PrevSubBands;
		float* Entropy;
		float* features;
		float* Periodicity;
		float* tempEntropy;
		float* tempPeriodicity;

}SubBandFeatures;

SubBandFeatures* initialSubBandFeatures(int nFFT,int fs,int stepsize);
void ComputeSubBandFeatures(SubBandFeatures* SubBandFeature, float* in);
void destroySubBandFeatures(SubBandFeatures** SubBandFeature);

#endif /* SUBBANDFEATURES_H_ */
