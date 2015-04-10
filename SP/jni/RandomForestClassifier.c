/*
 * RandomForestClassifier.c
 *
 *  Created on: May 1, 2014
 *      Author: Balaji Prasanna
 */

#include "RandomForestClassifier.h"
#include "Traindata.h"

int ntrees = 5,nclass=3;

void RandomForestClassifier(float* subBandFeature)
{
	   int tree_output[ntrees];
	   int i,CurrentNode,CutVar,max=0;
	   int vote1=0,vote2=0,vote3=0;                       // vote counters for 3 noise classes

	   /* tree# [4][numnodes] = {
	    * 						{cutvar}
	    * 						{cutvalue}
	    * 						{childnode}
	    * 						{nodelabel}
	    * 						};
	    */
//========================================================================================================================
	   //Traversing Tree1
	   CurrentNode = 0;
	   while (tree1[2][CurrentNode]!=0)
	       {
	           CutVar = tree1[0][CurrentNode];
	           if (subBandFeature[CutVar] < tree1[1][CurrentNode])
	           CurrentNode = tree1[2][CurrentNode]-1;
	           else CurrentNode = tree1[2][CurrentNode];
	               }
	               tree_output[0] = tree1[3][CurrentNode];
//========================================================================================================================
	               //Traversing Tree2
	               CurrentNode = 0;
	               while (tree5[2][CurrentNode]!=0)
	                   {
	                       CutVar = tree5[0][CurrentNode];
	                       if (subBandFeature[CutVar] < tree5[1][CurrentNode])
	                       CurrentNode = tree5[2][CurrentNode]-1;
	                       else CurrentNode = tree5[2][CurrentNode];
	                           }
	                           tree_output[1] = tree5[3][CurrentNode];
//========================================================================================================================
	                           //Traversing Tree3
	                           CurrentNode = 0;
	                           while (tree3[2][CurrentNode]!=0)
	                               {
	                                   CutVar = tree3[0][CurrentNode];
	                                   if (subBandFeature[CutVar] < tree3[1][CurrentNode])
	                                   CurrentNode = tree3[2][CurrentNode]-1;
	                                   else CurrentNode = tree3[2][CurrentNode];
	                                       }
	                                       tree_output[2] = tree3[3][CurrentNode];
//========================================================================================================================
	                                       //Traversing Tree4
	                                       CurrentNode = 0;
	                                       while (tree4[2][CurrentNode]!=0)
	                                           {
	                                               CutVar = tree4[0][CurrentNode];
	                                               if (subBandFeature[CutVar] < tree4[1][CurrentNode])
	                                               CurrentNode = tree4[2][CurrentNode]-1;
	                                               else CurrentNode = tree4[2][CurrentNode];
	                                                   }
	                                                   tree_output[3] = tree4[3][CurrentNode];
//========================================================================================================================
	                                                   //Traversing Tree5
	                                                   CurrentNode = 0;
	                                                   while (tree5[2][CurrentNode]!=0)
	                                                       {
	                                                           CutVar = tree5[0][CurrentNode];
	                                                           if (subBandFeature[CutVar] < tree5[1][CurrentNode])
	                                                           CurrentNode = tree5[2][CurrentNode]-1;
	                                                           else CurrentNode = tree5[2][CurrentNode];
	                                                               }
	                                                               tree_output[4] = tree5[3][CurrentNode];
//========================================================================================================================

for (i=0;i<ntrees;i++)
{
	if(tree_output[i] == 1) vote1++;
	if(tree_output[i] == 2) vote2++;
	if(tree_output[i] == 3) vote3++;
}

max = vote1; classdecision = 0;
if (max < vote2) {max = vote2; classdecision = 1;}
if (max < vote3) {max = vote3; classdecision = 2;}

}
