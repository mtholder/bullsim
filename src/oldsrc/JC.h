#ifndef _JC#define _JC#include "model.h"class JC : public Model	{	double beta;public:	JC();	JC(double pinv);	JC(int ncats,double gammaAlpha);	JC(double pinv,int ncats,double gammaAlpha);		void UpdatePMatrix(double **,double);	//void UpdatePmatGamma(double);	void InitializeParameters();	int GetEncodingType() {	return EncodingType(DNANoGap);}	int GetNStates()	{return 4;}};#endif