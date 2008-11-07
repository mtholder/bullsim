#ifndef _K2P#define _K2P#include "model.h"class K2P : public Model	{	double beta;	void CalculateBeta();		static double defKappa;public:enum params {kappa=0};	K2P();	K2P(double k);	K2P(double k, double pinv);	K2P(double k, int ncats,double gammaAlpha);	K2P(double k, double pinv,int ncats,double gammaAlpha);	~K2P();	void UpdatePMatrix(double **,double);	//void UpdatePmatGamma(double);	void InitializeParameters();	static void set_default(double k)	{		defKappa=k;		}	int GetEncodingType() {	return EncodingType(DNANoGap);}	int GetNStates()	{return 4;}	void ParameterHasChanged(Parameter *p)	{if(p==gammashape)	CalculateRates();											else	if(p==param[0])	CalculateBeta();											}};#endif