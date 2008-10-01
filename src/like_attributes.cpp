// This file is part of BULL, a program for phylogenetic simulations
// most of the code was written by Mark T. Holder.

//	This program is for internal use by the lab of Dr. Tandy Warnow only.
//	Do not redistribute the code.  It is distributed in the hope that it will be useful,
//	but WITHOUT ANY WARRANTY; without even the implied warranty of
//	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
//
//	Some of the code is from publically available source by Paul Lewis, Ziheng Yang, 
//	John Huelsenbeck, David Swofford , and others  (as noted in the code).
//	In fact the main structure of the program was created by modifying Paul Lewis' 
//	basiccmdline.cpp from his NCL
//
//	This code was used in Mark's dissertation, some changes were made in order to 
//	get it to compile on gcc.  It is possible that this porting introduced bugs (very
//	little debugging has been done on UNIX platforms).	I would suggest checking 
//	the simulator by generating data on trees with short branches, etc.
 

// This file is part of BULL, a program for phylogenetic simulations
// most of the code was written by Mark T. Holder.

//	This program is for internal use by the lab of Dr. Tandy Warnow only.
//	Do not redistribute the code.  It is distributed in the hope that it will be useful,
//	but WITHOUT ANY WARRANTY; without even the implied warranty of
//	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
//
//	Some of the code is from publically available source by Paul Lewis, Ziheng Yang, 
//	John Huelsenbeck, David Swofford , and others  (as noted in the code).
//	In fact the main structure of the program was created by modifying Paul Lewis' 
//	basiccmdline.cpp from his NCL
//
//	This code was used in Mark's dissertation, some changes were made in order to 
//	get it to compile on gcc.  It is possible that this porting introduced bugs (very
//	little debugging has been done on UNIX platforms).	I would suggest checking 
//	the simulator by generating data on trees with short branches, etc.
 

#include "like_attributes.hpp"
#include "bull.hpp"


using namespace bull;


int LikeAttr::currNChar(1);
int LikeAttr::currNStates(4);
int LikeAttr::currShortPerChar(1);
int LikeAttr::currNRateCats(1);
int LikeAttr::currNStatesInLastShort(1);
bool LikeAttr::currModelDirty(true);
double LikeAttr::Multiplier(1.53);

#if defined (CHAR_BY_CHAR) && CHAR_BY_CHAR
	int LikeAttr::currCharIndex(0);
#endif

LikeAttr::LikeAttr()//Constructor for internal node's attributes
{	model=NULL;
	Pmat=NULL;
	blen=NULL;
	blenMod=NULL;
}

void LikeAttr::SetBLen(double d)
{	if (blen)
		blen->SetCurrent(d);
	else {
		const int r = par(CUR)|par(MIN);
		blen= new BoundedParameter(d, 0.0, 0.0, r, bull::DEFAULT_START_EDGE_LEN);
	}
}	


LikeAttr *LikeAttr::Copy()
{	LikeAttr *tempLA;
	if (model)
		tempLA=new LikeAttr(model);
	else
		tempLA=new LikeAttr();
	if (blen)
		tempLA->blen=new BoundedParameter(blen->val,0.0,0.0,blen->GetSetting(),
											bull::DEFAULT_START_EDGE_LEN);
	//TEMPORARY COPY Leaves blen mod uncopied  MEMORY LEAK old blen still allocated I don't know why anymore!!!
	return tempLA;
}



TreeSimAttr::TreeSimAttr(int nc,Model *m)
{	model=m;
	rates=new double[nc];
	nChar=nc;
	hasRateHet=m->HasRateHet();
	freq=m->GetStateFreqs();
	nStates=m->GetNStates();
	ownsBrLens=false;

}


Model* LikeAttr::GetModel() {
return model;
}
LikeAttr::~LikeAttr()
{	//assert(SharePmatMemory());
	delete blen;
	delete blenMod;
}
void TreeSimAttr::GenerateRandomSequence(int nc,short *dest)
{	assert(nc <= nChar);
	for (int i=0; i < nc; i++)
		*dest++=GetRandomIndexFromFreqs(nStates,freq);
}

void TreeSimAttr::GenerateNewSetOfRates(int nc)
{	double pnv=model->GetPInv();
	if (pnv>0.0)
		{if(model->GetNRateCats()>1)
			{double sp=model->GetShapeParam();
			double *ra=rates;
			for (int i=0; i < nc; i++)
				*ra++=(RandomNumber()<pnv ? 0.0 : (rndgamma(sp))/sp);
			}
		else
			{double *ra=rates;
			for (int i=0; i < nc; i++)
				*ra++=(RandomNumber()<pnv ? 0.0 : 1.0);
			}
		}
	else
		{assert(model->GetNRateCats()>1);
		double sp=model->GetShapeParam();
		double *ra=rates;
		for (int i=0; i < nc; i++)
			*ra++=(rndgamma(sp))/sp;
		}
}

LikeAttr::LikeAttr(Model *m)//Constructor for internal node's attributes
{	model=m;
	Pmat=m->GetPmat();
	blen=NULL;
	blenMod=NULL;
}


SimNodeLikeAttr::SimNodeLikeAttr (int nc,Model *m)
	: LikeAttr(m)
{	CharsInShorts=new short[nc];
}

void SimNodeLikeAttr::SimulateSeq(int nc,short *ancSeq,double *rates)
{	//rate heterogeneity
	short* cis=CharsInShorts;
	double transformedBranch=(blen->val)/(1.0-model->GetPInv());
	if (model->GetNRateCats()>1)
		{for (int i=0; i < nc; i++)//rates could be anything
			{if(*rates>0.0)//don't iterate rates here cause you'll need it 
				{model->UpdatePRow(*Pmat,(*rates)*transformedBranch,*ancSeq);
				*cis++=GetRandomIndexFromFreqs(currNStates,(*Pmat)[*ancSeq++]);
				}
			else
				*cis++=*ancSeq++;
			rates++;
			}
		}
	else
		{model->UpdatePMatrix(*Pmat,transformedBranch); //get whole matrix, the sites will either be rate 1 or 0
		for (int i=0; i < nc; i++)
			{if(*rates++>0.0)
					*cis++=GetRandomIndexFromFreqs(currNStates,(*Pmat)[*ancSeq++]);
			else
					*cis++=*ancSeq++;
			}
		}
			
}	

void SimNodeLikeAttr::SimulateSeq(int nc,short *ancSeq)
{	//no rate heterogeneity
	short* cis=CharsInShorts;
	model->UpdatePMatrix(*Pmat,blen->val/(1.0-model->GetPInv()));
	for (int i=0; i < nc; i++)
		*cis++=GetRandomIndexFromFreqs(currNStates,(*Pmat)[*ancSeq++]);
}





