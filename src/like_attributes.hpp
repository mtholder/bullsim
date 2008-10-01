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
 

#ifndef LIKEATTRIBUTES
#define LIKEATTRIBUTES


#if defined HAVE_CONFIG_H && HAVE_CONFIG_H
#	include <config.h>
#endif 

#define BITSPERSHORT 16
#include "model.hpp"
#include "ssrf_codon_sub_mod.hpp"
#include "basic_bull.hpp"
#include "ncl/nxsdefs.h"
namespace bull {
enum mem {SHAREPMAT=1};
class	LikeAttr
{
	protected:
	//int memoryMode;
	Model *model;
	double ***Pmat;
	BoundedParameter *blen;
	PositiveParameter *blenMod;
	public :
	static double Multiplier,currNTax;
	static int currNChar,currNStates,currShortPerChar,currNRateCats,currNStatesInLastShort;

#	if defined (CHAR_BY_CHAR) && CHAR_BY_CHAR
		static int currCharIndex;
#	endif
	
	static bool currModelDirty;
	LikeAttr(); //Constructor for internal node's attributes
	virtual LikeAttr *Copy();
	void SetBLen(double d); 
	virtual ~LikeAttr();  //DANGER deletes blen and blenMod either add owns flag or make sure mem isn't deleted twice another way

	double* GetBLenPtr()	{return &(blen->val); }
	Model* GetModel();
	LikeAttr(Model *m); //Constructor for internal node's attributes
	void SetBLenPtr(BoundedParameter* b)	{/*delete blen; */ blen=b; }
	BoundedParameter* GetBLenParameterPtr() {return blen; }
	void DetachBLenPtr()	{blen=NULL; }
};

class SimNodeLikeAttr : public LikeAttr {
	//like a TerminalNodeLikeAttr, but it owns its Characters
	public:
	short *CharsInShorts;
	SimNodeLikeAttr (int,Model *m); //Constructor for Simulation node's attributes
	virtual ~SimNodeLikeAttr () {delete [] CharsInShorts; }
	void SimulateSeq(int nc,short *ancSeq,double *rates);
	void SimulateSeq(int nc,short *ancSeq);
		
};



class TreeSimAttr {
	public:
	Model *model;
	double **freq;
	int nChar,nStates;
	double *rates;
	bool hasRateHet,ownsBrLens;
	
	TreeSimAttr(int nc,Model *m);
	~TreeSimAttr()	{delete [] rates; }
	double **GetStateFreqs() {return freq; }
	int GetNChar() {return nChar; }
	int GetNStates() {return nStates; }
	Model* GetModel()	{return model; }
	void GenerateRandomSequence(int nc,short *dest);
	void GenerateNewSetOfRates(int nc);
	void SetOwnsBrLen(bool i) {ownsBrLens=i; } 
	bool GetOwnsBrLen() {return ownsBrLens; } 
	
};

}// namespace bull 


#endif
