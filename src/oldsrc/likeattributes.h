#ifndef LIKEATTRIBUTES#define LIKEATTRIBUTES#define BITSPERSHORT 16#include "models.h"#include "basicbulldefs.h"#include "nexusdefs.h"#include "AdvancedBullDefs.h"/* Like attributes are owned by a node and store the information to calculate likelihoods	for that node.  all pointer variable are pointers to pointers to allow for multiple model/partitoins	*/struct VirtualLikeAttExcep: public MTHException{public: VirtualLikeAttExcep(const char *p) : MTHException(p) {}};enum mem {SHAREPMAT=1};class   LikeAttr	{	protected:	//int memoryMode;	Model *model;	double ***Pmat;	BoundedParameter *blen;	PositiveParameter *blenMod;	public :	static double Multiplier,currNTax;	static int currNChar,currNStates,currShortPerChar,currNRateCats,currNStatesInLastShort;#ifdef CHARBYCHAR	static int currCharIndex;#endif	static bool currModelDirty;	LikeAttr();//Constructor for internal node's attributes	LikeAttr(Model *m);//Constructor for internal node's attributes	~LikeAttr();  //DANGER deletes blen and blenMod either add owns flag or make sure mem isn't deleted twice another way	void UpdatePmat(double x);	//bool SharePmatMemory();	virtual LikeAttr *Copy();	inline double ***getPmat();	Model* GetModel();	BoundedParameter* GetBLenParameterPtr()	{return blen;}	double* GetBLenPtr()	{return &(blen->val);}	void DetachBLenPtr()	{blen=NULL;}	void SetBLenPtr(BoundedParameter* b)	{delete blen; blen=b;}	void SetBLen(double d);		void MultBLen(double d);	void SetBLenModPtr(PositiveParameter* b)	{blenMod=b;}	PositiveParameter* GetBLenModPtr()	{return blenMod;}		//Virtuals that should be called by InternalNodeLikeAttr Only	virtual double *GetCondLikeArray() 		{throw VirtualLikeAttExcep("Calling virtual GetCondLikeArray"); return NULL;} //Should be called by InternalNodeLikeAttr only	virtual void SetThisCondLikeBasedOnAncArg(double *) {throw VirtualLikeAttExcep("Calling virtual SetThisCondLikeToLikeBelowBasedOnArg(double *)");} //Should be called by one of the InternalLikeAttr	virtual void SetThisCondLikeToLikeBelowBasedOnArg(double **){throw VirtualLikeAttExcep("Calling virtual SetThisCondLikeToLikeBelowBasedOnArg(double **)");} //Should be called by one of the InternalLikeAttr	//Virtuals that should be called by TerminalNodeLikeAttr Only	virtual void FillDoubleArrayWithCondLikeOfChars(double *cl) {throw VirtualLikeAttExcep("Calling virtual FillDoubleArrayWithCondLikeOfChars");}	virtual void RecodeDueToModelSizeChange()	{throw VirtualLikeAttExcep("Calling virtual RecodeDueToModelSizeChanges");}		//Virtuals that should be called by a NodeLikeAttrs Only	//virtual void SetArgToCondLikeOfThisSubTree(double *,bool) {throw VirtualLikeAttExcep("Calling virtual SetArgToCondLikeOfThisSubTree");} //Should be called by one of the NodeLikeAttr	//virtual void MultArgByCondLikeOfThisSubTree(double *,bool) {throw VirtualLikeAttExcep("Calling virtual SetArgToCondLikeOfThisSubTree");}	virtual void SetArgToCondLikeOfThisSubTree(double *) {throw VirtualLikeAttExcep("Calling virtual SetArgToCondLikeOfThisSubTree");} //Should be called by one of the NodeLikeAttr	virtual void MultArgByCondLikeOfThisSubTree(double *) {throw VirtualLikeAttExcep("Calling virtual MultArgByCondLikeOfThisSubTree");} //Should be called by one of the NodeLikeAttr	};class InternalNodeLikeAttr : public LikeAttr {	long sizeofcondlikeArray;		public:		double *condlike;		InternalNodeLikeAttr(int n,Model *m);//Constructor for internal node's attributes	InternalNodeLikeAttr(int n,Model *m,double *prealloc);	void Reset();	~InternalNodeLikeAttr();	double *GetCondLikeArray() {return condlike;}	//void SetArgToCondLikeOfThisSubTree(double *,bool);	void SetArgToCondLikeOfThisSubTree(double *);	//void MultArgByCondLikeOfThisSubTree(double *,bool);	void MultArgByCondLikeOfThisSubTree(double *);	void SetThisCondLikeBasedOnAncArg(double *);	void SetThisCondLikeToLikeBelowBasedOnArg(double *);	void SetThisCondLikeToLikeBelowBasedOnArg(double **);	};class TerminalNodeLikeAttr : public LikeAttr {	public:#ifdef ELIMINATEALLZEROS	int maxNObsStates;#endif	short *CharsInShorts;	int nObsStates,*OStates;	double *OStatesCL;	TerminalNodeLikeAttr(short *,Model *m);//Constructor for internal node's attributes	~TerminalNodeLikeAttr()	{delete [] OStates;delete [] OStatesCL;}	//void SetArgToCondLikeOfThisSubTree(double *,bool);	void SetArgToCondLikeOfThisSubTree(double *);	//void MultArgByCondLikeOfThisSubTree(double *,bool);	void MultArgByCondLikeOfThisSubTree(double *);	void FillDoubleArrayWithCondLikeOfChars(double *cl);		void RecodeDueToModelSizeChange();	void SetArgToCondLikeOfThisSubTreeWithOutSharedMatrix(double *cl);	void SetArgToCondLikeOfThisSubTreeWithZeroBranchLength(double *cl);	};class TreeLikeAttr 	{	double ***Pmat;	public:	Model *model;	double currentLikelihood;	double **freq;	double *likeBelow,*likeAbove,*prevl,*thisl;	int *isConstant;	bool ownsBrLens,needsToBeRescored;	int nChar,nStates,nRateCat,nColPerChar;	short *nTimesThisPatternOccurs,*constChars;	vector<PositiveParameter *> bLenModifiers;	vector<bool> ownsBLenModifier;int nspc,nsilc;/*these are owned at the Tree level because they are probably shared by all on the nodes in a tree (the use is most likely different rates for different characters, the exception to this isin a relaxed clock version in which there could be several bLenModifiers.Note that bLenModifiers may be shared by more than one Partition/Model so there is a vector of ownership flagsthis prevents double deleting, BUT DOES NOT assure that all other partition/models that use this TreeLikeAtt's bLenModifiersare notified that the modifiers are being deleted.*/	TreeLikeAttr(int n,Model *m,int *isConst,short *nReps,short *cc);	~TreeLikeAttr();		void UpdatePmat(double x);	//bool SharePmatMemory();	inline double ***getPmat();	Model* GetModel();	void SetLikelihood(double x)	{currentLikelihood=x;}	int GetNChar() {return nChar;}	int GetNStates() {return nStates;}	int GetNRateCats() {return nRateCat;}	int GetNShortsPerChar() {return nspc;}	int GetNStatesInLastShort() {return nsilc;}	bool GetOwnsBrLen() {return ownsBrLens;} //Should be called by TreeLikeAttr	double **GetStateFreqs() {return freq;}	int *GetIsConstant()	{return isConstant;}	short *GetConstPats()	{return constChars;}	short *GetNReps()	{return nTimesThisPatternOccurs;}	void SwapThisLikeArray()	{		double *temp;		temp=prevl;		prevl=thisl;		thisl=prevl;		}	double *GetThisLikeArray(bool swapThenRet) {		if(swapThenRet)	SwapThisLikeArray();		return thisl;		} //Should be called by TreeLikeAttr	void SetOwnsBrLen(bool i) {ownsBrLens=i;} //Should be called by TreeLikeAttr	void AddBrLenMod(PositiveParameter *blenMod,bool i)	{bLenModifiers.push_back(blenMod); ownsBLenModifier.push_back(i);}	PositiveParameter *GetBrLenMod(int blmi)	{return bLenModifiers[blmi];}	void NotifyModelOfChangedParameter(Parameter *p)		{needsToBeRescored=true; model->ParameterHasChanged(p);}	void NotifyModelOfChangedParameter(FreqParamGroup *f)	{needsToBeRescored=true; model->ParameterHasChanged(f);}	void NotifyModelsThatFreqParamChangesAreDone(FreqParamGroup *f)	{needsToBeRescored=true; model->FreqParamChangesShouldSumToOne(f);}#ifdef DZRATES	double GetOneBranchLike(vector<IntSet *> *,PositiveParameter *blen,PositiveParameter *blenMod);#else	double GetOneBranchLike(PositiveParameter *blen,PositiveParameter *blenMod);#endif	void AlertSharedMemory()	{model->AlertSharedMemory();}	double GetMultiplier()	{return model->GetMultiplier();}};class SetOfLikeAttr 	{	protected:	int nmodparts,nparts;//number of partiotions+models, and number of partitions	int *nmods;	LikeAttr*** likeAttribs;		public:		SetOfLikeAttr(int i);		~SetOfLikeAttr();		Model *GetModel(int p,int m) {return likeAttribs[p][m]->GetModel();}	void SetNumModsInPart(int partnum,int modnum);	void AddLikeAtt(int p,int m,LikeAttr *);	int GetNParts()	{return nparts;}	int GetNModels(int i) {return nmods[i];}	LikeAttr *GetLikeAtt(int p,int m)	{return likeAttribs[p][m];}	SetOfLikeAttr *Copy();		};class SetOfTreeLikeAttr  	{#ifdef DZRATES	vector<IntSet *> *rateSets;#endif	friend class Tree;	int nmodparts,nparts;//number of partiotions+models, and number of partitions	int *nmods;	TreeLikeAttr*** likeAttribs;		FreqParamGroup **modelMixingParam;	double **prevLArr,**thisLArr;	int *nChar,*nStates,*nRateCats,*packingType,*packingIndex;		public:	double *prevL,*thisL;#ifdef CODONHACK	double *thisLMult;#endif	SetOfTreeLikeAttr(int i);	~SetOfTreeLikeAttr();		void SetNumModsInPart(int partnum,int modnum);	void AddLikeAtt(int p,int m,TreeLikeAttr *);	int GetNParts()	{return nparts;}	int GetNModels(int i) {return nmods[i];}	Model *GetModel(int p,int m) {return likeAttribs[p][m]->GetModel();}	TreeLikeAttr *GetLikeAtt(int p,int m)	{return likeAttribs[p][m];}		void SetModelMixingParam(int i,FreqParamGroup *mmp)	{		assert(modelMixingParam && i<nparts);		modelMixingParam[i]=mmp;		}	void InitializeLikeAttrStatics(int partnum,int modn);	void FinishedAddingModels(int i,int pi);	void InitializeModelMixingParams();	int GetNChar(int p,int m) {return likeAttribs[p][m]->GetNChar();}	int GetNStates(int p,int m) {return likeAttribs[p][m]->GetNStates();}	int GetNRateCats(int p,int m) {return likeAttribs[p][m]->GetNRateCats();}	int GetNShortsPerChar(int p,int m) {return likeAttribs[p][m]->GetNShortsPerChar();}	int GetNStatesInLastShort(int p,int m) {return likeAttribs[p][m]->GetNStatesInLastShort();}	double GetLikelihoodFromCondLike(SetOfLikeAttr* rootSet,int partnum,bool OnlyScoreDirty);	double GetLikelihoodFromCondLike(LikeAttr * rootLA,TreeLikeAttr * treeLA);	void SetOwnerOfBrLengths(PartModIndex pmi,bool i)	{likeAttribs[pmi.partNum][pmi.modNum]->SetOwnsBrLen(i);}	bool GetOwnerOfBrLengths(PartModIndex pmi)	{return likeAttribs[pmi.partNum][pmi.modNum]->GetOwnsBrLen();}	void AddBrLenMod(PartModIndex pmi,PositiveParameter *blenMod,bool i) {likeAttribs[pmi.partNum][pmi.modNum]->AddBrLenMod(blenMod,i);}	PositiveParameter *GetBrLenMod(PartModIndex pmi,int blmi)	{return likeAttribs[pmi.partNum][pmi.modNum]->GetBrLenMod(blmi);}	double *GetLikeBelowArray(PartModIndex pmi)	{return likeAttribs[pmi.partNum][pmi.modNum]->likeBelow;}	double *GetLikeAboveArray(PartModIndex pmi)	{return likeAttribs[pmi.partNum][pmi.modNum]->likeAbove;}	void PrintModel(ostream &usrstream);#ifdef DZRATES	void AddRateSets(vector<IntSet *> *rs) {rateSets=rs;}	vector<IntSet *> *GetRateSets() {return rateSets;}#endif};#endif