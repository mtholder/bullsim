#ifndef LIKESETTINGS#define LIKESETTINGS#include <map>#include "basicbulldefs.h"#include "models.h"#include "MTHException.h"class LikeSettingException :public MTHException {public: LikeSettingException(const char *p) : MTHException(p) {}};/*		The purpose of these the settings classes is to provide an interface	with the program and allow checking that all of the user selections are legal without allocating memory 	to do the likelihood calculations.  Essentially these are the classes that would be directly affected	by a lset command.  Once a LScore (or search) are started.  These classes are used to set up the appropriate data	partitions, parameter alteration information, and the like attributes for the tree and nodes	LSettings constructor should be called right after the data is read into the file.  	IMPORTANT -- after being constructed all the taxa names are added to IncludedTaxa, the arguments should be the start and end of the 	data matrix.  Thus LSettings should be warned if any taxa or data are excluded	Right now I have it in the FinishCharactersBlock()*/class ModSettings{	bool stationaryModel,useCurrentBrLen;	double startingBrLen;	Model *mod;	PartModIndex ownerOfBlen,ownerOfBlenMod;//change ownerOfBlenMod to vector to allow relaxed clock starting points (multiple modifiers with one model)	int blenSetting,blenModSetting,numBrLenMod;	bool hasBrLenMod;	public :	ModSettings(Model *);	void SetUseCurrentBrLen(bool);	Model *GetModel();	bool StationaryModel();	PartModIndex GetBrLenOwner()	{return ownerOfBlen;}	PartModIndex GetBrLenModOwner()	{return ownerOfBlenMod;}	int GetBLenSetting()	{ return blenSetting;}	int GetBLenModSetting()	{ return blenModSetting;}	bool HasBrLenMod()	{return hasBrLenMod;}	int GetNumBrLenMod()	{return numBrLenMod;}};class PartSettings{	friend class LikePartSettings;	public :	bool AligningSeqs;	typedef pair<int,int> segment;	vector<segment> segmentsInPartition;		PartSettings(int start,int end);	int GetNSegments();	pair<int,int> GetSegment(int i);	bool DoingAlignment();	int GetRawDataSize();};class LikePartSettings{	public :	int nmodels;	vector<ModSettings *> modSettings;	int ModelMixingSetting;	FreqParamGroup *ModelMixingStartingVals;		LikePartSettings(PartSettings*);	void SetUseCurrentBrLen(bool);	int GetNModels();	Model *GetModel(int i);	bool AllStationaryModels();	int GetModelMixingParamSettings()	{return ModelMixingSetting;}	double GetStartingFreqnOfThisModel(int j);	PartModIndex GetBrLenOwner(int j)	{assert(j<nmodels);		return modSettings[j]->GetBrLenOwner();		}	PartModIndex GetBrLenModOwner(int j)	{assert(j<nmodels);		return modSettings[j]->GetBrLenModOwner();		}	int GetBLenSetting(int i)	{return modSettings[i]->GetBLenSetting();}	int GetBLenModSetting(int i)	{return modSettings[i]->GetBLenModSetting();}	bool HasBrLenMod(int j)	{return modSettings[j]->HasBrLenMod();}	int GetNumBrLenMod(int m)	{return modSettings[m]->GetNumBrLenMod();}};class DataSettings {	friend class LikeSettings;	friend class ParsSettings;	protected :	int nDataPartitions;	vector<nxsstring> IncludedTaxa;	bool userInputBrLen;	public:	vector<PartSettings *> dataPartSettings;	int GetNDataPartitions();	int GetRawDataSize(int i);	bool TreatGapsAsFifth()	{return false;}	void CopyAPartitionsRawData(short *destination,map<nxsstring,EncodedChars *>&rawData,int num,int taxnum);	int GetNIncludedTaxa()	{return IncludedTaxa.size();}	bool DoingAlignment();	bool DoingAlignment(int i);	nxsstring GetNameOfNthIncludedTaxon(int i)	{return IncludedTaxa[i];}	void AddTaxon(nxsstring i)	{IncludedTaxa.push_back(i);}};class LikeSettings {	bool summingOverAllAncStates;	public:	bool optimizingBrLens;	bool maximize;//true if ML false if just scoring one instantiation of parameters	int maxPasses;//max number of smoothing passes for each of the parameters sets	double delta;// smallest change in likelihood score that will keep optimization going	DataSettings *sharedDataSettings;		bool useRSBranchLengths;//Rogers Swofford approximation	bool sumOverAncStates;// whether using Felsensteins pruning vs calculating the like of one reconstruction	public:	vector<LikePartSettings *> partSettings;	LikeSettings(DataSettings *);//see note above	void SetUseCurrentBrLen(bool);	int GetNModelsInPart(int i);	Model *GetModel(int i,int j);	bool SummingOverAllAncStates();	bool AllStationaryModels();	int GetModelMixingParamSettings(int i)	{return partSettings[i]->GetModelMixingParamSettings();}	double GetStartingFreqnOfThisModel(int i,int j) {return partSettings[i]->GetStartingFreqnOfThisModel(j);}	int GetNSetsOfBranchLengths();	PartModIndex GetOwnerOfNthSetOfBrLengths(int ownint);	PartModIndex GetOwnerOfBrLens(int p,int  m);	int GetNSetsOfBranchLengthModifiers();	PartModIndex GetOwnerOfNthSetOfBrLengthModifiers(int ownint);	PartModIndex GetOwnerOfBrLenModifiers(int p,int  m);	int GetBLenSetting(int p,int m) 	{return partSettings[p]->GetBLenSetting(m);}	int GetBLenModSetting(int p,int m) 	{return partSettings[p]->GetBLenModSetting(m);}	bool HasBrLenMod(int p,int m)	{return partSettings[p]->HasBrLenMod(m);}	int GetNumBrLenMod(int p,int m)	{return partSettings[p]->GetNumBrLenMod(m);}};class ParsSettings {	ParsSettings(DataSettings *);//see note above	};#endif