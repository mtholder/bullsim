#ifndef BASICBULLDEFS#define BASICBULLDEFS#define DEBUGGING#define CHARBYCHAR#undef DZRATES#define DEFAULTSTARTINGBRLEN 0.1#define DEFAULTMAXPASSES	20#define DEFAULTDELTA	0.000001#define FTOL 0.00000003#define ZEPS 1.0e-10#define CGOLDEN 0.381966#define GOLD	1.618034#define GLIMIT 100.0#define SMALLDOUBLE 1.0e-10#define MAXDOUBLE 1.0e303#define MINPOSDOUBLE 1.0e-303#define OVERFLOWCHECK 1.0e293#define UNDERFLOWCHECK 1.0e-253#define MINBRLEN 1.0e-10#define MAXBRLEN 1.0e10class LikeUnderFlow{};struct PartModIndex	{	int partNum,modNum;	PartModIndex()	{partNum=0; modNum=0;}	PartModIndex(int p,int m)	{partNum=p; modNum=m;}	void Set(int p,int m) {partNum=p; modNum=m;}	bool operator<(const PartModIndex s)	const {		if(partNum==s.partNum)	return (modNum<s.modNum);		return (partNum<s.partNum);		}	bool operator!=(const PartModIndex s)	const {		return (partNum!=s.partNum)|(modNum!=s.modNum);		}	bool operator==(const PartModIndex s)	const {		return (partNum==s.partNum)&(modNum==s.modNum);		}};enum gencode {mito , nuc};#endif