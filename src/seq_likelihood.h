//	Copyright (C) 2008 Mark T. Holder
//
//	chimne_ssweep is free software; you can redistribute it and/or modify
//	it under the terms of the GNU General Public License as published by
//	the Free Software Foundation; either version 2 of the License, or
//	(at your option) any later version.
//
//	chimne_ssweep is distributed in the hope that it will be useful,
//	but WITHOUT ANY WARRANTY; without even the implied warranty of
//	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
//	GNU General Public License for more details.
//
//	You should have received a copy of the GNU General Public License
//	along with NCL; if not, write to the Free Software Foundation, Inc., 
//	59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

#if ! defined(SEQ_LIKELIHOOD_H)
#define SEQ_LIKELIHOOD_H

#ifdef __cplusplus
extern "C" 
{
#endif

#include "ncl/nxscdiscretematrix.h"
#include "cPhyProb/ccore/dsct_model.h"


typedef struct {
	PyObject_HEAD
	StateSetLookupStruct * sharedStateSetLookupStruct;
	LeafDataObj ** leafData;
	CLAObj ** clas;
	FullLAObj * treeLike;
	ASRVObj * asrv;
	DSCTModelObj ** model;
	PMatArrayObj ** pmats;
	unsigned nLeafData;
	unsigned nCLA;
	unsigned nPMat;
	unsigned nModels;
} LikeStructsBundle;


/* 
	Assumes that the NxsCDiscreteMatrix outlives the LikeStructsBundle (some
	pointers may be aliases.
*/
LikeStructsBundle newLikeStructs(const NxsCDiscreteMatrix matrix, unsigned nCLAs, unsigned nPMats, unsigned nRates);
LikeStructsBundle newPartitionedLikeStructs(const NxsCDiscreteMatrix matrix, unsigned nCLAs, unsigned nPMats,unsigned nSubsets);
void freeLikeStructFields(LikeStructsBundle *); 
void zeroLikeStructFields(LikeStructsBundle *); 


#ifdef __cplusplus
}
#endif
#endif
