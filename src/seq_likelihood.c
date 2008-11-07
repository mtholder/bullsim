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
#include <assert.h>
#include "seq_likelihood.h"
#define ALL_IN_ONE
#include "cPhyProb/ccore/dsct_model.c"

LikeStructsBundle newPartitionedLikeStructs(
  const NxsCDiscreteMatrix matrix, 
  unsigned nCLAs, 
  unsigned nPMats,
  unsigned nSubsets)
{
	unsigned arr_len = 2*matrix.nStates; 
	unsigned i, j;
	int ** p = 0L;
	LikeStructsBundle toReturn;
	zeroLikeStructFields(&toReturn);
	toReturn.nModels = nSubsets;

	if (!(matrix.matrix && matrix.stateList))
		goto errorExit;
		

	/*Allocate p - a copy of the stateList*/
	/* each of the "fundamental" states needs to spaces in the array" */
	for (i = matrix.nStates; i < matrix.nObservedStateSets; ++i) {
		arr_len += 1; //add one for the # of states element
		unsigned pos = matrix.stateListPos[i] ;
		arr_len += matrix.stateList[pos];
	}

	/*allocate a ragged two-D array with the memory contiguous.*/
	p = (int**)malloc(matrix.nObservedStateSets*sizeof(int*));
	if (p == 0L)
		goto errorExit;	
	p[0] = (int*)malloc(arr_len*sizeof(int));
	if (p[0] == 0L) 
		goto errorExit;

	int * curr_p = p[0];
	for (i = 0 ; i < matrix.nStates; ++i) {
		p[i] = curr_p;
		p[i][0] = 1;
		p[i][1] = i;
		curr_p += 2;
	}
	for (i = matrix.nStates; i < matrix.nObservedStateSets; ++i) {
		p[i] = curr_p;
		unsigned pos = matrix.stateListPos[i];
		p[i][0] = matrix.stateList[pos];
		curr_p++;
		unsigned j;
		for (j = 1; j <= p[i][0]; ++j) {
			*curr_p++ = matrix.stateList[pos + j];
		}
	}
	
	toReturn.sharedStateSetLookupStruct = sslookup_new(matrix.nStates, matrix.nObservedStateSets, p);
	if (toReturn.sharedStateSetLookupStruct == 0L)
		goto errorExit;
	p = 0L; /* p has given the pointer to toReturn.sharedStateSetLookupStruct so we set it to NULL to avoid double deletion */
	
	
	
	/* allocate the leafData and zero leafData pointer array*/
	toReturn.nLeafData = matrix.nTax;
	toReturn.leafData = (LeafDataObj **)malloc(matrix.nTax*sizeof(LeafDataObj*));	
	if (toReturn.leafData == 0L)
		goto errorExit;
	for (i = 0; i < matrix.nTax; ++i)
		toReturn.leafData[i] = NULL;
	/* allocate each leafData object and fill its data array with the same state codes as were used in the NCL NxsCDiscreteStateSet */
	for (i = 0; i < matrix.nTax; ++i) {
		toReturn.leafData[i] = leaf_data_new(matrix.nChar, 1, toReturn.sharedStateSetLookupStruct);
		if (!(toReturn.leafData[i]))
			goto errorExit;
		int * dest = toReturn.leafData[i]->ssind;
		assert(dest);
		NxsCDiscreteStateSet * sourceRow = matrix.matrix[i];
		if (!sourceRow)
			goto errorExit;
		for (j = 0; j < matrix.nChar; ++j)
			*dest++ = (int)(*sourceRow++);
	}
	
	toReturn.nCLA = nCLAs;
	toReturn.clas = (CLAObj **)malloc(nCLAs*sizeof(CLAObj *));	
	if (toReturn.clas == 0L)
		goto errorExit;
	for (i = 0; i < nCLAs; ++i)
		toReturn.clas[i] = NULL;
	for (i = 0; i < nCLAs; ++i) {
		toReturn.clas[i] = cla_ss_partitioned_new(matrix.nChar, matrix.nStates);
		if (!(toReturn.clas[i]))
			goto errorExit;
	}
	
	toReturn.treeLike = full_la_ss_partitioned_new(matrix.nChar, matrix.nStates);
	if (toReturn.treeLike == 0L)
		goto errorExit;
	
	for (i = 0; i < nSubsets; ++i) {
		(*toReturn.treeLike).state_categ_freqs[i][matrix.nStates] = 1.0;
	}

	toReturn.asrv = 0L;

	toReturn.model = (DSCTModelObj **)malloc(nSubsets*sizeof(DSCTModelObj *));	
	if (!toReturn.model)
		goto errorExit;
	for (i = 0; i < toReturn.nModels; ++i) {
		toReturn.model[i] = 0L;
	}
	for (i = 0; i < toReturn.nModels; ++i) {
		toReturn.model[i] = dsct_model_new(matrix.nStates);
		if (!toReturn.model[i])
			goto errorExit;
	}

	toReturn.nPMat = nPMats;
	toReturn.pmats = (PMatArrayObj **)malloc(nPMats*sizeof(PMatArrayObj *));	
	if (toReturn.pmats == 0L)
		goto errorExit;
	for (i = 0; i < nPMats; ++i)
		toReturn.pmats[i] = NULL;
	for (i = 0; i < nPMats; ++i) {
		toReturn.pmats[i] = cpmat_array_new(1, matrix.nStates);
		if (!(toReturn.pmats[i]))
			goto errorExit;
	}

	return toReturn;

	errorExit:
		PRINTF("In newLikeStructs errorExit\n");
		if ((toReturn.sharedStateSetLookupStruct == 0L) && p) {
			free(p[0]);
			free(p);
		}
		freeLikeStructFields(&toReturn);
		return toReturn;
}

LikeStructsBundle newLikeStructs(
  const NxsCDiscreteMatrix matrix, 
  unsigned nCLAs, 
  unsigned nPMats,
  unsigned nRates)
{
	unsigned arr_len = 2*matrix.nStates; 
	unsigned i, j;
	double rateCatFreq;
	int ** p = 0L;
	LikeStructsBundle toReturn;
	zeroLikeStructFields(&toReturn);
	toReturn.nModels = 1;

	if (!(matrix.matrix && matrix.stateList))
		goto errorExit;
		

	/*Allocate p - a copy of the stateList*/
	/* each of the "fundamental" states needs to spaces in the array" */
	for (i = matrix.nStates; i < matrix.nObservedStateSets; ++i) {
		arr_len += 1; //add one for the # of states element
		unsigned pos = matrix.stateListPos[i] ;
		arr_len += matrix.stateList[pos];
	}

	/*allocate a ragged two-D array with the memory contiguous.*/
	p = (int**)malloc(matrix.nObservedStateSets*sizeof(int*));
	if (p == 0L)
		goto errorExit;	
	p[0] = (int*)malloc(arr_len*sizeof(int));
	if (p[0] == 0L) 
		goto errorExit;

	int * curr_p = p[0];
	for (i = 0 ; i < matrix.nStates; ++i) {
		p[i] = curr_p;
		p[i][0] = 1;
		p[i][1] = i;
		curr_p += 2;
	}
	for (i = matrix.nStates; i < matrix.nObservedStateSets; ++i) {
		p[i] = curr_p;
		unsigned pos = matrix.stateListPos[i];
		p[i][0] = matrix.stateList[pos];
		curr_p++;
		unsigned j;
		for (j = 1; j <= p[i][0]; ++j) {
			*curr_p++ = matrix.stateList[pos + j];
		}
	}
	
	toReturn.sharedStateSetLookupStruct = sslookup_new(matrix.nStates, matrix.nObservedStateSets, p);
	if (toReturn.sharedStateSetLookupStruct == 0L)
		goto errorExit;
	p = 0L; /* p has given the pointer to toReturn.sharedStateSetLookupStruct so we set it to NULL to avoid double deletion */
	
	
	
	/* allocate the leafData and zero leafData pointer array*/
	toReturn.nLeafData = matrix.nTax;
	toReturn.leafData = (LeafDataObj **)malloc(matrix.nTax*sizeof(LeafDataObj*));	
	if (toReturn.leafData == 0L)
		goto errorExit;
	for (i = 0; i < matrix.nTax; ++i)
		toReturn.leafData[i] = NULL;
	/* allocate each leafData object and fill its data array with the same state codes as were used in the NCL NxsCDiscreteStateSet */
	for (i = 0; i < matrix.nTax; ++i) {
		toReturn.leafData[i] = leaf_data_new(matrix.nChar, nRates, toReturn.sharedStateSetLookupStruct);
		if (!(toReturn.leafData[i]))
			goto errorExit;
		int * dest = toReturn.leafData[i]->ssind;
		assert(dest);
		NxsCDiscreteStateSet * sourceRow = matrix.matrix[i];
		if (!sourceRow)
			goto errorExit;
		for (j = 0; j < matrix.nChar; ++j)
			*dest++ = (int)(*sourceRow++);
	}
	
	toReturn.nCLA = nCLAs;
	toReturn.clas = (CLAObj **)malloc(nCLAs*sizeof(CLAObj *));	
	if (toReturn.clas == 0L)
		goto errorExit;
	for (i = 0; i < nCLAs; ++i)
		toReturn.clas[i] = NULL;
	for (i = 0; i < nCLAs; ++i) {
		toReturn.clas[i] = cla_new(matrix.nChar, matrix.nStates, nRates);
		if (!(toReturn.clas[i]))
			goto errorExit;
	}
	
	toReturn.treeLike = full_la_new(matrix.nChar, matrix.nStates, nRates);
	if (toReturn.treeLike == 0L)
		goto errorExit;
	
	rateCatFreq = 1.0/((double) nRates);
	for (i = 0; i < nRates; ++i) {
		(*toReturn.treeLike).state_categ_freqs[i][matrix.nStates] = rateCatFreq;
	}

	toReturn.asrv = asrv_obj_new(nRates, 1, 0.5);
	if (toReturn.asrv == 0L)
		goto errorExit;

	toReturn.model = (DSCTModelObj **)malloc(1*sizeof(DSCTModelObj *));	
	if (!toReturn.model)
		goto errorExit;
	for (i = 0; i < toReturn.nModels; ++i) {
		toReturn.model[i] = 0L;
	}
	for (i = 0; i < toReturn.nModels; ++i) {
		toReturn.model[i] = dsct_model_new(matrix.nStates);
		if (!toReturn.model[i])
			goto errorExit;
	}

	toReturn.nPMat = nPMats;
	toReturn.pmats = (PMatArrayObj **)malloc(nPMats*sizeof(PMatArrayObj *));	
	if (toReturn.pmats == 0L)
		goto errorExit;
	for (i = 0; i < nPMats; ++i)
		toReturn.pmats[i] = NULL;
	for (i = 0; i < nPMats; ++i) {
		toReturn.pmats[i] = cpmat_array_new(nRates, matrix.nStates);
		if (!(toReturn.pmats[i]))
			goto errorExit;
	}

	return toReturn;

	errorExit:
		PRINTF("In newLikeStructs errorExit\n");
		if ((toReturn.sharedStateSetLookupStruct == 0L) && p) {
			free(p[0]);
			free(p);
		}
		freeLikeStructFields(&toReturn);
		return toReturn;
}


void zeroLikeStructFields(LikeStructsBundle * toFree)
{
	toFree->sharedStateSetLookupStruct = 0L;
	toFree->leafData = 0L;
	toFree->clas = 0L;
	toFree->treeLike = 0L;
	toFree->asrv = 0L;
	toFree->model = 0L;
	toFree->pmats = 0L;
}

// frees all of the memory that the LikeStructsBundle object points to, but does
// NOT free the LikeStructsBundle pointer itself.
void freeLikeStructFields(LikeStructsBundle * toFree)
{
	int i;
	if (toFree == 0L)
		return;
	if (toFree->sharedStateSetLookupStruct) {
		sslookup_dtor(toFree->sharedStateSetLookupStruct);
		toFree->sharedStateSetLookupStruct = 0L;
	}
	if (toFree->asrv) {
		asrv_obj_dtor(toFree->asrv);
		toFree->asrv = 0L;
	}

	if (toFree->leafData) {
		for (i = 0; i < toFree->nLeafData; ++i) {
			if (toFree->leafData[i]) {
				leaf_data_dtor(toFree->leafData[i]);
			}
		}
		free(toFree->leafData);
		toFree->leafData = 0L;
	}

	if (toFree->clas) {
		for (i = 0; i < toFree->nCLA; ++i) {
			if (toFree->clas[i]) {
				cla_dtor(toFree->clas[i]);
			}
		}
		free(toFree->clas);
		toFree->clas = 0L;
	}

	if (toFree->treeLike) {
		full_la_dtor(toFree->treeLike);
		toFree->treeLike = 0L;
	}
	if (toFree->model) {
		for (i = 0; i < toFree->nModels; ++i) {
			if (toFree->model[i]) {
				cdsctm_dtor(toFree->model[i]);
			}
		}
		free(toFree->model);
		toFree->model = 0L;
	}
	if (toFree->pmats) {
		for (i = 0; i < toFree->nPMat; ++i) {
			if (toFree->pmats[i]) {
				cpmat_array_dtor(toFree->pmats[i]);
			}
		}
		free(toFree->pmats);
		toFree->pmats = 0L;
	}
}


