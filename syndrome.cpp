
#include "syndrome.h"

/*
 *  Syndrome.h  v. 1.0.2	
 *  
 *  Created by Yongsoo Hwang on 4/09/13.
 *	Modified by Yongsoo Hwang on 2/16/15.
 *
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */
 

// build the commutativity relation table for a given operator, op.
// table[i] = syndrome between i-th operator and op 

void CommutativityTable(unsigned int** _gf_table, unsigned int* table, unsigned int op, unsigned int _size_gf_extension_ui)
{
	unsigned int opX, opZ;
	unsigned int targetX, targetZ;
	unsigned int tempA, tempB;
	
	unsigned int gf_size = 2 << (_gf_extension-1);
	unsigned int gf_size_1 = gf_size-1;
	unsigned int error_types = gf_size * gf_size;
	
	unsigned int* syndrome_poly; 
	syndrome_poly = (unsigned int*)malloc(sizeof(unsigned int)*_gf_extension);
	
	memset(table, 0, error_types*sizeof(unsigned int));

	opX = op>>_gf_extension; 
	opZ = op%gf_size;
	
	for(int j=0; j<error_types; ++j)
	{
		tempA = tempB = 0;
		
		memset(syndrome_poly, 0, _size_gf_extension_ui);
		
		targetX = j>>_gf_extension; 	
		targetZ = j%gf_size;

		if(opX && targetZ){ tempA = ((opX+targetZ-2)&gf_size_1)+1; }
		if(opZ && targetX){ tempB = ((opZ+targetX-2)&gf_size_1)+1; }

		for(int k=0; k<_gf_extension; ++k) syndrome_poly[k] =_gf_table[tempA][k] + _gf_table[tempB][k];
		for(int k=0; k<_gf_extension; ++k) table[j] = (table[j]<<_gf_extension) + syndrome_poly[k]&1; 
	}
	free(syndrome_poly);
}

void CommutativityTableB(bool** _gf_table, unsigned int* table, unsigned int op)
{
	unsigned int opX, opZ;
	unsigned int targetX, targetZ;
	unsigned int tempA, tempB;
	
	unsigned int gf_size = 2 << (_gf_extension-1);
	unsigned int gf_size_1 = gf_size-1;
	unsigned int error_types = gf_size * gf_size;
	
	bool* syndrome_poly; 
	syndrome_poly = (bool*)malloc(sizeof(bool)*_gf_extension);
	
	memset(table, 0, error_types*sizeof(unsigned int));

	opX = op>>_gf_extension; 
	opZ = op%gf_size;
	
	for(int j=0; j<error_types; ++j)
	{
		tempA = tempB = 0;
		
		memset(syndrome_poly, 0, sizeof(bool)*_gf_extension); 
		
		targetX = j>>_gf_extension; 	//j/_gf_size
		targetZ = j%gf_size;
		
		if(opX && targetZ) tempA = ((opX + targetZ - 2)&gf_size_1)+1;
		if(opZ && targetX) tempB = ((opZ + targetX - 2)&gf_size_1)+1;
		
		for(int k=0; k<_gf_extension; ++k) 
		{
			syndrome_poly[k] = _gf_table[tempA][k];
			if(_gf_table[tempB][k]) syndrome_poly[k]=!syndrome_poly[k];
		}
		for(int k=0; k<_gf_extension; ++k) table[j] = (table[j]<<_gf_extension) + syndrome_poly[k]; 
	}
	free(syndrome_poly);
}



// function to compute syndrome with an error operator of the length _row_weight
// usually called in updating msgQ

unsigned int ComputeSyndromeShorted(unsigned int** _gf_table, unsigned int* error, unsigned int row_idx, unsigned int _size_gf_extension_ui)
{
	unsigned int answer;
	unsigned int stabX, stabZ;
	unsigned int errorX, errorZ;
	unsigned int tempA, tempB;
	
	unsigned int gf_size = 2 << (_gf_extension-1);
	unsigned int gf_size_1 = gf_size-1;
	unsigned int error_types = gf_size * gf_size;
	
	unsigned int* syndrome_poly; 
	syndrome_poly = (unsigned int*)malloc(sizeof(unsigned int)*_gf_extension);
	
	memset(syndrome_poly, 0, _size_gf_extension_ui); 
	
	for(int j=0; j<_row_weight[row_idx]; ++j)
	{
		if(error[j])
		{
			tempA = tempB = 0;
			
			stabX = _row_op[row_idx][j]>>_gf_extension; ///_gf_size;
			stabZ = _row_op[row_idx][j]%gf_size;
			
			errorX = error[j]>>_gf_extension; 
			errorZ = error[j]%gf_size;
			
			if(stabX*errorZ) tempA = ((stabX+errorZ-2)%gf_size_1)+1;
			if(stabZ*errorX) tempB = ((stabZ+errorX-2)%gf_size_1)+1;
			
			for(int k=0; k<_gf_extension; ++k) syndrome_poly[k]+=_gf_table[tempA][k] + _gf_table[tempB][k];
		}
	}
	
	answer = 0;
	for(int k=0; k<_gf_extension; ++k) answer = (answer<<_gf_extension) + syndrome_poly[k]&1;
	
	free(syndrome_poly);
	
	return answer;
}



// function to compute syndrome with an error operator of the length N
void ComputeSyndromeExtended(unsigned int** _gf_table, unsigned int* syndrome, unsigned int* error, 
							 unsigned int _size_gf_extension_ui, unsigned int _M)
{
	unsigned int stabX, stabZ;
	unsigned int errorX, errorZ;
	unsigned int* syndrome_poly; 
	unsigned int tempA, tempB;
	
	unsigned int gf_size = 2 << (_gf_extension-1);
	unsigned int gf_size_1 = gf_size-1;
	
	unsigned int error_types = gf_size * gf_size;
	
	syndrome_poly = (unsigned int*)malloc(sizeof(unsigned int)*_gf_extension);
	
	for(int i=0; i<_M; ++i)
	{
		memset(syndrome_poly, 0, _size_gf_extension_ui); 
		
		for(int j=0; j<_row_weight[i]; ++j)
		{
			if(error[_row_id[i][j]])
			{
				tempA = tempB = 0;
			
				stabX = _row_op[i][j]>>_gf_extension; 		//_row_op[i][j]/_gf_size;
				stabZ = _row_op[i][j]%gf_size;
				
				errorX = error[_row_id[i][j]]>>_gf_extension; 
				errorZ = error[_row_id[i][j]]%gf_size;

				if(stabX*errorZ) tempA = ((stabX + errorZ - 2)%gf_size_1)+1;
				if(stabZ*errorX) tempB = ((stabZ + errorX - 2)%gf_size_1)+1;
				
				for(int k=0; k<_gf_extension; ++k) syndrome_poly[k] +=_gf_table[tempA][k] + _gf_table[tempB][k];	
			}
		}
		// To be honest, the following calculation is wrong !!! 
		// We have to find the right gf element compatible with the above polynomial expression.
		// But, this calculation is used for every case, therefore it can be used.

		for(int k=0; k<_gf_extension; ++k) syndrome[i]=(syndrome[i]<<_gf_extension) + (syndrome_poly[k]&1);
	}
	free(syndrome_poly);
}

// function to compute syndrome with an error operator of the length N
void ComputeSyndromeExtendedB(bool** _gf_table, unsigned int* syndrome, unsigned int* error, unsigned int _M)
{
	unsigned int stabX, stabZ;
	unsigned int errorX, errorZ;
	unsigned int tempA, tempB;
	
	unsigned int gf_size = 2 << (_gf_extension-1);
	unsigned int gf_size_1 = gf_size-1;
	
	unsigned int error_types = gf_size * gf_size;

	bool* syndrome_poly; 
	syndrome_poly = (bool*)malloc(sizeof(bool)*_gf_extension);
	
	for(int i=0; i<_M; ++i)
	{
		memset(syndrome_poly, false, sizeof(bool)*_gf_extension); 
		
		for(int j=0; j<_row_weight[i]; ++j)
		{
			if(error[_row_id[i][j]])
			{
				tempA = tempB = 0;
			
				stabX = _row_op[i][j]>>_gf_extension; 		
				stabZ = _row_op[i][j]%gf_size;
				
				errorX = error[_row_id[i][j]]>>_gf_extension; 
				errorZ = error[_row_id[i][j]]%gf_size;

				if(stabX*errorZ) tempA = ((stabX+errorZ-2)%gf_size_1)+1;
				if(stabZ*errorX) tempB = ((stabZ+errorX-2)%gf_size_1)+1;
				
				for(int k=0; k<_gf_extension; ++k) 
				{
					if(_gf_table[tempA][k]) syndrome_poly[k]=!syndrome_poly[k];
					if(_gf_table[tempB][k]) syndrome_poly[k]=!syndrome_poly[k];
				}
			}
		}
		
		for(int k=0; k<_gf_extension; ++k) syndrome[i]=(syndrome[i]<<_gf_extension) + syndrome_poly[k];
	}
	free(syndrome_poly);
}