
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

/*
 *  Syndrome.h  v. 1.0.1	
 *  
 *  Created by Yongsoo Hwang on 4/09/13.
 *
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */
 
 
extern unsigned int** _row_id;
extern unsigned int** _row_op;
extern unsigned int* _row_weight;

extern unsigned int _gf_extension;



// build the commutativity relation table for the given operator, op.
void CommutativityTable(unsigned int** _gf_table, unsigned int* table, unsigned int op, unsigned int _size_gf_extension_ui);
void CommutativityTableB(bool** _gf_table, unsigned int* table, unsigned int op);

// function to compute syndrome with an error operator of the length _row_weight
// usually called in updating msgQ
unsigned int ComputeSyndromeShorted(unsigned int** _gf_table, unsigned int* error, unsigned int row_idx, unsigned int _size_gf_extension_ui);



// function to compute syndrome with an error operator of the length N
void ComputeSyndromeExtended(unsigned int** _gf_table, unsigned int* syndrome, unsigned int* error, 
							 unsigned int _size_gf_extension_ui, unsigned int _M);
							 
void ComputeSyndromeExtendedB(bool** _gf_table, unsigned int* syndrome, unsigned int* error, unsigned int _M);
							 
