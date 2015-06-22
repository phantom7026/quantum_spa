
#include "galois_field.h"

/*
 *  GaloisField.h v. 1.0	
 *  
 *  Created by Yongsoo Hwang on 7/10/12.
 *
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */



void GetPrimitivePolynomial(unsigned int* primitive_poly, unsigned int gf_extension)
{
	// memory initialization
	memset(primitive_poly, 0, (gf_extension+1)*sizeof(unsigned int));

	// the characteristic of the Galois Field is 2, GF(2^m).	
	switch (gf_extension) 
	{
		case 1:
			// 1 + x = 0;
			primitive_poly[0] = primitive_poly[1] = 1;
			break;
		case 2:
			// 1 + x + X^2 = 0;
			primitive_poly[0] = primitive_poly[1] = primitive_poly[2] = 1;
			break;
		case 3: 
			// 1 + x + x^3 = 0
			primitive_poly[0] = primitive_poly[1] = primitive_poly[3] = 1; 
			break; 
		case 4: 
			// 1 + x + x^4 = 0
			primitive_poly[0] = primitive_poly[1] = primitive_poly[4] = 1; 
			break; 
		case 5: 
			// 1 + x^2 + x^5 = 0;
			primitive_poly[0] = primitive_poly[2] = primitive_poly[5] = 1; 
			break; 
		case 6: 
			// 1 + x + x^6 = 0;
			primitive_poly[0] = primitive_poly[1] = primitive_poly[6] = 1; 
			break; 
		case 7: 
			// 1 + x^3 + x^7 = 0;
			primitive_poly[0] = primitive_poly[3] = primitive_poly[7] = 1; 
			break; 
		case 8: 
			// 1 + x^2 + x^3 + x^4 + x^8 = 0
			primitive_poly[0] = primitive_poly[2] = primitive_poly[3] = primitive_poly[4] = primitive_poly[8] = 1; 
			break; 
		case 9: 
			// 1 + x^4 + x^9 = 0;
			primitive_poly[0] = primitive_poly[4] = primitive_poly[9] = 1; 
			break; 
		case 10: 
			// 1 + x^3 + x^10 = 0;
			primitive_poly[0] = primitive_poly[3] = primitive_poly[10] = 1; 
			break; 

	}
}


void MakeGfTable(unsigned int** gf_table, unsigned int* primitive_poly, unsigned int gf_extension)
{
	unsigned int gf_size = 2<<(gf_extension-1);
	
	for(int i=0; i<gf_size; ++i) 
		memset(gf_table[i], 0, gf_extension*sizeof(unsigned int));
	
	for(int i=1; i<gf_size; ++i)
	{
		if(i<=gf_extension)
		{
			gf_table[i][i-1] = 1;
		}
		else
		{
			for(int k=0; k<gf_extension; ++k)
			{
				if(gf_table[i-1][k]) gf_table[i][k+1] = 1;
			}
			if(gf_table[i][gf_extension])
			{
				for(int k=0; k<gf_extension; ++k)
					gf_table[i][k]+=primitive_poly[k];
			}
		}
		for(int j=0; j<gf_extension; ++j) gf_table[i][j]=gf_table[i][j]&1; //%2;
	}
}

void MakeGfTableB(bool** gf_table, bool* primitive_poly, unsigned int gf_extension)
{
	unsigned int gf_size = 2<<(gf_extension-1);

	for(int i=0; i<gf_size; ++i) 
		memset(gf_table[i], 0, gf_extension*sizeof(bool));
	
	for(int i=1; i<gf_size; ++i)
	{
		if(i<=gf_extension)
		{
			gf_table[i][i-1] = 1;
		}
		else
		{
			for(int k=0; k<gf_extension; ++k)
			{
				if(gf_table[i-1][k]) gf_table[i][k+1] = 1;
			}
			if(gf_table[i][gf_extension])
			{
				for(int k=0; k<gf_extension; ++k)
				{
					if(primitive_poly[k]) gf_table[i][k]=!gf_table[i][k];
					//gf_table[i][k]+=primitive_poly[k];
				}
			}
		}
		for(int j=0; j<gf_extension; ++j) gf_table[i][j]=gf_table[i][j]&1; //%2;
	}
}

void GetPrimitivePolynomialB(bool* primitive_poly, unsigned int gf_extension)
{
	// memory initialization
	memset(primitive_poly, 0, (gf_extension+1)*sizeof(bool));

	// the characteristic of the Galois Field is 2, GF(2^m).	
	switch (gf_extension) 
	{
		case 1:
			// 1 + x = 0;
			primitive_poly[0] = primitive_poly[1] = 1;
			break;
		case 2:
			// 1 + x + X^2 = 0;
			primitive_poly[0] = primitive_poly[1] = primitive_poly[2] = 1;
			break;
		case 3: 
			// 1 + x + x^3 = 0
			primitive_poly[0] = primitive_poly[1] = primitive_poly[3] = 1; 
			break; 
		case 4: 
			// 1 + x + x^4 = 0
			primitive_poly[0] = primitive_poly[1] = primitive_poly[4] = 1; 
			break; 
		case 5: 
			// 1 + x^2 + x^5 = 0;
			primitive_poly[0] = primitive_poly[2] = primitive_poly[5] = 1; 
			break; 
		case 6: 
			// 1 + x + x^6 = 0;
			primitive_poly[0] = primitive_poly[1] = primitive_poly[6] = 1; 
			break; 
		case 7: 
			// 1 + x^3 + x^7 = 0;
			primitive_poly[0] = primitive_poly[3] = primitive_poly[7] = 1; 
			break; 
		case 8: 
			// 1 + x^2 + x^3 + x^4 + x^8 = 0
			primitive_poly[0] = primitive_poly[2] = primitive_poly[3] = primitive_poly[4] = primitive_poly[8] = 1; 
			break; 
		case 9: 
			// 1 + x^4 + x^9 = 0;
			primitive_poly[0] = primitive_poly[4] = primitive_poly[9] = 1; 
			break; 
		case 10: 
			// 1 + x^3 + x^10 = 0;
			primitive_poly[0] = primitive_poly[3] = primitive_poly[10] = 1; 
			break; 

	}
}

