/* *  GaloisField.h v. 1.0	 *   *  Created by Yongsoo Hwang on 7/10/12. * *  Copyright 2011 __MyCompanyName__. All rights reserved. * */#include <string.h>void GetPrimitivePolynomialB(bool* primitive_poly, unsigned int gf_extension);void MakeGfTableB(bool** gf_table, bool* primitive_poly, unsigned int gf_extension);void MakeGfTable(unsigned int** gf_table, unsigned int* primitive_poly, unsigned int gf_extension);void GetPrimitivePolynomial(unsigned int* primitive_poly, unsigned int gf_extension);