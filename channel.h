
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <stdlib.h>
#include <limits.h>

//#include <random>       // new in c++11


// use random class
void depolarizing_channel_nonbinary(unsigned int* error, unsigned int error_weight, 
									unsigned int _error_types, unsigned int block_length);

// srand() & rand() + clock()
void depolarizing_channel_nonbinary_old(unsigned int* error, unsigned int error_weight,
                                    unsigned int _error_types, unsigned int block_length);
