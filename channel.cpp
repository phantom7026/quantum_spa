
#include "channel.h"

/*
void depolarizing_channel_nonbinary(unsigned int* error, unsigned int error_weight, unsigned int _error_types, unsigned int block_length)
{
	unsigned int qubit;
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::default_random_engine generator(seed);
	
	std::uniform_int_distribution<int> distribution_type(1, _error_types-1);
	std::uniform_int_distribution<int> distribution_qubit(0, block_length-1);
	//std::uniform_real_distribution<double> distribution_X_Z(0.0, 1.0);
	
	double rnd_dbl;
	int rnd_int, i;

	i=0;
	while((i++) < error_weight)
	{
		do{
			qubit = distribution_qubit(generator);
		}while(error[qubit]);
		
		//rnd_dbl = distribution_X_Z(generator);
		rnd_int = distribution_type(generator);

		//if(rnd_dbl<0.5) error[qubit] = gf_size * rnd_int;
		//else error[qubit] = rnd_int;
		error[qubit] = rnd_int;
	}
}
*/

void depolarizing_channel_nonbinary_old(unsigned int* error, unsigned int error_weight, unsigned int _error_types, unsigned int block_length)
{
    unsigned int qubit, idx, rnd_int;
    idx=0;
    
    while((idx++)<error_weight)
    {
        do{
            qubit = (rand() + clock())%(block_length-1);
        }while(error[qubit]);
        
        rnd_int = (rand() + clock())%(_error_types-1);
        error[qubit] = rnd_int;
    }
}

