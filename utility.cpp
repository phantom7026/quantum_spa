
#include "utility.h"
// used to find the most probable qudit-wise error

/*
unsigned int GetMax(double* array, int error_types)
{
	unsigned int answer=0;
	
	for(int i=1; i<error_types; ++i) if(array[answer] < array[i]) answer = i;
	return answer;
}
*/

// function for dividing the whole working coverage into that of each thread.
void GetParamsRange(const long long int size, int nprocs, int nprocs_2, int process_id, unsigned long long int* process_working)
{
	unsigned long long int temp1, temp2;
	temp1 = size >> nprocs_2;
	temp2 = size & (nprocs-1);
	
	process_working[0] = process_id*temp1 + (process_id>temp2?temp2:process_id); 
	process_working[1] = process_working[0]+temp1-1;
	
	if(temp2 > process_id) process_working[1]+=1;
}

void GetParamsRange_row(const int size, int nprocs, int process_id, unsigned int* process_working)
{
	unsigned int temp1, temp2;
	temp1 = size/nprocs;
	temp2 = size&(nprocs-1);
	
	process_working[0] = process_id*temp1+(process_id>temp2?temp2:process_id);
	process_working[1] = process_working[0] + temp1-1;
	
	if(temp2 > process_id) process_working[1]+=1;
}