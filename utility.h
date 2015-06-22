

// used to find the most probable qudit-wise error

//unsigned int GetMax(double* array, int error_types);


// function for dividing the whole working coverage into that of each thread.

void GetParamsRange(const long long int size, int nprocs, int nprocs_2, int process_id, unsigned long long int* process_working);

void GetParamsRange_row(const int size, int nprocs, int process_id, unsigned int* process_working);
