/*
 *  nonbinary spa.cpp v. 1.3.0	(random perturbation version)
 *  							
 *  stabilizer & msg -> reduction
 *
 *  Created by Yongsoo Hwang on 7/10/12.
 *	Modified by Yongsoo Hwang on 4/09/13.
 *	Modified by Yongsoo Hwang on 2/23/15. 
 *	- changed the parallel parts
 *	Modified by YHwang on 3/12/15.
 *
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */


// System Library 
#include <iostream>
#include <fstream>
#include <math.h>
#include <time.h>
#include <mpi.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>


// User Defined library
#include "galois_field.h"
#include "channel.h"
#include "syndrome.h"
#include "utility.h"


/* Definition of Global variables */

unsigned int* _row_weight;             // row weight of stabilizer generators
unsigned int** _row_op; 
unsigned int** _row_id;

// galois field
unsigned int _gf_extension;


/*
 * Main Function
 */

using namespace std;

inline unsigned int GetMax(double* array, int error_types)
{
	unsigned int answer=0;
	
	for(int i=1; i<error_types; ++i) if(array[answer] < array[i]) answer = i;
	return answer;
}


int main(int argc, char * argv[])
{
	/* variables for MPI */
	int process_id, nprocs;

	int process_root = 0;
	int tag = 100;
	
    MPI_Status status;
	/* MPI Settings */
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &process_id);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
	
	
	// argv[1]: simulation_option (0: naive SPA, 1: random perturbation, 2: enhanced decoding)
	// argv[2]: file name
	// argv[3]: # of simulations for each case
	// argv[4]: the characteristic of Galois Field, GF(2^m)
	
    if(argv[3]==NULL || argv[2]==NULL || argv[1]==NULL) 
    {
        if(process_id == process_root) 
        {
        	cout << "\n ===================================================================================================" << endl;
            cout << "\n ===================================================================================================" << endl;
            cout << "\t\t SPA decoding of Quantum LDPC codes " << endl;
            cout << " ===================================================================================================" << endl;
            cout << "\n USAGE : ";
            cout << " mpirun -np 'process N' 'execution program' 'data file' 'test cases' 'GF extension'" << endl;
            cout << "\n ===================================================================================================\n" << endl;
        }
        MPI_Finalize();
        exit(0);
    }
	/* variables for SPA */
	
	// code data
	string file = argv[1]; 
    
	// SPA parameters
	//unsigned int simulation_option = atoi(argv[1]);				// Applying Random Perturbation or Enhanced Feedback or Naive ?
	unsigned int simulation_cases = 10;            			 	// error cases
	unsigned int spa_iterations = 100;							// spa iterations
	unsigned int test_cases = atoi(argv[2]);					// simulations an each error case.
	
	unsigned int period = 10;
	double error_p = 0.01;
	double range = 1.0;
	
	double* error_prob; 
	error_prob = (double*)malloc(sizeof(double)*simulation_cases);
	
	unsigned int _M, _N;             			// spec of stabilizer generators

	unsigned int _max_row_weight;
	unsigned int _max_col_weight;


	// galois field setting
	_gf_extension = atoi(argv[3]);

	unsigned int _gf_base = 2;
	unsigned int _2_gf_extension = 2*_gf_extension;
	unsigned int _gf_size = _gf_base << (_gf_extension-1);
	unsigned int _gf_size_1 = _gf_size-1;
	unsigned int _error_types = _gf_size * _gf_size;
	unsigned int _error_types_1 = _error_types-1;
	
	
	///////////////////////////////////////////////////////////////////////////////////////////
	///																						///
	///								Read stabilizer generators								///
	///																						///
	///////////////////////////////////////////////////////////////////////////////////////////
	
	if(process_id == process_root)
	{
		ifstream input;
        
        try
        {
        	input.open(file.c_str());
           	input >> _M >> _N;
                        
            unsigned int* temp_counter_col; 
            temp_counter_col = (unsigned int*)malloc(sizeof(unsigned int)*_N);
            
            unsigned int temp_counter, temp_input;
            
            memset(temp_counter_col, 0, sizeof(unsigned int)*_N);
            
            for(int i=0; i<_M; ++i)
            {
                temp_counter = 0;
                
                for(int j=0; j<_N; ++j)
                {	
                    input >> temp_input;
                    
					if(temp_input > 0)
                    {
                        ++temp_counter;
                        ++temp_counter_col[j];
                    }
                }
                if(temp_counter>_max_row_weight) _max_row_weight = temp_counter;
            }
            
            for(int j=0; j<_N; ++j)
            {
                if(_max_col_weight<temp_counter_col[j]) _max_col_weight = temp_counter_col[j];
            }
            
            free(temp_counter_col);
            input.close();
        }
        catch(exception e)
        {
            cout << " Error during Read Data File ! " << endl;
        }
	}
	
	MPI_Bcast(&_M, 1, MPI_INT, process_root, MPI_COMM_WORLD);
	MPI_Bcast(&_N, 1, MPI_INT, process_root, MPI_COMM_WORLD);
	MPI_Bcast(&_max_row_weight, 1, MPI_INT, process_root, MPI_COMM_WORLD);
	MPI_Bcast(&_max_col_weight, 1, MPI_INT, process_root, MPI_COMM_WORLD);
	
	unsigned int** _col_id;
	unsigned int* _col_weight;

	_row_id = (unsigned int**)malloc(sizeof(unsigned int*)*_M);
	_row_op = (unsigned int**)malloc(sizeof(unsigned int*)*_M);
	_col_id = (unsigned int**)malloc(sizeof(unsigned int*)*_max_col_weight);
	_row_weight = (unsigned int*)malloc(sizeof(unsigned int)*_M);
	_col_weight = (unsigned int*)malloc(sizeof(unsigned int)*_N);
	
	for(int i=0; i<_M; ++i)
	{
		_row_id[i] = (unsigned int*)malloc(sizeof(unsigned int)*_max_row_weight);
		_row_op[i] = (unsigned int*)malloc(sizeof(unsigned int)*_max_row_weight);
	}
	for(int i=0; i<_max_col_weight; ++i)
	{
		_col_id[i] = (unsigned int*)malloc(sizeof(unsigned int)*_N);
	}
	
	//////////////////////////////////////////////////////////////////////////////////////
	// variables used to reduce computation
	
	unsigned int _size_gf_extension_ui;		// _gf_extension * sizeof(unsigned int);
	unsigned int _size_M_ui;				// _M * sizeof(unsigned int);
	unsigned int _size_2_gf_extension_d;
	unsigned int _size_max_row_weight_ui; 
	
	unsigned int _MRE, _NCE, _NE;
	unsigned int _MRE_d, _NCE_d, _NE_d, _RE_d;

	_MRE = _M*_max_row_weight<<_2_gf_extension;
	_NCE = _N*_max_col_weight<<_2_gf_extension;
	_NE = _N<<_2_gf_extension;

	_MRE_d = _MRE*sizeof(double);
	_NCE_d = _NCE*sizeof(double);
	_NE_d = _NE*sizeof(double);
	_RE_d = _M*_max_row_weight*sizeof(double);
	
	_size_max_row_weight_ui = _max_row_weight * sizeof(unsigned int);
	_size_gf_extension_ui = _gf_extension*sizeof(unsigned int);
	_size_M_ui = _M*sizeof(unsigned int);
	_size_2_gf_extension_d = sizeof(double)<<_2_gf_extension;

	// list_2_gf_extension[j] = j << _2_gf_extension;
	unsigned int* list_2_gf_extension;
	list_2_gf_extension = (unsigned int*)malloc(sizeof(unsigned int)*_N);
	
	for(int i=0; i<_N; ++i)
		list_2_gf_extension[i] = i<<_2_gf_extension;

	//////////////////////////////////////////////////////////////////////////////////////	
	
	ifstream input;
	try
    {
        input.open(file.c_str());
        
        input >> _M >> _N;															// reading the size of stabilizers	

        unsigned int temp_input, temp_idx;
            
        for(int i=0; i<_M; ++i) 
        {
            temp_idx = 0;
            
            for(int j=0; j<_N; ++j)
            {
                input >> temp_input;												// reading character
                
                if(temp_input)														// if the character is bigger than 1
                {																	// it is nontrivial stabilizer (Pauli operator)
                    _row_id[i][temp_idx] = j;
                    _row_op[i][temp_idx++] = temp_input;                    
                }
            }
            _row_weight[i] = temp_idx;
        }
            
        memset(_col_weight, 0, sizeof(unsigned int)*_N);
        
        for(int i=0; i<_M; ++i)
        {	
            for(int j=0; j<_row_weight[i]; ++j)
            {	
                temp_idx = _row_id[i][j];
                    
                _col_id[_col_weight[temp_idx]][temp_idx] = i;
                ++_col_weight[temp_idx];
            }
        }			
        input.close();
    }
    catch(exception e)
    {
        cout << " Error during Read Data File ! " << endl;
    }
	
	////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////
	
	// display the stabilizer generators.
	/*
	if(process_id == 0)
	{
		cout << " display of the stabilizer generators " << endl;
		for(unsigned int i=0; i<_M; ++i)
		{
			for(unsigned int j=0; j<_row_weight[i]; ++j)
			{
				cout << _row_id[i][j] << "(" << _row_op[i][j] << ") ";
			}
			cout << endl;
		}
		
	}
	*/
	
	
	// memory for messages: MSGQ and MSGC
	double msgQ[_max_col_weight][_NE]; 					//[_N<<_2_gf_extension];
	double pmsgQ[_max_col_weight][_NE]; 				//[_N<<_2_gf_extension]; 
	double msgC[_M][_max_row_weight<<_2_gf_extension]; 	// * _error_types];

	/*
	if(process_id == 0)
	{
		for(int i=0; i<_max_col_weight; ++i)
		{
			for(int j=0; j<_N; ++j)
			{
				printf("(");
				for(int k=0; k<_error_types; ++k)
				{
					printf("%.4f ", msgQ[i][j*_error_types + k]);
				}
				printf(")\t");
			}
			printf("\n");
		}
	}
	*/
	
	// memory for index used to update msgQ and msgC	
	unsigned int msgQ_targetR[_M][_max_row_weight];
	unsigned int msgC_targetC[_max_col_weight][_N];
    
    
	// msgQ_targetR, msgC_targetC setting..
	unsigned int target;
    
    // msgQ_targetR
    for(int i=0; i<_M; ++i)
    {
    	for(int j=0; j<_row_weight[i]; ++j)
        {
        	target = _row_id[i][j];
            for(int m=0; m<_col_weight[target]; ++m)
            {
            	if(_col_id[m][target] == i)
                {
                	msgQ_targetR[i][j] = m;
                    break;
                }
            }
        }
    }
    // msgC_targetC
    for(int j=0; j<_N; ++j)
    {
    	for(int i=0; i<_col_weight[j]; ++i)
        {
            target = _col_id[i][j];
            for(int m=0; m<_row_weight[target]; ++m)
            {
                if(_row_id[target][m] == j)
                {
                	msgC_targetC[i][j] = m;
                    break;
                }
            }
        }
    }
    
    
    
	double prior[_N][_error_types];									// memory for prior probability
	unsigned int real_syndrome[_M];		// +1 or 0					// memory for real syndrome from real error
	
	double nfactor; 												// memory for normalizing
	double answer[_error_types];
	double process_answer[_error_types];
	double temp[_error_types];


	// channel error probability	
	for(int i=0; i<simulation_cases; ++i) error_prob[i] = error_p*(i+1);

	unsigned int flag_success;				// +1 or 0
	unsigned int success_counter;
	unsigned int idx;
	
	unsigned long long int temp_number;
	
	unsigned int* incompatible; 		// +1 or 0		// variable to compare the real syndrome and the estimated syndrome
	incompatible = (unsigned int*)malloc(sizeof(unsigned int)*_M);
	
	
	unsigned int* estimated_error; 
	estimated_error = (unsigned int*)malloc(sizeof(unsigned int)*_N);
	
	unsigned int* estimated_syndrome; 
	estimated_syndrome = (unsigned int*)malloc(sizeof(unsigned int)*_M);
	
	unsigned int* error; 
	error = (unsigned int*)malloc(sizeof(unsigned int)*_N);


	unsigned int* temp_vector;
    unsigned int* real_vector;
    	        		
    temp_vector = (unsigned int*)malloc(sizeof(unsigned int)*_max_row_weight);
    real_vector = (unsigned int*)malloc(sizeof(unsigned int)*_max_row_weight);
    	        		
    unsigned int keep_going;						// 1 or 0
	unsigned int target_col, target_row;
	unsigned int error_weight, estimated_error_weight;
    
    double starting, ending;						// variables to check working time
      	
  	double prior_each;								// prob/_error_types;
	double sum;
	double temp_belief[_error_types];



	// variables for temporary index								
    unsigned int temp_idx, temp_idx_2, temp_idx_3, temp_weight;
    
/*
 * 	parallel working boundary setting.
 */
 	unsigned int process_bounds[2];
 	GetParamsRange_row(_M, nprocs, process_id, process_bounds);
 	
 	unsigned long long int* total_checking;
	total_checking = (unsigned long long int*)malloc(sizeof(unsigned long long int)*_M);
	
	for(int i=0; i<_M; ++i)
	{
		total_checking[i] = (unsigned long long int)pow(_error_types, (_row_weight[i]-1));
	}
	
 	
/*
 * 	Galois Field setting 
 */
 	bool* primitive_polynomialB;
 	primitive_polynomialB = (bool*)malloc(sizeof(bool)*(_gf_extension+1));
 	bool** _gf_tableB;
 	_gf_tableB = (bool**)malloc(sizeof(bool*)<<_gf_extension);
 	for(int i=0; i<_gf_size; ++i)
 	{
 		_gf_tableB[i] = (bool*)malloc(sizeof(bool)*_gf_extension);
 	}
 	
	
	GetPrimitivePolynomialB(primitive_polynomialB, _gf_extension);
	MakeGfTableB(_gf_tableB, primitive_polynomialB, _gf_extension);
	
	
	free(primitive_polynomialB);
	
	if(process_id == 0)
	{
		//system("clear");
		cout << "\n=========================================================================" << endl;
		cout << " Decoding Result of the code \"" << file.c_str() << "\"" << endl;
		cout << "=========================================================================" << endl;
		cout << " Decoding option : enhanced feedback " << endl;
		cout << "-------------------------------------------------------------------------" << endl;
    	starting = MPI_Wtime();
	}	
	///

	unsigned int stabX, stabZ;
	unsigned int errorX, errorZ;
	unsigned int tempA, tempB;
	unsigned int syndrome_answer;
	
	bool* syndrome_poly;
	syndrome_poly = (bool*)malloc(sizeof(bool)*_gf_extension);
	
	unsigned int* frustrated_checks; 		// the list of frustrated check nodes
	frustrated_checks = (unsigned int*)malloc(sizeof(unsigned int)*_M); 


//////////////////////////////////////////////////////////////////////////////////////////
//	    Main SPA working     
//////////////////////////////////////////////////////////////////////////////////////////

	//for(int qubit_idx=0; qubit_idx<1/*_N*/; ++qubit_idx)
	//for(int error_idx=1; error_idx<_error_types; ++error_idx)
	for(int sim_idx=0; sim_idx<simulation_cases; ++sim_idx)
	{
		success_counter = 0;
		keep_going = 1;

        error_weight = (unsigned int)(error_prob[sim_idx]*_N);
        
        // if error weight is nontrivial 
        if(error_weight > 0)
        {
        	// prior probability of the depolarizing channel
        	prior_each = error_prob[sim_idx]/(_error_types-1);
        	
        	for(int test_idx = 0; test_idx < test_cases; ++test_idx)
        	{
        		// memory initialization
        		memset(prior, 0, _NE_d); 			
        		memset(msgC, 0, _MRE_d); 
				memset(msgQ, 0, _NCE_d); 

				//feedback_iteration = 0;		// used for enhanced feedback
				
        		// initialization: prior according to the channel
				for(int j=0; j<_N; ++j)
				{
					prior[j][0] = 1-error_prob[sim_idx]; 
					for(int k=1; k<_error_types; ++k) prior[j][k] = prior_each;
				}
				
				// Initialization: msgC
				for(int i=0; i<_M; ++i)
				{
                	for(int j=0; j<_row_weight[i]; ++j)
                	{
	                    temp_idx = list_2_gf_extension[j]; 
	                    temp_idx_2 = _row_id[i][j];						// qubit index
	                    memcpy(&msgC[i][temp_idx], &prior[temp_idx_2][0], _size_2_gf_extension_d); 
					}
    	        }
				
				// error generation & its syndrome computation				
				memset(real_syndrome, 0, _size_M_ui);
				
				if(process_id == process_root)
				{
					memset(error, 0, sizeof(unsigned int)*_N);
					//error[qubit_idx] = error_idx;
					
					depolarizing_channel_nonbinary_old(error, error_weight, _error_types, _N);
					ComputeSyndromeExtendedB(_gf_tableB, real_syndrome, error, _M); 
                    
                    cout << " real error : " << endl;
                    for(unsigned int i=0; i<_N; ++i) cout << error[i] << " ";
                    cout << endl;
				}
				
				MPI_Bcast(real_syndrome, _M, MPI_INT, process_root, MPI_COMM_WORLD);

////////////////////////////////////////////////////////////////////////////////////////
// SPA iteration
////////////////////////////////////////////////////////////////////////////////////////

    	  		for(int spa_idx = 0; spa_idx < spa_iterations; ++spa_idx)
    	  		{
////////////////////////////////////////////////////////////////////////////////////////
// 								update msgQ
////////////////////////////////////////////////////////////////////////////////////////
					memset(pmsgQ, 0, _NCE_d); 
					if(process_id == 0) memset(msgQ, 0, _NCE_d);
					
					for(int i=process_bounds[0]; i<=process_bounds[1]; ++i)
					{
						for(int j=0; j<_row_weight[i]; ++j)
						{
							memset(answer, 0, _size_2_gf_extension_d); //sizeof(double)<<_2_gf_extension);
							
							for(unsigned long long int number = 0; number < total_checking[i]; ++number)
							{
								memset(temp_vector, 0, _size_max_row_weight_ui);
								temp_number = number;
								
								idx = 0;
								
								do{
									temp_vector[idx++] = temp_number&_error_types_1;
									temp_number>>=_2_gf_extension;
								}while(temp_number);
								
								idx = 0;
								
								for(int v_idx=0; v_idx<j; ++v_idx) 
									real_vector[v_idx] = temp_vector[idx++];
									
								for(int v_idx=j+1; v_idx<_row_weight[i]; ++v_idx)
									real_vector[v_idx] = temp_vector[idx++];
								
									
								for(int k=0; k<_error_types; ++k)
								{
									real_vector[j] = k;
									temp[k] = 1.0;
									
									//////////////////////////////////////////////////////
									// originally, the following part was in the function "ComputeSyndromeShorted"
				
									memset(syndrome_poly, 0, sizeof(bool)*_gf_extension);
									
									for(int j=0; j<_row_weight[i]; ++j)
									{
										if(real_vector[j])
										{
											tempA = tempB = 0;
							
											stabX = _row_op[i][j]>>_gf_extension; ///_gf_size;
											stabZ = _row_op[i][j]%_gf_size;
			
											errorX = real_vector[j]>>_gf_extension; 
											errorZ = real_vector[j]%_gf_size;
			
											if(stabX && errorZ) tempA = ((stabX + errorZ - 2)%(_gf_size_1)+1);
											if(stabZ && errorX) tempB = ((stabZ + errorX - 2)%(_gf_size_1)+1);
			
											for(int m=0; m<_gf_extension; ++m) 
											{
												if(_gf_tableB[tempA][m]) syndrome_poly[m]=!syndrome_poly[m];
												if(_gf_tableB[tempB][m]) syndrome_poly[m]=!syndrome_poly[m];
											}
										}
									}
				
									syndrome_answer = 0;
									for(int m=0; m<_gf_extension; ++m) 
										syndrome_answer = syndrome_answer*_gf_size + syndrome_poly[m];
										
									//////////////////////////////////////////////////////
									//if(ComputeSyndromeShorted(real_vector, i, _size_gf_extension_ui) == real_syndrome[i])
									if(syndrome_answer == real_syndrome[i])
									{
										for(int m=0; m<j; ++m) 
											temp[k]*=msgC[i][list_2_gf_extension[m] + real_vector[m]];	
											
										for(int m=j+1; m<_row_weight[i]; ++m)
											temp[k]*=msgC[i][list_2_gf_extension[m] + real_vector[m]];
											
										answer[k]+=temp[k];
									}
								}
							}
							
							target_col = _row_id[i][j];
							target_row = msgQ_targetR[i][j];
							
							sum = answer[0];
							for(int k=1; k<_error_types; ++k) sum+=answer[k];
							nfactor = 1.0/sum;
							
							temp_idx = list_2_gf_extension[target_col]; //target_col<<_2_gf_extension; //*_error_types;
							
							for(int k=0; k<_error_types; ++k) 
								pmsgQ[target_row][temp_idx++] = nfactor*answer[k];
						}
					}
					
					MPI_Barrier(MPI_COMM_WORLD);
					MPI_Reduce(pmsgQ, msgQ, _NCE, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); 
					
					if(process_id == 0)
					{	MPI_Send(msgQ, _NCE, MPI_DOUBLE, 1, tag, MPI_COMM_WORLD); 	}
					else if(process_id == 1)
					{	MPI_Recv(msgQ, _NCE, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD, &status);	}

					// display msg Q
					/*
					if(process_id == 0)
					{
						for(int i=0; i<_N; ++i)
						{
							for(int j=0; j<_col_weight[i]; ++j)
							{
								printf("[%d, %d] (", j, i);
								for(int k=0; k<_error_types; ++k)
								{
									printf("%.4f ", msgQ[j][i*_error_types + k]);
								}
								cout << ") ";
							}
							cout << endl;
						}
					}
					*/
///////////////////////////////////////////////////////////////////////////////////////
					
///////////////////////////////////////////////////////////////////////////////////////
//					update Belief & check stop condition
///////////////////////////////////////////////////////////////////////////////////////
					
					if(process_id == 0)
    	        	{
    	        		estimated_error_weight = 0;
						memset(estimated_error, 0, sizeof(unsigned int)*_N);
						
						for(int j=0; j<_N; ++j)
    	        		{
    	        			memcpy(temp_belief, prior[j], _size_2_gf_extension_d); //sizeof(double)<<_2_gf_extension); 
    	        			
    	        			temp_idx = list_2_gf_extension[j]; //j<<_2_gf_extension; // * _error_types;

    	        			for(int i=0; i<_col_weight[j]; ++i)
    	        			{
    	        				for(int k=0; k<_error_types; ++k)
    	        				{
    	        					temp_belief[k]*= msgQ[i][temp_idx+k];
    	        				}
    	        			}
    	        			
    	        			estimated_error[j] = GetMax(temp_belief, _error_types);
    	        			if(estimated_error[j]) ++estimated_error_weight;
						}
    	        		
// check stop condition
    	        		
    	        		//////////////////////////////////////////////////////////////////
    	        		// originally, this following part is in the function "ComputedSyndromeExtended"
    	        		//ComputeSyndromeExtended(_gf_table, estimated_syndrome, estimated_error);
    	        		
    	        		memset(estimated_syndrome, 0, _size_M_ui);
    	        		for(int i=0; i<_M; ++i)
    	        		{
    	        			memset(syndrome_poly, 0, sizeof(bool)*_gf_extension);
    	        			for(int j=0; j<_row_weight[i]; ++j)
    	        			{
    	        				if(estimated_error[_row_id[i][j]])
								{
									tempA = tempB = 0;
			
									stabX = _row_op[i][j]>>_gf_extension; 		//_row_op[i][j]/_gf_size;
									stabZ = _row_op[i][j]%_gf_size;
				
									errorX = estimated_error[_row_id[i][j]]>>_gf_extension; 
									errorZ = estimated_error[_row_id[i][j]]%_gf_size;

									if(stabX && errorZ) tempA = ((stabX + errorZ - 2)%(_gf_size_1)+1);
									if(stabZ && errorX) tempB = ((stabZ + errorX - 2)%(_gf_size_1)+1);
				
									for(int k=0; k<_gf_extension; ++k)
									{
										if(_gf_tableB[tempA][k]) syndrome_poly[k]=!syndrome_poly[k];
										if(_gf_tableB[tempB][k]) syndrome_poly[k]=!syndrome_poly[k];
									} 
								}
    	        			}
    	        			for(int k=0; k<_gf_extension; ++k) 
    	        				estimated_syndrome[i]=(estimated_syndrome[i]<<_gf_extension) + syndrome_poly[k];
    	        		}
    	        		
    	        		//////////////////////////////////////////////////////////////////
    	        		
    	        		flag_success = 1;
    	        		
    	        		for(int i=0; i<_M; ++i)
						{
							incompatible[i] = 0;
							if(estimated_syndrome[i]!=real_syndrome[i])
							{
								incompatible[i] = 1;
								flag_success = 0;
							}
						}
						
						if(flag_success)
						{
							 ++success_counter;
						}
    	        	}
 	        		
 	        		MPI_Barrier(MPI_COMM_WORLD);
    	        	MPI_Bcast(&flag_success, 1, MPI_INT, 0, MPI_COMM_WORLD);
					if(flag_success==1) break;
					
					///////////////////////////////////////////////////////////////////////////////////////////
					
					///////////////////////////////////////////////////////////////////////////////////////////
					//				update msg C
					///////////////////////////////////////////////////////////////////////////////////////////
					
					memset(msgC, 0, _MRE_d);
					
    	        	if(process_id == 1)
    	        	{
// update msgC
    	        		for(int j=0; j<_N; ++j)
						{
                        	temp_idx_2 = list_2_gf_extension[j]; 
							// temp_idx_2 = j*_error_types;
							
							for(int i=0; i<_col_weight[j]; ++i)
							{
								memcpy(temp, prior[j], _size_2_gf_extension_d);
							
								// from 0 to i
								for(int m=0; m<i; ++m)
									for(int k=0; k<_error_types; ++k) temp[k]*=msgQ[m][temp_idx_2+k];
									
								// from i+1 to col_weight[j]
								for(int m=i+1; m<_col_weight[j]; ++m)
									for(int k=0; k<_error_types; ++k) temp[k]*=msgQ[m][temp_idx_2+k];


								target_row = _col_id[i][j];
    	                        target_col = msgC_targetC[i][j];
                    
                    			sum = temp[0];
                    			for(int k=1; k<_error_types; ++k) sum+=temp[k];
                    			nfactor = 1.0/sum; 
                    			
	                            temp_idx_3 = list_2_gf_extension[target_col]; 
	                            //temp_idx_3 = target_col * _error_types;
    	                        
    	                        for(int k=0; k<_error_types; ++k) msgC[target_row][temp_idx_3++] = nfactor * temp[k];
    	                        //msgC[target_row][temp_idx_3+k] = nfactor * temp[k];
	                        }
						}
					}
					
// update prior	by random perturbation				
					if(process_id == 0)
					{
						if(spa_idx && spa_idx%period == 0)
						{
							// prior initialization
							for(int j=0; j<_N; ++j)
							{
								prior[j][0] = 1-error_prob[sim_idx]; 
								for(int k=1; k<_error_types; ++k) prior[j][k] = prior_each;
							}
							
							srand(time(NULL));
							for(int i=0; i<_M; ++i)
							{
								if(incompatible[i])
								{
                    		        for(int j=0; j<_row_weight[i]; ++j)
                    		        {
                    	                temp_idx = _row_id[i][j];
                     		            sum = 0;
                        		            
                    		            sum = temp[0] = prior[temp_idx][0];
											                        		            
                    		            for(int k=1; k<_error_types; ++k)
                        		        {
                        		           	temp[k] = (1+range*((rand()+clock())&255) /256 ) * prior[temp_idx][k];
                        	            	sum+=temp[k];
                        	            }
        									
        								nfactor = 1.0/sum;                		            
        									
										for(int k=0; k<_error_types; ++k) prior[temp_idx][k] = nfactor*temp[k];
                        	        }
                        	    }
							}
						}
    	        	}
    	        	
    	        	MPI_Barrier(MPI_COMM_WORLD);
    	        	MPI_Bcast(msgC, _MRE, MPI_DOUBLE, 1, MPI_COMM_WORLD);
    		        
    		        if(spa_idx && spa_idx%period==0) MPI_Bcast(prior, _NE, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    	        }
        	}
        }
        else
        {
        	success_counter = test_cases;
        }
        
        // checking the terminaiton condition
        if(process_id == process_root)
		{
			double ber = 1-(double)success_counter/test_cases;
			
			printf(" channel error %.3f, ", error_prob[sim_idx]);
			printf(" #success %5d, BER: %.4f \n", success_counter, ber);

			if(ber>=0.3) keep_going = 0;
		}
		
		MPI_Bcast(&keep_going, 1, MPI_INT, process_root, MPI_COMM_WORLD);
		if(!keep_going) break;
	}
	if(process_id == process_root)
	{
		ending = MPI_Wtime();
		cout << "\n working time is " << ending-starting << endl;
	}



///////////////////////////////////////////////////////////////////////////////////
// 		Memory Free
///////////////////////////////////////////////////////////////////////////////////
	{
	free(list_2_gf_extension);
	free(syndrome_poly);
	free(error);
	free(error_prob);
	free(frustrated_checks);
	free(estimated_error);
	free(estimated_syndrome);
	
	free(incompatible);
	 
	free(_row_weight);
	free(_col_weight);

	for(int i=0; i<_M; ++i)
	{
		free(_row_op[i]);
		free(_row_id[i]);
	}

	free(_row_op);
	free(_row_id);
	
	for(int i=0; i<_max_col_weight; ++i) free(_col_id[i]);
	free(_col_id);

	for(int i=0; i<_gf_size; ++i) free(_gf_tableB[i]);
	free(_gf_tableB);
	}
	
	MPI_Finalize();
}
