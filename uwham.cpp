#include <cassert>
#include <iostream>
#include <fstream>
#include <stdlib.h> //rand()
#include <sstream> // std::istringstream
#include <time.h> // random seed
#include <vector> //array.size
#include <math.h>  // floor
#include <stdio.h> // printf
#include <unistd.h>
#include <cmath>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

// g++ uwham.cpp -lm -lgsl -lgslcblas -o uwham.o

using std::vector;
using std::cout;
using std::getline;
using std::string;
using std::ifstream;
using std::istringstream;
using std::isfinite;
using std::size_t;

const double kB = 0.0019872041;
const double pi = 3.1415926;
const double cal2J = 4.184;
const int maxneq = 100; 

int nstates = 0; // how many walks or thermodynamic states
int totalndata = 0; // total number of observations
vector<int> nobs; // how many observations at each state

vector<double> density_of_states; // this is the density of states, totalndata
vector<double> partition_function;
vector<vector<double> > bias_matrix; // this is the bias matrix totalndata*nstates

int convex_newton_raphson_solver(int niteration, double max_error, double factor, int literation) {
	int iteration=0;
	double error=1.0e10;
	vector<double> oldfe;
	vector<double> newfe;

	// void print_out_matrix(gsl_matrix *pmatrix);
	double cal_partition_function(int index);
	double cal_density_of_states(int index);
	void normalize_partition_functions();	
	double cal_free_energy(int index);

	double alpha = -1.0/totalndata;			
	int nstatesm1 = nstates - 1;
	int signum;	
	gsl_vector * gradient = gsl_vector_alloc(nstates);	
	gsl_matrix * hessian = gsl_matrix_alloc(nstates, nstates);
	gsl_matrix * weight_matrix = gsl_matrix_alloc(totalndata, nstates);
	gsl_matrix * pimatrix = gsl_matrix_alloc(nstates, nstates);
	gsl_vector * onem = gsl_vector_alloc(nstates);
	gsl_vector * onen = gsl_vector_alloc(totalndata);
	gsl_matrix * tmpm1 = gsl_matrix_alloc(nstates, totalndata);
	gsl_matrix * tmpm2 = gsl_matrix_alloc(nstates, nstates);			
	gsl_vector * tmpv1 = gsl_vector_alloc(nstates);
	gsl_matrix * tmphessian = gsl_matrix_alloc(nstatesm1, nstatesm1);
	gsl_matrix * invhessian = gsl_matrix_alloc(nstatesm1, nstatesm1);		
	gsl_permutation * perm = gsl_permutation_alloc(nstatesm1);
	gsl_vector * tmpgradient = gsl_vector_alloc(nstatesm1);
	gsl_vector * product = gsl_vector_alloc(nstatesm1);		

	// initialize vectors and matrices
	gsl_vector_set_all(onem, 1.0);
	gsl_vector_set_all(onen, 1.0);
	gsl_matrix_set_all(pimatrix, 0.0);
	for (int i=0; i<nstates; i++) {
		double tmpelem = 1.0*nobs[i]/totalndata;
		gsl_matrix_set(pimatrix, i, i, tmpelem);
	}

	// print
	printf("\033[1;32miteration=%-12d\033[0m", iteration);			
	for (int i=0; i<nstates; i++) {
		double fe = cal_free_energy(i);
		printf("%12g ", fe);
		oldfe.push_back(fe);
		newfe.push_back(fe);
	}
	printf("\n");	


	normalize_partition_functions(); 
	// update density of states
	for (int i=0; i<totalndata; i++) {
		density_of_states[i] = cal_density_of_states(i);
	}

	while ((iteration < niteration) && (error > max_error)) {
		iteration++;
		printf("\033[1;32miteration=%-12d\033[0m", iteration);				

		// calculate the weight matrix
		for (int i=0; i<totalndata; i++) {
			for (int j=0; j<nstates; j++) {
				double tmpelem = density_of_states[i]*bias_matrix[j][i]/partition_function[j]*totalndata;
				gsl_matrix_set(weight_matrix, i, j, tmpelem);
			}
		}
		// calculate the gradient matrix			
		gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, pimatrix, weight_matrix, 0.0, tmpm1);
		gsl_blas_dgemv(CblasNoTrans, alpha, tmpm1, onen, 0.0, gradient);
		gsl_matrix_set_all(hessian, 0.0);
		for (int i=0; i< nstates; i++) {
			double tmpelem = -1.0*gsl_vector_get(gradient, i);
			gsl_matrix_set(hessian, i, i, tmpelem);
		}
		gsl_blas_dgemv(CblasNoTrans, 1.0, pimatrix, onem, 0.0, tmpv1);
		gsl_vector_add(gradient, tmpv1);
		// calculate the hessian matrix
		gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, tmpm1, weight_matrix, 0.0, tmpm2);			
		gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, alpha, tmpm2, pimatrix, 1.0, hessian);			
		// calculate partition function changes
		for (int i=1; i<nstates; i++) {
			for (int j=1; j<nstates; j++) {
				double tmpelem = gsl_matrix_get(hessian, i, j);
				gsl_matrix_set(tmphessian, (i-1), (j-1), tmpelem);
			}
		}
		gsl_linalg_LU_decomp(tmphessian, perm, &signum);
		gsl_linalg_LU_invert(tmphessian, perm, invhessian);

		for (int i=1; i<nstates; i++) {
			double tmpelem = gsl_vector_get(gradient, i);
			gsl_vector_set(tmpgradient, (i-1), tmpelem);
		}
		gsl_blas_dgemv(CblasNoTrans, 1.0, invhessian, tmpgradient, 0.0, product);

		// // update partition function
		if (iteration > literation) {
		  factor = ((factor*iteration)<1.0)?(factor*(iteration-literation)):1.0;
		}
		partition_function[0] = 1.0;
		for (int i=1; i<nstates; i++) {
		  partition_function[i] *= exp(-factor*gsl_vector_get(product, (i-1)));
		}
		
		// // update density of states
		for (int i=0; i<totalndata; i++) {
			density_of_states[i] = cal_density_of_states(i);
		}		

		error = 0.0;
		for (int i=0; i<nstates; i++) {
			double fe = cal_free_energy(i);
			printf("%12g ", fe);
			newfe[i] = fe;
			if (! isfinite(fe)) {
				printf("\n\033[1;31mCannot converge. Need more self consistent iteration steps!\033[0m\n"); //what color?
				return 0;
			}
			double tmperror = fabs((newfe[i] - oldfe[i])/oldfe[i]);
			if (tmperror > error) {
				error = tmperror;
			}
			oldfe[i] = newfe[i];			
		}
		printf("\tmax error=%12g", error);
		printf("\n");
	}
	return 1;
}

double cal_partition_function(int index) {
	double pf = 0.0;
	for (int j=0; j<totalndata; j++) {
		pf += density_of_states[j]*bias_matrix[index][j];
	}
	return pf;
}

void normalize_partition_functions() {
	for (int i=1; i<nstates; i++) {
		partition_function[i] /= partition_function[0];
	}
	partition_function[0] = 1.0;
}

double cal_density_of_states(int index) {
	double ds = 0.0;

	for (int i=0; i<nstates; i++) {
		ds += nobs[i]*bias_matrix[i][index]/partition_function[i];
	}
	ds = 1.0/ds;
	return ds;
}

double cal_free_energy(int index) {
	double fe=0.0;
	if (index != 0) 
		fe = -log(partition_function[index]/partition_function[0]);
	return fe;
}


void print_out_weights(vector<int>& printitems, int normalization) {
	FILE * write_file;
	vector<double> sum(nstates, 0.0);
	double tmpw;
	int stateid, datumid;

	write_file = fopen("weights.data", "w");		
	// -w 0 means printing out everything, only print out the essential weight.
	if ((printitems.size()==1) and (printitems[0]==0)) {
		stateid = 0;
		datumid = 0;
		for (int j=0; j<totalndata; j++) {
			tmpw = density_of_states[j]*bias_matrix[stateid][j];
			fprintf(write_file, "%12g\n", tmpw);			
			datumid++;
			if (datumid == nobs[stateid]) {
				stateid++;
				datumid = 0;
			}
		}
	}
	else {
		// normalize
		if (normalization == 1) {
			for (int j=0; j<totalndata; j++) {
				for (int i=0; i<nstates; i++) {
					tmpw = density_of_states[j]*bias_matrix[i][j];
					sum[i] += tmpw;
				}
			}
		}
		// print out
		for (int j=0; j<totalndata; j++) {
			for (int i=0; i<printitems.size(); i++) {
				if (normalization == 1) {
					tmpw = density_of_states[j]*bias_matrix[printitems[i]-1][j]/sum[printitems[i]-1]; // state start form 1 not 0
				}
				else {
					tmpw = density_of_states[j]*bias_matrix[printitems[i]-1][j]; // state start form 1 not 0
				}
				fprintf(write_file, "%12g ", tmpw);
			}
			fprintf(write_file, "\n");
		}
	}	
	fclose(write_file);
}

int self_consistent_iteration(int startiteration, int niteration, double max_error) {
	int iteration = startiteration;
	double error=1.0e10;
	vector<double> oldfe;
	vector<double> newfe;

	double cal_partition_function(int index);
	double cal_density_of_states(int index);
	double cal_free_energy(int index);
	void normalize_partition_functions();	    
	
	// printf("iteration=%-12d", iteration);
	// printf("\033[1;32miteration=%-12d\033[0m", iteration);		
	for (int i=0; i<nstates; i++) {
		double fe = cal_free_energy(i);
		// printf("%12g ", fe);
		oldfe.push_back(fe);
		newfe.push_back(fe);
	}
	// printf("\n");	

	while ((iteration <= niteration) && (error > max_error)) {
		// printf("iteration=%-12d", iteration);
		printf("\033[1;32miteration=%-12d\033[0m", iteration);				
		// update density of states
		for (int i=0; i<totalndata; i++) {
			density_of_states[i] = cal_density_of_states(i);
		}
		// update partition function
		for (int i=0; i<nstates; i++) {
			partition_function[i] = cal_partition_function(i);
		}
        normalize_partition_functions();         

		error = 0.0;
		for (int i=0; i<nstates; i++) {
			double fe = cal_free_energy(i);
			printf("%12g ", fe);
			newfe[i] = fe;
			double tmperror = fabs((newfe[i] - oldfe[i])/oldfe[i]);
			if (tmperror > error) {
				error = tmperror;
			}
			oldfe[i] = newfe[i];			
		}
		printf("\tmax error=%12g", error);
		printf("\n");
		iteration++;		
	}
	return 1;
}


void shift_bias_energy_matrix() {
    for (int j=0; j<totalndata; j++) {
		double base=1.0e100;		
        for (int i=0; i<nstates; i++) {
			if (bias_matrix[i][j] < base) {
				base = bias_matrix[i][j];
			}
		}
        for (int i=0; i<nstates; i++) {
			bias_matrix[i][j] -= base;
		}
	}
}
		

void convert_bias_energy_matrix_from_gromacs(int ncol) {
	int ndel;

	for (int i=0; i<ncol; i++) {
		bias_matrix.pop_back();
	}
	ndel= bias_matrix.size() - nstates;
	bias_matrix.erase(bias_matrix.begin(), bias_matrix.begin()+ndel);
}

void convert_energy_to_bias_matrix() {
    for (int j=0; j<totalndata; j++) {
        for (int i=0; i<nstates; i++) {
            bias_matrix[i][j] = exp(-bias_matrix[i][j]);
        }
    }
}

void read_matrix(char* filename, vector<vector<double> >& matrix, vector<int>& ndata) {
	/* this function read a matrix from file. 
	   the ith one dimensional vector is the ith volumn */
	ifstream read_file;
	string line;
	double dummyx;
	int tmpndata;
	char oldline0;
	
	read_file.open(filename);	
	assert(read_file.is_open());
	
	getline(read_file, line); // read the first line
	while ((line[0] == '#') || (line[0] == '@')) {
		getline(read_file, line); 
	}
	istringstream firststring(line); // break first line into fields
	while (firststring >> dummyx) {
		vector<double> dummyarray;							
		dummyarray.push_back(dummyx);
		matrix.push_back(dummyarray);
	}
	tmpndata = 1;
	oldline0 = line[0];
			
	while (getline(read_file, line)) {
		if ((line[0] != '#') and (line[0] != '@')) {			
			istringstream readstring(line);
			int k = 0;
			while (readstring >> dummyx) {
				matrix[k].push_back(dummyx);                
				k++;
			}
			tmpndata++;
		}
		else if ((oldline0 != '#') and (oldline0 != '@')) {
			ndata.push_back(tmpndata);
			tmpndata = 0;
		}
		oldline0 = line[0];
	}
		
	read_file.close();
	ndata.push_back(tmpndata);
}

void readlist(const char* liststring, vector<int>& listitems) {
	istringstream vstring(liststring);

	while(vstring.good()) {
		string substr;
		getline(vstring, substr, ',');
		size_t found = substr.find('-');
		if (found > 0) {
			int starti = atoi(substr.substr(0, found).c_str());
			int endi = atoi(substr.substr(found+1).c_str());
			for (int i=starti; i<=endi; i++) {
				listitems.push_back(i);
			}
		}
		else {
			listitems.push_back(atoi(substr.c_str()));
		}
	}
}

double others_2_kBT(int unit, double temper) {
	double factor;

	switch(unit) {
	case 0: factor = 1.0; break;
	case 1: factor = 1.0/(kB*temper); break;
	case 2: factor = 1.0/(kB*temper)/cal2J; break;
	default:
		abort();
	}
	return(factor);
}

double kBT_2_others(int unit, double temper) {
	double factor;

	switch(unit) {
	case 0: factor = 1.0; break;
	case 1: factor = kB*temper; break;
	case 2: factor = kB*temper*cal2J; break;
	default:
		abort();
	}
	return(factor);	
}

void change_unit_of_bias_energy_matrix(double factor) {
	for (int i=0; i<nstates; i++) {
		for (int j=0; j<totalndata; j++) {
			bias_matrix[i][j] *= factor;
		}
	}	
}

void final_print_out(int unit, double factor) {
	double fe;
	const char* unitname;
	
	switch(unit) {
	case 0: unitname = "kBT"; break;
	case 1: unitname = "kcal/mol"; break;
	case 2: unitname = "kJ/mol"; break;
	default:
		abort();
	}

	printf("\n");
	printf("\033[1;33mThe free energy differences in units of \033[0m\033[1;31m%s\033[0m:\n", unitname);
	for (int i=0; i<nstates; i++) {
		fe = cal_free_energy(i);
		fe *= factor;
		printf("%12g ", fe);
	}
	printf("\n");
}

void initial_print_out(int unit) {
	const char* unitname;
	
	switch(unit) {
	case 0: unitname = "kBT"; break;
	case 1: unitname = "kcal/mol"; break;
	case 2: unitname = "kJ/mol"; break;
	default:
		abort();
	}
	
	// print out information
	cout << "\033[1;33mStates:\033[0m";
	cout << string(15, ' ' );
	for (int i=0; i<nstates; i++) {
		printf("%12d ", (i+1));
	}
	printf("\n");
	cout << "\033[1;33m# of Observations:\033[0m";
	cout << string(4, ' ' );	
	for (int i=0; i<nstates; i++) {
		printf("%12d ", nobs[i]);
	}
	printf("\n");
	cout << "\033[1;33mUnit of the input data:\033[0m";
	printf("\033[1;31m%10s\033[0m\n\n", unitname);
}


/*****************************************************************************************************/
int main(int argc, char* argv[]) {
	// cout << "RepExchange.o temper nattempts totalcycle" << "\n";
	time_t start,end;
	time (&start);

	double dummyx;
	int dummym, dummyn;
	string dummys;
	char* potfile;
	ifstream read_file;
	string line;
	int niteration = RAND_MAX;
	int literation = 1;
	int nequilibrium = 1;
	int startneq;
	int method;
	double tolerance = 1.0e-25;
	double gamma = 1.0;
	int option_char;
	int gromacs = -1;
	const char* printlist="";
	vector<int> printitems;
	int normalization=0;	
	double temper = 300.0;
	int succtag;
	vector<double> eq_density_of_states; 
	vector<double> eq_partition_function;	
	int inunit = 0;
	int outunit = 0;
	double unit_factor_in, unit_factor_out;
	
	opterr = 0;	
	while ((option_char = getopt(argc, argv, "hzg:w:d:f:m:r:l:q:n:e:t:i:o:")) != -1) {
		switch (option_char) {
		case 'h':
			printf("%s -d potential_energy_input_file [-g gromacs_sim_type] -t temperature  -i input_unit -o output_unit -m iteration_method -q equilibrium_length -n iteration_number -e error_tolerance -w weights_print_list [-z]\n", argv[0]);
			return 0;
		case 'd': potfile = optarg; break;
		case 'g': gromacs = atoi(optarg); break;
		case 't': temper = atof(optarg); break;
		case 'm': method = atoi(optarg); break;
		case 'r': gamma = atof(optarg); break;
		case 'l': literation = atoi(optarg); break;
		case 'q': nequilibrium = atoi(optarg); break;			
		case 'n': niteration = atoi(optarg); break;
		case 'e': tolerance = atof(optarg); break;
		case 'w': printlist = optarg; break;
		case 'i': inunit = atoi(optarg); break;
		case 'o': outunit = atoi(optarg); break;
		case 'z': normalization = 1; break;
		case '?': 
			if ((optopt == 'd') || (optopt == 'm') || (optopt == 'r') || (optopt == 'l') || (optopt == 'q') || (optopt == 'n') || (optopt == 'e') || (optopt == 't') || (optopt == 'w') || (optopt == 'i') || (optopt == 'o') || (optopt == 'g'))
				fprintf (stderr, "Option -%c requires an argument.\n", optopt);
			else if (isprint (optopt))
				fprintf (stderr, "Unknown option `-%c'.\n", optopt);
			else
				fprintf (stderr, "Unknown option character `\\x%x'.\n", optopt);
			return 1;
		default:
			abort ();
		}
	}

	// print out the command line:
	cout << "#\033[1;33m#Command Line: \033[0m";	
	char* command;
	for(int i = 0; i < argc; i++) {
		printf("%s ", argv[i]);		
	}
	printf("\n");	

	// read the print list for weights
	if (printlist[0] != '\0') {
		readlist(printlist, printitems);
	}
	// read the bias
	read_matrix(potfile, bias_matrix, nobs);
	nstates = nobs.size();
	for (int i=0; i<nstates; i++) {
		totalndata += nobs[i];
	}
	if (gromacs >= 0) {
		convert_bias_energy_matrix_from_gromacs(gromacs);
		inunit = 2;
	}
	shift_bias_energy_matrix();	
	// get unit change constant
	unit_factor_in = others_2_kBT(inunit, temper);
	unit_factor_out = kBT_2_others(outunit, temper);
	// change the unit of the bias_energy_matrix 
	if (inunit != 0) {
		change_unit_of_bias_energy_matrix(unit_factor_in);		
	}
    convert_energy_to_bias_matrix();	
	
    // initialize density of states
	for (int i=0; i<totalndata; i++) {
		density_of_states.push_back(1.0/totalndata);
	}
	// initialize partition functions 
	for (int i=0; i<nstates; i++) {
		partition_function.push_back(1.0);			
	}
	// initial print out
	initial_print_out(inunit);
	
	startneq = 1;
	// run self_consistent_iteration once
equilibrium:
	cout << "\033[1;33mEquilibration:\033[0m\n";
	self_consistent_iteration(startneq, nequilibrium, tolerance);
	eq_density_of_states = density_of_states;
	eq_partition_function = partition_function;
	
	// iteration
	cout << "\033[1;33mIterations:\033[0m\n";	
	switch (method) {
	case 0: succtag = self_consistent_iteration(1, niteration, tolerance); break;
	case 1: succtag = convex_newton_raphson_solver(niteration, tolerance, gamma, literation); break;						
	}

	// if failed to converge
	if (succtag == 0) {
		startneq = nequilibrium + 1;
		if (nequilibrium == 0) {
			nequilibrium++;
		}
		else {
			nequilibrium *= 2;
		}
		if (nequilibrium < maxneq) {
			density_of_states = eq_density_of_states;
			partition_function = eq_partition_function;
			goto equilibrium;
		}
		else {
			printf("\033[1;31mFailed to converge!\033[0m\n");
		}
	}
	// final print out and weights
	if (succtag == 1) {
		final_print_out(outunit, unit_factor_out);
		if (printlist[0] != '\0') {
			print_out_weights(printitems, normalization);
		}
	}
	
	time (&end);
	double dif = difftime (end,start);
	printf ("# Elasped time is %.2lf seconds.\n", dif );
	
	return 0;
}
