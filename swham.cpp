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
#include <algorithm> // find

using std::vector;
using std::cout;
using std::getline;
using std::string;
using std::ifstream;
using std::istringstream;
using std::size_t;
using std::find;

const double kB = 0.0019872041; // unit kcal/mol
const double pi = 3.1415926;
const double cal2J = 4.184;

int nstates = 0; // how many walks or thermodynamic states
int npara = 0; // how many parameters
int ndata = 0; // how mnay data for one observation
long totalndata = 0; // total number of observations
vector<int> nobs; // how many observations at each state

vector<vector<double> > paralist; // state information
vector<vector<vector<double> > > datalist; // data points
vector<vector<int> > nearest_neighbors; // nearest neighbors of states

int pftype; // potential function type
int inunit; // input data unit
int outunit; // output data unit

double potential_energy(int datasid, int datanid, int stateid) {
	double potential;
	double unitfactor;
	
	switch (inunit) {
	case 0: if (pftype==0) {
			unitfactor = 1.0;
		}
		else {
			unitfactor = paralist[datasid][0]/paralist[stateid][0]; 
		}
		break;		
	case 1: unitfactor = 1.0/(kB*paralist[stateid][0]); break;
	case 2: unitfactor = 1.0/(cal2J*kB*paralist[stateid][0]); break;
	}
	
	switch (pftype) {
	case 0: potential = datalist[datasid][datanid][stateid]*unitfactor; break;
	case 1:	potential = datalist[datasid][datanid][0]*unitfactor; break;
	case 2:	potential = paralist[stateid][1]*datalist[datasid][datanid][0]*unitfactor; break;		
	case 3:	potential = (datalist[datasid][datanid][0] + paralist[stateid][1]*datalist[datasid][datanid][1])*unitfactor; break;
	case 4: {
		double deltax = fabs(datalist[datasid][datanid][0] - paralist[stateid][2]);
		if ((npara = 4) and (deltax > (0.5*paralist[stateid][3]))) {
			deltax = fabs(paralist[stateid][3] - deltax);
		}
		potential = 0.5*paralist[stateid][1]*deltax*deltax*unitfactor;
		break;}
	case 5: {
		double deltax = fabs(datalist[datasid][datanid][0] - paralist[stateid][3]);
		if ((npara = 6) and (deltax > (0.5*paralist[stateid][5]))) {
			deltax = fabs(paralist[stateid][5] - deltax);
		}
		double deltay = fabs(datalist[datasid][datanid][1] - paralist[stateid][4]);
		if ((npara = 7) and (deltay > (0.5*paralist[stateid][6]))) {
			deltay = fabs(paralist[stateid][6] - deltay);
		}	
		potential = 0.5*paralist[stateid][1]*deltax*deltax*unitfactor + 0.5*paralist[stateid][2]*deltay*deltay*unitfactor;
		break;}		
	}
	return potential;
}


void find_nearest_neighbors() {
	void oneD_nearest_neighbors(int para_index, int period_index);
	void twoD_nearest_neighbors(int para_index1, int para_index2, int period_index1, int period_index2);
	
	switch (pftype) {
	case 0: { // assign neighbors by the order of states, always treat as periodic
		vector<int> tmpneighbors;
		tmpneighbors.push_back(1);
		tmpneighbors.push_back(nstates-1);		
		nearest_neighbors.push_back(tmpneighbors);
		tmpneighbors.clear();
		for (int j=1; j<(nstates-1); j++) {
			tmpneighbors.push_back(j-1);			
			tmpneighbors.push_back(j+1);
			nearest_neighbors.push_back(tmpneighbors);
			tmpneighbors.clear();			
		}
		tmpneighbors.push_back(0);		
		tmpneighbors.push_back(nstates-2);
		nearest_neighbors.push_back(tmpneighbors);
		tmpneighbors.clear();		
		break;
	}
	case 1: oneD_nearest_neighbors(0, -1); break; // Temperature RE
	case 2: oneD_nearest_neighbors(1, -1); break; // Hamiltonian RE
	case 3: twoD_nearest_neighbors(0, 1, -1, -1); break; // 2D RE
	case 4: oneD_nearest_neighbors(2, 3); break; // 1D umbrella sampling, possible periodic
	case 5: twoD_nearest_neighbors(3, 4, 5, 6); break; // 2D umbrella sampling, possible periodic
	}
}

void oneD_nearest_neighbors(int para_index, int period_index) {
	double neg_diff, pos_diff, diff;
	int neg_stateid, pos_stateid;

	for (int i=0; i<nstates; i++) {
		vector<int> tmpneighbors;
		neg_diff = -1.0e100;
		pos_diff = 1.0e100;
		neg_stateid = -1;
		pos_stateid = -1;
		vector<double> testvalues;
		testvalues.push_back(paralist[i][para_index]);	   
		if ((period_index > 0) && (paralist.size() > period_index)) {
			testvalues.push_back(paralist[i][para_index]+paralist[i][period_index]);
			testvalues.push_back(paralist[i][para_index]-paralist[i][period_index]);			
		}
		for (int k=0; k<testvalues.size(); k++) {
			for (int j=0; j<nstates; j++) {
				if (i != j) {
					diff = paralist[j][para_index] - testvalues[k];
					if ( (diff < 0.0) && (diff > neg_diff)) {
						neg_diff = diff;
						neg_stateid = j;
					}
					if ( (diff > 0.0) && (diff < pos_diff)) {
						pos_diff = diff;
						pos_stateid = j;
					}				
				}
			}
		}
		if (neg_stateid >= 0) {
			tmpneighbors.push_back(neg_stateid);
		}
		if (pos_stateid >= 0) {
			tmpneighbors.push_back(pos_stateid);
		}
		nearest_neighbors.push_back(tmpneighbors);
	}
}

void twoD_nearest_neighbors(int para_index1, int para_index2, int period_index1, int period_index2) {
	double first_neg_diff, first_pos_diff, second_neg_diff, second_pos_diff, first_diff, second_diff;
	int first_neg_stateid, first_pos_stateid, second_neg_stateid, second_pos_stateid;

	for (int i=0; i<nstates; i++) {
		// prepare
		vector<int> tmpneighbors;
		first_neg_diff = -1.0e100;
		second_neg_diff = -1.0e100;		
		first_pos_diff = 1.0e100;
		second_pos_diff = 1.0e100;		
		first_neg_stateid = -1;
		second_neg_stateid = -1;		
		first_pos_stateid = -1;
		second_pos_stateid = -1;
		vector<double> testvalues1;
		testvalues1.push_back(paralist[i][para_index1]);
		if ((period_index1 > 0) && (paralist.size() > period_index1)) {
			testvalues1.push_back(paralist[i][para_index1]+paralist[i][period_index1]);
			testvalues1.push_back(paralist[i][para_index1]-paralist[i][period_index1]);			
		}
		vector<double> testvalues2;
		testvalues2.push_back(paralist[i][para_index2]);
		if ((period_index2 > 0) && (paralist.size() > period_index2)) {
			testvalues2.push_back(paralist[i][para_index2]+paralist[i][period_index2]);
			testvalues2.push_back(paralist[i][para_index2]-paralist[i][period_index2]);			
		}
		// search states
		for (int k1=0; k1<testvalues1.size(); k1++) {
			for (int k2=0; k2<testvalues2.size(); k2++) {			
				for (int j=0; j<nstates; j++) {
					if (i != j) {
						first_diff = paralist[j][para_index1] - testvalues1[k1];
						second_diff = paralist[j][para_index2] - testvalues2[k2];
						if ((second_diff == 0.0) && (first_diff < 0.0) && (first_diff > first_neg_diff)) {
							first_neg_diff = first_diff;
							first_neg_stateid = j;
						}
						if ((second_diff == 0.0) && (first_diff > 0.0) && (first_diff < first_pos_diff)) {
							first_pos_diff = first_diff;
							first_pos_stateid = j;
						}
						if ((first_diff == 0.0) && (second_diff < 0.0) && (second_diff > second_neg_diff)) {
							second_neg_diff = second_diff;
							second_neg_stateid = j;
						}
						if ((first_diff == 0.0) && (second_diff > 0.0) && (second_diff < second_pos_diff)) {
							second_pos_diff = second_diff;
							second_pos_stateid = j;
						}				
					}
				}
			}
		}
		// add to array
		if (first_neg_stateid >= 0) {
			tmpneighbors.push_back(first_neg_stateid);
		}
		if (first_pos_stateid >= 0) {
			tmpneighbors.push_back(first_pos_stateid);
		}
		if (second_neg_stateid >= 0) {
			tmpneighbors.push_back(second_neg_stateid);
		}
		if (second_pos_stateid >= 0) {
			tmpneighbors.push_back(second_pos_stateid);
		}
		nearest_neighbors.push_back(tmpneighbors);
	}
}

void read_datalist(const char* filename) {
	ifstream read_file;
	string line;
	double dummyx;
	char oldline0;
	vector<vector<double> > element;
	
	read_file.open(filename);	
	assert(read_file.is_open());
			
	while (getline(read_file, line)) {
		if ((line[0] == '#') or (line[0] == '@')) {
			if (!element.empty()) {
				datalist.push_back(element);
				nobs.push_back(element.size());
				element.clear();
			}
		}
		else {
			vector<double>  tmparray;		
			istringstream vstring(line);
			while (vstring >> dummyx) {
				tmparray.push_back(dummyx);
			}
			element.push_back(tmparray);			
		}
	}
	if (!element.empty()) {
		datalist.push_back(element);
		nobs.push_back(element.size());
	}
	read_file.close();
}

void read_paralist(const char* filename, vector<int>& printstates) {
	ifstream read_file;
	string line;
	double dummyx;
	int dummyn;
	vector<int> printtag;
	int sindex=0;

	read_file.open(filename);	
	assert(read_file.is_open());
	if (paralist.empty()) {
		for (int i=0; i<nstates; i++) {
			vector<double>  paraarray;
			paralist.push_back(paraarray);
		}
	}

	while (getline(read_file, line)) {
		vector<double>  paraarray;		
		istringstream vstring(line);
		while (vstring >> dummyx) {
			paralist[sindex].push_back(dummyx);
		}
		sindex++;
	}
	read_file.close();
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



class RandomWalker {
private:
	int state_id, datum_id;
	vector<double> datum;

public:
	void init_walker(int new_state_id);
	void change_id(int new_state_id, int new_datum_id);
	int get_state_id();
	int get_datum_id();
	void walking();
};

int RandomWalker::get_state_id() {
	return state_id;
}

int RandomWalker::get_datum_id() {
	return datum_id;
}

void RandomWalker::change_id(int new_state_id, int new_datum_id) {
	state_id = new_state_id;
	datum_id = new_datum_id;
}

void RandomWalker::init_walker(int new_state_id) {
	double randnum;

	state_id = new_state_id;
	randnum = (double)rand()/(double)(RAND_MAX + 1.0);
	datum_id = floor(randnum*datalist[state_id].size());
	assert(datum_id < datalist[state_id].size());
}

void RandomWalker::walking() {
	double randnum;

	randnum = (double)rand()/(double)(RAND_MAX + 1.0);
	datum_id = floor(randnum*datalist[state_id].size());
	assert(datum_id < datalist[state_id].size());	
}

class RepExSystem {
private:
	vector<RandomWalker> walker;
public:
	void init_RepExSystem();
	void printout(vector<int>& printstates, vector<int>& printitems);
	void update_database(int mwalker, int nwalker);
	void swap(int m, int n);
	void exchange_attempts(int nattempts);
	void MDupdate();
};

void RepExSystem::init_RepExSystem() {
	for (int i=0; i<nstates; i++) {
		RandomWalker tmpwalker;
		tmpwalker.init_walker(i);
		walker.push_back(tmpwalker);
	}
}

void RepExSystem::MDupdate() {
	for (int j=0; j<nstates; j++) {
		walker[j].walking();
	}
}

void RepExSystem::update_database(int mwalker, int nwalker) {
	int m_state_id, m_datum_id, n_state_id, n_datum_id;
	double mdatum, ndatum;
	
	m_state_id = walker[mwalker].get_state_id();
	m_datum_id = walker[mwalker].get_datum_id();
	n_state_id = walker[nwalker].get_state_id();
	n_datum_id = walker[nwalker].get_datum_id();

	for (int i=0; i<ndata; i++) {
		mdatum = datalist[m_state_id][m_datum_id][i];
		ndatum = datalist[n_state_id][n_datum_id][i];
		datalist[m_state_id][m_datum_id][i] = ndatum;
		datalist[n_state_id][n_datum_id][i] = mdatum;
	}
}

void RepExSystem::swap(int m, int n) {
	int new_m_state_id, new_m_datum_id, new_n_state_id, new_n_datum_id;

	new_m_state_id = walker[n].get_state_id();
	new_m_datum_id = walker[n].get_datum_id();
	new_n_state_id = walker[m].get_state_id();
	new_n_datum_id = walker[m].get_datum_id();

	walker[m].change_id(new_m_state_id, new_m_datum_id);
	walker[n].change_id(new_n_state_id, new_n_datum_id);
}


void RepExSystem::exchange_attempts(int nattempts) {
	double factor;
	int m, n, swappair, m_state_id, m_datum_id, n_state_id, n_datum_id;
	double randnum;
	
	if (nattempts < 0) {
		nattempts *= -nstates; 
	}
	for (int i=0; i<nattempts; i++) {
		randnum = (double)rand()/(double)(RAND_MAX + 1.0);
		m = floor(randnum*nstates);
		assert(m<nstates);
		randnum = (double)rand()/(double)(RAND_MAX + 1.0);
		n = floor(randnum*nstates);
		assert(n<nstates);
		while (n == m) {
			randnum = (double)rand()/(double)(RAND_MAX + 1.0);
			n = floor(randnum*nstates);
			assert(n<nstates);
		}
		swappair = 1;

		m_state_id = walker[m].get_state_id();
		m_datum_id = walker[m].get_datum_id();
		n_state_id = walker[n].get_state_id();
		n_datum_id = walker[n].get_datum_id();


		factor = exp(-(potential_energy(n_state_id, n_datum_id, m_state_id) +
					   potential_energy(m_state_id, m_datum_id, n_state_id) -
					   potential_energy(m_state_id, m_datum_id, m_state_id) -
					   potential_energy(n_state_id, n_datum_id, n_state_id)));			
		
		if (factor < 1) {
			randnum = (double)rand()/(double)(RAND_MAX + 1.0);
			if (randnum > factor) {
				swappair = 0;
			}
		}
		if (swappair == 1) {
			update_database(m, n);
			swap(m, n);
		}
	}
	
}
		

void RepExSystem::printout(vector<int>& printstates, vector<int>& printitems) {
	int walker_id[nstates];
	long index;
	
	for (int i=0; i<nstates; i++) {
		int i_state_id = walker[i].get_state_id();
		walker_id[i_state_id] = i;
	}

	for (int i=0; i<printstates.size(); i++) {
		int i_datum_id = walker[walker_id[printstates[i]-1]].get_datum_id();
		if (printitems.empty()) {
			// add .5 to positive numbers and subtract .5 from negative numbers before converting
			index = (long)(datalist[printstates[i]-1][i_datum_id][ndata-1] + 0.5);
			printf("%12d\t", index);
		}
		else {
			for (int j=0; j<printitems.size(); j++) {
				printf("%12g\t", datalist[printstates[i]-1][i_datum_id][printitems[j]-1]);
			}				
		}
	}
	printf("\n");
}

class SerTempSystem {
private:
	RandomWalker walker;
	vector<double> freeE;
	vector<double> sum_of_freeE;	
	vector<double> pi0;
	vector<long> pi1;
public:
	void init_SerTempSystem();
	void update_freeE(int state_id, int burnin, long icycle, double alpha, double beat);
	int global_jump();
	int local_jump(int njumps); 
	void MDupdate();
	void printfreeE(vector<int>& printstates, int nfreeE);
	double standard_deviation();
	void cal_DoS();
};

void SerTempSystem::init_SerTempSystem() {
	walker.init_walker(0);
	for (int i=0; i<nstates; i++) {
		freeE.push_back(0.0);
		sum_of_freeE.push_back(0.0);		
		pi0.push_back(nobs[i]*1.0/totalndata);
		pi1.push_back(1);
	}
	// pi1[0]++;
}

void SerTempSystem::cal_DoS() {
	FILE * write_file;	
	double dos; // density of states
	vector<double> pf; // partition function
	double basepe, biaspe; 

	// calculate the partition functions
	for (int i=0; i<nstates; i++) {
		double tmppf = exp(-(freeE[i] - freeE[0]));		
		pf.push_back(tmppf);
	}
	write_file = fopen("weights.data", "w");			
	// calculate density of states
	for (int i=0; i<nstates; i++) {
		for (int j=0; j<nobs[i]; j++) {
			dos = 0.0;
			basepe = potential_energy(i, j, i);
			for (int k=0; k<nstates; k++) {
				biaspe = potential_energy(i, j, k);
				dos += nobs[k]*exp(-(biaspe - basepe))/pf[k];
			}
			dos = 1.0/dos;
			fprintf(write_file, "%12g\n", dos);
		}
	}
	fclose(write_file);	
}

void SerTempSystem::MDupdate() {
	walker.walking();
}

void SerTempSystem::update_freeE(int state_id, int burnin, long icycle, double alpha, double beta) {
	double gamma;
	double base;
	double change;

	pi1[state_id]++;	
	if (icycle <= burnin) {
		gamma = pow(icycle, -alpha);
		change = gamma;
	}
	else {
		gamma = 1.0/(icycle - burnin + pow(burnin, alpha/beta));
		change = pow(gamma, beta);
		// change = gamma*sumpi1/pi1[state_id]; //from paper, but not working		
	}
	if (change > 1.0) { //make sure change is smaller than 1 kBT
		change = 1.0;
	}
	freeE[state_id] -= change;
	if (state_id == 0) {
		base = freeE[0];
		for (int i=0; i<nstates; i++) {
			freeE[i] -= base;
		}
	}
	walker.change_id(state_id, 0);
	if (icycle > burnin) {	
		for (int i=0; i<nstates; i++) {
			sum_of_freeE[i] += freeE[i];
		}
	}
}


int SerTempSystem::global_jump() {
	double tmpprob;
	int state_id, datum_id, new_state_id;
	vector<double> prob;
	double randnum;
	double sum=0.0;

	state_id = walker.get_state_id();
	datum_id = walker.get_datum_id();
	for (int i=0; i<nstates; i++) {
		// tmpprob = pi0[i]*exp(freeE[i])*exp(-potential_energy(state_id, datum_id, i));
		tmpprob = pi0[i]*exp(freeE[i]-potential_energy(state_id, datum_id, i));			
		sum += tmpprob;
		prob.push_back(tmpprob);
	}
	randnum = (double)rand()/(double)(RAND_MAX + 1.0)*sum;
	sum = 0.0;
	for (int i=0; i<nstates; i++) {
		sum += prob[i];
		if (sum > randnum) {
			new_state_id = i;
			break;
		}
	}
	return(new_state_id);
}

int SerTempSystem::local_jump(int njumps) {
	double factor;
	int state_id, datum_id, new_state_id, current_state_id;
	double randnum, prob_try_n2o, prob_try_o2n;
	int ijump, wjump;
	
	state_id = walker.get_state_id();
	datum_id = walker.get_datum_id();
	current_state_id = state_id;
	
	ijump = 0;
	while (ijump < njumps) {
		randnum = (double)rand()/(double)(RAND_MAX + 1.0);
		int tmpid = floor(randnum*nearest_neighbors[current_state_id].size());
		new_state_id = nearest_neighbors[current_state_id][tmpid];
		prob_try_o2n = 1.0/nearest_neighbors[current_state_id].size();
		if (find(nearest_neighbors[new_state_id].begin(), nearest_neighbors[new_state_id].end(), current_state_id)
			!= nearest_neighbors[new_state_id].end()) {
			prob_try_n2o = 1.0/nearest_neighbors[new_state_id].size();
		}
		else {
			prob_try_n2o = 0.0;
		}
		factor = pi0[new_state_id]/pi0[current_state_id]*exp((freeE[new_state_id]-potential_energy(state_id, datum_id, new_state_id))
				 - (freeE[current_state_id]-potential_energy(state_id, datum_id, current_state_id)))*prob_try_n2o/prob_try_o2n;
		// Metropolis
		wjump = 1;
		if (factor < 1) {
			randnum = (double)rand()/(double)(RAND_MAX + 1.0);
			if (randnum > factor) {
				wjump = 0;
			}
		}
		if (wjump == 1) {
			current_state_id = new_state_id;
		}
		ijump++;
	}
	// printf("%5d\n", current_state_id);
	return(current_state_id);
}

void SerTempSystem::printfreeE(vector<int>& printstates, int nfreeE) {
	double result;
	double unitfactor;
	double sdvalue;
	const char* unitname;	

	switch(outunit) {
	case 0: unitname = "kBT"; break;
	case 1: unitname = "kcal/mol"; break;
	case 2: unitname = "kJ/mol"; break;
	default:
		abort();
	}
	printf("\033[1;33mFree Energy (\033[0m\033[1;31m%s\033[0m\033[1;33m):\033[0m", unitname);			
	for (int i=0; i<printstates.size(); i++) {
		switch (outunit) {
		case 0: unitfactor = 1.0; break;
		case 1: unitfactor = kB*paralist[printstates[i]-1][0]; break;
		case 2: unitfactor = kB*paralist[printstates[i]-1][0]*cal2J; break;
		}
		result = sum_of_freeE[printstates[i]-1]/(nfreeE*1.0)*unitfactor;
		printf("%12g ", result);
		sum_of_freeE[i] = 0.0;
	}
	sdvalue = standard_deviation();
	for (int i=0; i<nstates; i++) {
		pi1[i] = 0;
	}
	printf("\033[1;33mSD:\033[0m");
	printf("%12g", sdvalue);
	printf("\n");		
}


double SerTempSystem::standard_deviation() {
	double result;
	double sum=0.0;
	double sum2=0.0;
	
	for (int i=0; i<nstates; i++) {
		sum += pi1[i];
	}
	for (int i=0; i<nstates; i++) {
		sum2 += pow((pi1[i]*1.0/sum - nobs[i]*1.0/totalndata), 2);
	}
	result = sqrt(sum2/(nstates*1.0));
	return(result);
}


/*****************************************************************************************************/
int main(int argc, char* argv[]) {
	time_t start,end;
	time (&start);

	const char* printlist="";
	const char* pstatelist="";
	vector<int> printitems;
	vector<int> printstates;
	const char* datafile="";
	const char* parafile="";
	vector<vector<int> > freeElist;

	double temper = -1.0; 
	long totalcycle = 1;
	int nattempts = 0;
	int equilibrium = 0;
	int printfreqG = 0;
	double alpha = 0.5;
	double beta = 0.75;
	
	RepExSystem myresim;
	SerTempSystem mystsim;	
	int option_char;

	int multiplyK = 0;
	int multiplyM = 0;
	int printwtag = 0;
	int multiply_factor = 1;		
	
	srand (time(NULL));
	
	opterr = 0;	
	while ((option_char = getopt(argc, argv, "hkmwf:d:u:t:i:o:q:n:x:s:p:g:a:b:")) != -1) {
		switch (option_char) {
		case 'h':
			printf("%s -d observation_data_file -f state_information_file -u potential_function_type \
-t temperature  -i input_unit -o output_unit -q equilibrium_length -n number_of_cycles \
-x number_of_attempts -s print_list_of_states -p print_list_of_properties -g free_energy_print_frequency [-k -m -w]\n", argv[0]);
			return 0;
		case 'd': datafile = optarg; break;			
		case 'f': parafile = optarg; break;
		case 'u': pftype = atoi(optarg); break;
		case 't': temper = atof(optarg); break;
		case 'i': inunit = atoi(optarg); break;
		case 'o': outunit = atoi(optarg); break;			
		case 'q': equilibrium = atoi(optarg); break;			
		case 'n': totalcycle = atol(optarg) ; break;
		case 'x': nattempts = atoi(optarg); break;
		case 's': pstatelist = optarg; break;			
		case 'p': printlist = optarg; break;
		case 'g': printfreqG = atoi(optarg); break;
		case 'a': alpha = atof(optarg); break;
		case 'b': beta = atof(optarg); break;
		case 'k': multiplyK = 1; break;
		case 'm': multiplyM = 1; break;
		case 'w': printwtag = 1; break;
		case '?': 
			if ((optopt == 'd') || (optopt == 'f') || (optopt == 'u') || (optopt == 't') || (optopt == 'i') \
				|| (optopt == 'o') || (optopt == 'e') || (optopt == 'n') || (optopt == 'x') || (optopt == 's') \
				|| (optopt == 'p') || (optopt == 'q') || (optopt == 'a') || (optopt == 'b')) 
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

	srand (time(NULL)); //random seed
	
	// print out the command line:
	cout << "#\033[1;33m#Command Line: \033[0m";	
	char* command;
	for(int i = 0; i < argc; i++) {
		printf("%s ", argv[i]);		
	}
	printf("\n");	
	
	// change numbers
	if (multiplyK == 1) {
		multiply_factor *= 1000;
	}
	if (multiplyM == 1) {
		multiply_factor *= 1000000;
	}
	if (multiply_factor > 1) {
		equilibrium *= multiply_factor;
		totalcycle *= multiply_factor;
		printfreqG *= multiply_factor;
	}

	// print state: information of which states will be printed out
	if (pstatelist[0] != '\0') {
		readlist(pstatelist, printstates);
	}
	
	// print list: which fields will be printed out for each data point
	if (printlist[0] != '\0') {
		readlist(printlist, printitems);
	}
	
	// read data list
	read_datalist(datafile);
	nstates = nobs.size();
	//// add index to each data point
	long dummym = 0;
	for (int i=0; i<nstates; i++) {
		for (int j=0; j<nobs[i]; j++) {
			datalist[i][j].push_back(dummym);
			dummym++;
		}
	}
	totalndata = dummym--;
	ndata = datalist[0][0].size(); 

	// state information
	// the temperature of each state is determined by variable temper
	if (temper > 0) {
		for (int i=0; i<nstates; i++) {
			vector<double>  paraarray;
			paraarray.push_back(temper);
			paralist.push_back(paraarray);
		}
	}
	// read state parameters		
	if (parafile[0] != '\0') {
		read_paralist(parafile, printstates);
		npara = paralist[0].size();
	}
	// decide nearest neighbors
	find_nearest_neighbors();
		
	// begin analysis
	if (printfreqG > 0) { // run Serial Tempering SWHAM
		if (equilibrium == 0) {
			equilibrium = printfreqG;
		}
		mystsim.init_SerTempSystem();

		long icycle = 0;
		while(icycle <= (totalcycle + equilibrium)) {
			int newid;
			mystsim.MDupdate();
			if (nattempts == 0) {
				newid = mystsim.global_jump();
			}
			else {
				newid = mystsim.local_jump(nattempts);				
			}
			mystsim.update_freeE(newid, equilibrium, icycle, alpha, beta);		
			if (icycle > equilibrium) {
				if ((icycle-equilibrium)%printfreqG==0) {
					mystsim.printfreeE(printstates, printfreqG);
				}
			}
			icycle++;		
		}
		if (printwtag == 1) {
			mystsim.cal_DoS();			
		}
	}
	else { // run Replica Exchange SWHAM
		myresim.init_RepExSystem();
		if (nattempts == 0) { // in case forget set -x
			nattempts = 1;
		}
		long icycle = 0;
		while(icycle < (totalcycle + equilibrium)) {
			myresim.MDupdate();
			myresim.exchange_attempts(nattempts);
			if (icycle > equilibrium) {
				myresim.printout(printstates, printitems);
			}
			icycle++;		
		}
	}
	
	// time
	time (&end);
	double dif = difftime (end,start);
	printf ("# Elasped time is %.2lf seconds.\n", dif );
	
	return 0;
}
