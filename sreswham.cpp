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

using std::vector;
using std::cout;
using std::getline;
using std::string;
using std::ifstream;
using std::istringstream;
using std::size_t;

const double kB = 0.0019872041; // unit kcal/mol
const double pi = 3.1415926;
const double cal2J = 4.184;

int nstates = 0; // how many walks or thermodynamic states
int npara = 0; // how many parameters
int ndata = 0; // how mnay data for one observation
long totalndata = 0; // total number of observations
vector<int> nobs; // how many observations at each state
int nmacro; // number of macrostate clusters
vector<vector<int> > nobsmacro; // how many observations in each macrostates

vector<vector<double> > paralist;
vector<vector<vector<vector<double> > > > datalist;
int pftype; // potential function type
int inunit; // input data unit
int outunit; // output data unit

double potential_energy(int datasid, int datamid, int datanid, int stateid) {
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
	case 0: potential = datalist[datasid][datamid][datanid][stateid]*unitfactor; break;  // uwham type input
	case 1:	potential = datalist[datasid][datamid][datanid][0]*unitfactor; break;  // temperature RE
	case 2:	potential = paralist[stateid][1]*datalist[datasid][datamid][datanid][0]*unitfactor; break; // Hamiltonian RE		
	case 3:	potential = (datalist[datasid][datamid][datanid][0] + paralist[stateid][1]*datalist[datasid][datamid][datanid][1])*unitfactor; break; // T and H RE
	}
	return potential;
}

void read_datalist(const char* filename) {
	ifstream read_file;
	string line;
	double dummyx;
	char oldline0;
	vector<vector<vector<double> > > element;
	int imacro;	
	
	read_file.open(filename);	
	assert(read_file.is_open());
			
	while (getline(read_file, line)) {
		if ((line[0] == '#') or (line[0] == '@')) {
			if (!element.empty()) {
				datalist.push_back(element);
				int sumnobs = 0;
				vector<int>  tmpnobs;
				for (int i=0; i<nmacro; i++) {
					tmpnobs.push_back(element[i].size());
					sumnobs += element[i].size();
				}
				nobs.push_back(sumnobs);
				nobsmacro.push_back(tmpnobs);
				element.clear();
			}
			for (int i=0; i<nmacro; i++) {
				vector<vector<double> > emptyarray;
				element.push_back(emptyarray);
			}
		}
		else {
			istringstream vstring(line);
			vector<double>  tmparray;
			while (vstring >> dummyx) {
				tmparray.push_back(dummyx);
			}
            // truncate toward zero, index of macrostate clusters starts from one, not zero 			
			imacro = (int)(tmparray.back()+0.5) - 1; 
			element[imacro].push_back(tmparray);			
		}
	}
	if (!element.empty()) {
		datalist.push_back(element);
		int sumnobs = 0;
		vector<int>  tmpnobs;
		for (int i=0; i<nmacro; i++) {
			tmpnobs.push_back(element[i].size());
			sumnobs += element[i].size();
		}
		nobs.push_back(sumnobs);
		nobsmacro.push_back(tmpnobs);		
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
	int state_id, macro_id, datum_id;
	vector<double> datum;

public:
	void init_walker(int new_state_id);
	void change_id(int new_state_id, int new_macro_id, int new_datum_id);
	int get_state_id();
	int get_macro_id();	
	int get_datum_id();
	void walking();
};

int RandomWalker::get_state_id() {
	return state_id;
}

int RandomWalker::get_macro_id() {
	return macro_id;
}

int RandomWalker::get_datum_id() {
	return datum_id;
}

void RandomWalker::change_id(int new_state_id, int new_macro_id, int new_datum_id) {
	state_id = new_state_id;
	macro_id = new_macro_id;  // this is not necessary
	datum_id = new_datum_id;
}

void RandomWalker::init_walker(int new_state_id) {
	double randnum;

	state_id = new_state_id;
	randnum = (double)rand()/(double)(RAND_MAX + 1.0);
	macro_id = floor(randnum*datalist[state_id].size());
	assert(macro_id < datalist[state_id].size());	
	randnum = (double)rand()/(double)(RAND_MAX + 1.0);
	datum_id = floor(randnum*datalist[state_id][macro_id].size());
	assert(datum_id < datalist[state_id][macro_id].size());
}

void RandomWalker::walking() {
	double randnum;
	long tmpid;

	randnum = (double)rand()/(double)(RAND_MAX + 1.0);	
	if ((int)(paralist[state_id].back() + 0.5) == 1) { // choose from the whole state
		tmpid = floor(randnum*nobs[state_id]);
		for (int i=0; i<nmacro; i++) {
			if ((tmpid+1) > nobsmacro[state_id][i]) {
				tmpid -= nobsmacro[state_id][i];
			}
			else {
				macro_id = i;
				datum_id = tmpid;
				break;
			}
		}
	}
	else { // choose from the same macrostate
		datum_id = floor(randnum*nobsmacro[state_id][macro_id]);
	}
	assert(datum_id < datalist[state_id][macro_id].size());	
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
	int m_state_id, m_macro_id, m_datum_id, n_state_id, n_macro_id, n_datum_id;
	double mdatum, ndatum;
	
	m_state_id = walker[mwalker].get_state_id();
	m_macro_id = walker[mwalker].get_macro_id();	
	m_datum_id = walker[mwalker].get_datum_id();
	n_state_id = walker[nwalker].get_state_id();
	n_macro_id = walker[nwalker].get_macro_id();		
	n_datum_id = walker[nwalker].get_datum_id();

	// copy
	vector<double> mtmpelem;
	vector<double> ntmpelem;		
	for (int i=0; i<ndata; i++) {
		mtmpelem.push_back(datalist[m_state_id][m_macro_id][m_datum_id][i]);
		ntmpelem.push_back(datalist[n_state_id][n_macro_id][n_datum_id][i]);		
	}
	// move the last
	for (int i=0; i<ndata; i++) {
		datalist[m_state_id][m_macro_id][m_datum_id][i] = datalist[m_state_id][m_macro_id][nobsmacro[m_state_id][m_macro_id]-1][i];
		datalist[n_state_id][n_macro_id][n_datum_id][i] = datalist[n_state_id][n_macro_id][nobsmacro[n_state_id][n_macro_id]-1][i];	
	}
	// delete the last
	nobsmacro[m_state_id][m_macro_id]--;
	datalist[m_state_id][m_macro_id].resize(nobsmacro[m_state_id][m_macro_id]);
	nobsmacro[n_state_id][n_macro_id]--;	
	datalist[n_state_id][n_macro_id].resize(nobsmacro[n_state_id][n_macro_id]);
	// append	
	datalist[n_state_id][m_macro_id].push_back(mtmpelem);
	nobsmacro[n_state_id][m_macro_id]++;
	datalist[m_state_id][n_macro_id].push_back(ntmpelem);
	nobsmacro[m_state_id][n_macro_id]++;
}

void RepExSystem::swap(int m, int n) {
	int new_m_state_id, new_m_macro_id, new_m_datum_id, new_n_state_id, new_n_macro_id, new_n_datum_id;

	new_m_state_id = walker[n].get_state_id();
	new_m_macro_id = walker[m].get_macro_id();
	new_m_datum_id = datalist[new_m_state_id][new_m_macro_id].size()-1;
	new_n_state_id = walker[m].get_state_id();
	new_n_macro_id = walker[n].get_macro_id();
	new_n_datum_id = datalist[new_n_state_id][new_n_macro_id].size()-1;

	walker[m].change_id(new_m_state_id, new_m_macro_id, new_m_datum_id);
	walker[n].change_id(new_n_state_id, new_n_macro_id, new_n_datum_id);
}


void RepExSystem::exchange_attempts(int nattempts) {
	double factor;
	int m, n, swappair, m_state_id, m_macro_id, m_datum_id, n_state_id, n_macro_id, n_datum_id;
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
		m_macro_id = walker[m].get_macro_id();
		m_datum_id = walker[m].get_datum_id();
		n_state_id = walker[n].get_state_id();
		n_macro_id = walker[n].get_macro_id();
		n_datum_id = walker[n].get_datum_id();

		factor = exp(-(potential_energy(n_state_id, n_macro_id, n_datum_id, m_state_id) +
					   potential_energy(m_state_id, m_macro_id, m_datum_id, n_state_id) -
					   potential_energy(m_state_id, m_macro_id, m_datum_id, m_state_id) -
					   potential_energy(n_state_id, n_macro_id, n_datum_id, n_state_id)));
		
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
		int i_macro_id = walker[walker_id[printstates[i]-1]].get_macro_id();
		int i_datum_id = walker[walker_id[printstates[i]-1]].get_datum_id();
		if (printitems.empty()) {
			// add .5 to positive numbers and subtract .5 from negative numbers before converting
			index = (long)(datalist[printstates[i]-1][i_macro_id][i_datum_id][ndata-1] + 0.5);
			printf("%12d\t", index);
		}
		else {
			for (int j=0; j<printitems.size(); j++) {
				if (printitems[j] == (ndata-2)) {
					printf("%12d\t", (int)(datalist[printstates[i]-1][i_macro_id][i_datum_id][printitems[j]-1] + 0.5));
				}
				else {
					printf("%12g\t", datalist[printstates[i]-1][i_macro_id][i_datum_id][printitems[j]-1]);
				}
			}				
		}
	}
	printf("\n");
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

	double temper = -1.0; 
	long totalcycle = 1;
	int nattempts = 1;
	int equilibrium = 0;
	
	RepExSystem myresim;

	int option_char;

	int multiplyK = 0;
	int multiplyM = 0;
	int multiply_factor = 1;		
	
	srand (time(NULL));
	
	opterr = 0;	
	while ((option_char = getopt(argc, argv, "hkmf:d:a:u:t:i:o:q:n:x:s:p:")) != -1) {
		switch (option_char) {
		case 'h':
			printf("%s -d observation_data_file -f state_information_file -u potential_function_type \
-t temperature  -i input_unit -o output_unit -q equilibrium_length -n number_of_cycles \
-x number_of_exchange_attempts -s print_list_of_states -p print_list_of_properties [-k -m]\n", argv[0]);
			return 0;
		case 'd': datafile = optarg; break;			
		case 'f': parafile = optarg; break;
		case 'a': nmacro = atoi(optarg); break;			
		case 'u': pftype = atoi(optarg); break;
		case 't': temper = atof(optarg); break;
		case 'i': inunit = atoi(optarg); break;
		case 'o': outunit = atoi(optarg); break;			
		case 'q': equilibrium = atoi(optarg); break;			
		case 'n': totalcycle = atol(optarg) ; break;
		case 'x': nattempts = atoi(optarg); break;
		case 's': pstatelist = optarg; break;			
		case 'p': printlist = optarg; break;
		case 'k': multiplyK = 1; break;
		case 'm': multiplyM = 1; break;
		case '?': 
			if ((optopt == 'd') || (optopt == 'f') || (optopt == 'u') || (optopt == 't') || (optopt == 'i') \
				|| (optopt == 'o') || (optopt == 'e') || (optopt == 'n') || (optopt == 'x') || (optopt == 's') \
				|| (optopt == 'p') || (optopt == 'q') || (optopt == 'a')) 
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
		for (int j=0; j<nmacro; j++) {
			for (int k=0; k<nobsmacro[i][j]; k++) {			
				datalist[i][j][k].push_back(dummym);
				dummym++;
			}
		}
	}
	ndata = datalist[0][0][0].size();
	
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

	// begin analysis

	myresim.init_RepExSystem();

	long icycle = 0;
	while(icycle < (totalcycle + equilibrium)) {
		myresim.MDupdate();
		myresim.exchange_attempts(nattempts);
		if (icycle > equilibrium) {
			myresim.printout(printstates, printitems);
		}
		icycle++;		
	}

	// time
	time (&end);
	double dif = difftime (end,start);
	printf ("# Elasped time is %.2lf seconds.\n", dif );
	
	return 0;
}
