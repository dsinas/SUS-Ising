/********************************************************************************
 *
 * Usage: A program that computes the probability distribution of order parameter
 * of the 2D Ising model using Successive Umbrella Sampling.
 *
 * This program is made freely available with the understanding that every copy
 * of this file must include src files and headers and that it comes without any
 * WITHOUT ANY WARRANTY.
 *
 ********************************************************************************/

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <ctime>
#include <climits>
#include <iomanip>
#include <getopt.h>

#include "CLASS_ISING.h"

using namespace std;

// function to print a help message for the program to run
void print_usage() {
	cerr << "Purpose:\t Perform Successive Umbrella Sampling to measure orderparameter distribution of Ising2D model" << endl;
	cerr << "Usage:\t\t SUS [options]" << endl;
	cerr << "Available options (default in []):" << endl;
	cerr << "\t\t-h\t print this help message" << endl;
	cerr << "\t\t-L arg\t set size of the lattice L - default: [20]" << endl;
	cerr << "\t\t-J arg\t set Ising coupling J - default: [0.5)]" << endl;
}

/** MAIN FUNCTION ************************************************************/

int main(int argc, char* argv[]) {

    int Lx = 20, Ly = Lx; // the default lattice size

    double ln_f = 0.; // change for large systems at low temperature (large coupling J)

    double temperature = 1.; // set constant: Ising coupling changes
    double J_const = 0.5; // the default value for Ising coupling

    double sD = 0.5; // initial composiotion of the membrane: 50-50
    double sU = 1.-sD;

    srand((unsigned) time(NULL));
    seed((unsigned long) (rand()*102011));

    int file_counter = 0;
    
    char opt;
	while ((opt = getopt(argc, argv, "hL:J:")) != -1) {
		switch (opt) {
		case 'h':
			print_usage();
			return 0;
		case 'L':
			Lx = atoi(optarg);
			Ly = Lx;
			break;
		case 'J':
			J_const = atof(optarg);
			break;
		default:
			print_usage();
			return 1;
		}
	}
	
	cout << "\nSystem runs at" << endl;
	cout << "L = " << Lx << endl;
    cout << "J = " << J_const << "\n" << endl;
	cout << setfill('=') << setw(60) << "\n" << endl;
	cout << setfill(' ');
	
    IsingClass m_system(Lx, Ly);
    m_system.set_temperature(temperature);
    m_system.set_J_const(J_const);
    m_system.set_spin_composition(sD,sU);

    m_system.initialization();

    m_system.set_file_counter(file_counter);
    m_system.set_ln_f(ln_f);

    m_system.get_data_files();

}
