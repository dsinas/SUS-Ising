/*******************************************************************************
* 
* CLASS_ISING.h: A class for implementing Successive Umbrella Sampling
* to measure the probability distribution of the order parameter.
*
* This program is made freely available with the understanding that every 
* copy of this file must include this header and that it comes without any 
* WITHOUT ANY WARRANTY.
*******************************************************************************/

#ifndef CLASS_ISING_H
#define CLASS_ISING_H

#include "CLASS_MC.h"

using namespace std;

/** definition of IsingClass **********************************************/
class IsingClass : public MC_SystemClass {

protected:

    vector<double> ln_P;
    vector<int> op_vector;
    vector<double> ave_energy, ave_energy2;
    vector<double> ave_isingE, ave_isingE2;
    int n_bins;
    int w_win;
    int min_OP, max_OP;

    int f_counter;
    double ln_f;
    
    char filename[50];

public:

    IsingClass(const int _Lx, const int _Ly):
            MC_SystemClass(_Lx, _Ly)
    {
        n_bins = (int)(m_grid.get_volume()+1);
    }

    virtual ~IsingClass() {}

    void set_parameters() {

        ln_P.resize(n_bins);
        ave_energy.resize(n_bins);
        ave_energy2.resize(n_bins);
        ave_isingE.resize(n_bins);
        ave_isingE2.resize(n_bins);
        op_vector.resize(n_bins);
        min_OP = -m_grid.get_volume();
        max_OP = +m_grid.get_volume();
        w_win = 2;

        for (int bin=0;bin<n_bins;bin++)
            op_vector[bin] = min_OP + bin * w_win;

//         m_energy = get_energy();
//         m_orderParameter = get_order_parameter();

    }

    void get_data_files() {
        set_parameters();
        successiveUmbrellaSampling();

        print_prob_data();
        // print_ave_energy_data();
        // print_ave_isingE_data();
    }

    /** performing Successive Umbrella Sampling ******************************/
    void successiveUmbrellaSampling() {
        int Nmax = 20000; // # sweeps given experimentally [100k]
        double log_ratio = 0.;
        ln_P[0] = 0.;
        ave_energy[0] = m_energy;
        ave_energy2[0] = m_energy * m_energy;

        ave_isingE[0] = ising_E;
        ave_isingE2[0] = ising_E * ising_E;
        unsigned int nSteps = Nmax * m_grid.get_volume(); // 
        
        cout << setprecision(6) << fixed;
        
        for (int bin=1;bin<n_bins;bin++) {
            unsigned long int n_MC_move = 0;
            unsigned long int left_hist = 0, right_hist = 0;
            int min_win = min_OP + (bin-1) * w_win;
            int max_win = min_win + w_win; // (double)(bin+1) * w_bin;
            int mid_win = min_win + 1;
            int temp_orderParameter;
            cout<<internal<<"bin = "<<bin<<" \t";

            /** performing Successive Umbrella Sampling for the current window **/
            while (left_hist < nSteps || right_hist < nSteps) {
                int cell_index = (int) (drandom()*(double)m_grid.get_volume());

                int initial_state = m_grid.get_cellSpin(cell_index);

                temp_orderParameter = m_orderParameter - 2 * initial_state;

                if (min_win <= temp_orderParameter && temp_orderParameter <= max_win) {
                    double factor;
                    if (m_orderParameter<mid_win && temp_orderParameter>mid_win)
                        factor = ln_f;
                    else if (m_orderParameter>mid_win && temp_orderParameter<mid_win)
                        factor = -ln_f;
                    else
                        factor = 0.;

                    MC_move(cell_index, factor);

                    n_MC_move++;

                    if (n_MC_move==ULONG_MAX) {
                        cout<<"n_MC_move = "<<n_MC_move<<" exceeds long limit."<<endl;
                        abort();
                    }

                    ave_energy[bin] += m_energy;
                    ave_energy2[bin] += m_energy * m_energy;

                    ave_isingE[bin] += ising_E;
                    ave_isingE2[bin] += ising_E * ising_E;
                }

                if (m_orderParameter <= mid_win)
                    left_hist++;

                if (left_hist==ULONG_MAX) {
                    cout<<"left_hist = "<<left_hist<<" exceeds long limit."<<endl;
                    abort();
                }

                if (m_orderParameter >= mid_win)
                    right_hist++;

                if (right_hist==ULONG_MAX) {
                    cout<<"right_hist = "<<right_hist<<" exceeds long limit."<<endl;
                    abort();
                }

                /** do comment in long simulation **/
                // if ((left_hist!=0 && left_hist%1000000==0) || (right_hist!=0 && right_hist%1000000==0))
                    // cout<<left_hist<<"\t"<<right_hist<<endl;
            }

            /** preparing starting point for the next window *****************/
            while (m_orderParameter != max_win) {
                int cell_index = (int) (drandom()*(double)m_grid.get_volume());
                int initial_state = m_grid.get_cellSpin(cell_index);
                temp_orderParameter = m_orderParameter - 2 * initial_state;

                if (min_win <= temp_orderParameter && temp_orderParameter <= max_win) {
                    double factor;
                    if (m_orderParameter<mid_win && temp_orderParameter>mid_win)
                        factor = ln_f;
                    else if (m_orderParameter>mid_win && temp_orderParameter<mid_win)
                        factor = -ln_f;
                    else
                        factor = 0.;
                    MC_move(cell_index, factor);
                }
            }

            double ratio = (double)right_hist/(double)left_hist;
            log_ratio = (double)(log(ratio));
            ln_P[bin] = ln_P[bin-1] + log_ratio - ln_f;

            ln_f = ln_P[bin-1]-ln_P[bin];

            cout << "ln_f = " << ln_f << " \t" << right << "ratio = " << ratio << endl;

            ave_energy[bin] /= (double)(n_MC_move);
            ave_energy2[bin] /= (double)(n_MC_move);
            ave_isingE[bin] /= (double)(n_MC_move);
            ave_isingE2[bin] /= (double)(n_MC_move);

            /** printing the spin and height lattices at the middle bin **/
            if ( bin == n_bins/2 ) {
                print_spin_lattice();
            }

	    print_spin_lattice(bin);

        }
    }

    /** printing the probability and energy moments to files *****************/
    void print_prob_data () {
        sprintf(filename,"prob_data_%d_j%.2f.dat",m_lx,J_const);
        ofstream prob_file(filename);
        for (int bin=0;bin<n_bins;bin++)
            prob_file<<op_vector[bin]<<"\t"<<ln_P[bin]<<endl;
        prob_file.close();
    }

    void print_ave_energy_data() {
		sprintf(filename,"aveE_data_%d_j%.2f.dat",m_lx,J_const);
        ofstream aveE_file(filename);
        for (int bin=0;bin<n_bins;bin++)
            aveE_file<<ave_energy[bin]<<"\t"<<ave_energy2[bin]<<endl;
        aveE_file.close();
    }

    void print_ave_isingE_data() {
		sprintf(filename,"isingE_data_%d_j%.2f.dat",m_lx,J_const);
        ofstream ising_file(filename);
        for (int bin=0;bin<n_bins;bin++)
            ising_file<<ave_isingE[bin]<<"\t"<<ave_isingE2[bin]<<endl;
        ising_file.close();
    }

    /** printing spin lattice ************************************************/
    void print_spin_lattice() {
		sprintf(filename,"spin_%d_j%.2f.dat",m_lx,J_const);
        ofstream spin_file(filename);
        for (int i=0;i<m_lx;i++) {
            for (int j=0;j<m_ly;j++) {
                int cell_index = m_grid.get_cell_index(i,j);
                int cell_State = m_grid.get_cellSpin(cell_index);
                spin_file<<cell_State<<"\t";
            }
            spin_file<<endl;
        }
        spin_file.close();
    }

    /** **********************************************************************/
    void print_spin_lattice(int fn) {
        sprintf(filename,"spins/spin_%d_j%.2f_%d.dat",m_lx,J_const,fn);
        ofstream spin_file(filename);
        for (int i=0;i<m_lx;i++) {
            for (int j=0;j<m_ly;j++) {
                int cell_index = m_grid.get_cell_index(i,j);
                int cell_State = m_grid.get_cellSpin(cell_index);
                spin_file<<cell_State<<"\t";
            }
            spin_file<<endl;
        }
        spin_file.close();
    }

    /** **********************************************************************/
    void set_file_counter(int _c) {
        f_counter = _c;
    }

    void set_ln_f(double _f) {
        ln_f = _f;
    }

};

#endif
