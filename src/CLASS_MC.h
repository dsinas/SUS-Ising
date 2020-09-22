/*******************************************************************************
* 
* CLASS_MC.h: A class for implementing Monte Carlo move on 2D Ising model.
*
* This program is made freely available with the understanding that every 
* copy of this file must include this header and that it comes without any 
* WITHOUT ANY WARRANTY.
*******************************************************************************/

#ifndef CLASS_MC_H
#define CLASS_MC_H

#include "RndGen.h"
#include "CLASS_GRID.h"

using namespace std;

const double Boltzmann_const = 1.;
const int SpinD = -1;
const int SpinU = +1;

class MC_SystemClass {

protected:
    LatticeClass m_grid;
    int m_lx, m_ly;
    double m_beta, m_temperature;
    double J_const;
    double m_energy;
    int m_orderParameter;

    double ising_E, delta_ising_E;

    double spinD_frac;
    double spinU_frac;
    int initial_OP;

public:
    MC_SystemClass(const int _Lx, const int _Ly):
            m_grid(_Lx,_Ly), m_lx(_Lx), m_ly(_Ly), J_const(0.)

    {
    }

    ~MC_SystemClass() {}

    /* ********************************************************************* */
    void initialization() {

        initialize_lattices();

        m_orderParameter = get_order_parameter();
        m_energy = get_energy();

    }

    /** initializing lattices ************************************************/
    void initialize_lattices() {

        /** initializing the lattice states **/
        for (int i=0;i<m_grid.get_volume();i++) {
//             int spin;
//             if (drandom()<spinD_frac)
//                 spin = SpinD;
//             else
//                 spin = SpinU;

            m_grid.set_cellSpin(i,SpinD); // for SUS, starting configuration is all spins down
        }

}

    /** set order parameter corresponding to initial compositions ************/
    void set_order_parameter() {

        while (get_order_parameter()<initial_OP) {
            int idx = (int) (drandom()*(double)m_grid.get_volume());
            int cell_state = m_grid.get_cellSpin(idx);
            if (cell_state<0)
                m_grid.set_cellSpin(idx,SpinU);
        }

        while (get_order_parameter()>initial_OP) {
            int idx = (int) (drandom()*(double)m_grid.get_volume());
            int cell_state = m_grid.get_cellSpin(idx);
            if (cell_state>0)
                m_grid.set_cellSpin(idx,SpinD);
        }

    }

    /** **********************************************************************/
    int get_order_parameter() {
        int result(0);
        for (int i=0;i<m_grid.get_volume();i++)
            result += m_grid.get_cellSpin(i);
        return result;
    }

    /** **********************************************************************/
    /** calculating the total energy of the system ***************************/
    double get_energy() {

        double E0(0.);
        int spin(0), s_sum(0);

        for (int idx=0;idx<m_grid.get_volume();idx++) {

            /** energy contribution of spin-spin interaction **/
            spin = m_grid.get_cellSpin(idx);
            s_sum = m_grid.get_sum_neighbors_spin(idx);
            E0 += (double)(spin*s_sum);

        }

        ising_E = -0.5 * E0; // initial value for Ising energy term

        E0 *= -0.5 * J_const;

        return E0;

    }

    /** **********************************************************************/
    /** Monte Carlo move                                                     */
    /** **********************************************************************/
    void MC_move(int cell_index, double factor) {

        double Pacc, delta_E(0.);

        // for SUS, get the cell index from argument + factor
/*
	int cell_index;
        do {
             cell_index = (int) (drandom()*(double)m_grid.get_volume());
        } while (fixed_number(cell_index));
*/

        int spin = m_grid.get_cellSpin(cell_index);
        int s_sum = m_grid.get_sum_neighbors_spin(cell_index);

        delta_ising_E = 2.*spin*s_sum;
        delta_E = J_const*delta_ising_E;

        Pacc = exp(-m_beta*delta_E + factor);
        if (Pacc>=1. || drandom()<Pacc) {
            m_grid.flip_spin(cell_index);
            m_orderParameter -= 2*spin;
            m_energy += delta_E;
            ising_E += delta_ising_E;
        }

    }

    /** **********************************************************************/
    bool fixed_number(int idx) {
        double frac = 0.1;

        int minOP = -(int)(frac * (double)m_grid.get_volume());
        int maxOP = +(int)(frac * (double)m_grid.get_volume());
        minOP += initial_OP;
        maxOP += initial_OP;
        int cell_state = m_grid.get_cellSpin(idx);
        int temp_orderParameter = m_orderParameter - 2*cell_state;
        if (minOP <= temp_orderParameter && temp_orderParameter <= maxOP)
            return false;
        else
            return true;
    }

    /** performing KAWASAKI Monte Carlo ************************************/
    void Kawasaki() {

        double Pacc, delta_E(0.);
        int idx[2], spin[2], s_sum[2];

        do {
            idx[0] = (int) (drandom()*(double)m_grid.get_volume());
            idx[1] = (int) (drandom()*(double)m_grid.get_volume());

            spin[0] = m_grid.get_cellSpin(idx[0]);
            spin[1] = m_grid.get_cellSpin(idx[1]);

            s_sum[0] = m_grid.get_sum_neighbors_spin(idx[0]);
            s_sum[1] = m_grid.get_sum_neighbors_spin(idx[1]);
        }
        while (spin[0]==spin[1] || s_sum[0]==s_sum[1]);

	delta_E = J_const*(spin[1]-spin[0])*(s_sum[1]-s_sum[0]);

	Pacc = exp(-m_beta*delta_E);
        if (Pacc>=1 || drandom()<Pacc) {
            m_grid.swap_spins(idx[0],idx[1]);
            m_energy += delta_E;
        }

    }

    /** **********************************************************************/
    /** setting parameters ***************************************************/
    void set_beta(const double t_beta) {
        m_beta = t_beta;
    }

    void set_temperature(const double t_temperature) {
        m_temperature = t_temperature;
        m_beta = 1./(Boltzmann_const*t_temperature);
    }

    void set_J_const(const double t_J) {
        J_const = t_J;
    }

    /** **********************************************************************/
    void set_spin_composition(double _sd, double _su) {
        spinD_frac = _sd;
        spinU_frac = _su;
        initial_OP = (int)((spinD_frac*(double)SpinD + spinU_frac*(double)SpinU)*(double)m_grid.get_volume());
    }

};

#endif
