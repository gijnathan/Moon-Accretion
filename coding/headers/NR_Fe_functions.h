// NR_Fe_functions.h

#ifndef NR_FE_FUNCTIONS_H
#define NR_FE_FUNCTIONS_H


using namespace std;

#include <iostream>
#include <iomanip>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <random>
#include <fstream>


// string line; 

// INPUT PARAMETER DECLARATIONS (read from file: 3isotope_input.txt, later to be sent in as input variables in nbdiff)
double DFe;// = 0.0; // phase partition coeff D for Fe
double TotEM;// = 0.0; //total Equilibrating Mantle (EM)(silicate)
double TotEC;// = 0.0; //total Equilibrating Core (EC)(metal)
double capDelta57;// = 0.0; // cap Delta 57 _ sil-met
double t54;// = 0.0; //total Fe54 (moles) defined in input file
double t56;// = 0.0; // total Fe56 (moles) defined	 in input file
double t57;// = 0.0; // total Fe56 (moles) defined	 in input file
int counter = 0; 
int counter_index=12; 

extern const int n; // these are defined in the main Nbdiff file
extern const int m; // these are defined in the main Nbdiff file 
// IRMM constants
const double Fe54_IRMM = 0.05845;
const double Fe56_IRMM = 0.91754;
const double Fe57_IRMM = 0.021191;
// other var
const double theta56 = 0.678; // CC in mathematica // ratio of delta56Fe to delta57Fe
const double R56 = Fe56_IRMM/Fe54_IRMM; // Fe56/Fe54 IRMM ratio
const double R57 = Fe57_IRMM/Fe54_IRMM; // Fe57/Fe54 IRMM ratio

const double R56_54_IRMM = Fe56_IRMM / Fe54_IRMM; //using fe56
const double R57_54_IRMM = Fe57_IRMM / Fe54_IRMM; //using fe57

const double Fe_56_Total_IRMM = Fe54_IRMM + Fe56_IRMM; // For normalizing when not including all isotopes of Fe...using fe56
const double Fe_57_Total_IRMM = Fe54_IRMM + Fe57_IRMM; // For normalizing when not including all isotopes of Fe..using fe57
const double Fe_Total_IRMM = Fe54_IRMM + Fe56_IRMM + Fe57_IRMM; // For normalizing when including 3 isotopes of Fe


// ofstream outfile0;

// declare functions
// used for epsilon: 
inline double get_m54(double s54,
					  double TotEM, double TotEC, double DFe, double capDelta57, double t54, double t56, double t57); // gets m54, given s54; (m54 = t54 - s54)
inline double get_m57(double s57,
					  double TotEM, double TotEC, double DFe, double capDelta57, double t54, double t56, double t57); // gets m57, given s57; (m57 = t57 - s57)
inline double get_epsilon(double s54, double m54, double s57, double m57, 
						  double TotEM, double TotEC, double DFe, double capDelta57, double t54, double t56, double t57); // calculate epsilon given s54,m54,s57,m57; (using equ (5))
inline double use_epsilon(double s54, double s57, 
						  double TotEM, double TotEC, double DFe, double capDelta57, double t54, double t56, double t57); // calculates values to send to get_epsilon, requiring only s54, and s57 (calculates m54 and m57 from those values)
inline double NR_epsilon(double s54, 
						 double TotEM, double TotEC, double DFe, double capDelta57, double t54, double t56, double t57); // this will take in s54 and return s57

// used for delta:
inline double get_m56(double s56, 
					  double TotEM, double TotEC, double DFe, double capDelta57, double t54, double t56, double t57); // gets m56, given s56; (m56 = t56 - s56)
inline double get_s56(double s54, 
					  double TotEM, double TotEC, double DFe, double capDelta57, double t54, double t56, double t57); // gets s56, given s54 (and s57 calculated within); (rearranged equation (4))
inline double get_delta(double s54, double s56, double s57,
						double TotEM, double TotEC, double DFe, double capDelta57, double t54, double t56, double t57); // this uses equation (6) which requires inputs of s54, s56, s57
inline double use_delta(double s54, 
						double TotEM, double TotEC, double DFe, double capDelta57, double t54, double t56, double t57); // if s54 and s57 are sent to this, it will use g_s56 to get s56 and can calculate delta
inline void NR_Fe(double (&Fe_sil_isotope)[3], double (&Fe_met_isotope)[3],
                  double equil_moles[], double equil_moles_H2O, double equil_moles_S_sil, double equil_moles_C_sil,
                  double equil_moles_oxygen, double equil_moles_H, double equil_moles_S_met, double equil_moles_C_met,
                  double P, double T); // 



inline double get_m54(double s54, 
					  double TotEM, double TotEC, double DFe, double capDelta57, double t54, double t56, double t57)
{	
	double m54 = t54 - s54; // equ (1)
	return m54;
}

inline double get_m57(double s57,
					  double TotEM, double TotEC, double DFe, double capDelta57, double t54, double t56, double t57)
{	
	double m57 = t57 - s57; // equ (3)
	return m57;
}

inline double get_epsilon(double s54, double m54, double s57, double m57, 
						  double TotEM, double TotEC, double DFe, double capDelta57, double t54, double t56, double t57)
{	
	
	double epsilon; 
	epsilon = capDelta57 - (s57/s54 - m57/m54)*(1000.0/R57);// equ (5)
	
	// cout << "() capDelta57: " << capDelta57 << endl;
	// cout << "() s57: " << s57 << endl;
	// cout << "() s54: " << s54 << endl;
	// cout << "() m57: " << m57 << endl;
	// cout << "() m54: " << m54 << endl;
	// cout << "() R57: " << R57 << endl;


	return epsilon;
}

inline double use_epsilon(double s54, double s57, 
						  double TotEM, double TotEC, double DFe, double capDelta57, double t54, double t56, double t57)
{
	double epsilon;

	double m54 = get_m54(s54,
						 TotEM, TotEC, DFe, capDelta57, t54, t56, t57); //get_m54(s54)
	double m57 = get_m57(s57, 
					     TotEM, TotEC, DFe, capDelta57, t54, t56, t57); //get_m57(s57)

	// cout << "a capDelta57: " << capDelta57 << endl;
	epsilon = get_epsilon(s54, m54, s57, m57, 
						  TotEM, TotEC, DFe, capDelta57, t54, t56, t57); //get_epsilon()
	// cout << "b capDelta57: " << capDelta57 << endl;
	
	return epsilon;

}

inline double NR_epsilon(double s54, 
						 double TotEM, double TotEC, double DFe, double capDelta57, double t54, double t56, double t57) // this will take in s54 and return s57
{
	// cout << " NR_epsilon() -----------------------------------------" << endl;
	double s57_guess1; 
	double s57_guess2;
	double s57_guess_new;
	double s57_guess_here;
	double epsilon_guess_new; 
	double epsilon_guess1;
	double epsilon_guess2;
	double m; 
	double b; 

	// cout << "NR_epsilon START" << endl;
	// cout << "s54: " << s54 << endl; 
	
	s57_guess1 = 0.000001 * t57; //just a little above zero to avoid nans 
	s57_guess2 = 0.999999 * t57; //just a little below total to avoid nans

	// cout << "t57: " << t57 << endl;

	// cout << "s57_guess1: " << s57_guess1 << endl;
	// cout << "s57_guess2: " << s57_guess2 << endl;

	// cout << "d capDelta57: " << capDelta57 << endl;
	
	epsilon_guess1 = use_epsilon(s54, s57_guess1, 
								 TotEM, TotEC, DFe, capDelta57, t54, t56, t57);
	epsilon_guess2 = use_epsilon(s54, s57_guess2, 
								 TotEM, TotEC, DFe, capDelta57, t54, t56, t57);
	
	// cout << "e capDelta57: " << capDelta57 << endl;
	
	// cout << "epsilon_guess1: " << epsilon_guess1 << endl;
	// cout << "epsilon_guess2: " << epsilon_guess2 << endl;

	epsilon_guess_new = 1.0; 
	int max_iteration = 100; 
	int iteration = 0; 
	// cout << "1 s57_guess_new: " << s57_guess_new << endl;
	
	while ((abs(epsilon_guess_new) > 1.0e-8) && (iteration < max_iteration))
	{
		// cout << "2 s57_guess_new: " << s57_guess_new << endl;

		//m = (y2 - y1)/(x2 - x1)
		m = (epsilon_guess1 - epsilon_guess2)/(s57_guess1 - s57_guess2);
		//b = y1 - mx1
		b = epsilon_guess1 - m * s57_guess1; 
		s57_guess_new = -b/m;
		// cout << "3 s57_guess_new: " << s57_guess_new << endl;
	 	
	 	// cout << "f capDelta57: " << capDelta57 << endl;
	
		epsilon_guess_new = use_epsilon(s54, s57_guess_new, 
										TotEM, TotEC, DFe, capDelta57, t54, t56, t57);
		
		// cout << "g capDelta57: " << capDelta57 << endl;
	
		if (epsilon_guess_new * epsilon_guess1 > 0.0)
		{
			s57_guess1 = s57_guess_new; 
			epsilon_guess1 = epsilon_guess_new; 
			iteration+=1;
		}
		else 
		{
			s57_guess2 = s57_guess_new; 
			epsilon_guess2 = epsilon_guess_new; 
			iteration+=1;
		}
		
	}
	// cout << "s57_guess_new: " << s57_guess_new << endl;
	// cout << "NR_epsilon END" << endl;
	
	return s57_guess_new; 
}

inline double get_m56(double s56,
					  double TotEM, double TotEC, double DFe, double capDelta57, double t54, double t56, double t57)
{	
	double m56 = t56 - s56; // equ (2)
	return m56;
}

inline double get_s56(double s54, 
					  double TotEM, double TotEC, double DFe, double capDelta57, double t54, double t56, double t57)
{
	// cout << "uio capDelta57: " << capDelta57 << endl;
	double m54 = get_m54(s54,
						 TotEM, TotEC, DFe, capDelta57, t54, t56, t57);
	// cout << "1" << endl;
	double s57 = NR_epsilon(s54, 
							TotEM, TotEC, DFe, capDelta57, t54, t56, t57); // this will take in s54 and return s57
	double m57 = get_m57(s57, 
						 TotEM, TotEC, DFe, capDelta57, t54, t56, t57); // gets m57, given s57; (m57 = t57 - s57)
	// s56 = (m54 + m56 + m57) * (TotEM/TotEC) / DFe - s54 - s57; // need to further rearrange equation (4) from this point...
	double s56 = ((m54 + m57 + t56) * (TotEM/(TotEC*DFe))  - s54 - s57) * (1.0 + (TotEC*DFe)/(TotEC * DFe + TotEM)); // rearranged equation (4)
	return s56; 
}

inline double get_delta(double s54, double s56, double s57, 
						double TotEM, double TotEC, double DFe, double capDelta57, double t54, double t56, double t57)
{	
	double delta; 
	delta = theta56 - (s56/(s54 * R56) - 1.0)/(s57/(s54 * R57) - 1.0);// equ (6)
	return delta;
}

inline double use_delta(double s54, 
						double TotEM, double TotEC, double DFe, double capDelta57, double t54, double t56, double t57)
{
	double delta; 
	// s54 = input
	// cout << "2" << endl;
	double s57 = NR_epsilon(s54, 
							TotEM, TotEC, DFe, capDelta57, t54, t56, t57);
	double s56 = get_s56(s54, 
						 TotEM, TotEC, DFe, capDelta57, t54, t56, t57); 
	// cout << "1 s54: " << s54 << endl;
	// cout << "1 s57: " << s57 << endl;
	// cout << "1 s56: " << s56 << endl;

	delta = get_delta(s54, s56, s57, 
					  TotEM, TotEC, DFe, capDelta57, t54, t56, t57);
	// cout << "2 s54: " << s54 << endl;
	// cout << "2 delta: " << delta << endl;
	// cout << "2 s57: " << s57 << endl;
	// cout << "2 s56: " << s56 << endl; 
	return delta; 	
}

inline void NR_Fe(double (&Fe_sil_isotope)[3], double (&Fe_met_isotope)[3],
                  double equil_moles[], double equil_moles_H2O, double equil_moles_S_sil, double equil_moles_C_sil,
                  double equil_moles_oxygen, double equil_moles_H, double equil_moles_S_met, double equil_moles_C_met,
                  double P, double T) //(13 total variables: 7 input, 6 output pass by reference (the 6 isotopes))
{	
	// ofstream outfile1;
	// outfile1.open ("outfile1.txt");
 //  	outfile1.setf(ios::showpoint);
	// for (int i=0;i<100000;i++){
	// 	double S_54_guess = i*(t54*0.01*(1.-1.E-9))/1000.+0.01*t54*1.E-9;
	// 	double delta = use_delta(S_54_guess);
	// 	outfile1 << S_54_guess << "\t" << delta << endl;
	// }
	// outfile1.close();

	double mseq_FeO = equil_moles[n - m];
	double mseq_Fe = equil_moles[n];

	// 7 inputs:
	double TotEM;
	double TotEC;
	double DFe;
	double capDelta57;
	double t54;
	double t56;
	double t57;

	counter+=1; 
	cout << "counter: " << counter << endl;

	
	TotEM = 0; //input 1
	for (int i = 0; i < n; i++)
	{
	  TotEM += equil_moles[i];
	}
	TotEM += equil_moles_H2O;
	TotEM += equil_moles_S_sil;
	TotEM += equil_moles_C_sil;
	
	TotEC = 0; //input 2
	for (int i = n; i < (n + m); i++)
	{
	  TotEC += equil_moles[i];
	}
	TotEC += equil_moles_oxygen;
	TotEC += equil_moles_H;
	TotEC += equil_moles_S_met;
	TotEC += equil_moles_C_met;

	t54 = Fe_sil_isotope[0] * mseq_FeO + Fe_met_isotope[0] * mseq_Fe; //input 5 total 54Fe in moles;
	t56 = Fe_sil_isotope[1] * mseq_FeO + Fe_met_isotope[1] * mseq_Fe; //input 6 total 56Fe in moles;
	t57 = Fe_sil_isotope[2] * mseq_FeO + Fe_met_isotope[2] * mseq_Fe; //input 7 total 57Fe in moles;

	// cout << "sum: " << Fe_met_isotope[0] + Fe_met_isotope[1] + Fe_met_isotope[2] << endl;

	DFe = (mseq_Fe / TotEC) / (mseq_FeO / TotEM); //input 3
	// cout << "mseq_Fe: " << mseq_Fe << endl;
	// cout << "TotEC: " << TotEC << endl;
	// cout << "mseq_FeO: " << mseq_FeO << endl;
	// cout << "TotEM: " << TotEM << endl;
	
	capDelta57 = 0.1; //input 4 this needs to have a set of if statements to handle delta_version determining what the value of capDelta57 is 
	
	
	// cout << "DFe: " << DFe << endl;
	// cout << "TotEM: " << TotEM << endl;
	// cout << "TotEC: " << TotEC << endl;
	// cout << "capDelta57: " << capDelta57 << endl;
	// cout << "t54: " << t54 << endl;
	// cout << "t56: " << t56 << endl;
	// cout << "t57: " << t57 << endl;
				  
	// variables for Newton-Raphson guess loop:
	double s54_guess1; 
	double s54_guess2;
	double s54_guess_new;
	double s54_guess_here;

	double s57_guess1; 
	double s57_guess2;

	double delta_here; 
	double delta_guess_new; 
	double delta_guess1;
	double delta_guess2;
	double m; 
	double b; 

	// cout << "---------------------------------------------" << endl;
	// cout << "\t Fe_sil_isotope[0]: " << Fe_sil_isotope[0] << endl;
	// cout << "\t Fe_sil_isotope[1]: " << Fe_sil_isotope[1] << endl;
	// cout << "\t Fe_sil_isotope[2]: " << Fe_sil_isotope[2] << endl;
	// cout << "\t Fe_met_isotope[0]: " << Fe_met_isotope[0] << endl;
	// cout << "\t Fe_met_isotope[1]: " << Fe_met_isotope[1] << endl;
	// cout << "\t Fe_met_isotope[2]: " << Fe_met_isotope[2] << endl;
	// cout << "---------------------------------------------" << endl;

	// s54_guess1 = 0.000001 * t54; //just a little above zero to avoid nans 
	// s54_guess2 = 0.999999 * t54; //just a little below total to avoid nans 

	// a better way of choosing the starting s54guesses 
	// finding the first pair with opposite signs, given increments of 1/100 * t54
	double smallDelta = 1e-9; 
	const int numIncrements = 100;
	double increment = t54 / static_cast<double>(numIncrements);

	double s54first = 0.0;
	double s54secnd = 0.0; 
	bool foundPair = false; 
	for (int i = 0; i <= numIncrements; i++) {
	    double S_54_guess1 = i * increment;
	    double S_54_guess2 = (i + 1) * increment;
	    if (i == 0) {
        S_54_guess1 += smallDelta;
        S_54_guess2 += smallDelta;
	    } else if (i == numIncrements) {
	        S_54_guess1 -= smallDelta;
	        S_54_guess2 -= smallDelta;
	    }
	    // use s54guess1,s54guess2
	    double delta1 = use_delta(S_54_guess1, 
							      TotEM, TotEC, DFe, capDelta57, t54, t56, t57);
	    double delta2 = use_delta(S_54_guess2,
	    						  TotEM, TotEC, DFe, capDelta57, t54, t56, t57);
	    // Check if delta1 and delta2 have opposite signs
	    if ((delta1 > 0 && delta2 < 0) || (delta1 < 0 && delta2 > 0)) {
	        s54first = S_54_guess1;
	        s54secnd = S_54_guess2;
	        foundPair = true;
	        break;  // Exit the loop since we've found the pair
	    }	
	}


	s54_guess1 = s54first; //just a little above zero to avoid nans 
	s54_guess2 = s54secnd; //just a little below total to avoid nans 

	
	delta_guess1 = use_delta(s54_guess1, 
							 TotEM, TotEC, DFe, capDelta57, t54, t56, t57);
	delta_guess2 = use_delta(s54_guess2, 
							 TotEM, TotEC, DFe, capDelta57, t54, t56, t57);

	// cout << "delta_guess1: " << delta_guess1 << endl;
	// cout << "delta_guess2: " << delta_guess2 << endl;

	// outfile0 << "A\t" << "s54_guess1" << "\t" << "delta_guess1" << "\tB\t" << "s54_guess2" << "\t" << "delta_guess2" << endl;
	// outfile0 << "A\t" << s54_guess1 << "\t" << delta_guess1 << "\tB\t" << s54_guess2 << "\t" << delta_guess2 << endl;
	// outfile0 << "B\t" << s54_guess2 << "\t" << delta_guess2 << endl;	
	// outfile0 << std::fixed << std::setprecision(15);

	// outfile0 << std::left << "iteration" << " " << "s54_guess1" << " "
 //    << std::setw(19) << "s54_guess2" << " "
 //    << std::setw(19) << "s54_guess_new" << " "
 //    << " 0 "
 //    << std::setw(19) << "delta_guess1" << " "
	// << std::setw(19) << "delta_guess2" << " "
	// << std::setw(19) << "delta_guess_new" 
 //    << std::endl;

	delta_guess_new = 1.0; 
	int max_iteration = 10000000; 
	int iteration = 0; 
	while ((abs(delta_guess_new) > 1.0e-9) && (iteration < max_iteration))
	{	

		m = (delta_guess1 - delta_guess2)/(s54_guess1 - s54_guess2);
		b = delta_guess1 - m * s54_guess1; 
		s54_guess_new = -b/m; 
		delta_guess_new = use_delta(s54_guess_new, 
							 		TotEM, TotEC, DFe, capDelta57, t54, t56, t57);
		
		if (delta_guess_new * delta_guess1 > 0.0) // same sign for delta_guess_new and delta_guess1
		{
			// outfile0 << std::fixed << std::setprecision(15);

			// outfile0 << std::left << iteration << " " << std::setw(19) << s54_guess1 << " "
	  //       << std::setw(19) << s54_guess2 << " "
	  //       << std::setw(19) << s54_guess_new << " "
	  //       << " A "
	  //       << std::setw(19) << delta_guess1 << " "
			// << std::setw(19) << delta_guess2 << " "
			// << std::setw(19) << delta_guess_new 
	  //       << std::endl;
		
			s54_guess1 = s54_guess_new; 
			delta_guess1 = delta_guess_new; 
			iteration += 1;
		}
		else // opposite signs for delta_guess_new and delta_guess1
		{
			// outfile0 << std::fixed << std::setprecision(15);

			// outfile0 << std::left << iteration << " " << std::setw(19) << s54_guess1 << " "
	  //       << std::setw(19) << s54_guess2 << " "
	  //       << std::setw(19) << s54_guess_new << " "
	  //       << " B "
	  //       << std::setw(19) << delta_guess1 << " "
			// << std::setw(19) << delta_guess2 << " "
			// << std::setw(19) << delta_guess_new 
	  //       << std::endl;
			
			s54_guess2 = s54_guess_new; 
			delta_guess2 = delta_guess_new; 
			iteration += 1;
			
		}

		
	}


	// Fe54 
	double s54_output = s54_guess_new; 
	double m54_output = get_m54(s54_guess_new, 
							    TotEM, TotEC, DFe, capDelta57, t54, t56, t57);
	
	// Fe57 obtained from that as well
	// cout << "3" << endl;
	double s57_output = NR_epsilon(s54_guess_new, 
								   TotEM, TotEC, DFe, capDelta57, t54, t56, t57);
	double m57_output = get_m57(s57_output, 
								TotEM, TotEC, DFe, capDelta57, t54, t56, t57);

	// Fe56 obtained from that as well
	double s56_output = get_s56(s54_guess_new, 
								TotEM, TotEC, DFe, capDelta57, t54, t56, t57);
	double m56_output = get_m56(s56_output,
								TotEM, TotEC, DFe, capDelta57, t54, t56, t57); 

 //    cout << "****************************************" << endl;

 //    cout << "s54_output: " << s54_output << endl;
 //    cout << "m54_output: " << m54_output << endl;
 //    cout << "s56_output: " << s56_output << endl;
 //    cout << "m56_output: " << m56_output << endl;
 //    cout << "s57_output: " << s57_output << endl;
 //    cout << "m57_output: " << m57_output << endl;
    

	// cout << "****************************************" << endl;

	if (counter==counter_index) cout << "s56_output: " << s56_output << endl;
	if (counter==counter_index) cout << "s54_output: " << s54_output << endl;
	if (counter==counter_index) cout << "s57_output: " << s57_output << endl;

	double mantle_d56 = ((s56_output / s54_output) / (Fe56_IRMM / Fe54_IRMM) - 1.0) * 1000.0 ;
	double mantle_d57 = ((s57_output / s54_output) / (Fe57_IRMM / Fe54_IRMM) - 1.0) * 1000.0 ;

	if (counter==counter_index) cout << "mantle_d56: " << mantle_d56 << endl;
	if (counter==counter_index) cout << "mantle_d57: " << mantle_d57 << endl; 
	if (counter==counter_index) cout << endl;
	double core_d56 = ((m56_output / m54_output) / (Fe56_IRMM / Fe54_IRMM) - 1.0) * 1000.0 ;
	double core_d57 = ((m57_output / m54_output) / (Fe57_IRMM / Fe54_IRMM) - 1.0) * 1000.0 ;
	// cout << "core_d56: " << core_d56 << endl;
	// cout << "core_d57: " << core_d57 << endl; 

	// cout << "****************************************" << endl;

	// UPDATE Fe_sil/met_isotope arrays: 
	Fe_sil_isotope[0] = s54_output; 
	Fe_sil_isotope[1] = s56_output; 
	Fe_sil_isotope[2] = s57_output; 
	Fe_met_isotope[0] = m54_output; 
	Fe_met_isotope[1] = m56_output; 
	Fe_met_isotope[2] = m57_output; 

	if (isnan(mantle_d56)){
	cout << "--------------------NAN--------------------" << endl;
	cout << "mantle_d56 = nan" << endl;
	cout << "here are the input variables: " << endl << endl;
	cout << "\t "  <<  DFe << " DFe" << endl;
	cout << "\t "  <<  TotEM << " TotEM" << endl;
	cout << "\t "  <<  TotEC << " TotEC" << endl;
	cout << "\t "  <<  capDelta57 << " capDelta57" << endl;
	cout << "\t "  <<  t54 << " t54" << endl;
	cout << "\t "  <<  t56 << " t56" << endl;
	cout << "\t "  <<  t57 << " t57" << endl;
	cout << "----------------------------------------" << endl << endl;
	}

	// else {
	// cout << "mantle_d56 is defined" << endl;
	// cout << "here are the input variables: " << endl << endl;
	// cout << "\t "  <<  DFe << " DFe" << endl;
	// cout << "\t "  <<  TotEM << " TotEM" << endl;
	// cout << "\t "  <<  TotEC << " TotEC" << endl;
	// cout << "\t "  <<  capDelta57 << " capDelta57" << endl;
	// cout << "\t "  <<  t54 << " t54" << endl;
	// cout << "\t "  <<  t56 << " t56" << endl;
	// cout << "\t "  <<  t57 << " t57" << endl;

	// }


	// cout << "---------------------------------------------" << endl;
	// cout << "\t Fe_sil_isotope[0]: " << Fe_sil_isotope[0] << endl;
	// cout << "\t Fe_sil_isotope[1]: " << Fe_sil_isotope[1] << endl;
	// cout << "\t Fe_sil_isotope[2]: " << Fe_sil_isotope[2] << endl;
	// cout << "\t Fe_met_isotope[0]: " << Fe_met_isotope[0] << endl;
	// cout << "\t Fe_met_isotope[1]: " << Fe_met_isotope[1] << endl;
	// cout << "\t Fe_met_isotope[2]: " << Fe_met_isotope[2] << endl;
	// cout << "---------------------------------------------" << endl;

	return; //s54_guess_new; 
}



#endif // NR_FUNCTIONS_H