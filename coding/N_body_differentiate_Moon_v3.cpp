// N_body_differentiate.cpp
// Multistage differentiation of the accreting bodies in N-body simulation
// Performs either a one-off calculation or a least squares regression (depending on value of "fit" in file input.txt)
// 5 INPUT FILES REQUIRED: input.txt; distance_data.txt; NB_data.txt; final_planets.txt, GI.txt
// Two output files (one for silicate/s and one for metal/cores) * - but see below under v-11
// Volatile elements can be omitted by specifying n=14 (giving m=9)
// Mass of bodies can be adjusted by parameter "mass_adjust" in order to end with a planet of 1 earth-mass
// All bodies differentiate at time=0
// Computes O'Brien (2006) and Grand Tack (GT) models
// Number of starting compositions can be varied (maximum number is specified by  const "max_comps"
// Propagation of Kd uncertainties to final Earth mantle & core concentrations added in version 3
// Output of P, T and fO2 in output file "core.txt" added in version v-4
// Output of embryo and planetesimal masses in version v-5
// Propagated uncertainties on log Kd (Si) in version v-6
// Impactor mantle added only after metal-silicate equilibrium in version v-7 - eliminates need for disequilibrium
// Calculates fraction of silicate involved in equilibration using Eq. (6) of Deguen et al (2011, EPSL) in version v-8
// Version 9 includes W partitioning with dependence on silicate melt composition (Wade et al 2012)
// Version 10 includes Mo partitioning with dependence on silicate melt composition (Wade et al 2012)
// Version 11 includes S in oxidised composition and assumes that S partitions completely into the Fe melt
// Version 11 uses a modified input file; outputs compositions of core-forming metal increments to 3rd output file
// Version 12 - Mo IS TEMPORARILY EXCLUDED FROM CHI SQ CALCULATION;
// **v-12 is based on COMPOSITIONAL FeO GRADIENT**
// Version 12; 23th October 2012; D. Rubie
// Version 14 - replaces version 12 which had too many modifications (such as collisional erosion for "Mercury")
//    - earth planet# and mass_adjust parameter values moved from input.txt to NB_data.txt file
//    - FeO gradient redefined to give small (non-zero) value of xfe_met_3 end of gradient to avoid mass balance crashing
//    - maximum xfe_met value input from input.txt as "xfe_met_1" (replaces "xfe_met_0" which was formerly set as a constant to 0.999)
// Version 15 - includes water in compositions from beyond the snow line
// Version 15 - 5 March 2013; D. Rubie
// Version 16 - Array sizes defined using constants in order to handle large NB files (lines 55-57)
// Version 16 - 27 March 2013; D. Rubie
// Version 17 - includes choice of three composition-distance models
// version 17 - 2 April 2013; D. Rubie
// Version 18 - defines three distances for composition-distance profile:
//             dist_1 within which 99.9% Fe as metal and Si in metal
//             dist_2 within which Fe is partly oxidised
//             dist_3 within which Fe is completely oxidised and beyond which bodies are H2O bearing
// Version 18 - 9 July 2013, D. Rubie
// Version 18 - 2 October 2013: includes H in core composition
// Version 19 - revised peridotite melting curves, based on Zhang & Herzberg, Tronnes & Frost and Liebske & Frost
// Version 19 - includes new version of REACT and implements composition redefinition when xFe_liq = 0 in 2nd react call
// Oxygen Kd is corrected in version 19 - now gives correct values
// REACT convergence problems solved in v.19
// CORE INCREMENT CALCULATIONS CURRENTLY ONLY VALID WITH 100% METAL EQUILIBRATING
// Version 19 - 30 November 2013, Authors: D. Rubie, S. Jacobson, A. Morbidelli
// version 20 - includes new gradient model (steep gradient between dist_1 and dist_2 and shallow gradient between dist_2 and dist_3
// Arrays values etc re-set to zero on each func call (17/02/2014)
// "Switch" section deleted in order to use new input files in which new index # is always that of the largest body (10 March 2014)
// Version 21 includes sulphur - in this version, maximum S content of core-forming liquid is given by (Fe,Ni)S
// Version 21 includes option of specifying a fixed value of p_factor for embryos and only refining the value for planetesimals
// Version 21 - 2 May 2014. Authors: D. Rubie, S. Jacobson, A. Morbidelli
// Version 23 develops S and volatile elements using reverse gradient model (cf v-22 which uses CV component - unsuccessfully)
// Vera's S-saturation model incorporated - requires slightly lower equlibration pressures than the O'Neill model - see function S_conc
// FeS + HSE segregation up to the start of late veneer; determines time of last GI for embryos
// Propagated compositional uncertainties are output in mantle file
// Version 23 - 20th January 2015.
// VERSION 23A (25th August 2015) outputs two additional files listing evolution of mantle and core compositions of all embyros (excluding the final planets)
// 2nd September 2015 - Definitions of metal mole fractions changed (lines 3407-3410 and 3518-3521) to include "o_new" term
// VERSION 25 - 26 November 2015 - includes the effect of S on the metal-silicate partitioning of HSE
// VERSION 26  - 4 December 2015 - includes tungsten
// 14 December - Mg/Si gradient incorporated based on forsterite fractionation (enhancement of depletion). The gradient is currently linear from 0.7 to 3.0 AU
// and is defined by fo_en_max value that is read from the input file. By setting its value to 1.0, there is no gradient.
// 20 January 2016 - corrections made to treatment of sulfur (S-contents no longer limited to Xs=0.5 and S in equilibrating silicate partitions entirely into metal)
// 22 January 2016 - line 3499 in react corrected (definition of lithophile elements previously included all elements up to Pt)
// VERSION 27 - 12 February 2016 - S partitioning between metal and silicate is implemented (Boujibar et al but without effect of oxygen in metal)
// 15 February 2016 - Multistage FeS segregation occurs only in events that follow giant impacts, with the following exception:
//       In the case of both multistage and single stage FeS segregation, FeS segregation continues in all events after final GI until onset of late veneer
// VERSION 28 - 9 August 2016 - New react call, fits time_magma, matte_P_factor, and dist_0_S with a gradient S model using the HSEs as a constraint
// VERSION 29 - 14 March 2016 - fits dist_0_C with a gradient C model using Mo and potentially C and W as constraints
// VERSION 30 - 4 February 2018 - fit_cond is used to indicate which free parameters are being fit and the accordingly fitted elements (0 -> MSEs, 1 -> Na, 2 -> HSEs, S and C)
//        fit_cond -> 0 fits the MSEs by adjusting the free parameters in the PARAM_MSE.txt file
//        fit_cond -> 1 fits the NA by adjusting the parameter dist_0[3]
//        fit_cond -> 2 fits the HSEs, S, and C (not W or Mo) by adjusting the free parameters in the PARAM_HSE.txt file
// VERSION 31 Early 2020 - adding fit_mars (G. Nathan)
//      fit_mars -> 0 does not fit mars
//      fit_mars -> 1 fits only mars and not earth analogue
//      fit_mars -> 2 fits mars and earth
//      fit_mars -> 3 fits analogue "planet" to the moon (H. Nguyen)
// VERSION 32 November 2020 - adding ability to track and fractionate isotopes during core formation (G. Nathan)
// VERSION 33 May 2023 - more variation on iron isotope fractionation and adding Vanadium isotope fractionation
// VERSION 33b_3isotope July 2023 - tracking three isotopes of iron
// VERSION 35 - this version uses the Newton-Raphson approach to solving the Fe isotope composition following equilibration 
// VERSION 36 - this version uses the newly redone 3isotope ANALYTIC solution to solving the Fe isotope composition following equilibration 

// VERSION MOON v1 - lunar composition evolution from post-impact disk. Starting composition is Bulk Silicate Earth
// 9 June 2025
// Modified 'elements.h' to include moon_mantle estimates in wt% and fractional uncertainty for all elements. Elements that have values of -1.0 have unavailable data (H. Nguyen)
// Added additional fit_mars == 3 condition to fit for ONLY THE MOON. Parameters would fit the forming body to the Moon's mantle abundances estimated by Warren (2005) and Delano (1986). 
// Fixed an issue with the way the code read INPUT_Moon_ v1.txt when it read volatile elements in. When a body containing volatile elements existed within 1 AU, it made all of the material disappear. 
// Initially thought third column in COMPOSITIONAL DATA in INPUT_Moon_ v1.txt was ONLY the mass adjustment, but for volatile elements it was actually dist_0. 
// Fixed by making dist_0 for Na closer to the Sun in INPUT_Moon_ v1 

// 13 June 2025
// Version 1:
// Modified variables in code to relate to the Moon rather than Earth.
// Added var called "moon_mass" and assigned value of 7.34767E22 kg (as looked up on Google)
// Replaced all instances of "earth_mass" with "moon_mass"
// Replaced all instances of "earth_planet" with "moon_body" so that code reads in index of moon analog as moon_body instead of earth number
// Replaced all mentions of Earth's mantle/ Earth's core with Moon's mantle/ Moon's core in log.txt


// 26 June 2025
// Version 2:
// Added new moon parameters to control moon oxidation gradient:
// moon_xfe_met_1
// moon_xfe_met_2
// moon_xsi_met_1
// moon_xsi_met_0
// moon_disk_start (Earth radii)
// moon_disk_end (Earth radii)
// Updated INPUT_Moon file to version 2 and to accept these parameters (INPUT_Moon_v2.txt)
// Chanegd name of fit_mars to fit_planet

// 27 June 2025
// Removed moon_xsi variables (Si in metal is just 0)
// Added moon_dist_1 and moon_dist_2 to control gradient
// Added moon variables to Fe_metalfraction and write_outfile_headers
// Added dist_comp == 5 (it MUST be 5 in order to run Moon ver)
// NOTE: (it MUST be 5 in order to run Moon ver)
// Added new moon variables to input_parameter.h, which is read by p_inputs 
// Added new chi_sq punishment for (fit_cond == 0 && dist_comp == 5)
// New boundaries punish chi square when moon variables are unrealistic (when xfe > 1 or xfe < 0 or  dist 1 > dist 2)
// NEED TO UPDATE PARAM_MSE (June 27, 2025)
// Edited outfile // OUTPUT FINAL PLANET COMPOSITION DETAILS TO FILES to output distance in terms of Earth radii

// 2 July 2025
// Changed old distance_data reading stream to not be restricted to start at index of 4
// now can start at any index

// 7 July 2025
// Changed vers constant to Nbdiff-Moon-v2

// 10 July 2025
// Edited fit_mse in elements.h to fit for only these elements: Fe, Si, Ni, Co, Cr, and V

// 21 July 2025
// Fixed issue with chi sq contributions not summing up to the total chi sq. 
// Problem was due to chi sq contributions being hard coded to be divided by 4 while the total chi
// sq was divided by the degrees of freedom. FIxed this issue by having the chi sq contributions be 
// divided by dof instead of 4

// Moon ver. 3
// 24 July 2025
// Added new fit conditions and new parameters to control moon gradient

// Authors: D. Rubie, S. Jacobson, A. Morbidelli


using namespace std;
//#include "stdafx.h"
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <fstream>
#include <math.h>
#include <string>
#include "headers/elements.h"
#include "headers/input_parameters.h"
#include "headers/nr.h"
#include "headers/amoeba.h"
#include "headers/amotry.h"
// #include "headers/NR_Fe_functions.h" // Needed for version #35, no longer needed 


const string vers = "Nbdiff-Moon-v2";               // Version Moon #2
const string prog_name = "N_body_differentiate";

const double earth_mass = 5.97219E24;     // Mass of Earth in kg
const double moon_mass = 7.34767E22;      // Mass of Moon in kg
const double pi = 3.141593;
const double g_const = 6.673e-11;         // Gravitational constant
const double R = 8.314;                // Gas constant

const int n = 18;          // TOTAL NUMBER OF ELEMENTS USED (excluding oxygen and S) - must not exceed the number listed in file "elements.cpp"
// With n=12 excludes HSEs; n=17 includes HSE's and W
const int m = n - 4;       // total number of siderophile elements used (excluding oxygen and S) - "4" is currently the number of lithophile elements
const int nm = n + m;
const int nm3 = n + m + 3;
const int n_hse = 4;       // number of HSEs

// CONSTANTS THAT DEFINE ARRAY SIZES:
const int n_bodies = 6900; // Maximum number of embryos + planetesimals in simulation
const int n_imps = 1500;      // Maximum number of impacts suffered by any one body
const int n_emb = 405;     // Maximum number of embyros in simulation

const char tab = '\t';

// CONSTANTS FOR ISOTOPES: //fractionate
// Fe isotopic abundances
// IRMM-014 standard values, from Table 8 in Taylor et al 1992
const double Fe54_IRMM = 0.05845; 
const double Fe56_IRMM = 0.91754; 
const double Fe57_IRMM = 0.021191;
// const double Fe58_IRMM = 0.002819

const double R56 = Fe56_IRMM/Fe54_IRMM; //(*Ratio 56/54 IRMM-014 isotopic standard*)
const double R57 = Fe57_IRMM/Fe54_IRMM; //(*Ratio 57/54 IRMM-014 isotopic standard*)

// const double R56_54_IRMM = Fe56_IRMM / Fe54_IRMM; //using fe56
// const double R57_54_IRMM = Fe57_IRMM / Fe54_IRMM; //using fe57

// const double Fe_56_Total_IRMM = Fe54_IRMM + Fe56_IRMM; // For normalizing when not including all isotopes of Fe...using fe56
// const double Fe_57_Total_IRMM = Fe54_IRMM + Fe57_IRMM; // For normalizing when not including all isotopes of Fe..using fe57
// const double Fe_Total_IRMM = Fe54_IRMM + Fe56_IRMM + Fe57_IRMM; // For normalizing when including 3 isotopes of Fe


// CONSTANTS DEFINED IN FILE "elements.cpp":
// atomic weights and oxide molecular weights for the elements (defined in file "elements.cpp")
extern const double AW[];
extern const double MW[];
extern const double AW_O;              // Oxygen atomic weight
extern const double a[];               // Log10 Kd parameters (log Kd = a + b/T + cP/T) and uncertainties "sig.."
extern const double sig_a[];
extern const double lowp_a[];          // P < 5 GPa
extern const double lowp_sig_a[];      // P < 5 GPa
extern const double b[];
extern const double sig_b[];
extern const double lowp_b[];          // P < 5 GPa
extern const double lowp_sig_b[];      // P < 5 GPa
extern const double c[];
extern const double sig_c[];
extern const double lowp_c[];          // P < 5 GPa
extern const double lowp_sig_c[];      // P < 5 GPa
extern const double d[];
extern const double sig_d[];
extern const int val[];                // valence of each element

extern const double sulf_a[];          // Log10 Kd parameters (log Kd = a + b/T + cP/T) for FeS-silicate partitioning and uncertainties
extern const double sulf_sig_a[];
extern const double sulf_b[];
extern const double sulf_sig_b[];
extern const double sulf_c[];
extern const double sulf_sig_c[];

extern const int sil_c[];              // Value specifies concentration level in mantle output file (0=wt%, 1=ppm, 2=ppb)
extern const int met_c[];              // ditto for core output files
extern const char fit_mse[];               // Values ('y' or 'n') specifies whether or not the element is included in chi square calculation
extern const char fit_na[];               //
extern const char fit_hse[];               //
extern const char el_type[];           // Specifies type of element
extern const double EM[];              // Earth mantle oxide/element concentrations, variably as wt%, ppm or ppb
extern const double p_mantle[];           // Earth Mantle concentrations ALL AS WT%
extern const double sig_p_mantle[];       // Earth Fractional uncertainties on above values
extern const double m_mantle[];           // Mars Mantle concentrations as WT%
extern const double sig_m_mantle[];       // Mars sigma uncertainty on above value
extern const double moon_mantle[];       // Moon Mantle concentrations ALL AS WT%
extern const double sig_moon_mantle[];    // Moon fractional uncertainty on above values
extern const double MW_H2O;
extern const double MW_S2O;
extern const double AW_H;              // Atomic weight of H
extern const double AW_S;              // Atomic weight of S
extern const double AW_C;             // Atomic weight of C


// FUNCTION PROTOTYPE DECLARATIONS:
double func(Vec_I_DP &x);              // Function that models accretion and core formation and returns chi sq
// Following converts moles to wt% and outputs to files:
void weight_comp(double moles[], double moles_O, double moles_H2O, double moles_H, double moles_S_mantle, double moles_S_core, double moles_C_mantle, double moles_C_core, const double AW[],
                 const double MW[]);
void weight_comp_em(double moles[], double moles_O, double moles_H2O, double moles_H, double moles_S_mantle, double moles_S_core, double moles_C_mantle, double moles_C_core, const double AW[],
                    const double MW[]);
void weight_comp_imp(double moles[], double moles_O, double moles_H2O, double moles_H, double moles_S_mantle, double moles_S_core, double moles_C_mantle, double moles_C_core, const double AW[],
                     const double MW[]);
void moles_comp_source(double moles[], double moles_O, double moles_H2O, double moles_H, double moles_S_mantle, double moles_S_core, double moles_C_mantle, double moles_C_core,
                       double total_moles[], double total_moles_O, double total_moles_H2O, double total_moles_H, double total_moles_S_mantle, double total_moles_S_core, double total_moles_C_mantle, double total_moles_C_core);
void incr_comp(double moles[], double moles_oxygen, double moles_H, double moles_S_core, double moles_C_core, double sum_core, const double AW[]);
double pressure(double object_mass);            // calculates P(CMB)
double temperature(double P);       // calculates liquidus T for peridotite
std::string format(double f, int n);
// Calculates fraction of targets mentle that equilibrates with impactor's core:

void distrib_coeff(double (& kd)[n], double P, double T, double Xs);                                  // Kd values (Si, Ni, Co and trace elements)
double oxygen_kd_Frost(double xfeo_mw, double P, double T); // determines oxygen Kd (Frost et al 2010)
double oxygen_kd_Rubie(double moles[], double moles_oxygen, double moles_S_met,
                       double moles_S_sil, double moles_C_met, double moles_C_sil, double pkd[],
                       double feo_new, double & nio_new, double & sio2_new,
                       double & fe_new, double & ni_new, double & si_new, double & o_new,
                       double & xfeo_mw, double & xsi_liq, double & xfe_liq, double & xo_liq, double & xni_liq);
void react(double (&new_moles)[nm3], double moles[], double moles_oxygen, double moles_S_met, double moles_S_sil, double moles_C_met, double moles_C_sil, double P, double T, double pkd[]);   // mass balance/partitioning calculation - see Rubie et al (2011, suppl. info)
double S_partitioning_blanchard(double moles[], double moles_oxygen, double moles_S_met, double moles_S_sil, double moles_C_met, double moles_C_sil, double P, double T);       // calculated Ds and distribution of S between metal and silicate
double C_partitioning(double moles[], double moles_oxygen, double moles_S_met, double moles_S_sil, double moles_C_met, double moles_C_sil, double P, double T);       // calculated Dc and distribution of C between metal and silicate
void WMo_partitioning(double (& moles)[nm3], double moles_oxygen, double moles_S_met, double moles_S_sil, double moles_C_met, double moles_C_sil, double P, double T);
void FeS_HSE(double h_moles[], double h_moles_H2O, double h_moles_S, double h_moles_FeS, double (& pmoles_HSE)[2 * n_hse], double P);      // calculates HSE partitioning between FeS and silicate mantle
void write_headings();                                                  // writes headings in planet output files
void write_headings_em();                                               // writes headings in embyro output files
double S_conc(double moles[], double moles_H2O, double moles_S_mantle, double moles_C_mantle, double P);
double mass_f(double moles[], double moles_H2O, double moles_O, double moles_H, double moles_S_mantle, double moles_S_core, double moles_C_mantle, double moles_C_core, const double AW[],
              const double MW[]);  // Calculates body mass from number of accumulated moles

void write_outfile_header(ofstream &outfile, string simulation_no, int dist_comp, double p_factor, double em_equil, double pf_em, double dist_1,
                          double dist_2, double dist_3, double erf_width_1, double erf_width_2, double xfe_met_1, double xfe_met_2, double xfe_met_3, double xsi_met_0, double wtci_H2O,
                          double mass_adjust, double xmet_pl, double xmet_em, double x_mantle, double em_xsil, double t_overstep, const char el_type[],
                          char el_sym[][2], double dist_0[], double dist_0_S, double dist_0_C, double matte_HSE, double matte_P_factor, double pf, 
                          double moon_xfe_met_1, double moon_xfe_met_2, double moon_dist_1, double moon_dist_2, double moon_disk_start, double moon_disk_end); // Moon_version 
void update_contributions(int no_1, int no_2, double mass_1, double mass_2);
double Fe_metalfraction(int dist_comp, double dist, double dist_1, double dist_2, double dist_3, double xfe_met_1,  double xfe_met_2, double xfe_met_3, double erf_width_1, double erf_width_2, 
                        double moon_xfe_met_1, double moon_xfe_met_2, double moon_dist_1, double moon_dist_2, double moon_disk_start, double moon_disk_end); // Moon_version 
double Si_metalfraction(int dist_comp, double dist, double dist_1, double dist_2, double dist_3, double xsi_met_0, double erf_width_1, double erf_width_2);
void starting_molar_proportions(double dist,
                                double (&moles)[nm], double &moles_H2O, double &moles_S_sil, double &moles_C_sil,
                                double &moles_oxygen, double &moles_H, double &moles_S_met, double &moles_C_met);
double get_mass_factor(double moles[], double moles_H2O, double moles_S_sil, double moles_C_sil,
                       double moles_oxygen, double moles_H, double moles_S_met, double moles_C_met);
void undifferentiated_composition(double mass, double mass_factor,
                                  double (&moles)[nm], double &moles_H2O, double &moles_S_sil, double &moles_C_sil,
                                  double &moles_oxygen, double &moles_H, double &moles_S_met, double &moles_C_met);
void differentiated_composition(double mass, double mass_factor,
                                double (&moles)[nm], double &moles_H2O, double &moles_S_sil, double &moles_C_sil,
                                double &moles_oxygen, double &moles_H, double &moles_S_met, double &moles_C_met);
void metal_silicate_equilibration(double mass,
                                  double (&mseq_moles)[nm], double &mseq_moles_H2O, double &mseq_moles_S_mantle, double &mseq_moles_C_mantle,
                                  double &mseq_moles_O, double &mseq_moles_H, double &mseq_moles_S_core, double &mseq_moles_C_core,
                                  double &P, double &T, double &fo2, double Fe_isotope, int &zero_xfe_liq);
void initiate_bulk_composition(double mass, int e_no, int coll_no,
                               double moles[], double moles_H2O, double moles_S_sil, double moles_C_sil,
                               double moles_oxygen, double moles_H, double moles_S_met, double moles_C_met,
                               double P, double T, double fo2);
void initiate_isotopic_composition(int e_no, int coll_no, double Fe_sil_isotope[], double Fe_met_isotope[]); //fractionate //3isotope
void update_bulk_composition(int e_no, int coll_no,
                             double equil_moles[], double equil_moles_H2O, double equil_moles_S_sil, double equil_moles_C_sil,
                             double equil_moles_oxygen, double equil_moles_H, double equil_moles_S_met, double equil_moles_C_met,
                             double noneq_moles[], double noneq_moles_H2O, double noneq_moles_S_sil, double noneq_moles_C_sil,
                             double noneq_moles_oxygen, double noneq_moles_H, double noneq_moles_S_met, double noneq_moles_C_met,
                             double P, double T, double fo2);
void accretion_no_equilibration(double (&moles_target)[nm], double & moles_H2O_target, double & moles_S_mantle_target, double & moles_C_mantle_target,
                                double & moles_O_target, double & moles_H_target, double & moles_S_target, double & moles_C_target,
                                double moles_projectile[], double moles_H2O_projectile, double moles_S_mantle_projectile, double moles_C_mantle_projectile,
                                double moles_O_projectile, double moles_H_projectile, double moles_S_projectile, double moles_C_projectile);
double silicate_fraction(double mass_1, double mass_2,
                         double no_1_moles[], double no_1_moles_H2O, double no_1_moles_S_mantle, double no_1_moles_C_mantle,
                         double no_1_moles_oxygen, double no_1_moles_H, double no_1_moles_S_core, double no_1_moles_C_core,
                         double no_2_moles[], double no_2_moles_H2O, double no_2_moles_S_mantle, double no_2_moles_C_mantle,
                         double no_2_moles_oxygen, double no_2_moles_H, double no_2_moles_S_core, double no_2_moles_C_core);
void record_silicate_fraction(double mass_2, double xsil);
void identify_molar_reservoirs(double xsil, double xmet,
                               double no_1_moles[], double no_1_moles_H2O, double no_1_moles_S_mantle, double no_1_moles_C_mantle,
                               double no_1_moles_oxygen, double no_1_moles_H, double no_1_moles_S_core, double no_1_moles_C_core,
                               double no_1_Fe_sil_isotope[3], double no_1_Fe_met_isotope[3],
                               double no_2_moles[], double no_2_moles_H2O, double no_2_moles_S_mantle, double no_2_moles_C_mantle,
                               double no_2_moles_oxygen, double no_2_moles_H, double no_2_moles_S_core, double no_2_moles_C_core,
                               double no_2_Fe_sil_isotope[3], double no_2_Fe_met_isotope[3],
                               double (&equil_moles)[nm], double & equil_moles_H2O, double & equil_moles_S_mantle, double & equil_moles_C_mantle,
                               double & equil_moles_oxygen, double & equil_moles_H, double & equil_moles_S_core, double & equil_moles_C_core,
                               double (&equil_Fe_sil_isotope)[3], double (&equil_Fe_met_isotope)[3],
                               double (&noneq_moles)[nm], double & noneq_moles_H2O, double & noneq_moles_S_mantle, double & noneq_moles_C_mantle,
                               double & noneq_moles_oxygen, double & noneq_moles_H, double & noneq_moles_S_core, double & noneq_moles_C_core,
                               double (&noneq_Fe_sil_isotope)[3], double (&noneq_Fe_met_isotope)[3]);
void update_sources(int e_no, int coll_no, double xsil, double xmet,
                    double no_1_moles[], double no_1_moles_H2O, double no_1_moles_S_mantle, double no_1_moles_C_mantle,
                    double no_1_moles_oxygen, double no_1_moles_H, double no_1_moles_S_core, double no_1_moles_C_core,
                    double no_2_moles[], double no_2_moles_H2O, double no_2_moles_S_mantle, double no_2_moles_C_mantle,
                    double no_2_moles_oxygen, double no_2_moles_H, double no_2_moles_S_core, double no_2_moles_C_core,
                    double equil_moles[], double equil_moles_H2O, double equil_moles_S_mantle, double equil_moles_C_mantle,
                    double equil_moles_oxygen, double equil_moles_H, double equil_moles_S_core, double equil_moles_C_core,
                    double noneq_moles[], double noneq_moles_H2O, double noneq_moles_S_mantle, double noneq_moles_C_mantle,
                    double noneq_moles_oxygen, double noneq_moles_H, double noneq_moles_S_core, double noneq_moles_C_core);
void update_sources_hse(int e_no, int coll_no, double h_moles[], double h_moles_S, double h_moles_FeS, double pmoles_HSE[]);
void starting_isotopic_proportions(double (&Fe_sil_isotope)[3], double (&Fe_met_isotope)[3]);//GN fractionate //3isotope
double Liu_Partition(double T, double P); //GN fractionate
void Fe_fractionate(double (&Fe_sil_isotope)[3], double (&Fe_met_isotope)[3], //3isotope
                    double equil_moles[], double equil_moles_H2O, double equil_moles_S_sil, double equil_moles_C_sil,
                    double equil_moles_oxygen, double equil_moles_H, double equil_moles_S_met, double equil_moles_C_met,
                    double P, double T); //fractionate
void update_Fe_composition(int e_no, int coll_no,
                           double equil_moles[], double equil_Fe_sil_isotope[], double equil_Fe_met_isotope[],
                           double noneq_moles[], double noneq_Fe_sil_isotope[], double noneq_Fe_met_isotope[]); //fractionate //3isotope

ofstream outfile0, outfile1, outfile2, outfile3, outfile4, outfile5, outfile6, outfile7, outfile8, outfile9, outfile10, outfile12, outfile13, outfile14; // outfile11, //noMarsNumber

// ofstream testfile;


// INPUT PARAMETER DELARATIONS (read from file: input.txt)
int dist_comp = 0;                        // Choice of composition-distance model (0=constant; 1=step; 2=gradient)
int fit = 0;                              // Value of 0 -> one-off calculation; value of 1 -> least squares refinement
int fit_cond = 0;
int gi_out = 0;                           // Value of 0 -> full calculation; value of 1 returns only GI.txt to be edited by user
int em_equil = 0;                         // Value of 0 - pressure factor for embryo impacts = p_factor; value of 1- pressure factor for embryo impacts = pf_em
double p_factor(0);                       // pressure factor for planetesimals; also for embryos when em_equil=0
double pf_em(0);                       // pressure factor for embryos when em_equil=1
double pf(0);                          // pressure factor
double em_xsil = 1.0;                        // Factor by which xsil is changed for embryo-embryo impacts, relative to Deguen model
double x_mantle = 0.0;                    // Fraction of mantle that equilibrates with metal; value of zero - Deguen model applies; value <zero - Deguen model is overridden
int comps(0);                          // Number of starting compositions used

// Parameters that describe composition as a function of distance (AU) (* indicates a fitting parameter)
double xfe_met_1 = 0.0;
double xfe_met_2 = 0.0;
double xfe_met_3 = 0.0;
double xfe_met = 0.0;                     // Fraction of total Fe initially in metal
double xsi_met_0 = 0.0;                   // *Maximum fraction of total Si initially in metal
double xsi_met = 0.0;                     // Fraction of total Si initially in metal
double dist_1 = 0.0;                   // *Distance (AU) at which xsi_met becomes zero
double dist_2 = 0.0;                   // *Distance (AU) at which xfe_met becomes zero
double dist_3 = 0.0;                   // Distance (AU) beyond which compositions are H2O-bearing
double dist_min = 0.0;                    // Minimum distance in "distance_data" file
double erf_width_1 = 1.0;
double erf_width_2 = 1.0;

double mass_adjust(1.0);                  // Mass of all bodies are adjusted by this factor in order to end with one lunar mass
double xmet_pl = 1.0;                     // fraction of impacting planetesimal's core that equilibrates
double xmet_em = 1.0;                     // fraction of impacting embryo's core that equilibrates
double t_overstep(0);                     // Temperature adjustment relative to melting curve
int earth_planet(0);                      // Final Earth planet number for chi square calculation
// int mars_planet(0);                       // Final Mars planet number for chi square calculation //noMarsNumber
int moon_body(0);                            // Final Moon body number for chi square calculation

char el_sym[n][2];                        // Element symbols
double wtci[n] = {0.0};                   // Wt% CI element concentration (n = number of elements)
double en_dep[n] = {0.0};                  // Depletion/enhancement factors of n elements rel.to CI
double dist_0[n] = {0.0};                  // Only for voltile elements to define comp. gradient
char fit_elem[n] = {'n'};

// S variables:
char el_sym_S;                         // Element symbol for S
string el_sym_H2O;
double wtci_S = 0.0;                   // Wt% CI element concentration for S
double dist_0_S = 0.0;                    // Distance at which zero S is condensed

// C variables:
char el_sym_C;                         // Element symbol for C
double wtci_C = 0.0;                   // Wt% CI element concentration for C
double dist_0_C = 0.0;                    // Distance at which zero C is condensed
double Dc = 0.0;


// Isotopic Fractionation delta variables: //fractionate
int delta_version_Fe = 0;  // different versions of delta for isotopic fractionation. possible values: 1,2,3. See function: Fe_fractionate()
double delta_amount = 0.0; // cap-delta silicate-metal hard entry for use with delta_version_Fe=5 
double nishahar_version = 0; // choosing which version of Ni and Shahar 2022 data version

// Moon variables (June 27 2025)
double moon_xfe_met_1 = 0.0;
double moon_xfe_met_2 = 0.0;
double moon_dist_1 = 0.0; // (Earth radii)
double moon_dist_2 = 0.0; // (Earth radii)
double moon_disk_start = 0.0; // (Earth radii)
double moon_disk_end = 0.0; // (Earth radii)


double wtci_H2O = 0.0;                    // Wt% CI element concentration for H2O

int b_status[n_bodies] = {0};    // Tracks differentiation status of bodies (no previous record = 0; defined and differentiated = 1; previous collision history = 2; no longer exists = 3)
int e_coll[n_bodies] = {0};                // collision number for each body (only used for embryos)

double e_pl_mass[n_emb] = {0.0}; // Mass of each embryo that is derived from planetesimals
double e_em_mass[n_emb] = {0.0}; // Mass of each embryo that is derived from embryos
double e_SN_mass[n_emb] = {0.0}; // Mass of each embryo that is derived from bodies from beyond snow line
int GI[n_emb] = {0};             // Value is changed to 1 when there is a giant impact
int last_GI[n_emb] = {0};        // Value is changed to 1 when previous event was a giant impact - used to determine timing of FeS segregation

double e_time[n_emb][n_imps] = {0.0, 0.0};        // Time of accretion+differentiation event for each embryo (e = embryo)
double e_mass[n_emb][n_imps] = {0.0, 0.0};        // Mass of accreting embryos following each differentiation event
double e_P[n_emb][n_imps] = {0.0, 0.0};           // Equilibration pressure of each differentiation event
double e_T[n_emb][n_imps] = {0.0, 0.0};           // Equilibration temperature of each differentiation event
double e_fO2[n_emb][n_imps] = {0.0, 0.0};            // Equilibration fO2 of each differentiation event
double e_moles[n_emb][n + m][n_imps] = {0.0, 0.0, 0.0};   // Composition of acccreting embryos following each differentiation event
double e_moles_H2O[n_emb][n_imps] = {0.0, 0.0};      // Moles of H2O in mantle
double e_moles_S_mantle[n_emb][n_imps] = {0.0, 0.0}; // Moles of S in mantle
double e_moles_C_mantle[n_emb][n_imps] = {0.0, 0.0}; // Moles of C in mantle
double e_moles_oxygen[n_emb][n_imps] = {0.0, 0.0};   // Moles of O in core (indexes are body# and impact# respectively)
double e_moles_H[n_emb][n_imps] = {0.0, 0.0};        // Moles of H in core
double e_moles_S_core[n_emb][n_imps] = {0.0, 0.0};   // Moles of S in core
double e_moles_C_core[n_emb][n_imps] = {0.0, 0.0};   // Moles of C in core

double e_equil[n_emb][n_imps] = {0.0, 0.0};                       // Degree of equilibration of targets mantle
double i_equil[n_emb][n_imps] = {0.0, 0.0};                       // Degree of equilibration of impactors core

double i_mass[n_emb][n_imps] = {0.0, 0.0};        // Mass of accreted impactors following each differentiation event
double i_moles[n_emb][n + m][n_imps] = {0.0, 0.0, 0.0};   // Composition of accreted impactors following each differentiation event
double i_moles_H2O[n_emb][n_imps] = {0.0, 0.0};      // Moles of H2O in mantle
double i_moles_S_mantle[n_emb][n_imps] = {0.0, 0.0}; // Moles of S in mantle
double i_moles_C_mantle[n_emb][n_imps] = {0.0, 0.0}; // Moles of C in mantle
double i_moles_oxygen[n_emb][n_imps] = {0.0, 0.0};   // Moles of O in core (indexes are body# and impact# respectively)
double i_moles_H[n_emb][n_imps] = {0.0, 0.0};        // Moles of H in core
double i_moles_S_core[n_emb][n_imps] = {0.0, 0.0};   // Moles of S in core
double i_moles_C_core[n_emb][n_imps] = {0.0, 0.0};   // Moles of C in core


double s_mass[n_emb][n_imps] = {0.0, 0.0};        // Mass of accreted impactors following each differentiation event
double s_moles[n_emb][n + m][n_imps] = {0.0, 0.0, 0.0};   // Composition of accreted impactors following each differentiation event
double s_moles_H2O[n_emb][n_imps] = {0.0, 0.0};      // Moles of H2O in mantle
double s_moles_S_mantle[n_emb][n_imps] = {0.0, 0.0}; // Moles of S in mantle
double s_moles_C_mantle[n_emb][n_imps] = {0.0, 0.0}; // Moles of C in mantle
double s_moles_oxygen[n_emb][n_imps] = {0.0, 0.0};   // Moles of O in core (indexes are body# and impact# respectively)
double s_moles_H[n_emb][n_imps] = {0.0, 0.0};        // Moles of H in core
double s_moles_S_core[n_emb][n_imps] = {0.0, 0.0};   // Moles of S in core
double s_moles_C_core[n_emb][n_imps] = {0.0, 0.0};   // Moles of C in core

double incr_moles[n_emb][n + m][n_imps] = {0.0, 0.0, 0.0};   // Composition of core-forming increments for each differentiation event
double incr_moles_oxygen[n_emb][n_imps] = {0.0, 0.0};   // Ditto for O (indexes are body# and impact# respectively)
double incr_moles_H[n_emb][n_imps] = {0.0, 0.0};     // Ditto for H
double incr_moles_S[n_emb][n_imps] = {0.0, 0.0};     // Ditto for S
double incr_moles_C[n_emb][n_imps] = {0.0, 0.0};     // Ditto for C

// ISOTOPE TRACKING //fractionate
// double e_Fe_sil_isotope[n_emb][3][n_imps] = {0.0, 0.0, 0.0};   // Iron isotopic composition of acccreting embryos in silicate following each differentiation event
// dimensions: [number of embryos][3: 54Fe, 56Fe, 57Fe][number of impacts] //GN fractionate
 double e_Fe_sil_isotope[n_bodies][3][n_imps] = {0.0, 0.0, 0.0};   // Iron isotopic composition of acccreting embryos in silicate following each differentiation event

// double e_Fe_met_isotope[n_emb][3][n_imps] = {0.0, 0.0, 0.0};   // Iron isotopic composition of acccreting embryos in metal following each differentiation event
// dimensions: [number of embryos][3: 54Fe, 56Fe, 57Fe][number of impacts] //GN fractionate
double e_Fe_met_isotope[n_bodies][3][n_imps] = {0.0, 0.0, 0.0};   // Iron isotopic composition of acccreting embryos in metal following each differentiation event



double matte_P_factor = 0.0;        // Equilibration P is multiplied by this factor to give effective pressure of mantle S saturation
// therefore P_FeS = matte_P_factor*P = matte_P_factor*pf*pcmb and ks = matte_P_factor*pf
double time_magma = 0.0;            // Timescale of magma ocean crystallization
int matte_HSE = 0;                  // No HSE extraction by matte when value is 0; HSE extraction when value is 1

int fit_planet = 0;                   // 0 -> fit earth, 1 -> fit mars, 2 -> fit both, 3 -> fit moon

// Values of the following 3 parameters are read from "NB_data" input file (body mass parameters are declared below):
double max_pl = 0.0;                   // Maximum planetesimal mass (input from 2nd line of "NB_data.txt)
double t(0);                           // time of impact
int no_1(0), no_2(0);                     // Impactor body numbers

// Values of the following parameters are read from "distance_data.txt" input file
int body = 0;
double dist[n_bodies] = {0};

// Output file formatting parameters:
int norm_conc(0);                      // When value is 1, final errors and earth mantle composition are output for each planet
int tm_col(0);                             // value determines if time column is required in output files (function: write_headers)

// Declare vectors x and y, and matrix p:    // Vectors and arrays of amoeba fitting parameters
double chi_sq = 1000.0;
int chi_sq_multiplier = 0;

int count(0), count1(0);                  // Count number of iterations in "oxygen_kd" and "react" functions respectively
int write_file = 0;                       // When equal to 1, "func" will write final results to outfiles
int write_errors = 0;                     // When equal to 1, propagated errors are calculated

// Error propagation variables:
double sig_si = 0.0;                      // fraction of sigKD(Si) uncertainty added/subtracted to/from Log KD(Si); value defined in "input.txt" file
int error = 0;                         // When equal to +1 or -1 progates upper and lower uncertainties on log Kd parameters to mantle and core compositions
// PROPOGATED UNCERTAINTY LIMITS on final compositions:
double weightplus[10][n + m] = {0.0}, weightminus[10][n + m] = {0.0};

// The following parameters are used for estimate uncertainties
double wtminus[n_emb][n + m] = {0.0, 0.0};
double wtplus[n_emb][n + m] = {0.0, 0.0};
double wtbest[n_emb][n + m] = {0.0, 0.0};
double err[n_emb][n + m] = {0.0, 0.0};
double wtbest_O[n_emb] = {0.0};
double wtbest_H2O[n_emb] = {0.0};
double wtbest_S_mantle[n_emb] = {0.0};
double wtbest_S_core[n_emb] = {0.0};
double wtbest_C_mantle[n_emb] = {0.0};
double wtbest_C_core[n_emb] = {0.0};
double wtbest_H[n_emb] = {0.0};

// Last giant impact times and late veneer parameters
double t_last_GI[n_emb] = {0.0};
double em_tmax[n_emb] = {0.0};    // Time of final impact event
double t0_LV[n_emb] = {0.0};         // Time of late veneer initiation
double t_matte[n_emb] = {0.0};    // Time beyond which sulfide segregation can occur (specified in GI.txt file for each embryo)
int em_max = 0;                     // maximum embryo index #
int em_no[n_emb] = {0};            // embryo index numbers
int planet_no[10] = {0};        // Final planet index numbers
int number_embryos = 0;
int em_target[n_emb] = {0};        // index numbers of targets with which embryos merge (only actual bodies)
double em_merge_time[n_emb] = {0.0}; // time at which an embryo merges with a larger body


string title, line;

// Parameters used to specify average degrees of silicate equilibration (as well as max and min values)
int n_pl = 0;
int n0_pl = 0;
int n_em = 0;
double xsil_pl = 0.0;
double xsil_em = 0.0;
double xsil_pl_max = 0.0;
double xsil_pl_min = 0.999;
double xsil_em_max = 0.0;
double xsil_em_min = 0.999;

int mseq_call_no = 0;
int react_call_no = 0;
int no_react_call_no = 0;
int fractionate_call_no = 0;
int no_fractionate_call_no = 0;
//______________________________________________________________________________________________________________________________________________
// MAIN PROGRAM -- READS PRIMARY INPUT VARIABLES AND CALLS AMOEBA WHICH MINIMUMS THE CHI SQ FROM func()
int main(void)
{
// Ouput vers to screen to make sure using the correct one.
  cout << prog_name << ": " << vers << endl;

// testfile.open("test.txt");

//__________________________________________________________________________________________________________
// INPUT STARTING COMPOSITIONS AND MODEL PARAMETERS (Read INPUT_Moon_v3.txt)
// Primary input parameters always read, but not necessarily used depending on dist_comp and fit
  ifstream infile ("INPUT_Moon_v3.txt");
  if (infile.is_open())
  {
    getline (infile, line);
    infile >> fit_cond;           getline (infile, line);
    infile >> dist_comp;          getline (infile, line);
    infile >> fit;                getline (infile, line);
    infile >> gi_out;             getline (infile, line);
    infile >> sig_si;             getline (infile, line);
    infile >> em_equil;           getline (infile, line);
    infile >> p_factor;           getline (infile, line);
    infile >> pf_em;              getline (infile, line);
    infile >> dist_1;             getline (infile, line);
    infile >> dist_2;             getline (infile, line);
    infile >> dist_3;             getline (infile, line);
    infile >> xsi_met_0;          getline (infile, line);
    infile >> xfe_met_1;          getline (infile, line);
    infile >> xfe_met_2;          getline (infile, line);
    infile >> xfe_met_3;          getline (infile, line);
    infile >> xmet_pl;            getline (infile, line);
    infile >> xmet_em;            getline (infile, line);
    infile >> em_xsil;            getline (infile, line);
    infile >> x_mantle;           getline (infile, line);
    infile >> t_overstep;         getline (infile, line);
    infile >> matte_P_factor;     getline (infile, line);
    infile >> matte_HSE;          getline (infile, line);
    infile >> time_magma;         getline (infile, line);
    infile >> fit_planet;           getline (infile, line);
    infile >> Dc;                 getline (infile, line);
    infile >> delta_version_Fe;   getline (infile, line);
    infile >> delta_amount;       getline (infile, line);
    infile >> nishahar_version;   getline (infile, line);
    infile >> moon_xfe_met_1;     getline (infile, line); // NEW From Moon_v2
    infile >> moon_xfe_met_2;     getline (infile, line);
    infile >> moon_dist_1;     getline (infile, line);
    infile >> moon_dist_2;     getline (infile, line);
    infile >> moon_disk_start;    getline (infile, line);
    infile >> moon_disk_end;      getline (infile, line);
    getline (infile, line);       getline (infile, line);
    infile >> el_sym_H2O;
    infile >> wtci_H2O;           getline(infile, line);
    getline (infile, line);       getline (infile, line);

//______________________________________________________________________
// NEXT ALL ELEMENTS (el_symb[]) ARE READ IN TO SET THE CI COMPOSITIONS: IN NBDiff Moon version, this composition is the bulk silicate Earth
// the CI wt% of all elements are read in: wtci[]
// refractory elements have enhancement/depletion factors en_dep[]
// and volatile elements have a zero concentration distance dist_0[]
// Furthermore for volatile elements: interior to dist_0[], it is zero wt%,
// exterior to dist_3, it is wtci[], and in between is a gradient between the two
// FIRST S, WHICH IS CONSIDERED HERE TO PARTITION 100% INTO METAL
    infile >> el_sym_S; // element symbol
    infile >> wtci_S; // CI wt%
    infile >> dist_0_S; // zero wt% dist
    getline (infile, line);
    infile >> el_sym_C; // element symbol
    infile >> wtci_C; // CI wt%
    infile >> dist_0_C; // zero wt% dist
    getline (infile, line);
// NOW ALL OTHER ELEMENTS.
    for (int i = 0; i < n; i++)
    { 
      for (int j = 0; j < 2; j++)
      {
        if ((el_sym[i][0] == 'V') || (el_sym[i][0] == 'W')) continue; //for element symbols with only 1 character
        infile >> el_sym[i][j]; // element symbol
      }
      
      infile >> wtci[i]; // CI wt%
      if (el_type[i] == 'v' || el_type[i] == 'H') {
        infile >> dist_0[i]; // for volatiles, zero wt% dist
        
      }
      else {
        infile >> en_dep[i]; // for refractories, enhancement factor 
      }
      
    }
    getline (infile, line);
  
  }
  infile.close();

//__________________________________________________________________________________________________________
// GET INITIAL HELIOCENTRIC DISTANCE DATA FOR STARTING BODIES (Read distance_data.txt)
  ifstream infile2("distance_data.txt");

// // old version: h++ gives initial value of 4; 1-3 correspond to Jupiter, Saturn and Sun
// // h++ gives initial value of 4; 1-3 correspond to Jupiter, Saturn and Sun
//   int h(3);
//   if (infile2.is_open()) getline (infile2, line);
//   infile2 >> dist_min;
// // Loop through every line of infile2 (distance_data.txt)
// // Each line is a body index and it initial heliocentric distance
//   do
//   {
//     h++;
//     infile2 >> body >> dist[h];
//     getline (infile2, line);
//   }
//   while (!infile2.eof() );
//   infile2.close();

// NEW VERSION: assign body index to 'h'  — June 30, 2025
  int h;
  if (infile2.is_open()) getline(infile2, line);
  infile2 >> dist_min;
  do {
    infile2 >> h >> dist[h];
    getline(infile2, line);
  } while (!infile2.eof());
  infile2.close();

//__________________________________________________________________________________________________________
// GET ACCRETION HISTORY FOR ALL EMBRYOS AND DETERMINE THE LAST GIANT IMPACT (GI) TIMES FOR EACH EMBRYO (Read NB_data.txt and create GI.txt)
  string sim_no;
  double mass1(0.0), mass2(0.0), t(0.0);
  double max_pl(6);
  int no1(0), no2(0);

// Sets indices for embryos
  for (int i = 4; i < n_emb; i++) em_no[i - 4] = i;

  ifstream infile3("NB_data.txt");
  if (infile3.is_open()) getline (infile3, sim_no);
// reads planet# for earth analogue
  infile3 >> moon_body;         getline (infile3, line);
// reads planet # for mars analogue DO NOT LOOK FOR MARS_PLANET IN NB_data //noMarsNumber
  // infile3 >> mars_planet;         getline (infile3, line);
// Adjusts masses to give earth one Earth mass
  infile3 >> mass_adjust;          getline (infile3, line);
  infile3 >> max_pl;                getline (infile3, line);
// Loop through each accretion event in of inputfile3 (NB_data.txt)
// storing the event and finding the last giant impact of each embryo
  while (!infile3.eof() )
  {
// Read first 5 parameters from each line of data file
    infile3 >> t >> no1 >> mass1 >> no2 >> mass2;   getline (infile3, line);
// Finds maximum embryo number em_max
    if (no1 > em_max) em_max = no1;
// Records latest giant impact on embryo no1
    if (mass2 > max_pl * mass_adjust) t_last_GI[no1] = t;
// Records latest accretion event on embryo no1
    if (t > em_tmax[no1]) em_tmax[no1] = t;
// Records target and time when projectile embryos strike an embryo target (used in the output for EMBRYO_*.txt files)
    if (mass2 > max_pl * mass_adjust)
    {
      for (int i = 0; i < n_emb; i++)
        if (no2 == em_no[i])
        {
          em_target[i] = no1;
          em_merge_time[i] = t;
        }
    }
  }
  infile3.close();

//__________________________________________________________________________________________________________
// CREATE GI.txt file IF IT DOES NOT EXIST
// GI.txt contains a table of final giant impact times, final impact times, start of late veneer and start of hadean matte times
// The first two columns of times are inputs from NB_data.txt and so cannot be altered by the user
// However, the second two columns of times are input parameters which the user can alter if so desired
// To do so, run the code with the gi_out variable in the INPUT_*.txt file set to 1, then alter the GI.txt file as desired
// Otherwise, run the code with gi_out variable in the INPUT_*.txt file set to 0, then rely on the code only to manipulate t0_late veneer and t0_matte
  if (gi_out == 1) remove("GI.txt");
  ifstream infile0("GI.txt");
  if (!infile0.is_open() || gi_out == 1)
  {
    ofstream outfile;
    outfile.open("GI.txt");
    outfile.setf(ios::showpoint);
    outfile << sim_no << endl;
// Default time beyond which the late veneer begins
    double default_LV = 5.0e8;
// Default time beyond which FeS exsolution and segregation can occur
    double default_t_matte = 0.0e8;
    outfile << "t0_late veneer = 5.0e8 years is a default value that results in no late veneer" << endl;
    outfile << "t0_matte = 0.0 years is a default value that results in FeS segregation THROUGHOUT accretion (whenever S saturation occurs)" << endl;
    outfile << "All times in years" << endl;
    outfile << "Em#        t_final_GI          t_final_imp         t0_late veneer      t0_matte" << endl;
    for (int i = 4; i < em_max + 1; i++)
    {
      if (em_tmax[i] > 0.0)
      {
        outfile << fixed << setw(3) << i << scientific << setw(20) << t_last_GI[i] << setw(20) << em_tmax[i] << setw(20) << default_LV << setw(20) << default_t_matte << endl;
      }
    }
    outfile.close();
    if (gi_out == 1) exit(0);
  }
  infile0.close();

//__________________________________________________________________________________________________________
// OPEN OUTPUT FILES

  outfile0.open ("ATTEMPTS.txt", std::ios_base::app);
  outfile0.setf(ios::showpoint);

  outfile1.open ("MANTLE.txt");
  outfile1.setf(ios::showpoint);
  outfile1 << "C++ code: " << prog_name << ".cpp. Version: " << vers << endl;
  if (fit == 1) outfile1 << "Results of least squares fit" << endl;
  else if (fit == 0) outfile1 << "Results of one-off calculation" << endl;
  outfile1 << "MANTLE EVOLUTION" << endl;
  outfile1 << "Input file: INPUT_Moon_v3.txt" << endl;

  outfile2.open ("CORE.txt");
  outfile2.setf(ios::showpoint);
  outfile2 << "C++ code: " << prog_name << ".cpp. Version: " << vers << endl;
  if (fit == 1) outfile2 << "Results of least squares fit" << endl;
  else if (fit == 0) outfile2 << "Results of one-off calculation" << endl;
  outfile2 << "CORE EVOLUTION" << endl;
  outfile2 << "Input file: INPUT_Moon_v3.txt" << endl;

  outfile3.open ("CORE_INCRE.txt");
  outfile3.setf(ios::showpoint);
  outfile3 << "C++ code: " << prog_name << ".cpp. Version: " << vers << endl;
  if (fit == 1) outfile3 << "Results of least squares fit" << endl;
  else if (fit == 0) outfile3 << "Results of one-off calculation" << endl;
  outfile3 << "CORE INCREMENTS" << endl;
  outfile3 << "Input file: INPUT_Moon_v3.txt" << endl;

  outfile4.open ("EMBRYO_MANTLE.txt");
  outfile4.setf(ios::showpoint);
  outfile4 << "C++ code: " << prog_name << ".cpp. Version: " << vers << endl;
  if (fit == 1) outfile4 << "Results of least squares fit" << endl;
  else if (fit == 0) outfile4 << "Results of one-off calculation" << endl;
  outfile4 << "EMBRYO MANTLE EVOLUTION" << endl;
  outfile4 << "Input file: INPUT_Moon_v3.txt" << endl;

  outfile5.open ("EMBRYO_CORE.txt");
  outfile5.setf(ios::showpoint);
  outfile5 << "C++ code: " << prog_name << ".cpp. Version: " << vers << endl;
  if (fit == 1) outfile5 << "Results of least squares fit" << endl;
  else if (fit == 0) outfile5 << "Results of one-off calculation" << endl;
  outfile5 << "EMBRYO CORE EVOLUTION" << endl;
  outfile5 << "Input file: INPUT_Moon_v3.txt" << endl;

  outfile6.open ("IMPACTOR_MANTLE.txt");
  outfile6.setf(ios::showpoint);
  outfile6 << "C++ code: " << prog_name << ".cpp. Version: " << vers << endl;
  if (fit == 1) outfile6 << "Results of least squares fit" << endl;
  else if (fit == 0) outfile6 << "Results of one-off calculation" << endl;
  outfile6 << "EMBRYO MANTLE EVOLUTION" << endl;
  outfile6 << "Input file: INPUT_Moon_v3.txt" << endl;

  outfile7.open ("IMPACTOR_CORE.txt");
  outfile7.setf(ios::showpoint);
  outfile7 << "C++ code: " << prog_name << ".cpp. Version: " << vers << endl;
  if (fit == 1) outfile7 << "Results of least squares fit" << endl;
  else if (fit == 0) outfile7 << "Results of one-off calculation" << endl;
  outfile7 << "EMBRYO CORE EVOLUTION" << endl;
  outfile7 << "Input file: INPUT_Moon_v3.txt" << endl;

  outfile8.open ("SOURCE_MANTLE.txt");
  outfile8.setf(ios::showpoint);
  outfile8 << "C++ code: " << prog_name << ".cpp. Version: " << vers << endl;
  if (fit == 1) outfile8 << "Results of least squares fit" << endl;
  else if (fit == 0) outfile8 << "Results of one-off calculation" << endl;
  outfile8 << "EMBRYO CORE EVOLUTION" << endl;
  outfile8 << "Input file: INPUT_Moon_v3.txt" << endl;

  outfile9.open ("SOURCE_CORE.txt");
  outfile9.setf(ios::showpoint);
  outfile9 << "C++ code: " << prog_name << ".cpp. Version: " << vers << endl;
  if (fit == 1) outfile9 << "Results of least squares fit" << endl;
  else if (fit == 0) outfile9 << "Results of one-off calculation" << endl;
  outfile9 << "EMBRYO CORE EVOLUTION" << endl;
  outfile9 << "Input file: INPUT_Moon_v3.txt" << endl;

  outfile10.open ("FINAL.txt");
  outfile10.setf(ios::showpoint);
  outfile10 << "C++ code: " << prog_name << ".cpp. Version: " << vers << endl;
  if (fit == 1) outfile10 << "Results of least squares fit" << endl;
  else if (fit == 0) outfile10 << "Results of one-off calculation" << endl;
  outfile10 << "EMBRYO CORE EVOLUTION" << endl;
  outfile10 << "Input file: INPUT_Moon_v3.txt" << endl;

  //noMarsNumber
  // outfile11.open ("FINAL_MARS.txt");
  // outfile11.setf(ios::showpoint);
  // outfile11 << "C++ code: " << prog_name << ".cpp. Version: " << vers << endl;
  // if (fit == 1) outfile11 << "Results of least squares fit" << endl;
  // else if (fit == 0) outfile11 << "Results of one-off calculation" << endl;
  // outfile11 << "MARS EMBRYO CORE EVOLUTION" << endl;
  // outfile11 << "Input file: INPUT_Moon_v3.txt" << endl;

  outfile12.open ("COMPOSITION_Fe.txt");
  outfile12.setf(ios::showpoint);
  outfile12 << "C++ code: " << prog_name << ".cpp. Version: " << vers << endl;
  if (fit == 1) outfile12 << "Results of least squares fit" << endl;
  else if (fit == 0) outfile12 << "Results of one-off calculation" << endl;
  outfile12 << "ISOTOPIC COMPOSITION OF FE" << endl;
  outfile12 << "Input file: INPUT_Moon_v3.txt" << endl;

  outfile13.open ("BODY_28.txt");
  outfile13.setf(ios::showpoint);
  // outfile13 << "INFORMATION FOR BODY 28" << endl;
  // outfile13 << "Input file: INPUT_Moon_v3.txt" << endl;

  outfile14.open ("testfile.txt");
  outfile14.setf(ios::showpoint);

//__________________________________________________________________________________________________________
//
// initialize NP and MP variables
  int NP;
  int MP;
// FIT CONDITION = 0, Using one-off calculation
  if (fit_cond == 0)
  {
    for (int i = 0; i < n; i++) fit_elem[i] = fit_mse[i];


    if (dist_comp == 1) NP = 2;
    else if (dist_comp == 2) NP = 6;
    else if (dist_comp == 3) NP = 5;
    else if (dist_comp == 4) NP = 6;
    else if (dist_comp == 5) NP = 5;
    MP = NP + 1;

// ASSIGN THE FITTED PARAMETERS
// Must match the list in input_parameters.h
    Vec_DP x(NP);
    if (dist_comp == 1)  // CONSTANT
    {
      x[0] = p_factor;
      x[1] = xfe_met_2;
    }
    else if (dist_comp == 2)  // STEP
    {
      x[0] = p_factor;
      x[1] = dist_1;
      x[2] = dist_2;
      x[3] = xsi_met_0;
      x[4] = xfe_met_2;
      x[5] = xfe_met_3;
    }
    else if (dist_comp == 3)  // GRADIENT
    {
      x[0] = p_factor;
      x[1] = dist_1;
      x[2] = dist_2;
      x[3] = xsi_met_0;
      x[4] = xfe_met_2;
    }
    else if (dist_comp == 4)  // ERF
    {
      x[0] = p_factor;
      x[1] = dist_1;
      x[2] = dist_2;
      x[3] = xsi_met_0;
      x[4] = erf_width_1;
      x[5] = erf_width_2;
    }
    else if (dist_comp == 5)  // MOON GRADIENT (Moon ver.)
    {
      x[0] = p_factor;
      x[1] = moon_xfe_met_1;
      x[2] = moon_xfe_met_2;
      x[3] = moon_dist_1;
      x[4] = moon_dist_2;
    }

//__________________________________________________________________________________________________________
// CODE PERFORMS A ONE-OFF CALCULATION
    if (fit == 0)
    {
      // Output the chi sq and the 1 sigma uncertainty from Kd
      cout << "ONE-OFF CALCULATION" << endl;
      cout << "CHI-SQ VALUES WITH Kd's +/- 1 SIGMA" << endl;
      write_file = 0;   // Make sure the output files are NOT written yet
      write_errors = 1; // Records error data including when error is -1, +1, and 0
      error = -1;
      cout << endl << "MINUS ERROR CALCULATION:" << endl;
      cout << "CHI SQ = " << func(x) << endl;
      error = +1;
      cout << endl << "PLUS ERROR CALCULATION:" << endl;
      cout << "CHI SQ = " << func(x) << endl;
      error = 0;
      write_file = 1;   // Make sure the output files ARE written now
      cout << endl << "BEST FIT CALCULATION:" << endl;
      chi_sq = func(x);
      cout << "CHI SQ = " << chi_sq << endl;
    }
// CODE PERFORMS A LEAST SQUARES REGRESSION
    else if (fit == 1)
    {
      Vec_DP y(MP);
      Mat_DP p(MP, NP);
      // Number of Amoeba evaluations
      int nfunc;
//__________________________________________________________________________________________________________
// FOR FITTING CALCULATIONS, READ PARAM.txt for initial parameter guesses
// reads MP guesses of NP fit parameters (MP=NP+1) for AMOEBA
      ifstream paramfile("PARAM_MSE.txt");
      if (paramfile.is_open())
      {
        for (int i = 0; i < NP; i++)
        {
          for (int j = 0; j < MP; j++) {
            paramfile >> p[j][i];
          }
          getline(paramfile, line);
        }
      }
      paramfile.close();

// NP is the number of fitting parameters and MP is the number of initial guess arrays
// p is the matrix of initial guesses read from the input file (PARAM.txt)
// so x is an array of one set of initial guesses
// and y is the corresponding chi sq
      write_file = 0;   // Make sure the output files are NOT written yet
      write_errors = 0; // Does not record error data
      error = 0;        // Uses the nominal Kd values

// REPORT THE INITIAL SIMPLEX
      cout << endl << "INITIAL MATRIX OF MINIMIZATION PARAMETERS AND CHI SQ VALUES:" << endl;
      for (int i = 0; i < MP; i++)
      {
        for (int j = 0; j < NP; j++) x[j] = p[i][j];
        y[i] = func(x);
      }

// MINIMIZATION CALCULATION USING AMOEBA
      cout << endl << "MINIMIZATION" << endl;
      cout << "PARAMETER GUESSES AND CHI SQ VALUES:" << endl;
      const DP FTOL = 1.0e-3;
      NR::amoeba2(p, y, FTOL, func, nfunc);

// REPORT THE TOTAL NUMBER OF EVALUATIONS
      cout << endl << "Number of function evaluations: " << nfunc << endl;

// REPORT THE FINAL SIMPLEX
      cout << endl << "PARAMETERS AT VERTICES OF FINAL SIMPLEX AND CHI SQ VALUES:" << endl;
      for (int i = 0; i < MP; i++)
      {
        for (int j = 0; j < NP; j++) cout << fixed << setprecision(5) << p[i][j] << tab;
        cout << "chi sq = " << fixed << setprecision(5) << y[i] << endl;
      }

// CALCULATE THE FINAL AVERAGE VALUES IN THE CENTER OF THE SIMPLEX
      cout << endl << "PARAMETERS AT CENTER OF FINAL SIMPLEX AND CHI SQ VALUE:" << endl;
      for (int j = 0; j < NP; j++)
      {
        DP sum = 0.0;
        for (int i = 0; i < MP; i++) sum += p[i][j];
        x[j] = sum / MP;
      }
      DP yfinal = func(x);

// SELECT PARAMETERS THAT GIVE LOWEST CHI SQ AS FINAL
      for (int i = 0; i < MP; i++)
      {
        if (y[i] < yfinal)
        {
          for (int j = 0; j < NP; j++) x[j] = p[i][j];
          yfinal = y[i];
        }
      }
      cout << endl << "FINAL BEST FIT PARAMETERS AND CHI SQ FROM FIT:" << endl;
      for (int j = 0; j < NP; j++) cout << fixed << setprecision(5) << x[j] << tab;
      cout << "chi sq = " << fixed << setprecision(5) << yfinal << endl;
      for (int j = 0; j < NP; j++) outfile0 << fixed << setprecision(6) << x[j] << tab;
      outfile0 << fixed << setprecision(6) << yfinal << endl;

// Output the final chi sq and the 1 sigma uncertainty from Kd
      cout << endl << "CHI-SQ VALUES WITH Kd's +/- 1 SIGMA" << endl;
      write_file = 0;   // Make sure the output files are NOT written yet
      write_errors = 1; // Records error data including when error is -1, +1, and 0
      error = -1;
      cout << endl << "MINUS ERROR CALCULATION:" << endl;
      chi_sq = func(x);
      cout << "CHI SQ = " << chi_sq << endl;
      error = +1;
      cout << endl << "PLUS ERROR CALCULATION:" << endl;
      chi_sq = func(x);
      cout << "CHI SQ = " << chi_sq << endl;
      error = 0;
      write_file = 1;   // Make sure the output files ARE written now
      cout << endl << "BEST FIT CALCULATION:" << endl;
      chi_sq = func(x);
      cout << "CHI SQ = " << chi_sq << endl;
    }
  }
  else if (fit_cond == 1)
  {
    for (int i = 0; i < n; i++) fit_elem[i] = fit_na[i];

    NP = 1;
    MP = NP + 1;

// ASSIGN THE FITTED PARAMETERS
// Must match the list in input_parameters.h
    Vec_DP x(NP);
    x[0] = dist_0[3];

    if (fit == 0)
    {
      // Output the chi sq and the 1 sigma uncertainty from Kd
      cout << "ONE-OFF CALCULATION" << endl;
      cout << "CHI-SQ VALUES WITH Kd's +/- 1 SIGMA" << endl;
      write_file = 0;   // Make sure the output files are NOT written yet
      write_errors = 1; // Records error data including when error is -1, +1, and 0
      error = -1;
      cout << endl << "MINUS ERROR CALCULATION:" << endl;
      cout << "CHI SQ = " << func(x) << endl;
      error = +1;
      cout << endl << "PLUS ERROR CALCULATION:" << endl;
      cout << "CHI SQ = " << func(x) << endl;
      error = 0;
      write_file = 1;   // Make sure the output files ARE written now
      cout << endl << "BEST FIT CALCULATION:" << endl;
      chi_sq = func(x);
      cout << "CHI SQ = " << chi_sq << endl;
    }
// CODE PERFORMS A LEAST SQUARES REGRESSION
    else if (fit == 1)
    {
      Vec_DP y(MP);
      Mat_DP p(MP, NP);
      // Number of Amoeba evaluations
      int nfunc;
//__________________________________________________________________________________________________________
// MP guesses of NP fit parameters (MP=NP+1) for AMOEBA
      p[0][0] = -5;
      p[1][0] = 5;

// NP is the number of fitting parameters and MP is the number of initial guess arrays
// p is the matrix of initial guesses read from the input file (PARAM.txt)
// so x is an array of one set of initial guesses
// and y is the corresponding chi sq
      write_file = 0;   // Make sure the output files are NOT written yet
      write_errors = 0; // Does not record error data
      error = 0;        // Uses the nominal Kd values

// REPORT THE INITIAL SIMPLEX
      cout << endl << "INITIAL MATRIX OF MINIMIZATION PARAMETERS AND CHI SQ VALUES:" << endl;
      for (int i = 0; i < MP; i++)
      {
        for (int j = 0; j < NP; j++) x[j] = p[i][j];
        y[i] = func(x);
      }

// MINIMIZATION CALCULATION USING AMOEBA
      cout << endl << "MINIMIZATION" << endl;
      cout << "PARAMETER GUESSES AND CHI SQ VALUES:" << endl;
      const DP FTOL = 1.0e-2;
      NR::amoeba(p, y, FTOL, func, nfunc);

// REPORT THE TOTAL NUMBER OF EVALUATIONS
      cout << endl << "Number of function evaluations: " << nfunc << endl;

// REPORT THE FINAL SIMPLEX
      cout << endl << "PARAMETERS AT VERTICES OF FINAL SIMPLEX AND CHI SQ VALUES:" << endl;
      for (int i = 0; i < MP; i++)
      {
        for (int j = 0; j < NP; j++) cout << fixed << setprecision(5) << p[i][j] << tab;
        cout << "chi sq = " << fixed << setprecision(5) << y[i] << endl;
      }

// CALCULATE THE FINAL AVERAGE VALUES IN THE CENTER OF THE SIMPLEX
      cout << endl << "PARAMETERS AT CENTER OF FINAL SIMPLEX AND CHI SQ VALUE:" << endl;
      for (int j = 0; j < NP; j++)
      {
        DP sum = 0.0;
        for (int i = 0; i < MP; i++) sum += p[i][j];
        x[j] = sum / MP;
      }
      DP yfinal = func(x);

// SELECT PARAMETERS THAT GIVE LOWEST CHI SQ AS FINAL
      for (int i = 0; i < MP; i++)
      {
        if (y[i] < yfinal)
        {
          for (int j = 0; j < NP; j++) x[j] = p[i][j];
          yfinal = y[i];
        }
      }
      cout << endl << "FINAL BEST FIT PARAMETERS AND CHI SQ FROM FIT:" << endl;
      for (int j = 0; j < NP; j++) cout << fixed << setprecision(5) << x[j] << tab;
      cout << "chi sq = " << fixed << setprecision(5) << yfinal << endl;
      for (int j = 0; j < NP; j++) outfile0 << fixed << setprecision(6) << x[j] << tab;
      outfile0 << fixed << setprecision(6) << yfinal << endl;

// Output the final chi sq and the 1 sigma uncertainty from Kd
      cout << endl << "CHI-SQ VALUES WITH Kd's +/- 1 SIGMA" << endl;
      write_file = 0;   // Make sure the output files are NOT written yet
      write_errors = 1; // Records error data including when error is -1, +1, and 0
      error = -1;
      cout << endl << "MINUS ERROR CALCULATION:" << endl;
      chi_sq = func(x);
      cout << "CHI SQ = " << chi_sq << endl;
      error = +1;
      cout << endl << "PLUS ERROR CALCULATION:" << endl;
      chi_sq = func(x);
      cout << "CHI SQ = " << chi_sq << endl;
      error = 0;
      write_file = 1;   // Make sure the output files ARE written now
      cout << endl << "BEST FIT CALCULATION:" << endl;
      chi_sq = func(x);
      cout << "CHI SQ = " << chi_sq << endl;
    }
  }
  else if (fit_cond == 2)
  {
    for (int i = 0; i < n; i++) fit_elem[i] = fit_hse[i];

    NP = 4;
    MP = NP + 1;

// ASSIGN THE FITTED PARAMETERS
// Must match the list in input_parameters.h
    Vec_DP x(NP);
    x[0] = dist_0_S;
    x[1] = dist_0_C;
    x[2] = matte_P_factor;
    x[3] = time_magma;

    if (fit == 0)
    {
      // Output the chi sq and the 1 sigma uncertainty from Kd
      cout << "ONE-OFF CALCULATION" << endl;
      cout << "CHI-SQ VALUES WITH Kd's +/- 1 SIGMA" << endl;
      write_file = 0;   // Make sure the output files are NOT written yet
      write_errors = 1; // Records error data including when error is -1, +1, and 0
      error = -1;
      cout << endl << "MINUS ERROR CALCULATION:" << endl;
      cout << "CHI SQ = " << func(x) << endl;
      error = +1;
      cout << endl << "PLUS ERROR CALCULATION:" << endl;
      cout << "CHI SQ = " << func(x) << endl;
      error = 0;
      write_file = 1;   // Make sure the output files ARE written now
      cout << endl << "BEST FIT CALCULATION:" << endl;
      chi_sq = func(x);
      cout << "CHI SQ = " << chi_sq << endl;
    }
    // CODE PERFORMS A LEAST SQUARES REGRESSION
    else if (fit == 1)
    {
      Vec_DP y(MP);
      Mat_DP p(MP, NP);
      // Number of Amoeba evaluations
      int nfunc;
//__________________________________________________________________________________________________________
// FOR FITTING CALCULATIONS, READ PARAM.txt for initial parameter guesses
// reads MP guesses of NP fit parameters (MP=NP+1) for AMOEBA
      ifstream paramfile("PARAM_HSE.txt");
      if (paramfile.is_open())
      {
        for (int i = 0; i < NP; i++)
        {
          for (int j = 0; j < MP; j++) {
            paramfile >> p[j][i];
          }
          getline(paramfile, line);
        }
      }
      paramfile.close();

// NP is the number of fitting parameters and MP is the number of initial guess arrays
// p is the matrix of initial guesses read from the input file (PARAM.txt)
// so x is an array of one set of initial guesses
// and y is the corresponding chi sq
      write_file = 0;   // Make sure the output files are NOT written yet
      write_errors = 0; // Does not record error data
      error = 0;        // Uses the nominal Kd values

// REPORT THE INITIAL SIMPLEX
      cout << endl << "INITIAL MATRIX OF MINIMIZATION PARAMETERS AND CHI SQ VALUES:" << endl;
      for (int i = 0; i < MP; i++)
      {
        for (int j = 0; j < NP; j++) x[j] = p[i][j];
        y[i] = func(x);
      }

// MINIMIZATION CALCULATION USING AMOEBA
      cout << endl << "MINIMIZATION" << endl;
      cout << "PARAMETER GUESSES AND CHI SQ VALUES:" << endl;
      const DP FTOL = 1.0e-3;
      NR::amoeba2(p, y, FTOL, func, nfunc);

// REPORT THE TOTAL NUMBER OF EVALUATIONS
      cout << endl << "Number of function evaluations: " << nfunc << endl;

// REPORT THE FINAL SIMPLEX
      cout << endl << "PARAMETERS AT VERTICES OF FINAL SIMPLEX AND CHI SQ VALUES:" << endl;
      for (int i = 0; i < MP; i++)
      {
        for (int j = 0; j < NP; j++) cout << fixed << setprecision(5) << p[i][j] << tab;
        cout << "chi sq = " << fixed << setprecision(5) << y[i] << endl;
      }

// CALCULATE THE FINAL AVERAGE VALUES IN THE CENTER OF THE SIMPLEX
      cout << endl << "PARAMETERS AT CENTER OF FINAL SIMPLEX AND CHI SQ VALUE:" << endl;
      for (int j = 0; j < NP; j++)
      {
        DP sum = 0.0;
        for (int i = 0; i < MP; i++) sum += p[i][j];
        x[j] = sum / MP;
      }
      DP yfinal = func(x);

// SELECT PARAMETERS THAT GIVE LOWEST CHI SQ AS FINAL
      for (int i = 0; i < MP; i++)
      {
        if (y[i] < yfinal)
        {
          for (int j = 0; j < NP; j++) x[j] = p[i][j];
          yfinal = y[i];
        }
      }
      cout << endl << "FINAL BEST FIT PARAMETERS AND CHI SQ FROM FIT:" << endl;
      for (int j = 0; j < NP; j++) cout << fixed << setprecision(5) << x[j] << tab;
      cout << "chi sq = " << fixed << setprecision(5) << yfinal << endl;
      for (int j = 0; j < NP; j++) outfile0 << fixed << setprecision(6) << x[j] << tab;
      outfile0 << fixed << setprecision(6) << yfinal << endl;

// Output the final chi sq and the 1 sigma uncertainty from Kd
      cout << endl << "CHI-SQ VALUES WITH Kd's +/- 1 SIGMA" << endl;
      write_file = 0;   // Make sure the output files are NOT written yet
      write_errors = 1; // Records error data including when error is -1, +1, and 0
      error = -1;
      cout << endl << "MINUS ERROR CALCULATION:" << endl;
      chi_sq = func(x);
      cout << "CHI SQ = " << chi_sq << endl;
      error = +1;
      cout << endl << "PLUS ERROR CALCULATION:" << endl;
      chi_sq = func(x);
      cout << "CHI SQ = " << chi_sq << endl;
      error = 0;
      write_file = 1;   // Make sure the output files ARE written now
      cout << endl << "BEST FIT CALCULATION:" << endl;
      chi_sq = func(x);
      cout << "CHI SQ = " << chi_sq << endl;
    }
  }
  // Moon ver.
  // FIT CONDITION = 3, fit conditions for Moon (fit_moon_1): Fe, Si, lunar core mass
  if (fit_cond == 3)
  {
    for (int i = 0; i < n; i++) fit_elem[i] = fit_moon_1[i];


    if (dist_comp == 1) NP = 2;
    else if (dist_comp == 2) NP = 6;
    else if (dist_comp == 3) NP = 5;
    else if (dist_comp == 4) NP = 6;
    else if (dist_comp == 5) NP = 5;
    MP = NP + 1;

// ASSIGN THE FITTED PARAMETERS
// Must match the list in input_parameters.h
    Vec_DP x(NP);
    if (dist_comp == 1)  // CONSTANT
    {
      x[0] = p_factor;
      x[1] = xfe_met_2;
    }
    else if (dist_comp == 2)  // STEP
    {
      x[0] = p_factor;
      x[1] = dist_1;
      x[2] = dist_2;
      x[3] = xsi_met_0;
      x[4] = xfe_met_2;
      x[5] = xfe_met_3;
    }
    else if (dist_comp == 3)  // GRADIENT
    {
      x[0] = p_factor;
      x[1] = dist_1;
      x[2] = dist_2;
      x[3] = xsi_met_0;
      x[4] = xfe_met_2;
    }
    else if (dist_comp == 4)  // ERF
    {
      x[0] = p_factor;
      x[1] = dist_1;
      x[2] = dist_2;
      x[3] = xsi_met_0;
      x[4] = erf_width_1;
      x[5] = erf_width_2;
    }
    else if (dist_comp == 5)  // MOON GRADIENT (Moon ver.)
    {
      x[0] = p_factor;
      x[1] = moon_xfe_met_1;
      x[2] = moon_xfe_met_2;
      x[3] = moon_dist_1;
      x[4] = moon_dist_2;
    }

//__________________________________________________________________________________________________________
// CODE PERFORMS A ONE-OFF CALCULATION
    if (fit == 0)
    {
      // Output the chi sq and the 1 sigma uncertainty from Kd
      cout << "ONE-OFF CALCULATION" << endl;
      cout << "CHI-SQ VALUES WITH Kd's +/- 1 SIGMA" << endl;
      write_file = 0;   // Make sure the output files are NOT written yet
      write_errors = 1; // Records error data including when error is -1, +1, and 0
      error = -1;
      cout << endl << "MINUS ERROR CALCULATION:" << endl;
      cout << "CHI SQ = " << func(x) << endl;
      error = +1;
      cout << endl << "PLUS ERROR CALCULATION:" << endl;
      cout << "CHI SQ = " << func(x) << endl;
      error = 0;
      write_file = 1;   // Make sure the output files ARE written now
      cout << endl << "BEST FIT CALCULATION:" << endl;
      chi_sq = func(x);
      cout << "CHI SQ = " << chi_sq << endl;
    }
// CODE PERFORMS A LEAST SQUARES REGRESSION
    else if (fit == 1)
    {
      Vec_DP y(MP);
      Mat_DP p(MP, NP);
      // Number of Amoeba evaluations
      int nfunc;
//__________________________________________________________________________________________________________
// FOR FITTING CALCULATIONS, READ PARAM.txt for initial parameter guesses
// reads MP guesses of NP fit parameters (MP=NP+1) for AMOEBA
      ifstream paramfile("PARAM_MSE.txt");
      if (paramfile.is_open())
      {
        for (int i = 0; i < NP; i++)
        {
          for (int j = 0; j < MP; j++) {
            paramfile >> p[j][i];
          }
          getline(paramfile, line);
        }
      }
      paramfile.close();

// NP is the number of fitting parameters and MP is the number of initial guess arrays
// p is the matrix of initial guesses read from the input file (PARAM.txt)
// so x is an array of one set of initial guesses
// and y is the corresponding chi sq
      write_file = 0;   // Make sure the output files are NOT written yet
      write_errors = 0; // Does not record error data
      error = 0;        // Uses the nominal Kd values

// REPORT THE INITIAL SIMPLEX
      cout << endl << "INITIAL MATRIX OF MINIMIZATION PARAMETERS AND CHI SQ VALUES:" << endl;
      for (int i = 0; i < MP; i++)
      {
        for (int j = 0; j < NP; j++) x[j] = p[i][j];
        y[i] = func(x);
      }

// MINIMIZATION CALCULATION USING AMOEBA
      cout << endl << "MINIMIZATION" << endl;
      cout << "PARAMETER GUESSES AND CHI SQ VALUES:" << endl;
      const DP FTOL = 1.0e-3;
      NR::amoeba2(p, y, FTOL, func, nfunc);

// REPORT THE TOTAL NUMBER OF EVALUATIONS
      cout << endl << "Number of function evaluations: " << nfunc << endl;

// REPORT THE FINAL SIMPLEX
      cout << endl << "PARAMETERS AT VERTICES OF FINAL SIMPLEX AND CHI SQ VALUES:" << endl;
      for (int i = 0; i < MP; i++)
      {
        for (int j = 0; j < NP; j++) cout << fixed << setprecision(5) << p[i][j] << tab;
        cout << "chi sq = " << fixed << setprecision(5) << y[i] << endl;
      }

// CALCULATE THE FINAL AVERAGE VALUES IN THE CENTER OF THE SIMPLEX
      cout << endl << "PARAMETERS AT CENTER OF FINAL SIMPLEX AND CHI SQ VALUE:" << endl;
      for (int j = 0; j < NP; j++)
      {
        DP sum = 0.0;
        for (int i = 0; i < MP; i++) sum += p[i][j];
        x[j] = sum / MP;
      }
      DP yfinal = func(x);

// SELECT PARAMETERS THAT GIVE LOWEST CHI SQ AS FINAL
      for (int i = 0; i < MP; i++)
      {
        if (y[i] < yfinal)
        {
          for (int j = 0; j < NP; j++) x[j] = p[i][j];
          yfinal = y[i];
        }
      }
      cout << endl << "FINAL BEST FIT PARAMETERS AND CHI SQ FROM FIT:" << endl;
      for (int j = 0; j < NP; j++) cout << fixed << setprecision(5) << x[j] << tab;
      cout << "chi sq = " << fixed << setprecision(5) << yfinal << endl;
      for (int j = 0; j < NP; j++) outfile0 << fixed << setprecision(6) << x[j] << tab;
      outfile0 << fixed << setprecision(6) << yfinal << endl;

// Output the final chi sq and the 1 sigma uncertainty from Kd
      cout << endl << "CHI-SQ VALUES WITH Kd's +/- 1 SIGMA" << endl;
      write_file = 0;   // Make sure the output files are NOT written yet
      write_errors = 1; // Records error data including when error is -1, +1, and 0
      error = -1;
      cout << endl << "MINUS ERROR CALCULATION:" << endl;
      chi_sq = func(x);
      cout << "CHI SQ = " << chi_sq << endl;
      error = +1;
      cout << endl << "PLUS ERROR CALCULATION:" << endl;
      chi_sq = func(x);
      cout << "CHI SQ = " << chi_sq << endl;
      error = 0;
      write_file = 1;   // Make sure the output files ARE written now
      cout << endl << "BEST FIT CALCULATION:" << endl;
      chi_sq = func(x);
      cout << "CHI SQ = " << chi_sq << endl;
    }
  }

// Close output files
  outfile0.close(); outfile1.close(); outfile2.close(); outfile3.close(); outfile4.close(); outfile5.close(); outfile6.close(); outfile7.close(); outfile8.close(); outfile9.close(); outfile10.close(); //outfile11.close(), //noMarsNumber
  outfile12.close(); outfile13.close(); outfile14.close();// testfile.close();

//   system("osascript -e 'display dialog \"Finished!\" buttons \"OK\"'");

  return 0;
}

//______________________________________________________________________________________________________________________________________________
// FUNCTION FOR CORE FORMATION CALCULATION AND CHI SQUARED DETERMINATION CALLED BY AMOEBA
DP func(Vec_I_DP &x)
{
  int NP = x.size();
  /*
     GUIDE TO GLOBAL VARIABLES THAT MUST BE SET:

     EQUILIBRATION STATE:
     double e_P[i][j]                    pressure of metal-silicate equilibration in the ith embryo during the jth accretion event
     double e_T[i][j]                    temperature of metal-silicate equilibration in the ith embryo during the jth accretion event
     double e_fO2[i][j]                  oxygen fugacity of metal-silicate equilibration in the ith embryo during the jth accretion event

     MANTLE COMPOSITION:
     double e_moles[i][x][j]             number of moles of the xth silicate species where x varies from 0 to n-1 in the ith embryo during the jth accretion event
     double e_moles_H2O[i][j]            number of moles of water in the ith embryo during the jth accretion event
     double e_moles_S_mantle[i][j]       number of moles of sulfur in the ith embryo during the jth accretion event

     CORE COMPOSITION:
     double e_moles[i][x][j]             number of moles of the xth metal species where x varies from n to n+m-1 in the ith embryo during the jth accretion event
     double e_moles_oxygen[i][j]         number of moles of oxygen in the ith embryo during the jth accretion event
     double e_moles_H[i][j]              number of moles of hydrogen in the ith embryo during the jth accretion event
     double e_moles_S_core[i][j]         number of moles of sulfur in the ith embryo during the jth accretion event

     CORE FORMING LIQUID COMPOSITION:
     double incr_moles[i][x][j]          number of moles of the xth metal species where x varies from n to n+m-1 in the ith embryo during the jth accretion event
     double incr_moles_oxygen[i][j]      number of moles of oxygen in the ith embryo during the jth accretion event
     double incr_moles_H[i][j]           number of moles of hydrogen in the ith embryo during the jth accretion event
     double incr_moles_S[i][j]           number of moles of sulfur in the ith embryo during the jth accretion event
   */

// LOCAL VARIABLES USED TO PASS GLOBAL VARIABLES TO SUBFUNCTIONS OF weigh_comp and react
// They lack two dimensions (embryo and accretion event) compared to the equivalent global variable
  double P(0.0), T(0.0), fo2(0.0); // Pressure, temperature and fO2 of metal silicate equilibration event
  int zero_xfe_liq;
  double t(0.0);
  //int no_1(0), no_2(0);
  double mass_1(0.0), mass_2(0.0);

  double moles[n + m + 3] = {0.0}; // Similar to e_moles, moles of each species in silicate liquid (0 -- n-1) and metal liquid (n -- n+m-1)
  // but there are three extra elements to pass the oxygen in metal, fO2, and the presence of metal, the lack of which is indicated by 0.0
  double moles_H2O(0.0);     // Similar to e_moles_H2O, water in silicate liquid
  double moles_S_mantle(0.0); // Similar to e_moles_S_mantle, sulfur in silicate liquid
  double moles_C_mantle(0.0); // Similar to e_moles_C_mantle, sulfur in silicate liquid
  double moles_oxygen(0.0);  // Similar to e_moles_oxygen, oxygen in metal liquid
  double moles_H(0.0);       // Similar to e_moles_H, hydrogen in metal liquid
  double moles_S_core(0.0);  // Similar to e_moles_S_core, sulfur in metal liquid
  double moles_C_core(0.0);  // Similar to e_moles_C_core, sulfur in metal liquid

  double mass_factor = 0.0;           // Relation between no. of moles and mass = 1.0 kg for a given composition

  double no_1_moles[n + m] = {0.0};
  double no_1_moles_H2O = 0.0;
  double no_1_moles_S_mantle = 0.0;
  double no_1_moles_C_mantle = 0.0;
  double no_1_moles_oxygen = 0.0;
  double no_1_moles_H = 0.0;
  double no_1_moles_S_core = 0.0;
  double no_1_moles_C_core = 0.0;

  double no_1_Fe_sil_isotope[3] = {0.0};
  double no_1_Fe_met_isotope[3] = {0.0};

  double no_2_moles[n + m] = {0.0};
  double no_2_moles_H2O = 0.0;
  double no_2_moles_S_mantle = 0.0;
  double no_2_moles_C_mantle = 0.0;
  double no_2_moles_oxygen = 0.0;
  double no_2_moles_H = 0.0;
  double no_2_moles_S_core = 0.0;
  double no_2_moles_C_core = 0.0;
  int no_2_met = 1;          // Tracks whether there is metal present in the projectile; default assumes there is metal present (Value of "1" means that metal is present; when not, value is changed to "0")

  double no_2_Fe_sil_isotope[3] = {0.0};
  double no_2_Fe_met_isotope[3] = {0.0};

  double equil_moles[n + m] = {0.0};
  double equil_moles_H2O = 0.0;
  double equil_moles_S_mantle = 0.0;
  double equil_moles_C_mantle = 0.0;
  double equil_moles_oxygen = 0.0;
  double equil_moles_H = 0.0;
  double equil_moles_S_core = 0.0;
  double equil_moles_C_core = 0.0;

  double equil_Fe_sil_isotope[3] = {0.0};
  double equil_Fe_met_isotope[3] = {0.0};

  double noneq_moles[n + m] = {0.0};
  double noneq_moles_H2O = 0.0;
  double noneq_moles_S_mantle = 0.0;
  double noneq_moles_C_mantle = 0.0;
  double noneq_moles_oxygen = 0.0;
  double noneq_moles_H = 0.0;
  double noneq_moles_S_core = 0.0;
  double noneq_moles_C_core = 0.0;

  double noneq_Fe_sil_isotope[3] = {0.0};
  double noneq_Fe_met_isotope[3] = {0.0};

  double no1_fe = 0.0;
  double no2_fe = 0.0;

// OTHER VARIABLES
  double xsil(0.0);          // Fraction of embryo's mantle that equilibrates with impactor's core

//__________________________________________________________________________________________________________
// Re-set following global variables to zero for each instance of func
  mass_factor = 0.0;         // Scaling factor that is a unit Earth mass divided by the total number of starting moles
  // used to scale the starting molar composition to the mass of the body

  // Re-set these counters and value trackers
  n_pl = 0;
  n0_pl = 0;
  n_em = 0;
  xsil_pl = 0.0;
  xsil_em = 0.0;
  xsil_pl_max = 0.0;
  xsil_pl_min = 0.999;
  xsil_em_max = 0.0;
  xsil_em_min = 0.999;

  chi_sq_multiplier = 0;

// Re-set variables associated with every body
  for (int i = 0; i < n_bodies; i++) {
    e_coll[i] = 0;
    b_status[i] = 0;
  }

// Re-set variables associated with every embryo
  for (int i = 0; i < n_emb; i++)
  {
    e_SN_mass[i] = 0.0;
    e_em_mass[i] = 0.0;
    e_pl_mass[i] = 0.0;
    GI[i] = 0.0;
    last_GI[i] = 0.0;

// Re-set variables associated with every accretion event on every embryo
    for (int j = 0; j < n_imps; j++)
    {
      e_time[i][j] = 0.0;
      e_mass[i][j] = 0.0;
      e_P[i][j] = 0.0;
      e_T[i][j] = 0.0;
      e_fO2[i][j] = 0.0;
      e_moles_oxygen[i][j] = 0.0;
      e_moles_S_core[i][j] = 0.0;
      e_moles_S_mantle[i][j] = 0.0;
      e_moles_C_core[i][j] = 0.0;
      e_moles_C_mantle[i][j] = 0.0;
      e_moles_H2O[i][j] = 0.0;
      e_moles_H[i][j] = 0.0;
      i_moles_oxygen[i][j] = 0.0;
      i_moles_S_core[i][j] = 0.0;
      i_moles_S_mantle[i][j] = 0.0;
      i_moles_C_core[i][j] = 0.0;
      i_moles_C_mantle[i][j] = 0.0;
      i_moles_H2O[i][j] = 0.0;
      i_moles_H[i][j] = 0.0;
      s_moles_oxygen[i][j] = 0.0;
      s_moles_S_core[i][j] = 0.0;
      s_moles_S_mantle[i][j] = 0.0;
      s_moles_C_core[i][j] = 0.0;
      s_moles_C_mantle[i][j] = 0.0;
      s_moles_H2O[i][j] = 0.0;
      s_moles_H[i][j] = 0.0;
      incr_moles_oxygen[i][j] = 0.0;
      incr_moles_H[i][j] = 0.0;
      incr_moles_S[i][j] = 0.0;
      incr_moles_C[i][j] = 0.0;
// Re-set variables associated with every species, every accretion event and every embryo
      for (int k = 0; k < n + m; k++)
      {
        e_moles[i][k][j] = 0;
        i_moles[i][k][j] = 0;
        s_moles[i][k][j] = 0;
        incr_moles[i][k][j] = 0;
      }
    }
  }

//__________________________________________________________________________________________________________
// Define input parameters for core formation function "func" using input_parameters.h
// *****THIS MUST BE REVISED AS NECESSARY WHEN FITTING PARAMETERS ARE CHANGED
  if (fit_cond == 0)
  {
    double p_input[12];
    fit_mse_parameters(p_input, dist_comp, fit, p_factor, dist_1, dist_2, xsi_met_0, xfe_met_2, xfe_met_3, erf_width_1, erf_width_2, 
                       moon_xfe_met_1, moon_xfe_met_2, moon_dist_1, moon_dist_2, x); // Moon ver.
    p_factor = p_input[0];
    dist_1 = p_input[1];
    dist_2 = p_input[2];
    xsi_met_0 = p_input[3];
    xfe_met_2 = p_input[4];
    xfe_met_3 = p_input[5];
    erf_width_1 = p_input[6];
    erf_width_2 = p_input[7];
    moon_xfe_met_1 = p_input[8];
    moon_xfe_met_2 = p_input[9];
    moon_dist_1 = p_input[10];
    moon_dist_2 = p_input[11];
  }
  else if (fit_cond == 1)
  {
    dist_0[3] = x[0];
  }
  else if (fit_cond == 2)
  {
    double p_input[4];
    fit_hse_parameters(p_input, dist_0_S, dist_0_C, matte_P_factor, time_magma, x);
    dist_0_S = p_input[0];
    dist_0_C = p_input[1];
    matte_P_factor = p_input[2];
    time_magma = p_input[3];
  }
  else if (fit_cond == 3) // Moon ver.
  {
        double p_input[12];
    fit_mse_parameters(p_input, dist_comp, fit, p_factor, dist_1, dist_2, xsi_met_0, xfe_met_2, xfe_met_3, erf_width_1, erf_width_2, 
                       moon_xfe_met_1, moon_xfe_met_2, moon_dist_1, moon_dist_2, x); // Moon ver.
    p_factor = p_input[0];
    dist_1 = p_input[1];
    dist_2 = p_input[2];
    xsi_met_0 = p_input[3];
    xfe_met_2 = p_input[4];
    xfe_met_3 = p_input[5];
    erf_width_1 = p_input[6];
    erf_width_2 = p_input[7];
    moon_xfe_met_1 = p_input[8];
    moon_xfe_met_2 = p_input[9];
    moon_dist_1 = p_input[10];
    moon_dist_2 = p_input[11];
  }
//__________________________________________________________________________________________________________
// IF UNREALISTIC VALUES ARE PROPOSED BY AMOEBA, CHI SQ IS ARBITARILY INCREASED BY 10% WITHOUT PERFORMING A CALCULATION
// *****THIS MUST BE REVISED AS NECESSARY WHEN FITTING PARAMETERS ARE CHANGED
  if (fit == 1)
  {
    cout << "fit_cond: " << fit_cond << " dist_comp: " << dist_comp << endl;
    if ((fit_cond == 0 && dist_comp == 1 && (x[0] < 0.0 || x[0] > 1.0 || x[1] < 0.0 || x[1] > 1.0)) ||
        (fit_cond == 0 && dist_comp == 2 && (x[0] < 0.0 || x[0] > 1.0 || x[1] > x[2] || x[2] > dist_3 || x[3] < 0.0 || x[3] > 1.0 || x[4] < x[5] || x[4] > 1.0 || x[5] < 0.0)) ||
        (fit_cond == 0 && dist_comp == 3 && (x[0] < 0.0 || x[0] > 1.0 || x[1] > x[2] || x[2] > dist_3 || x[3] < 0.0 || x[3] > 1.0 || x[4] < 0.0 || x[4] > 1.0)) ||
        (fit_cond == 0 && dist_comp == 4 && (x[0] < 0.0 || x[0] > 1.0 || x[1] > x[2] || x[2] > dist_3 || x[3] < 0.0 || x[3] > 1.0 || x[4] < 0.0 || x[5] < 0.0)) ||
        (fit_cond == 0 && dist_comp == 5 && (x[0] < 0.0 || x[0] > 1.0 || x[1] < 0.0 || x[1] > 1.0 || x[2] < 0.0 || x[2] > 1.0 || x[3] > x[4])) || // Moon_ver
        (fit_cond == 2 && (x[0] > dist_3 || x[1] > dist_3 || x[2] < 0.0 || x[2] > 1.0 || x[3] < 0.0)) ||
        (fit_cond == 3 && dist_comp == 5 && (x[0] < 0.0 || x[0] > 1.0 || x[1] < 0.0 || x[1] > 1.0 || x[2] < 0.0 || x[2] > 1.0 || x[3] > x[4]))) // Moon_ver

    {
// In the case of unrealistic values, chi sq is arbitarily increased by 10%
      chi_sq = chi_sq * 1.1;

// Write the input parameters and bad chi sq to screen
      for (int i = 0; i < NP; i++) cout << fixed << setprecision(5) << x[i] << tab;
      cout << "chi sq*= " << fixed << setprecision(5) << chi_sq << endl;

      return chi_sq;
    }
    else
    {
// Write the input parameters to screen and to outfile0 (ATTEMPTS.txt)
      for (int i = 0; i < NP; i++)
      {
        cout << fixed << setprecision(5) << x[i] << tab;
      }
    }
  }

//__________________________________________________________________________________________________________
// OPEN infile3 ("NB_data.txt") THAT LISTS N-BODY SIMULATION HISTORY
// In order to be read properly, there must be a carriage return at the end of the final line
  string simulation_no;
  ifstream infile3("NB_data.txt"); // declaration of 3rd input file
  if (infile3.is_open()) getline (infile3, simulation_no); // reads simulation # as string
  infile3 >> moon_body; getline (infile3, line);
  // infile3 >> mars_planet;   getline (infile3, line); //noMarsNumber
  infile3 >> mass_adjust;          getline (infile3, line);// Adjusts masses to give earth one Earth mass
  infile3 >> max_pl;               getline (infile3, line);


//__________________________________________________________________________________________________________
// Write final (i.e. optimized) simulation details to output files:
  if (write_file == 1)
  {
// outfile1 is MANTLE.txt
    write_outfile_header(outfile1, simulation_no, dist_comp, p_factor, em_equil, pf_em, dist_1, dist_2, dist_3, erf_width_1, erf_width_2, xfe_met_1, xfe_met_2, xfe_met_3,
                         xsi_met_0, wtci_H2O, mass_adjust, xmet_pl, xmet_em, x_mantle, em_xsil, t_overstep, el_type, el_sym, dist_0, dist_0_S, dist_0_C, matte_HSE,
                         matte_P_factor, pf,
                         moon_xfe_met_1, moon_xfe_met_2, moon_dist_1, moon_dist_2, moon_disk_start, moon_disk_end); // Moon_version 
// outfile2 is CORE.txt
    write_outfile_header(outfile2, simulation_no, dist_comp, p_factor, em_equil, pf_em, dist_1, dist_2, dist_3, erf_width_1, erf_width_2, xfe_met_1, xfe_met_2, xfe_met_3,
                         xsi_met_0, wtci_H2O, mass_adjust, xmet_pl, xmet_em, x_mantle, em_xsil, t_overstep, el_type, el_sym, dist_0, dist_0_S, dist_0_C, matte_HSE,
                         matte_P_factor, pf,
                         moon_xfe_met_1, moon_xfe_met_2, moon_dist_1, moon_dist_2, moon_disk_start, moon_disk_end); // Moon_version 

// outfile3 is CORE_INCRE.txt
    write_outfile_header(outfile3, simulation_no, dist_comp, p_factor, em_equil, pf_em, dist_1, dist_2, dist_3, erf_width_1, erf_width_2, xfe_met_1, xfe_met_2, xfe_met_3,
                         xsi_met_0, wtci_H2O, mass_adjust, xmet_pl, xmet_em, x_mantle, em_xsil, t_overstep, el_type, el_sym, dist_0, dist_0_S, dist_0_C, matte_HSE,
                         matte_P_factor, pf,
                         moon_xfe_met_1, moon_xfe_met_2, moon_dist_1, moon_dist_2, moon_disk_start, moon_disk_end); // Moon_version 

// outfile4 is EMBRYO_MANTLE.txt
    write_outfile_header(outfile4, simulation_no, dist_comp, p_factor, em_equil, pf_em, dist_1, dist_2, dist_3, erf_width_1, erf_width_2, xfe_met_1, xfe_met_2, xfe_met_3,
                         xsi_met_0, wtci_H2O, mass_adjust, xmet_pl, xmet_em, x_mantle, em_xsil, t_overstep, el_type, el_sym, dist_0, dist_0_S, dist_0_C, matte_HSE,
                         matte_P_factor, pf,
                         moon_xfe_met_1, moon_xfe_met_2, moon_dist_1, moon_dist_2, moon_disk_start, moon_disk_end); // Moon_version 

// outfile5 is EMBRYO_CORE.txt
    write_outfile_header(outfile5, simulation_no, dist_comp, p_factor, em_equil, pf_em, dist_1, dist_2, dist_3, erf_width_1, erf_width_2, xfe_met_1, xfe_met_2, xfe_met_3,
                         xsi_met_0, wtci_H2O, mass_adjust, xmet_pl, xmet_em, x_mantle, em_xsil, t_overstep, el_type, el_sym, dist_0, dist_0_S, dist_0_C, matte_HSE,
                         matte_P_factor, pf,
                         moon_xfe_met_1, moon_xfe_met_2, moon_dist_1, moon_dist_2, moon_disk_start, moon_disk_end); // Moon_version 

// outfile7 is IMPACTOR_MANTLE.txt
    write_outfile_header(outfile6, simulation_no, dist_comp, p_factor, em_equil, pf_em, dist_1, dist_2, dist_3, erf_width_1, erf_width_2, xfe_met_1, xfe_met_2, xfe_met_3,
                         xsi_met_0, wtci_H2O, mass_adjust, xmet_pl, xmet_em, x_mantle, em_xsil, t_overstep, el_type, el_sym, dist_0, dist_0_S, dist_0_C, matte_HSE,
                         matte_P_factor, pf,
                         moon_xfe_met_1, moon_xfe_met_2, moon_dist_1, moon_dist_2, moon_disk_start, moon_disk_end); // Moon_version 

// outfile7 is IMPACTOR_CORE.txt
    write_outfile_header(outfile7, simulation_no, dist_comp, p_factor, em_equil, pf_em, dist_1, dist_2, dist_3, erf_width_1, erf_width_2, xfe_met_1, xfe_met_2, xfe_met_3,
                         xsi_met_0, wtci_H2O, mass_adjust, xmet_pl, xmet_em, x_mantle, em_xsil, t_overstep, el_type, el_sym, dist_0, dist_0_S, dist_0_C, matte_HSE,
                         matte_P_factor, pf,
                         moon_xfe_met_1, moon_xfe_met_2, moon_dist_1, moon_dist_2, moon_disk_start, moon_disk_end); // Moon_version 

// outfile8 is SOURCE_MANTLE.txt
    write_outfile_header(outfile8, simulation_no, dist_comp, p_factor, em_equil, pf_em, dist_1, dist_2, dist_3, erf_width_1, erf_width_2, xfe_met_1, xfe_met_2, xfe_met_3,
                         xsi_met_0, wtci_H2O, mass_adjust, xmet_pl, xmet_em, x_mantle, em_xsil, t_overstep, el_type, el_sym, dist_0, dist_0_S, dist_0_C, matte_HSE,
                         matte_P_factor, pf,
                         moon_xfe_met_1, moon_xfe_met_2, moon_dist_1, moon_dist_2, moon_disk_start, moon_disk_end); // Moon_version 

// outfile9 is SOURCE_CORE.txt
    write_outfile_header(outfile9, simulation_no, dist_comp, p_factor, em_equil, pf_em, dist_1, dist_2, dist_3, erf_width_1, erf_width_2, xfe_met_1, xfe_met_2, xfe_met_3,
                         xsi_met_0, wtci_H2O, mass_adjust, xmet_pl, xmet_em, x_mantle, em_xsil, t_overstep, el_type, el_sym, dist_0, dist_0_S, dist_0_C, matte_HSE,
                         matte_P_factor, pf,
                         moon_xfe_met_1, moon_xfe_met_2, moon_dist_1, moon_dist_2, moon_disk_start, moon_disk_end); // Moon_version 

// outfile10 is FINAL.txt
    write_outfile_header(outfile10, simulation_no, dist_comp, p_factor, em_equil, pf_em, dist_1, dist_2, dist_3, erf_width_1, erf_width_2, xfe_met_1, xfe_met_2, xfe_met_3,
                         xsi_met_0, wtci_H2O, mass_adjust, xmet_pl, xmet_em, x_mantle, em_xsil, t_overstep, el_type, el_sym, dist_0, dist_0_S, dist_0_C, matte_HSE,
                         matte_P_factor, pf,
                         moon_xfe_met_1, moon_xfe_met_2, moon_dist_1, moon_dist_2, moon_disk_start, moon_disk_end); // Moon_version 

// outfile11 is FINAL_MARS.txt - not used in this version!!!! noMarsNumber
    // write_outfile_header(outfile11, simulation_no, dist_comp, p_factor, em_equil, pf_em, dist_1, dist_2, dist_3, erf_width_1, erf_width_2, xfe_met_1, xfe_met_2, xfe_met_3,
    //                      xsi_met_0, wtci_H2O, mass_adjust, xmet_pl, xmet_em, x_mantle, em_xsil, t_overstep, el_type, el_sym, dist_0, dist_0_S, dist_0_C, matte_HSE,
    //                      matte_P_factor, pf);
// outfile12 is COMPOSITION_Fe.txt
    write_outfile_header(outfile12, simulation_no, dist_comp, p_factor, em_equil, pf_em, dist_1, dist_2, dist_3, erf_width_1, erf_width_2, xfe_met_1, xfe_met_2, xfe_met_3,
                         xsi_met_0, wtci_H2O, mass_adjust, xmet_pl, xmet_em, x_mantle, em_xsil, t_overstep, el_type, el_sym, dist_0, dist_0_S, dist_0_C, matte_HSE,
                         matte_P_factor, pf,
                         moon_xfe_met_1, moon_xfe_met_2, moon_dist_1, moon_dist_2, moon_disk_start, moon_disk_end); // Moon_version 
// outfile13 is BODY_28.txt
    // write_outfile_header(outfile13, simulation_no, dist_comp, p_factor, em_equil, pf_em, dist_1, dist_2, dist_3, erf_width_1, erf_width_2, xfe_met_1, xfe_met_2, xfe_met_3,
    //                      xsi_met_0, wtci_H2O, mass_adjust, xmet_pl, xmet_em, x_mantle, em_xsil, t_overstep, el_type, el_sym, dist_0, dist_0_S, dist_0_C, matte_HSE,
    //                      matte_P_factor, pf);
    // DONT WRITE A HEADER FOR OUTFILE13

  }

//__________________________________________________________________________________________________________
// READ infile5 ("GI.txt") THAT LISTS GIANT IMPACT HISTORY
  string sim_title;
  ifstream infile5("GI.txt");
// Check to make sure GI.txt exists
  if (infile5.is_open())
    getline (infile5, sim_title); // Reads simulation # as string
  else
  {
    cout << endl << "GI.txt file is missing." << endl;
    exit(1);
  }

// Ensures the correct GI.txt file
  if (sim_title != simulation_no)
  {
    cout << "title: " << sim_title << " " << simulation_no << endl;
    cout << endl << "File 'GI.txt' is for the wrong simulation." << endl;
    exit(1);
  }
  getline (infile5, line);
  getline (infile5, line);
  getline (infile5, line);
  getline (infile5, line);

// Read late veneer parameters from GI.txt file
  int e_index;
  while (!infile5.eof())
  {
    infile5 >> e_index;
    infile5 >> t_last_GI[e_index] >> em_tmax[e_index] >> t0_LV[e_index] >> t_matte[e_index];
// Currently, only the Earth-like planet has an imposed late veneer
// The rest of the embryos have late veneers set in GI.txt and if unchanged by the user are very late, i.e. 5.0e8 yrs
// The late veneer begines at time_magma (Myr) after the last giant impact
//    if (e_index == earth_planet) {
//    if (t_last_GI[e_index] > 0.0)
//    {
    t0_LV[e_index] = t_last_GI[e_index] + time_magma * 1.e6;
//    }

    getline (infile5, line);
  }
  number_embryos = e_index;

//__________________________________________________________________________________________________________
// GIANT LOOP THROUGH THE NB_data.txt FILE, MERGING/DIFFERENTIATING AFTER EACH ACCRETION EVENT
// EACH ITERATION OF THE LOOP, THERE ARE DIFFERENT COLLIDING BODIES #1 and #2
  while (!infile3.eof())     // Continue until end of file is reached
  {
    infile3 >> t >> no_1 >> mass_1 >> no_2 >> mass_2; // read first 5 parameters from each line of data file
    getline (infile3, line);
    // if (no_1 == 28 && write_file == 1) cout << endl << "___________________________________________________________________" << endl << " 28 hit by " << no_2 << endl << endl;
    // if (no_1 == 27 && write_file == 1) cout << endl << "___________________________________________________________________" << endl << " 27 hit by " << no_2 << endl << endl;
    // if (no_1 == 37 && write_file == 1) cout << endl << "___________________________________________________________________" << endl << " 37 hit by " << no_2 << endl << endl;
    // if (no_1 == 403 && write_file == 1) cout << endl << "___________________________________________________________________" << endl << " 403 hit by " << no_2 << endl << endl;

    // HfW:
    // for body 1 and body 2--> check their radioactive_time_counter updated to time t from whichever time they were before
    // array for length of how many times each body is impacted 
    // each body has some hf that is decaying and needs its counter re-updated after each impact that resets the hf/w ratio in the body

  
        
// Body #1 should always be more massive than body #2
    if (mass_1 < mass_2)
    {
      cout << "Time = " << t << tab << "PROGRAM TERMINATED: FIRST BODY IS SMALLER THAN THE SECOND!" << endl;
      exit(1);
    }
// Masses are adjusted to end up with the embryo indexed moon_body with the same mass as Moon
    mass_1 = mass_1 * mass_adjust;
    mass_2 = mass_2 * mass_adjust;

//______________________________________________________________________
// IN THE CASE THAT THE COMPOSITION OF no_1 IS NOT DEFINED
        
    // if (no_1 == 12 && write_file == 1 && e_time[no_1][e_coll[no_1]]/1e6 < 0.37 && e_time[no_1][e_coll[no_1]]/1e6 > 0.337  ) cout << endl << "the time is: "<< e_time[no_1][e_coll[no_1]]/1e6 <<  " and this is e_fe_sil SPOT 1aaaaa BEFORE bstatusno1=0 = " << no_1_Fe_sil_isotope[1] <<  " and e_Fe_sil_isotope[no_1][1][e_coll[no_1]]:" << e_Fe_sil_isotope[no_1][1][e_coll[no_1]] << endl << "and here bstatus[12] = " << b_status[no_1] << endl;
     
   if (b_status[no_1] == 0)
   {
    // Initialize the accretion event counter
      e_coll[no_1] = 0;
      e_mass[no_1][e_coll[no_1]] = mass_1;
// GET INITIAL COMPOSITION OF BODY #1
// Determine the fraction of Fe and Si of the body that is initially metal due to heliocentric distance
      xfe_met = Fe_metalfraction(dist_comp, dist[no_1], dist_1, dist_2, dist_3, xfe_met_1, xfe_met_2, xfe_met_3, erf_width_1, erf_width_2,
                                 moon_xfe_met_1, moon_xfe_met_2, moon_dist_1,moon_dist_2, moon_disk_start, moon_disk_end); // Moon_version 

      xsi_met = Si_metalfraction(dist_comp, dist[no_1], dist_1, dist_2, dist_3, xsi_met_0, erf_width_1, erf_width_2);

      starting_molar_proportions(dist[no_1],
                                 no_1_moles, no_1_moles_H2O, no_1_moles_S_mantle, no_1_moles_C_mantle,
                                 no_1_moles_oxygen, no_1_moles_H, no_1_moles_S_core, no_1_moles_C_core);
      
      starting_isotopic_proportions(no_1_Fe_sil_isotope, no_1_Fe_met_isotope); //GN
      
      mass_factor = get_mass_factor(no_1_moles, no_1_moles_H2O, no_1_moles_S_mantle, no_1_moles_C_mantle,
                                    no_1_moles_oxygen, no_1_moles_H, no_1_moles_S_core, no_1_moles_C_core);
// DIFFERENTIATE BODY #1
// FULLY OXIDIZED COMPOSITION - therefore no metal, no differentiation, all S in silicate
      if (xfe_met < 1e-20)
      { 
         // if (no_1 == 12 && write_file == 1 && e_time[no_1][e_coll[no_1]]/1e6 < 0.37 && e_time[no_1][e_coll[no_1]]/1e6 > 0.337  ) cout << endl << "the time is: "<< e_time[no_1][e_coll[no_1]]/1e6 <<  " and this is e_fe_sil SPOT 1 8a = " << no_1_Fe_sil_isotope[1] <<  " and e_Fe_sil_isotope[no_1][1][e_coll[no_1]]:" << e_Fe_sil_isotope[no_1][1][e_coll[no_1]] << endl;
      
// Calculate the metal-free composition
        undifferentiated_composition(mass_1, mass_factor,
                                     no_1_moles, no_1_moles_H2O, no_1_moles_S_mantle, no_1_moles_C_mantle,
                                     no_1_moles_oxygen, no_1_moles_H, no_1_moles_S_core, no_1_moles_C_core);
// No need to equilibrate so set equilibration conditions to zero
        P = T = fo2 = 0;
      }
// PARTIALLY OXIDIZED COMPOSITION - therefore differentiated
      // if (no_1 == 12 && write_file == 1 && e_time[no_1][e_coll[no_1]]/1e6 < 0.37 && e_time[no_1][e_coll[no_1]]/1e6 > 0.337  ) cout << endl << "the time is: "<< e_time[no_1][e_coll[no_1]]/1e6 <<  " and this is e_fe_sil SPOT 1 7a = " << no_1_Fe_sil_isotope[1] <<  " and e_Fe_sil_isotope[no_1][1][e_coll[no_1]]:" << e_Fe_sil_isotope[no_1][1][e_coll[no_1]] << endl;
      else
      {
        
        differentiated_composition(mass_1, mass_factor,
                                   no_1_moles, no_1_moles_H2O, no_1_moles_S_mantle, no_1_moles_C_mantle,
                                   no_1_moles_oxygen, no_1_moles_H, no_1_moles_S_core, no_1_moles_C_core);        
// Calculate the differentiated composition
        metal_silicate_equilibration(mass_1,
                                     no_1_moles, no_1_moles_H2O, no_1_moles_S_mantle, no_1_moles_C_mantle,
                                     no_1_moles_oxygen, no_1_moles_H, no_1_moles_S_core, no_1_moles_C_core,
                                     P, T, fo2, no1_fe, zero_xfe_liq);
       
        
        Fe_fractionate(no_1_Fe_sil_isotope, no_1_Fe_met_isotope,
        			   no_1_moles, no_1_moles_H2O, no_1_moles_S_mantle, no_1_moles_C_mantle,
        			   no_1_moles_oxygen, no_1_moles_H, no_1_moles_S_core, no_1_moles_C_core,
        			   P, T);//fractionate
        
      }
// Store the composition in the array of mantle and core compositions as a function of accretion event (e.g. e_moles[][][])
        
      
      initiate_bulk_composition(mass_1, no_1, e_coll[no_1],
                                no_1_moles, no_1_moles_H2O, no_1_moles_S_mantle, no_1_moles_C_mantle,
                                no_1_moles_oxygen, no_1_moles_H, no_1_moles_S_core, no_1_moles_C_core,
                                P, T, fo2);
      
      initiate_isotopic_composition(no_1, e_coll[no_1], no_1_Fe_sil_isotope, no_1_Fe_met_isotope);

      // cout << "Fe no_1_moles: " << no_1_moles[4] << endl;

      

// Keep track of mass accreted from H2O-bearing compositions
      if (dist[no_1] >= dist_3) e_SN_mass[no_1] = mass_1;
// Composition of #1 is now defined but it has no collisional history as yet
      b_status[no_1] = 1;
      // if (no_1 == 12 && write_file == 1 && e_time[no_1][e_coll[no_1]]/1e6 < 0.37 && e_time[no_1][e_coll[no_1]]/1e6 > 0.337  ) cout << endl << "the time is: "<< e_time[no_1][e_coll[no_1]]/1e6 <<  " and this is SPOT 8 and no_1_Fe_sil_isotope[1] = " << no_1_Fe_sil_isotope[1] << endl;
    
      // if (no_1 == 12 && write_file == 1 && e_time[no_1][e_coll[no_1]]/1e6 < 0.37 && e_time[no_1][e_coll[no_1]]/1e6 > 0.337  ) cout << endl << "the time is: "<< e_time[no_1][e_coll[no_1]]/1e6 <<  " and this is e_fe_sil SPOT 1 3a = " << no_1_Fe_sil_isotope[1] <<  " and e_Fe_sil_isotope[no_1][1][e_coll[no_1]]:" << e_Fe_sil_isotope[no_1][1][e_coll[no_1]] << endl;

   }
 
    // if (no_1 == 12 && write_file == 1 && e_time[no_1][e_coll[no_1]]/1e6 < 0.37 && e_time[no_1][e_coll[no_1]]/1e6 > 0.337  ) cout << endl << "the time is: "<< e_time[no_1][e_coll[no_1]]/1e6 <<  " and this is e_fe_sil SPOT 1aaa = " << no_1_Fe_sil_isotope[1] <<  " and e_Fe_sil_isotope[no_1][1][e_coll[no_1]]:" << e_Fe_sil_isotope[no_1][1][e_coll[no_1]] << endl;

// Define composition of no_1 when its differentiation/accretion history already exists (i.e. when no_1 is an embryo) //3isotope
      // cout << "1 :" << no_1_Fe_sil_isotope[1] << endl;
      // cout << "2 :" << no_1_Fe_sil_isotope[1] << endl;
      // // cout << endl << "the time is: "<< e_time[no_1][e_coll[no_1]]/1e6 <<  " and this is e_fe_sil SPOT g1= " << no_1_Fe_sil_isotope[1] <<  " and e_Fe_sil_isotope[no_1][1][e_coll[no_1]]:" << e_Fe_sil_isotope[no_1][1][e_coll[no_1]] << endl;
      // cout << "3 :" << no_1_Fe_sil_isotope[1] << endl;
    else if (b_status[no_1] == 2)
    { 
      
      // if (no_1 == 12 && write_file == 1) cout << endl << "the time is: "<< e_time[no_1][e_coll[no_1]]/1e6 <<  " and this is e_fe_sil SPOT 1aa = " << no_1_Fe_sil_isotope[1] <<  " and e_Fe_sil_isotope[no_1][1][e_coll[no_1]]:" << e_Fe_sil_isotope[no_1][1][e_coll[no_1]] << endl;

      for (int i = 0; i < (n + m); i++) no_1_moles[i] = e_moles[no_1][i][e_coll[no_1]];
      no_1_moles_H2O = e_moles_H2O[no_1][e_coll[no_1]];
      no_1_moles_S_mantle = e_moles_S_mantle[no_1][e_coll[no_1]];
      no_1_moles_C_mantle = e_moles_C_mantle[no_1][e_coll[no_1]];
      no_1_moles_oxygen = e_moles_oxygen[no_1][e_coll[no_1]];
      no_1_moles_H = e_moles_H[no_1][e_coll[no_1]];
      no_1_moles_S_core = e_moles_S_core[no_1][e_coll[no_1]];
      no_1_moles_C_core = e_moles_C_core[no_1][e_coll[no_1]];

     
      no_1_Fe_met_isotope[0] = e_Fe_met_isotope[no_1][0][e_coll[no_1]];
      no_1_Fe_met_isotope[1] = e_Fe_met_isotope[no_1][1][e_coll[no_1]];
      no_1_Fe_met_isotope[2] = e_Fe_met_isotope[no_1][2][e_coll[no_1]]; //3isotope
      no_1_Fe_sil_isotope[0] = e_Fe_sil_isotope[no_1][0][e_coll[no_1]];
      no_1_Fe_sil_isotope[1] = e_Fe_sil_isotope[no_1][1][e_coll[no_1]];
      no_1_Fe_sil_isotope[2] = e_Fe_sil_isotope[no_1][2][e_coll[no_1]]; //3isotope
     
      
    }
//______________________________________________________________________
// IN THE CASE THAT THE COMPOSITION OF no_2 IS NOT DEFINED


    if (b_status[no_2] == 0)
    {
// Collision number for no_2 when it is an embryo
      if (mass_2 > max_pl * mass_adjust && no_2 <= em_max )
      {
        e_coll[no_2] = 0;
        e_mass[no_2][e_coll[no_2]] = mass_2;
      }
// GET INITIAL COMPOSITION OF BODY #2
// Determine the fraction of Fe and Si of the body that is initially metal due to heliocentric distance
      xfe_met = Fe_metalfraction(dist_comp, dist[no_2], dist_1, dist_2, dist_3, xfe_met_1, xfe_met_2, xfe_met_3, erf_width_1, erf_width_2,
                                 moon_xfe_met_1, moon_xfe_met_2, moon_dist_1,moon_dist_2, moon_disk_start, moon_disk_end); // Moon_version 

      xsi_met = Si_metalfraction(dist_comp, dist[no_2], dist_1, dist_2, dist_3, xsi_met_0, erf_width_1, erf_width_2);
      starting_molar_proportions(dist[no_2],
                                 no_2_moles, no_2_moles_H2O, no_2_moles_S_mantle, no_2_moles_C_mantle,
                                 no_2_moles_oxygen, no_2_moles_H, no_2_moles_S_core, no_2_moles_C_core);
     
      starting_isotopic_proportions(no_2_Fe_sil_isotope, no_2_Fe_met_isotope);
     

      mass_factor = get_mass_factor(no_2_moles, no_2_moles_H2O, no_2_moles_S_mantle, no_2_moles_C_mantle,
                                    no_2_moles_oxygen, no_2_moles_H, no_2_moles_S_core, no_2_moles_C_core);
      
// DIFFERENTIATE BODY #2
// FULLY OXIDIZED COMPOSITION - therefore no metal, no differentiation, all S in silicate
      if (xfe_met < 1e-20)
      {
        
// Calculate the metal-free composition
        undifferentiated_composition(mass_2, mass_factor,
                                     no_2_moles, no_2_moles_H2O, no_2_moles_S_mantle, no_2_moles_C_mantle,
                                     no_2_moles_oxygen, no_2_moles_H, no_2_moles_S_core, no_2_moles_C_core);
// No need to equilibrate so set equilibration conditions to zero
        P = T = fo2 = 0;
        
      }
// PARTIALLY OXIDIZED COMPOSITION - therefore differentiated - have some iron metal in body
      else
      {
        
        differentiated_composition(mass_2, mass_factor,
                                   no_2_moles, no_2_moles_H2O, no_2_moles_S_mantle, no_2_moles_C_mantle,
                                   no_2_moles_oxygen, no_2_moles_H, no_2_moles_S_core, no_2_moles_C_core);

// Calculate the differentiated composition
        metal_silicate_equilibration(mass_2,
                                     no_2_moles, no_2_moles_H2O, no_2_moles_S_mantle, no_2_moles_C_mantle,
                                     no_2_moles_oxygen, no_2_moles_H, no_2_moles_S_core, no_2_moles_C_core,
                                     P, T, fo2, no2_fe, zero_xfe_liq);
        

        Fe_fractionate(no_2_Fe_sil_isotope, no_2_Fe_met_isotope,
                       no_2_moles, no_2_moles_H2O, no_2_moles_S_mantle, no_2_moles_C_mantle,
                       no_2_moles_oxygen, no_2_moles_H, no_2_moles_S_core, no_2_moles_C_core,
                       P, T);//fractionate
        // if (no_2 == 49 && write_file == 1) cout << "no_2 = " << no_2 << " and no_2_Fe_sil_isotope[0] = " << no_2_Fe_sil_isotope[0] << endl;
        
        // NR_Fe(no_2_Fe_sil_isotope, no_2_Fe_met_isotope,
        //       no_2_moles, no_2_moles_H2O, no_2_moles_S_mantle, no_2_moles_C_mantle,
        //       no_2_moles_oxygen, no_2_moles_H, no_2_moles_S_core, no_2_moles_C_core,
        //       P, T);
        // if (no_2 == 49 && write_file == 1) cout << "no_2 = " << no_2 << " and no_2_Fe_sil_isotope[0] = " << no_2_Fe_sil_isotope[0] << endl;
        
        

      }
// Store the composition in the array of mantle and core compositions as a function of accretion event, only for embryos
      if (mass_2 > max_pl * mass_adjust && no_2 <= em_max )
      {
        
        initiate_bulk_composition(mass_2, no_2, e_coll[no_2],
                                  no_2_moles, no_2_moles_H2O, no_2_moles_S_mantle, no_2_moles_C_mantle,
                                  no_2_moles_oxygen, no_2_moles_H, no_2_moles_S_core, no_2_moles_C_core,
                                  P, T, fo2);
        
        
        initiate_isotopic_composition(no_2, e_coll[no_2], no_2_Fe_sil_isotope, no_2_Fe_met_isotope);

        // cout << "Fe no_1_moles undiff: " << no_1_moles[4] << endl;

        
        

      }
// Keep track of mass accreted from H2O-bearing compositions
      if (dist[no_2] >= dist_3) e_SN_mass[no_2] = mass_2;
// Composition of #2 is now defined but it has no collisional history
      b_status[no_2] = 1;
    }
// Define composition of no_2 when its differentiation/accretion history already exists (i.e. when no_2 is an embryo)
    else if (b_status[no_2] == 2)
    {

      for (int i = 0; i < (n + m); i++) no_2_moles[i] = e_moles[no_2][i][e_coll[no_2]];
      no_2_moles_H2O = e_moles_H2O[no_2][e_coll[no_2]];
      no_2_moles_S_mantle = e_moles_S_mantle[no_2][e_coll[no_2]];
      no_2_moles_C_mantle = e_moles_C_mantle[no_2][e_coll[no_2]];
      no_2_moles_oxygen = e_moles_oxygen[no_2][e_coll[no_2]];
      no_2_moles_H = e_moles_H[no_2][e_coll[no_2]];
      no_2_moles_S_core = e_moles_S_core[no_2][e_coll[no_2]];
      no_2_moles_C_core = e_moles_C_core[no_2][e_coll[no_2]];

      
      no_2_Fe_met_isotope[0] = e_Fe_met_isotope[no_2][0][e_coll[no_2]];
      no_2_Fe_met_isotope[1] = e_Fe_met_isotope[no_2][1][e_coll[no_2]];
      no_2_Fe_met_isotope[2] = e_Fe_met_isotope[no_2][2][e_coll[no_2]]; //3isotope
      no_2_Fe_sil_isotope[0] = e_Fe_sil_isotope[no_2][0][e_coll[no_2]];
      no_2_Fe_sil_isotope[1] = e_Fe_sil_isotope[no_2][1][e_coll[no_2]];
      no_2_Fe_sil_isotope[2] = e_Fe_sil_isotope[no_2][2][e_coll[no_2]]; //3isotope
      
      
      
    }

// METAL ANND SILICATE COMPOSITIONS OF BOTH (NEW OR OLD, DIFFERENTIATED OR OXIDIZED/UNDIFFERENTIATED) BODIES ARE NOW DEFINED
//______________________________________________________________________________________________________________________________________________
// COLLISION, ACCRETION AND METAL-SILICATE EQUILIBRATION EVENTS

// Update collison numbers, masses, impact times, contributions from embryos and planetesimals and bodies from beyond SL:
// Collision numbers
    (e_coll[no_1])++;
    (e_coll[no_2])++;
// Time of collision
    e_time[no_1][e_coll[no_1]] = t;
// Total mass after collision
    e_mass[no_1][e_coll[no_1]] = mass_1 + mass_2;

// DAUPHAS INFORMATION
    i_mass[no_1][e_coll[no_1]] = mass_2;
    for (int i = 0; i < (n + m); i++) i_moles[no_1][i][e_coll[no_1]] = no_2_moles[i];
    i_moles_H2O[no_1][e_coll[no_1]] = no_2_moles_H2O;
    i_moles_S_mantle[no_1][e_coll[no_1]] = no_2_moles_S_mantle;
    i_moles_C_mantle[no_1][e_coll[no_1]] = no_2_moles_C_mantle;
    i_moles_oxygen[no_1][e_coll[no_1]] = no_2_moles_oxygen;
    i_moles_H[no_1][e_coll[no_1]] = no_2_moles_H;
    i_moles_S_core[no_1][e_coll[no_1]] = no_2_moles_S_core;
    i_moles_C_core[no_1][e_coll[no_1]] = no_2_moles_C_core;

// Track respective contributions of embryos and planetesimals respectively, with respect to total mass accreted
// E.g.: collison numbers, masses, impact times, contributions from embryos and planetesimals and bodies from beyond SL:
    update_contributions(no_1, no_2, mass_1, mass_2);
// Track whether the collision is a giant impact
    if (mass_2 > max_pl * mass_adjust) GI[no_1] = 1;
    else GI[no_1] = 0;
// Track whether the projectile has a core
    if ((no_2_moles[n] > 1e-20)) no_2_met = 1;
    else no_2_met = 0;

// DEFINE P-T EQUILIBRATION CONDITIONS:
// Different Pf for embryos option (pressure factor = Pf)
    if (mass_2 > max_pl * mass_adjust && em_equil == 1) pf = pf_em;
    else pf = p_factor;
    P = pf * pressure(mass_1 + mass_2);
    T = temperature(P) + t_overstep;
    fo2 = 0.0;

// SPECIFY EXTENT OF IMPACTOR'S CORE EQUILIBRATION (currently two possibilities - DEPENDING ON MASS OF IMPACTOR):
// xmet is the fraction of impactor's core that equilibrates
    double xmet(0.0);
    if (mass_2 <= max_pl * mass_adjust) xmet = xmet_pl;
    else if (mass_2 > max_pl * mass_adjust) xmet = xmet_em;

// Check that there is sufficient oxygen to oxidize all metal in case of late veneer
    /*    double required_moles_oxygen = 0.0;
        for (int i = n; i < (n + m); i++) required_moles_oxygen += no_2_moles[i] * val[i - m] / 2.0;
        double available_moles_oxygen = no_1_moles_H2O + no_2_moles_H2O + no_2_moles_oxygen;*/

// OUTLINE OF DIFFERENTIATION STEPS
// I. Primary Differentiation.
//    Does the accretion event take place during the late veneer?
//       N. Is there interacting projectile metal? (if no, projectile is from outer solar system and (currently) it is also completely oxidized)
//          N. Is the projectile completely oxidized?
//               Y. Add the projectile to the target's mantle.
//               N. Then the projectile's core is added to the target's core, and the projectile's mantle is added to the target's mantle.
//          Y. Is there enough water in the interacting mantle to oxidize the projectile's metal?
//             Y. Add the projectile to the target's mantle; change the water budget to reflect this oxidation.
//             N. Is there sulfur?
//                N. Use react to determine what the composition of the new metal and oxide liquids; then add them to the mantle and core.
//                Y. Same as above, but iterate the react solution until the S is also equilibrated.
//       Y. Is there enough water in the interacting mantle to oxidize the projectile's metal?
//          Y. Add the projectile to the target's mantle; change the water budget to reflect this oxidation.
//          N. Do the above anyways, but report to the screen that there is not enough water. (&punish chi sq)
// II. Hadean Matte
//    Are the conditions for the Hadean Matte met?
//       N. Do not have a Hadean matte.
//       Y. Do have a Hadean matte.
//__________________________________________________________________________________________________________
// PRIMARY DIFFERENTIATION
// MAIN STAGE ACCRETION (PRE-LATE VENEER)
    if (t < t0_LV[no_1] ) //|| (t >= t0_LV[no_1] && required_moles_oxygen > available_moles_oxygen && earth_planet != no_1)
    {

//______________________________________________________________________
// IS THERE INTERACTING PROJECTILE METAL?
// PROJECTILE HAS NO METAL OR NONE OF THE METAL INTERACTS WITH THE TARGET'S MANTLE
// First condition is if there is no metal-silicate equilibration and the second condition is if the projectile has no metal
      if ((xmet == 0.0) || (no_2_met == 0))
      {

// Count such collisions
        n0_pl++;
// Error check
        if (xmet != 0.0)
        {
          for (int i = n; i < (n + m); i++) if (no_2_moles[i] > 1e-20) cout << "Projectile with metal combined as if there was no metal." << tab << react_call_no << endl;
          if (no_2_moles_oxygen > 1e-20) cout << "Projectile with metal combined as if there was no metal." << endl;
          if (no_2_moles_H > 1e-20) cout << "Projectile with metal combined as if there was no metal." << endl;
          if (no_2_moles_S_core > 1e-20) cout << "Projectile with metal combined as if there was no metal." << endl;
          if (no_2_moles_C_core > 1e-20) cout << "Projectile with metal combined as if there was no metal." << endl;
        }
// With no core, there is no metal-silicate equilibrating material
        xsil = 0.0;
        xmet = 0.0;
// DAUPHAS INFORMATION
        e_equil[no_1][e_coll[no_1]] = xsil;
        i_equil[no_1][e_coll[no_1]] = xmet;
// Switch from no_1 and no_2 reservoirs to equilibrating and nonequilibrating reservoirs
        
       
        identify_molar_reservoirs(xsil, xmet,
                                  no_1_moles, no_1_moles_H2O, no_1_moles_S_mantle, no_1_moles_C_mantle,
                                  no_1_moles_oxygen, no_1_moles_H, no_1_moles_S_core, no_1_moles_C_core,
                                  no_1_Fe_sil_isotope, no_1_Fe_met_isotope,
                                  no_2_moles, no_2_moles_H2O, no_2_moles_S_mantle, no_2_moles_C_mantle,
                                  no_2_moles_oxygen, no_2_moles_H, no_2_moles_S_core, no_2_moles_C_core,
                                  no_2_Fe_sil_isotope, no_2_Fe_met_isotope,
                                  equil_moles, equil_moles_H2O, equil_moles_S_mantle, equil_moles_C_mantle,
                                  equil_moles_oxygen, equil_moles_H, equil_moles_S_core, equil_moles_C_core,
                                  equil_Fe_sil_isotope, equil_Fe_met_isotope,
                                  noneq_moles, noneq_moles_H2O, noneq_moles_S_mantle, noneq_moles_C_mantle,
                                  noneq_moles_oxygen, noneq_moles_H, noneq_moles_S_core, noneq_moles_C_core,
                                  noneq_Fe_sil_isotope, noneq_Fe_met_isotope);
        


// No need to equilibrate so set equilibration conditions to zero
        P = 0.0;
        T = 0.0;
        fo2 = 0.0;
// Store the composition in the array of mantle and core compositions as a function of accretion event (e.g. e_moles[][][])
        update_bulk_composition(no_1, e_coll[no_1],
                                equil_moles, equil_moles_H2O, equil_moles_S_mantle, equil_moles_C_mantle,
                                equil_moles_oxygen, equil_moles_H, equil_moles_S_core, equil_moles_C_core,
                                noneq_moles, noneq_moles_H2O, noneq_moles_S_mantle, noneq_moles_C_mantle,
                                noneq_moles_oxygen, noneq_moles_H, noneq_moles_S_core, noneq_moles_C_core,
                                P, T, fo2);
        
        
        update_Fe_composition(no_1, e_coll[no_1],
                              equil_moles, equil_Fe_sil_isotope, equil_Fe_met_isotope,
                              noneq_moles, noneq_Fe_sil_isotope, noneq_Fe_met_isotope);
        
        update_sources(no_1, e_coll[no_1], xsil, xmet,
                       no_1_moles, no_1_moles_H2O, no_1_moles_S_mantle, no_1_moles_C_mantle,
                       no_1_moles_oxygen, no_1_moles_H, no_1_moles_S_core, no_1_moles_C_core,
                       no_2_moles, no_2_moles_H2O, no_2_moles_S_mantle, no_2_moles_C_mantle,
                       no_2_moles_oxygen, no_2_moles_H, no_2_moles_S_core, no_2_moles_C_core,
                       equil_moles, equil_moles_H2O, equil_moles_S_mantle, equil_moles_C_mantle,
                       equil_moles_oxygen, equil_moles_H, equil_moles_S_core, equil_moles_C_core,
                       noneq_moles, noneq_moles_H2O, noneq_moles_S_mantle, noneq_moles_C_mantle,
                       noneq_moles_oxygen, noneq_moles_H, noneq_moles_S_core, noneq_moles_C_core);

// There is no core-forming liquid metal
        for (int i = n; i < (n + m); i++) incr_moles[no_1][i][e_coll[no_1]] = 0.0;
        incr_moles_oxygen[no_1][e_coll[no_1]]  = 0.0;
        incr_moles_H[no_1][e_coll[no_1]] = 0.0;
        incr_moles_S[no_1][e_coll[no_1]] = 0.0;
        incr_moles_C[no_1][e_coll[no_1]] = 0.0;
      }
//______________________________________________________________________
// PROJECTILE HAS METAL AND THE METAL INTERACTS WITH THE TARGET'S MANTLE
// EITHER COMPLETE OR PARTIAL METAL-SILICATE EQUILIBRATION:
// Define composition of equilibrating system (the array "moles" and "moles_oxygen" are input for the mass balance function "react")
// xmet = fraction of impactor's core that equilibrates and xsil = fraction of proto-planet's mantle that equilibrates
      else
      {
// Determine xsil (fraction of silicate mantle that equilibrates) using Eq(6) of Deguen et al (2011, EPSL 310, 303-313):
        if (x_mantle == 0.0) xsil = silicate_fraction(mass_1, mass_2,
                                      no_1_moles, no_1_moles_H2O, no_1_moles_S_mantle, no_1_moles_C_mantle,
                                      no_1_moles_oxygen, no_1_moles_H, no_1_moles_S_core, no_1_moles_C_core,
                                      no_2_moles, no_2_moles_H2O, no_2_moles_S_mantle, no_2_moles_C_mantle,
                                      no_2_moles_oxygen, no_2_moles_H, no_2_moles_S_core, no_2_moles_C_core);
// Deguen model is overridden FOR ALL IMPACTS
        else if (x_mantle > 0.0) xsil = x_mantle;
// Changes value of xsil by the factor em_xsil for embryo impacts
        if (mass_2 > max_pl * mass_adjust) xsil = em_xsil * xsil;
        if (xsil > 1.0 ) xsil = 1.0;
// DAUPHAS INFORMATION
        e_equil[no_1][e_coll[no_1]] = xsil;
        i_equil[no_1][e_coll[no_1]] = xmet;
//Store max, min and av values of xsil:
        if (error == 0) record_silicate_fraction(mass_2, xsil);
// Switch from no_1 and no_2 reservoirs to equilibrating and nonequilibrating reservoirs
       
      
        identify_molar_reservoirs(xsil, xmet,
                                  no_1_moles, no_1_moles_H2O, no_1_moles_S_mantle, no_1_moles_C_mantle,
                                  no_1_moles_oxygen, no_1_moles_H, no_1_moles_S_core, no_1_moles_C_core,
                                  no_1_Fe_sil_isotope, no_1_Fe_met_isotope,
                                  no_2_moles, no_2_moles_H2O, no_2_moles_S_mantle, no_2_moles_C_mantle,
                                  no_2_moles_oxygen, no_2_moles_H, no_2_moles_S_core, no_2_moles_C_core,
                                  no_2_Fe_sil_isotope, no_2_Fe_met_isotope,
                                  equil_moles, equil_moles_H2O, equil_moles_S_mantle, equil_moles_C_mantle,
                                  equil_moles_oxygen, equil_moles_H, equil_moles_S_core, equil_moles_C_core,
                                  equil_Fe_sil_isotope, equil_Fe_met_isotope,
                                  noneq_moles, noneq_moles_H2O, noneq_moles_S_mantle, noneq_moles_C_mantle,
                                  noneq_moles_oxygen, noneq_moles_H, noneq_moles_S_core, noneq_moles_C_core,
                                  noneq_Fe_sil_isotope, noneq_Fe_met_isotope);
        
        


//__________________________________
// IS THERE ENOUGH WATER TO OXIDIZE INTERACTING PROJECTILE CORE MASS?
// Determine moles of oxygen required to completely oxidize the equilibrating fraction of impactor's core:
        double required_moles_oxygen = 0.0;
        for (int i = n; i < (n + m); i++) required_moles_oxygen += equil_moles[i] * val[i - m] / 2.0;
        required_moles_oxygen += equil_moles_H * 0.5;
// Determine moles of oxygen available in the equilibrating melt:
        double available_moles_oxygen = equil_moles_oxygen + equil_moles_H2O;
// THERE IS ENOUGH WATER TO OXIDIZE INTERACTING PROJECTILE CORE MASS AND THE PROJECTILE IS A PLANETESIMAL
        if (available_moles_oxygen >= required_moles_oxygen && mass_2 <= max_pl * mass_adjust )
        {
// For reference, below are the indexes for elements in elements.h
// Element index: 0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22
// Element name: Al Mg Ca Na Fe Si Ni Co Nb Ta  V Cr Pt Pd Ru Ir  W Mo Mn Ga Zn In Cu
// Element index:            17 18 19 20 21 22 23 24 25 26 27 28 29
// Metals start at Fe at index n then continue to n+m where m = n - 4
// If n = 17, then only up to W is included
// All metals are oxidized and added to the silicate
          
          
          
          equil_Fe_sil_isotope[0] = (equil_Fe_sil_isotope[0] * equil_moles[n - m] + equil_Fe_met_isotope[0] * equil_moles[n]) / (equil_moles[n - m] + equil_moles[n]);
          equil_Fe_sil_isotope[1] = (equil_Fe_sil_isotope[1] * equil_moles[n - m] + equil_Fe_met_isotope[1] * equil_moles[n]) / (equil_moles[n - m] + equil_moles[n]);
          equil_Fe_sil_isotope[2] = (equil_Fe_sil_isotope[2] * equil_moles[n - m] + equil_Fe_met_isotope[2] * equil_moles[n]) / (equil_moles[n - m] + equil_moles[n]); //3isotope

          
      
          for (int i = (n - m); i < n; i++) equil_moles[i] = equil_moles[i] + equil_moles[i + m];
          for (int i = n; i < (n + m); i++) equil_moles[i] = 0.0;
          if (required_moles_oxygen <= equil_moles_oxygen) equil_moles_H2O = equil_moles_H2O;
          else equil_moles_H2O = equil_moles_H2O - required_moles_oxygen + equil_moles_oxygen;
          equil_moles_S_mantle = equil_moles_S_mantle + equil_moles_S_core;
          equil_moles_C_mantle = equil_moles_C_mantle + equil_moles_C_core;
          equil_moles_oxygen = 0.0;
          equil_moles_H = 0.0;
          equil_moles_S_core = 0.0;
          equil_moles_C_core = 0.0;
          

// No need to equilibrate so set equilibration conditions to zero
          P = 0.0;
          T = 0.0;
          fo2 = 0.0;
// Store the composition in the array of mantle and core compositions as a function of accretion event (e.g. e_moles[][][])
          update_bulk_composition(no_1, e_coll[no_1],
                                  equil_moles, equil_moles_H2O, equil_moles_S_mantle, equil_moles_C_mantle,
                                  equil_moles_oxygen, equil_moles_H, equil_moles_S_core, equil_moles_C_core,
                                  noneq_moles, noneq_moles_H2O, noneq_moles_S_mantle, noneq_moles_C_mantle,
                                  noneq_moles_oxygen, noneq_moles_H, noneq_moles_S_core, noneq_moles_C_core,
                                  P, T, fo2);
          

          
          update_Fe_composition(no_1, e_coll[no_1],
                                equil_moles, equil_Fe_sil_isotope, equil_Fe_met_isotope,
                                noneq_moles, noneq_Fe_sil_isotope, noneq_Fe_met_isotope);
         
          
          update_sources(no_1, e_coll[no_1], xsil, xmet,
                         no_1_moles, no_1_moles_H2O, no_1_moles_S_mantle, no_1_moles_C_mantle,
                         no_1_moles_oxygen, no_1_moles_H, no_1_moles_S_core, no_1_moles_C_core,
                         no_2_moles, no_2_moles_H2O, no_2_moles_S_mantle, no_2_moles_C_mantle,
                         no_2_moles_oxygen, no_2_moles_H, no_2_moles_S_core, no_2_moles_C_core,
                         equil_moles, equil_moles_H2O, equil_moles_S_mantle, equil_moles_C_mantle,
                         equil_moles_oxygen, equil_moles_H, equil_moles_S_core, equil_moles_C_core,
                         noneq_moles, noneq_moles_H2O, noneq_moles_S_mantle, noneq_moles_C_mantle,
                         noneq_moles_oxygen, noneq_moles_H, noneq_moles_S_core, noneq_moles_C_core);
        }
//__________________________________
// THERE IS NOT ENOUGH WATER TO OXIDIZE INTERACTING PROJECTILE CORE MASS
        else
        {
// Calculate the differentiated composition
          if (no_1 == 28 && no_2 == 37 && write_file == 1) cout << xsil << " " << equil_moles[n] / equil_moles[n - m] << endl;
          metal_silicate_equilibration(mass_1,
                                       equil_moles, equil_moles_H2O, equil_moles_S_mantle, equil_moles_C_mantle,
                                       equil_moles_oxygen, equil_moles_H, equil_moles_S_core, equil_moles_C_core,
                                       P, T, fo2, no1_fe, zero_xfe_liq);

          // if (no_1 == 12 && write_file == 1) cout << endl << "the time is: "<< e_time[no_1][e_coll[no_1]]/1e6 <<  " and this is e_fe_sil SPOT wow14 = " << no_1_Fe_sil_isotope[1] <<  " and e_Fe_sil_isotope[no_1][1][e_coll[no_1]]:" << e_Fe_sil_isotope[no_1][1][e_coll[no_1]] << endl;



          if (zero_xfe_liq == 1)
          {
// Store the composition in the array of mantle and core compositions as a function of accretion event (e.g. e_moles[][][])
            

            
        
            Fe_fractionate(equil_Fe_sil_isotope, equil_Fe_met_isotope,
                           equil_moles, equil_moles_H2O, equil_moles_S_mantle, equil_moles_C_mantle,
                           equil_moles_oxygen, equil_moles_H, equil_moles_S_core, equil_moles_C_core,
                           P, T);//fractionate
            // if (no_1 == 49 && write_file == 1) cout << "no_1 = " << no_1 << " and no_1_Fe_sil_isotope[0] = " << no_1_Fe_sil_isotope[0] << endl;
        
	        // NR_Fe(equil_Fe_sil_isotope, equil_Fe_met_isotope,
         //          equil_moles, equil_moles_H2O, equil_moles_S_mantle, equil_moles_C_mantle,
         //          equil_moles_oxygen, equil_moles_H, equil_moles_S_core, equil_moles_C_core,
         //          P, T);
	        // if (no_1 == 49 && write_file == 1) cout << "no_1 = " << no_1 << " and no_1_Fe_sil_isotope[0] = " << no_1_Fe_sil_isotope[0] << endl;
	        
            

            update_bulk_composition(no_1, e_coll[no_1],
                                    equil_moles, equil_moles_H2O, equil_moles_S_mantle, equil_moles_C_mantle,
                                    equil_moles_oxygen, equil_moles_H, equil_moles_S_core, equil_moles_C_core,
                                    noneq_moles, noneq_moles_H2O, noneq_moles_S_mantle, noneq_moles_C_mantle,
                                    noneq_moles_oxygen, noneq_moles_H, noneq_moles_S_core, noneq_moles_C_core,
                                    P, T, fo2);
    
            
            update_Fe_composition(no_1, e_coll[no_1],
                                  equil_moles, equil_Fe_sil_isotope, equil_Fe_met_isotope,
                                  noneq_moles, noneq_Fe_sil_isotope, noneq_Fe_met_isotope);
            
            
          
            update_sources(no_1, e_coll[no_1], xsil, xmet,
                           no_1_moles, no_1_moles_H2O, no_1_moles_S_mantle, no_1_moles_C_mantle,
                           no_1_moles_oxygen, no_1_moles_H, no_1_moles_S_core, no_1_moles_C_core,
                           no_2_moles, no_2_moles_H2O, no_2_moles_S_mantle, no_2_moles_C_mantle,
                           no_2_moles_oxygen, no_2_moles_H, no_2_moles_S_core, no_2_moles_C_core,
                           equil_moles, equil_moles_H2O, equil_moles_S_mantle, equil_moles_C_mantle,
                           equil_moles_oxygen, equil_moles_H, equil_moles_S_core, equil_moles_C_core,
                           noneq_moles, noneq_moles_H2O, noneq_moles_S_mantle, noneq_moles_C_mantle,
                           noneq_moles_oxygen, noneq_moles_H, noneq_moles_S_core, noneq_moles_C_core);
                       
          }
          
          else
          { 
      
            equil_Fe_sil_isotope[0] = (equil_Fe_sil_isotope[0] * equil_moles[n - m] + equil_Fe_met_isotope[0] * equil_moles[n]) / (equil_moles[n - m] + equil_moles[n]);
            equil_Fe_sil_isotope[1] = (equil_Fe_sil_isotope[1] * equil_moles[n - m] + equil_Fe_met_isotope[1] * equil_moles[n]) / (equil_moles[n - m] + equil_moles[n]);
            equil_Fe_sil_isotope[2] = (equil_Fe_sil_isotope[2] * equil_moles[n - m] + equil_Fe_met_isotope[2] * equil_moles[n]) / (equil_moles[n - m] + equil_moles[n]); //3isotope
            
            
      
            // All metals are oxidized and added to the silicate
            for (int i = (n - m); i < n; i++) equil_moles[i] = equil_moles[i] + equil_moles[i + m];
            for (int i = n; i < (n + m); i++) equil_moles[i] = 0.0;
            equil_moles_H2O = 0.0;
            equil_moles_S_mantle = equil_moles_S_mantle + equil_moles_S_core;
            equil_moles_C_mantle = equil_moles_C_mantle + equil_moles_C_core;
            equil_moles_oxygen = 0.0;
            equil_moles_H = 0.0;
            equil_moles_S_core = 0.0;
            equil_moles_C_core = 0.0;

            // No need to equilibrate so set equilibration conditions to zero
            P = 0.0;
            T = 0.0;
            fo2 = 0.0;
            // Store the composition in the array of mantle and core compositions as a function of accretion event (e.g. e_moles[][][])
            update_bulk_composition(no_1, e_coll[no_1],
                                    equil_moles, equil_moles_H2O, equil_moles_S_mantle, equil_moles_C_mantle,
                                    equil_moles_oxygen, equil_moles_H, equil_moles_S_core, equil_moles_C_core,
                                    noneq_moles, noneq_moles_H2O, noneq_moles_S_mantle, noneq_moles_C_mantle,
                                    noneq_moles_oxygen, noneq_moles_H, noneq_moles_S_core, noneq_moles_C_core,
                                    P, T, fo2);
            
          
            update_Fe_composition(no_1, e_coll[no_1],
                                  equil_moles, equil_Fe_sil_isotope, equil_Fe_met_isotope,
                                  noneq_moles, noneq_Fe_sil_isotope, noneq_Fe_met_isotope);
            
            update_sources(no_1, e_coll[no_1], xsil, xmet,
                           no_1_moles, no_1_moles_H2O, no_1_moles_S_mantle, no_1_moles_C_mantle,
                           no_1_moles_oxygen, no_1_moles_H, no_1_moles_S_core, no_1_moles_C_core,
                           no_2_moles, no_2_moles_H2O, no_2_moles_S_mantle, no_2_moles_C_mantle,
                           no_2_moles_oxygen, no_2_moles_H, no_2_moles_S_core, no_2_moles_C_core,
                           equil_moles, equil_moles_H2O, equil_moles_S_mantle, equil_moles_C_mantle,
                           equil_moles_oxygen, equil_moles_H, equil_moles_S_core, equil_moles_C_core,
                           noneq_moles, noneq_moles_H2O, noneq_moles_S_mantle, noneq_moles_C_mantle,
                           noneq_moles_oxygen, noneq_moles_H, noneq_moles_S_core, noneq_moles_C_core);
            
          }
        }
      }
         
//__________________________________________________________________________________________________________
// HADEAN MATTE - FeS segregation to core
// Now determine if the mantle is supersaturated with respect to S and, if so, by how much. In the case of superstaturation, FeS liquid also segregates to the core.
// The parameter S_frac_matte is the actual S concentration in mantle divided by the mantle S saturation concentration; there is supersaturation when this parameters is >1.0
      for (int i = 0; i < n; i++) moles[i] = e_moles[no_1][i][e_coll[no_1]];
      moles_H2O = e_moles_H2O[no_1][e_coll[no_1]];
      moles_S_mantle = e_moles_S_mantle[no_1][e_coll[no_1]];
      moles_C_mantle = e_moles_C_mantle[no_1][e_coll[no_1]];
      P = pf * pressure(mass_1 + mass_2);
      double S_frac_matte = S_conc(moles, moles_H2O, moles_S_mantle, moles_C_mantle, P);
// Mantle is supersaturated in S and FeS liquid (matte) (+ HSEs) segregates to the core
      if ((t >= t_matte[no_1] && last_GI[no_1] == 1 && S_frac_matte > 1.0 && moles[4] > 0.0 && matte_HSE == 1) ||
          (t >= t_matte[no_1] && t > t_last_GI[no_1] && S_frac_matte > 1.0 && moles[4] > 0.0 && matte_HSE == 1))
      {
        // if (no_1 == 12 && write_file == 1) cout << endl << "the time is: "<< e_time[no_1][e_coll[no_1]]/1e6 <<  " and this is e_fe_sil SPOT wow11 = " << no_1_Fe_sil_isotope[1] <<  " and e_Fe_sil_isotope[no_1][1][e_coll[no_1]]:" << e_Fe_sil_isotope[no_1][1][e_coll[no_1]] << endl;

// Moles of exsolved FeS liquid
        double h_moles_FeS = 0.0;
// Moles of S that remain dissolved in mantle
        double h_moles_S = 0.0;

// Moles FeS that segregate
        h_moles_FeS = moles_S_mantle * (1.0 - 1.0 / S_frac_matte);
// Check to make sure there is enough FeO to create that much FeS
        if (moles[4] < h_moles_FeS) h_moles_FeS = moles[4];
// Moles S remaining in mantle
        h_moles_S = moles_S_mantle - h_moles_FeS;

// Now calculate partitioning of HSEs between FeS and silicate mantle and determine new mantle and core concentrations
        double h_moles[n] = {0.0};
        for (int i = 0; i < n; i++) h_moles[i] = e_moles[no_1][i][e_coll[no_1]];
        double h_moles_H2O = e_moles_H2O[no_1][e_coll[no_1]];

        double pmoles_HSE[2 * n_hse] = {0.0}; // pmoles_HSE give HSE concentrations after partitioning between FeS and silicate

        FeS_HSE(h_moles, h_moles_H2O, h_moles_S, h_moles_FeS, pmoles_HSE, P);

        // New concentrations after extraction of FeS to the core
        e_moles_S_core[no_1][e_coll[no_1]] += h_moles_FeS;
        e_moles_S_mantle[no_1][e_coll[no_1]] = h_moles_S;
        e_moles[no_1][4][e_coll[no_1]] -= h_moles_FeS; // New FeO concentration in mantle
        e_moles[no_1][n][e_coll[no_1]] += h_moles_FeS; // New Fe content of core

        if (e_fO2[no_1][e_coll[no_1]] == 0.0) e_fO2[no_1][e_coll[no_1]] = 99.0;

        // Now retrieve pmoles_HSE and calculate new mantle and core HSE concentrations
        int hse_no = 0;
        for (int i = 0; i < n; i++) // mantle concentrations
        {
          if (el_type[i] == 'h' || el_type[i] == 'H')
          {
            e_moles[no_1][i][e_coll[no_1]] = e_moles[no_1][i][e_coll[no_1]] - pmoles_HSE[hse_no]; // Subtract the moles of HSEs contained in FeS from the mantle:
            hse_no++;
          }
        }

        // Add the moles of HSEs contained in FeS to the core:
        hse_no = 0;
        for (int i = n; i < n + m; i++)
        {
          if (el_type[i - m] == 'h' || el_type[i] == 'H')
          {
            e_moles[no_1][i][e_coll[no_1]] = e_moles[no_1][i][e_coll[no_1]] + pmoles_HSE[hse_no];
            hse_no++;
          }
        }

        update_sources_hse(no_1, e_coll[no_1], h_moles, h_moles_S, h_moles_FeS, pmoles_HSE);
      }
      // if (no_1 == 12 && write_file == 1) cout << endl << "the time is: "<< e_time[no_1][e_coll[no_1]]/1e6 <<  " and this is e_fe_sil SPOT wow11 = " << no_1_Fe_sil_isotope[1] <<  " and e_Fe_sil_isotope[no_1][1][e_coll[no_1]]:" << e_Fe_sil_isotope[no_1][1][e_coll[no_1]] << endl;


      // Store compositions of metal that are being added to proto-core (composition of core-forming increments)
      for (int i = n; i < (n + m); i++) incr_moles[no_1][i][e_coll[no_1]] = e_moles[no_1][i][e_coll[no_1]] - e_moles[no_1][i][e_coll[no_1] - 1];
      incr_moles_oxygen[no_1][e_coll[no_1]]  = e_moles_oxygen[no_1][e_coll[no_1]] - e_moles_oxygen[no_1][e_coll[no_1] - 1];
      incr_moles_H[no_1][e_coll[no_1]] = e_moles_H[no_1][e_coll[no_1]] - e_moles_H[no_1][e_coll[no_1] - 1];
      incr_moles_S[no_1][e_coll[no_1]] = e_moles_S_core[no_1][e_coll[no_1]] - e_moles_S_core[no_1][e_coll[no_1] - 1];
      incr_moles_C[no_1][e_coll[no_1]] = e_moles_C_core[no_1][e_coll[no_1]] - e_moles_C_core[no_1][e_coll[no_1] - 1];
    }
//__________________________________________________________________________________________________________
// LATE VENEER - ANY ACCRETED METAL IS STRANDED AND FULLY OXIDIZED BY WATER (if available)
// For reference, below are the indexes for elements in elements.h
// Element index: 0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16
// Element name: Al Mg Ca Na Fe Si Ni Co Nb Ta  V Cr Pt Pd Ru Ir  W Mo Mn Ga Zn In Cu
// Element index:            17 18 19 20 21 22 23 24 25 26 27 28 29
// Metals start at Fe at index n then continue to n+m where m = n - 4
// If n = 17, then only up to W is included
// All metals are oxidized and added to the silicate
    else         // "else" means t>=t0_LV for body no_1
    {
      // if (no_1 == 12 && write_file == 1) cout << endl << "the time is: "<< e_time[no_1][e_coll[no_1]]/1e6 <<  " and this is e_fe_sil SPOT wow10 = " << no_1_Fe_sil_isotope[1] <<  " and e_Fe_sil_isotope[no_1][1][e_coll[no_1]]:" << e_Fe_sil_isotope[no_1][1][e_coll[no_1]] << endl;

// Check that there is sufficient water to oxidize all metal:
      double required_moles_oxygen = 0.0;
      for (int i = n; i < (n + m); i++) required_moles_oxygen += no_2_moles[i] * val[i - m] / 2.0;
      double available_moles_oxygen = no_1_moles_H2O + no_2_moles_H2O + no_2_moles_oxygen;
// Determine moles of oxygen available in the equilibrating melt:
      if (required_moles_oxygen > available_moles_oxygen)
      {
        //if (earth_planet == no_1) chi_sq_multiplier++; //punish chi-squ with chi_sq_multiplier if there is not enough water

// Determine xsil (fraction of silicate mantle that equilibrates) using Eq(6) of Deguen et al (2011, EPSL 310, 303-313):
        if (x_mantle == 0.0) xsil = silicate_fraction(mass_1, mass_2,
                                      no_1_moles, no_1_moles_H2O, no_1_moles_S_mantle, no_1_moles_C_mantle,
                                      no_1_moles_oxygen, no_1_moles_H, no_1_moles_S_core, no_1_moles_C_core,
                                      no_2_moles, no_2_moles_H2O, no_2_moles_S_mantle, no_2_moles_C_mantle,
                                      no_2_moles_oxygen, no_2_moles_H, no_2_moles_S_core, no_2_moles_C_core);
// Deguen model is overridden FOR ALL IMPACTS
        else if (x_mantle > 0.0) xsil = x_mantle;
// Changes value of xsil by the factor em_xsil for embryo impacts
        if (mass_2 > max_pl * mass_adjust) xsil = em_xsil * xsil;
        if (xsil > 1.0 ) xsil = 1.0;
// DAUPHAS INFORMATION
        e_equil[no_1][e_coll[no_1]] = xsil;
        i_equil[no_1][e_coll[no_1]] = xmet;
//Store max, min and av values of xsil:
        if (error == 0) record_silicate_fraction(mass_2, xsil);
// Switch from no_1 and no_2 reservoirs to equilibrating and nonequilibrating reservoirs
        
       
        identify_molar_reservoirs(xsil, xmet,
                                  no_1_moles, no_1_moles_H2O, no_1_moles_S_mantle, no_1_moles_C_mantle,
                                  no_1_moles_oxygen, no_1_moles_H, no_1_moles_S_core, no_1_moles_C_core,
                                  no_1_Fe_sil_isotope, no_1_Fe_met_isotope,
                                  no_2_moles, no_2_moles_H2O, no_2_moles_S_mantle, no_2_moles_C_mantle,
                                  no_2_moles_oxygen, no_2_moles_H, no_2_moles_S_core, no_2_moles_C_core,
                                  no_2_Fe_sil_isotope, no_2_Fe_met_isotope,
                                  equil_moles, equil_moles_H2O, equil_moles_S_mantle, equil_moles_C_mantle,
                                  equil_moles_oxygen, equil_moles_H, equil_moles_S_core, equil_moles_C_core,
                                  equil_Fe_sil_isotope, equil_Fe_met_isotope,
                                  noneq_moles, noneq_moles_H2O, noneq_moles_S_mantle, noneq_moles_C_mantle,
                                  noneq_moles_oxygen, noneq_moles_H, noneq_moles_S_core, noneq_moles_C_core,
                                  noneq_Fe_sil_isotope, noneq_Fe_met_isotope);
        
// Calculate the differentiated composition
        metal_silicate_equilibration(mass_1,
                                     equil_moles, equil_moles_H2O, equil_moles_S_mantle, equil_moles_C_mantle,
                                     equil_moles_oxygen, equil_moles_H, equil_moles_S_core, equil_moles_C_core,
                                     P, T, fo2, no1_fe, zero_xfe_liq);
        if (zero_xfe_liq == 1)
        {
         
         
        
          Fe_fractionate(equil_Fe_sil_isotope, equil_Fe_met_isotope,
                         equil_moles, equil_moles_H2O, equil_moles_S_mantle, equil_moles_C_mantle,
                         equil_moles_oxygen, equil_moles_H, equil_moles_S_core, equil_moles_C_core,
                         P, T);//fractionate

          // if (no_1 == 49 && write_file == 1) cout << "no_1 = " << no_1 << " and no_1_Fe_sil_isotope[0] = " << no_1_Fe_sil_isotope[0] << endl;
        
          // NR_Fe(equil_Fe_sil_isotope, equil_Fe_met_isotope,
          //       equil_moles, equil_moles_H2O, equil_moles_S_mantle, equil_moles_C_mantle,
          //       equil_moles_oxygen, equil_moles_H, equil_moles_S_core, equil_moles_C_core,
          //       P, T);
          // if (no_1 == 49 && write_file == 1) cout << "no_1 = " << no_1 << " and no_1_Fe_sil_isotope[0] = " << no_1_Fe_sil_isotope[0] << endl;
        
         
// Store the composition in the array of mantle and core compositions as a function of accretion event (e.g. e_moles[][][])
          update_bulk_composition(no_1, e_coll[no_1],
                                  equil_moles, equil_moles_H2O, equil_moles_S_mantle, equil_moles_C_mantle,
                                  equil_moles_oxygen, equil_moles_H, equil_moles_S_core, equil_moles_C_core,
                                  noneq_moles, noneq_moles_H2O, noneq_moles_S_mantle, noneq_moles_C_mantle,
                                  noneq_moles_oxygen, noneq_moles_H, noneq_moles_S_core, noneq_moles_C_core,
                                  P, T, fo2);
          
          update_Fe_composition(no_1, e_coll[no_1],
                                equil_moles, equil_Fe_sil_isotope, equil_Fe_met_isotope,
                                noneq_moles, noneq_Fe_sil_isotope, noneq_Fe_met_isotope);
          


          update_sources(no_1, e_coll[no_1], xsil, xmet,
                         no_1_moles, no_1_moles_H2O, no_1_moles_S_mantle, no_1_moles_C_mantle,
                         no_1_moles_oxygen, no_1_moles_H, no_1_moles_S_core, no_1_moles_C_core,
                         no_2_moles, no_2_moles_H2O, no_2_moles_S_mantle, no_2_moles_C_mantle,
                         no_2_moles_oxygen, no_2_moles_H, no_2_moles_S_core, no_2_moles_C_core,
                         equil_moles, equil_moles_H2O, equil_moles_S_mantle, equil_moles_C_mantle,
                         equil_moles_oxygen, equil_moles_H, equil_moles_S_core, equil_moles_C_core,
                         noneq_moles, noneq_moles_H2O, noneq_moles_S_mantle, noneq_moles_C_mantle,
                         noneq_moles_oxygen, noneq_moles_H, noneq_moles_S_core, noneq_moles_C_core);
          

        }
        else
        {
         
          equil_Fe_sil_isotope[0] = (equil_Fe_sil_isotope[0] * equil_moles[n - m] + equil_Fe_met_isotope[0] * equil_moles[n]) / (equil_moles[n - m] + equil_moles[n]);
          equil_Fe_sil_isotope[1] = (equil_Fe_sil_isotope[1] * equil_moles[n - m] + equil_Fe_met_isotope[1] * equil_moles[n]) / (equil_moles[n - m] + equil_moles[n]);
          equil_Fe_sil_isotope[2] = (equil_Fe_sil_isotope[2] * equil_moles[n - m] + equil_Fe_met_isotope[2] * equil_moles[n]) / (equil_moles[n - m] + equil_moles[n]); //3isotope
          
          

          // All metals are oxidized and added to the silicate
          // if (write_file == 1) cout << "we're here no.2, no fractionation occuring" << endl;
          for (int i = (n - m); i < n; i++) equil_moles[i] = equil_moles[i] + equil_moles[i + m];
          for (int i = n; i < (n + m); i++) equil_moles[i] = 0.0;
          equil_moles_H2O = 0.0;
          equil_moles_S_mantle = equil_moles_S_mantle + equil_moles_S_core;
          equil_moles_C_mantle = equil_moles_C_mantle + equil_moles_C_core;
          equil_moles_oxygen = 0.0;
          equil_moles_H = 0.0;
          equil_moles_S_core = 0.0;
          equil_moles_C_core = 0.0;

          // No need to equilibrate so set equilibration conditions to zero
          P = 0.0;
          T = 0.0;
          fo2 = 0.0;
          // Store the composition in the array of mantle and core compositions as a function of accretion event (e.g. e_moles[][][])
          update_bulk_composition(no_1, e_coll[no_1],
                                  equil_moles, equil_moles_H2O, equil_moles_S_mantle, equil_moles_C_mantle,
                                  equil_moles_oxygen, equil_moles_H, equil_moles_S_core, equil_moles_C_core,
                                  noneq_moles, noneq_moles_H2O, noneq_moles_S_mantle, noneq_moles_C_mantle,
                                  noneq_moles_oxygen, noneq_moles_H, noneq_moles_S_core, noneq_moles_C_core,
                                  P, T, fo2);
          
          
          update_Fe_composition(no_1, e_coll[no_1],
                                equil_moles, equil_Fe_sil_isotope, equil_Fe_met_isotope,
                                noneq_moles, noneq_Fe_sil_isotope, noneq_Fe_met_isotope);
          
          
          update_sources(no_1, e_coll[no_1], xsil, xmet,
                         no_1_moles, no_1_moles_H2O, no_1_moles_S_mantle, no_1_moles_C_mantle,
                         no_1_moles_oxygen, no_1_moles_H, no_1_moles_S_core, no_1_moles_C_core,
                         no_2_moles, no_2_moles_H2O, no_2_moles_S_mantle, no_2_moles_C_mantle,
                         no_2_moles_oxygen, no_2_moles_H, no_2_moles_S_core, no_2_moles_C_core,
                         equil_moles, equil_moles_H2O, equil_moles_S_mantle, equil_moles_C_mantle,
                         equil_moles_oxygen, equil_moles_H, equil_moles_S_core, equil_moles_C_core,
                         noneq_moles, noneq_moles_H2O, noneq_moles_S_mantle, noneq_moles_C_mantle,
                         noneq_moles_oxygen, noneq_moles_H, noneq_moles_S_core, noneq_moles_C_core);
        }
        
        // Store compositions of metal that are being added to proto-core (composition of core-forming increments)
        for (int i = n; i < (n + m); i++) incr_moles[no_1][i][e_coll[no_1]] = e_moles[no_1][i][e_coll[no_1]] - e_moles[no_1][i][e_coll[no_1] - 1];
        incr_moles_oxygen[no_1][e_coll[no_1]]  = e_moles_oxygen[no_1][e_coll[no_1]] - e_moles_oxygen[no_1][e_coll[no_1] - 1];
        incr_moles_H[no_1][e_coll[no_1]] = e_moles_H[no_1][e_coll[no_1]] - e_moles_H[no_1][e_coll[no_1] - 1];
        incr_moles_S[no_1][e_coll[no_1]] = e_moles_S_core[no_1][e_coll[no_1]] - e_moles_S_core[no_1][e_coll[no_1] - 1];
        incr_moles_C[no_1][e_coll[no_1]] = e_moles_C_core[no_1][e_coll[no_1]] - e_moles_C_core[no_1][e_coll[no_1] - 1];

      }
      else //available moles oxygen > required moles oxygen; a special case of update_bulk_composition:
      {
        
      
        equil_Fe_sil_isotope[0] = (equil_Fe_sil_isotope[0] * equil_moles[n - m] + equil_Fe_met_isotope[0] * equil_moles[n]) / (equil_moles[n - m] + equil_moles[n]);
        equil_Fe_sil_isotope[1] = (equil_Fe_sil_isotope[1] * equil_moles[n - m] + equil_Fe_met_isotope[1] * equil_moles[n]) / (equil_moles[n - m] + equil_moles[n]);
        equil_Fe_sil_isotope[2] = (equil_Fe_sil_isotope[2] * equil_moles[n - m] + equil_Fe_met_isotope[2] * equil_moles[n]) / (equil_moles[n - m] + equil_moles[n]); //3isotope
      
        for (int i = 0; i < n - m; i++) e_moles[no_1][i][e_coll[no_1]] = no_1_moles[i] + no_2_moles[i];
        for (int i = n - m; i < n; i++) e_moles[no_1][i][e_coll[no_1]] = no_1_moles[i] + no_2_moles[i] + no_2_moles[i + m];
        if (required_moles_oxygen <= no_2_moles_oxygen) e_moles_H2O[no_1][e_coll[no_1]] = no_1_moles_H2O + no_2_moles_H2O;
        else e_moles_H2O[no_1][e_coll[no_1]] = no_1_moles_H2O + no_2_moles_H2O - required_moles_oxygen + no_2_moles_oxygen;
        e_moles_S_mantle[no_1][e_coll[no_1]] = no_1_moles_S_mantle + no_2_moles_S_mantle + no_2_moles_S_core;
        e_moles_C_mantle[no_1][e_coll[no_1]] = no_1_moles_C_mantle + no_2_moles_C_mantle + no_2_moles_C_core;

// Composition of core of no_1 remains unchanged during differentiation scenario III:
        for (int i = n; i < (n + m); i++) e_moles[no_1][i][e_coll[no_1]] = no_1_moles[i];
        e_moles_oxygen[no_1][e_coll[no_1]] = no_1_moles_oxygen;
        e_moles_H[no_1][e_coll[no_1]]  = no_1_moles_H;
        e_moles_S_core[no_1][e_coll[no_1]] = no_1_moles_S_core;
        e_moles_C_core[no_1][e_coll[no_1]] = no_1_moles_C_core;

// Store P, T, log(fO2) for each equilibration event
        e_P[no_1][e_coll[no_1]] = 0.0;
        e_T[no_1][e_coll[no_1]] = 0.0;
        e_fO2[no_1][e_coll[no_1]] = 0.0;

        
          
        update_Fe_composition(no_1, e_coll[no_1],
                              equil_moles, equil_Fe_sil_isotope, equil_Fe_met_isotope,
                              noneq_moles, noneq_Fe_sil_isotope, noneq_Fe_met_isotope);
        
        
          

        // update_source_composition (special case):
        for (int i = 0; i < n - m; i++) s_moles[no_1][i][e_coll[no_1]] = no_2_moles[i];
        for (int i = n - m; i < n; i++) s_moles[no_1][i][e_coll[no_1]] = no_2_moles[i] + no_2_moles[i + m];
        s_moles_H2O[no_1][e_coll[no_1]] = no_2_moles_H2O;
        s_moles_S_mantle[no_1][e_coll[no_1]] = no_2_moles_S_mantle + no_2_moles_S_core;
        s_moles_C_mantle[no_1][e_coll[no_1]] = no_2_moles_C_mantle + no_2_moles_C_core;

// There is no core-forming liquid metal
        for (int i = n; i < (n + m); i++) incr_moles[no_1][i][e_coll[no_1]] = 0.0;
        incr_moles_oxygen[no_1][e_coll[no_1]]  = 0.0;
        incr_moles_H[no_1][e_coll[no_1]] = 0.0;
        incr_moles_S[no_1][e_coll[no_1]] = 0.0;
        incr_moles_C[no_1][e_coll[no_1]] = 0.0;
      }
    }
// END of differentiation option III
    // if (no_1 == 12 && write_file == 1) cout << endl << "the time is: "<< e_time[no_1][e_coll[no_1]]/1e6 <<  " and this is e_fe_sil SPOT wow1 = " << no_1_Fe_sil_isotope[1] <<  " and e_Fe_sil_isotope[no_1][1][e_coll[no_1]]:" << e_Fe_sil_isotope[no_1][1][e_coll[no_1]] << " and bstatusno1: "<< b_status[no_1] << endl;

    b_status[no_1] = 2; // body no_1 now has a collisional history
    b_status[no_2] = 3; // body no_2 no longer exists
    // if (no_1 == 12 && write_file == 1) cout << endl << "the time is: "<< e_time[no_1][e_coll[no_1]]/1e6 <<  " and this is e_fe_sil SPOT wow2 = " << no_1_Fe_sil_isotope[1] <<  " and e_Fe_sil_isotope[no_1][1][e_coll[no_1]]:" << e_Fe_sil_isotope[no_1][1][e_coll[no_1]] << " and bstatusno1: "<< b_status[no_1] << endl;

    last_GI[no_1] = GI[no_1];
  }
// END of running through NB_data input file line by line
  infile3.close();

//______________________________________________________________________________________________________________________________________________
// CALCULATE UNCERTAINTIES IN FINAL PLANET AND MERGED EMBRYO COMPOSITIONS and CHI SQUARED (for planet specified to be earth-like)

  double weight[n + m + 1] = {0};
  double weight_H2O = 0.0;
  double weight_H = 0.0;
  double weight_S_mantle = 0.0;
  double weight_S_core = 0.0;
  double weight_C_mantle = 0.0;
  double weight_C_core = 0.0;

  int j(0), kk(0);
  for (kk = 1; kk < n_emb; kk++)
  {
    for (int i = 0; i < (n + m); i++)
    {
      if (b_status[kk] == 2) moles[i] = e_moles[kk][i][e_coll[kk]]; // Final terrestrial planets
      else if (b_status[kk] == 3) moles[i] = e_moles[kk][i][e_coll[kk] - 1]; // embryos that have been accreted to other proto-planets
    }
    moles_H2O = e_moles_H2O[kk][e_coll[kk]];
    moles_oxygen = e_moles_oxygen[kk][e_coll[kk]];
    moles_H = e_moles_H[kk][e_coll[kk]];
    moles_S_mantle = e_moles_S_mantle[kk][e_coll[kk]];
    moles_S_core = e_moles_S_core[kk][e_coll[kk]];
    moles_C_mantle = e_moles_C_mantle[kk][e_coll[kk]];
    moles_C_core = e_moles_C_core[kk][e_coll[kk]];

    double sum_sil = 0;
    for (int i = 0; i < n; i++)         // Determine final mantle composition
    {
      weight[i] = moles[i] * MW[i];
      sum_sil += weight[i];
    }
    weight_H2O = moles_H2O * MW_H2O;
    sum_sil += weight_H2O;
    weight_S_mantle = moles_S_mantle * AW_S;
    sum_sil += weight_S_mantle;
    weight_C_mantle = moles_C_mantle * AW_C;
    sum_sil += weight_C_mantle;

    double sum_met = 0;
    for (int i = n; i < (n + m); i++)      // Determine final core composition
    {
      weight[i] = moles[i] * AW[i - m];
      sum_met += weight[i];
    }
    weight[n + m] = moles_oxygen * AW_O;
    sum_met += weight[n + m];
    weight_H = moles_H * AW_H;
    sum_met += weight_H;
    weight_S_core = moles_S_core * AW_S;
    sum_met += weight_S_core;
    weight_C_core = moles_C_core * AW_C;
    sum_met += weight_C_core;

    if (sum_met > 0.)
    {
      for (int i = n; i < (n + m); i++)
        weight[i] = weight[i] * 100.0 / sum_met;
      weight[n + m] = weight[n + m] * 100.0 / sum_met; // oxygen concentration
      weight_H = weight_H * 100.0 / sum_met;
      weight_S_core = weight_S_core * 100.0 / sum_met;
      weight_C_core = weight_C_core * 100.0 / sum_met;
    }
    else
    {
      for (int i = n; i < (n + m); i++) weight[i] = 0.0;
      weight[n + m] = 0.0; // oxygen concentration
      weight_H = 0.0;
      weight_S_core = 0.0;
      weight_C_core = 0.0;
    }

    for (int i = 0; i < n; i++)
    {
      if (((el_sym[i][0] == 'N') && (el_sym[i][1] == 'a')) || (i > (n - m + 1)))
        weight[i] = weight[i] * (100.0 / sum_sil) * (AW[i] / MW[i]); // element concentrations
      else weight[i] = weight[i] * 100.0 / sum_sil; // oxide concentrations
    }
    weight_H2O = weight_H2O * 100.0 / sum_sil;
    weight_S_mantle = weight_S_mantle * 100.0 / sum_sil;
    weight_C_mantle = weight_C_mantle * 100.0 / sum_sil;

    if (fit == 0 && kk == moon_body && write_file == 1)
    {
      cout << setprecision(2) << weight_H2O << " wt% H2O in Moon's mantle" << endl;
      cout << setprecision(2) << weight_S_mantle << " wt% S in Moon's mantle" << endl;
      cout << setprecision(2) << weight_S_core << " wt% S in Moon's core" << endl;
      cout << setprecision(2) << weight_C_mantle << " wt% C in Moon's mantle" << endl;
      cout << setprecision(2) << weight_C_core << " wt% C in Moon's core" << endl;
    }

    // Calculate and store final uncertainty limit values:
    if (write_errors == 1)
    {
      if (error == -1)
      {
        for (int i = 0; i < n + m; i++)
          wtminus[kk][i] = weight[i];
      }
      else if (error == +1)
      {
        for (int i = 0; i < n + m; i++)
          wtplus[kk][i] = weight[i];
      }
      else if (error == 0)
      {
        for (int i = 0; i < n + m; i++)
          wtbest[kk][i] = weight[i];
      }
      for (int i = 0; i < n + m; i++)
      {
        err[kk][i] = (fabs(wtbest[kk][i] - wtminus[kk][i]) + fabs(wtbest[kk][i] - wtplus[kk][i])) / 2.0;
        if ((i == 3) || (i > 5 && i < n)) // Na and elements beyond Si as elemental concentrations, others as oxides
          err[kk][i] = err[kk][i] * AW[i] / MW[i];
      }
      wtbest_H2O[kk] = weight_H2O;
      wtbest_S_mantle[kk] = weight_S_mantle;
      wtbest_C_mantle[kk] = weight_C_mantle;
      wtbest_O[kk] = weight[n + m];
      wtbest_S_core[kk] = weight_S_core;
      wtbest_C_core[kk] = weight_C_core;
      wtbest_H[kk] = weight_H;
    }

    // CHI SQUARED CALCULATION:
    int num_cons = 0;

    // FIT ONLY EARTH: note: this works now, it is exactly like the original.
    if (fit_planet == 0)
      // CHI SQUARED CALCULATION INCLUDES EARTH PLANET
    {
      if ((write_file == 1) && (kk == moon_body))
        cout << endl << endl << "ELEMENTS USED IN CHI SQUARE CALCULATION:" << endl;
      if (kk == moon_body)
      {
        chi_sq = 0.0;
        for (int i = 0; i < n; i++)
        {
          double old_chisq(0); //added to see individual chi contribution
          if (fit_elem[i] == 'y')
          {
            old_chisq = chi_sq; //added to see individual chi contribution
            chi_sq += pow(((p_mantle[i] - weight[i]) / (p_mantle[i] * sig_p_mantle[i])), 2.0);
            num_cons += 1;
            if (write_file == 1) for (int j = 0; j < 2; j++) cout << el_sym[i][j];
            if (write_file == 1) cout << " contributes " << (chi_sq - old_chisq) / 4.0; //added to see individual chi contribution
            if (write_file == 1) cout << endl;
          }
        }

        double nb_ta(0); // Contribution of Nb/Ta ratio to chi squared
        int nb(0);
        double old_chisq(0); //added to see individual chi contribution
        for (int i = 0; i < n; i++)
        {
          if ((el_sym[i][0] == 'N') && (el_sym[i][1] == 'b'))
            nb = i;
        }
        if (fit_elem[nb] == 'y' && fit_elem[nb + 1] == 'y')
        {
          nb_ta = (weight[nb] / weight[nb + 1]);
          old_chisq = chi_sq; //added to see individual chi contribution
          chi_sq += pow(((14.0 - nb_ta) / 0.3), 2.0); // Mantle value is 14.0 +/-0.3 (Mnker et al)
          num_cons += 1;
          if (write_file == 1) cout << "Nb_Ta contributes " << (chi_sq - old_chisq) / 4.0; //added to see individual chi contribution
          if (write_file == 1) cout << endl; //added to see individual chi contribution
        }

        /*
        chi_sq += pow(((0.0715 - weight_H2O) / 0.01), 2.0);
        num_cons += 1;
        if (write_file == 1) cout << "H2O" << endl;
        */

        if (fit_cond == 2)
        {
          //chi_sq += pow(((0.0225 - weight_S_mantle) / 0.0050), 2.0);
          num_cons += 1;
          if (write_file == 1) cout << "S" << endl;

          //chi_sq += pow(((1.7 - weight_S_core) / 0.3), 2.0);
          num_cons += 1;
          if (write_file == 1) cout << "S in core" << endl;

          /*        chi_sq += pow(((0.0075 - weight_C_mantle) / 0.0020), 2.0);
                  num_cons += 1;
                  if (write_file == 1) cout << "C" << endl;

                  chi_sq += pow(((0.2 - weight_C_core) / 0.05), 2.0);
                  num_cons += 1;
                  if (write_file == 1) cout << "and C in core" << endl;*/
        }

        chi_sq = chi_sq * (1 + chi_sq_multiplier);
        int dof;
        if (num_cons - NP > 0) dof = num_cons - NP;
        else dof = 1;
        chi_sq = chi_sq / dof; // divide by dof if you want
        if (write_file == 1 || fit == 1) {
          cout << "chi sq = " << fixed << setprecision(2) << chi_sq << endl;
        }
        if (write_file == 1) cout << "red chi sq = " << fixed << setprecision(2) << chi_sq << tab << "#cons = " <<  num_cons << tab << "#param = " << NP << tab << "#dof = " << dof << endl;
      }
    }
    // FIT ONLY MARS: note: does not include fit_cond==2 or nb/ta ratio
    else if (fit_planet == 1) // THIS DOESNT WORK IN THIS noMarsNumber version!!!!!!!!!!!!!!
      // commented out includes mars, uncommented is the same as fit_planet==0
    {
      
      if ((write_file == 1) && (kk == moon_body))
        cout << endl << endl << "ELEMENTS USED IN CHI SQUARE CALCULATION:" << endl;
      if (kk == moon_body)
      {
        chi_sq = 0.0;
        for (int i = 0; i < n; i++)
        {
          double old_chisq(0); //added to see individual chi contribution
          if (fit_elem[i] == 'y')
          {
            old_chisq = chi_sq; //added to see individual chi contribution
            chi_sq += pow(((p_mantle[i] - weight[i]) / (p_mantle[i] * sig_p_mantle[i])), 2.0);
            num_cons += 1;
            if (write_file == 1) for (int j = 0; j < 2; j++) cout << el_sym[i][j];
            if (write_file == 1) cout << " contributes " << (chi_sq - old_chisq) / 4.0; //added to see individual chi contribution
            if (write_file == 1) cout << endl;
          }
        }

        double nb_ta(0); // Contribution of Nb/Ta ratio to chi squared
        int nb(0);
        double old_chisq(0); //added to see individual chi contribution
        for (int i = 0; i < n; i++)
        {
          if ((el_sym[i][0] == 'N') && (el_sym[i][1] == 'b'))
            nb = i;
        }
        if (fit_elem[nb] == 'y' && fit_elem[nb + 1] == 'y')
        {
          nb_ta = (weight[nb] / weight[nb + 1]);
          old_chisq = chi_sq; //added to see individual chi contribution
          chi_sq += pow(((14.0 - nb_ta) / 0.3), 2.0); // Mantle value is 14.0 +/-0.3 (Mnker et al)
          num_cons += 1;
          if (write_file == 1) cout << "Nb_Ta contributes " << (chi_sq - old_chisq) / 4.0; //added to see individual chi contribution
          if (write_file == 1) cout << endl; //added to see individual chi contribution
        }

        /*
        chi_sq += pow(((0.0715 - weight_H2O) / 0.01), 2.0);
        num_cons += 1;
        if (write_file == 1) cout << "H2O" << endl;
        */

        if (fit_cond == 2)
        {
          //chi_sq += pow(((0.0225 - weight_S_mantle) / 0.0050), 2.0);
          num_cons += 1;
          if (write_file == 1) cout << "S" << endl;

          //chi_sq += pow(((1.7 - weight_S_core) / 0.3), 2.0);
          num_cons += 1;
          if (write_file == 1) cout << "S in core" << endl;

          /*        chi_sq += pow(((0.0075 - weight_C_mantle) / 0.0020), 2.0);
                  num_cons += 1;
                  if (write_file == 1) cout << "C" << endl;

                  chi_sq += pow(((0.2 - weight_C_core) / 0.05), 2.0);
                  num_cons += 1;
                  if (write_file == 1) cout << "and C in core" << endl;*/
        }

        chi_sq = chi_sq * (1 + chi_sq_multiplier);
        int dof;
        if (num_cons - NP > 0) dof = num_cons - NP;
        else dof = 1;
        chi_sq = chi_sq / dof; // divide by dof if you want
        if (write_file == 1 || fit == 1) cout << "chi sq = " << fixed << setprecision(2) << chi_sq << endl;
        if (write_file == 1) cout << "red chi sq = " << fixed << setprecision(2) << chi_sq << tab << "#cons = " <<  num_cons << tab << "#param = " << NP << tab << "#dof = " << dof << endl;
      }
      // if ((write_file == 1) && (kk == mars_planet)) //noMarsNumber
      //   cout << endl << endl << "ELEMENTS USED IN CHI SQUARE CALCULATION:" << endl;
      // if (kk == mars_planet)
      // {

      //   chi_sq = 0.0;
      //   for (int i = 0; i < n; i++)
      //   {
      //     if (fit_elem[i] == 'y')
      //     {
      //       chi_sq += pow(((m_mantle[i] - weight[i]) / (m_mantle[i] * sig_m_mantle[i])), 2.0);
      //       num_cons += 1;
      //       if (write_file == 1) for (int j = 0; j < 2; j++) cout << el_sym[i][j];
      //       if (write_file == 1) cout << tab << fixed << setprecision(2) << pow(((m_mantle[i] - weight[i]) / (m_mantle[i] * sig_m_mantle[i])), 2.0);
      //       if (write_file == 1) cout << endl;
      //     }
      //   }

      //   // not including because i can't find a good value for it in yoshizaki or taylor
      //   // double nb_ta(0); // Contribution of Nb/Ta ratio to chi squared
      //   // int nb(0);
      //   // double old_chisq(0); //added to see individual chi contribution
      //   // for (int i = 0; i < n; i++)
      //   // {
      //   //   if ((el_sym[i][0] == 'N') && (el_sym[i][1] == 'b'))
      //   //     nb = i;
      //   // }
      //   // if (fit_elem[nb] == 'y' && fit_elem[nb + 1] == 'y')
      //   // {
      //   //   nb_ta = (weight[nb] / weight[nb + 1]);
      //   //   old_chisq = chi_sq; //added to see individual chi contribution
      //   //   chi_sq += pow(((14.0 - nb_ta) / 0.3), 2.0); // Mantle value is 14.0 +/-0.3 (Mnker et al)
      //   //   num_cons += 1;
      //   //   if (write_file == 1) cout << "Nb_Ta contributes " << (chi_sq - old_chisq)/4.0; //added to see individual chi contribution
      //   //   if (write_file == 1) cout << endl; //added to see individual chi contribution
      //   // }


      //   chi_sq = chi_sq * (1 + chi_sq_multiplier);
      //   int dof;
      //   if (num_cons - NP > 0) dof = num_cons - NP;
      //   else dof = 1;
      //   if (write_file == 1 || fit == 1) cout << "chi sq = " << fixed << setprecision(2) << chi_sq << endl;
      //   chi_sq = chi_sq / dof; // divide by dof if you want
      //   if (write_file == 1) cout << "red chi sq = " << fixed << setprecision(2) << chi_sq << tab << "#cons = " <<  num_cons << tab << "#param = " << NP << tab << "#dof = " << dof << endl;
    }
    // FIT BOTH EARTH AND MARS: no, it doesnt do this anymore, commented out. 
    else if (fit_planet == 2)// THIS DOESNT WORK IN THIS noMarsNumber version!!!!!!!!!!!!!!
      // commented out includes mars, uncommented is the same as fit_planet==0
    {

      if ((write_file == 1) && (kk == moon_body))
        cout << endl << endl << "ELEMENTS USED IN CHI SQUARE CALCULATION:" << endl;
      if (kk == moon_body)
      {
        chi_sq = 0.0;
        for (int i = 0; i < n; i++)
        {
          double old_chisq(0); //added to see individual chi contribution
          if (fit_elem[i] == 'y')
          {
            old_chisq = chi_sq; //added to see individual chi contribution
            chi_sq += pow(((p_mantle[i] - weight[i]) / (p_mantle[i] * sig_p_mantle[i])), 2.0);
            num_cons += 1;
            if (write_file == 1) for (int j = 0; j < 2; j++) cout << el_sym[i][j];
            if (write_file == 1) cout << " contributes " << (chi_sq - old_chisq) / 4.0; //added to see individual chi contribution
            if (write_file == 1) cout << endl;
          }
        }

        double nb_ta(0); // Contribution of Nb/Ta ratio to chi squared
        int nb(0);
        double old_chisq(0); //added to see individual chi contribution
        for (int i = 0; i < n; i++)
        {
          if ((el_sym[i][0] == 'N') && (el_sym[i][1] == 'b'))
            nb = i;
        }
        if (fit_elem[nb] == 'y' && fit_elem[nb + 1] == 'y')
        {
          nb_ta = (weight[nb] / weight[nb + 1]);
          old_chisq = chi_sq; //added to see individual chi contribution
          chi_sq += pow(((14.0 - nb_ta) / 0.3), 2.0); // Mantle value is 14.0 +/-0.3 (Mnker et al)
          num_cons += 1;
          if (write_file == 1) cout << "Nb_Ta contributes " << (chi_sq - old_chisq) / 4.0; //added to see individual chi contribution
          if (write_file == 1) cout << endl; //added to see individual chi contribution
        }

        /*
        chi_sq += pow(((0.0715 - weight_H2O) / 0.01), 2.0);
        num_cons += 1;
        if (write_file == 1) cout << "H2O" << endl;
        */

        if (fit_cond == 2)
        {
          //chi_sq += pow(((0.0225 - weight_S_mantle) / 0.0050), 2.0);
          num_cons += 1;
          if (write_file == 1) cout << "S" << endl;

          //chi_sq += pow(((1.7 - weight_S_core) / 0.3), 2.0);
          num_cons += 1;
          if (write_file == 1) cout << "S in core" << endl;

          /*        chi_sq += pow(((0.0075 - weight_C_mantle) / 0.0020), 2.0);
                  num_cons += 1;
                  if (write_file == 1) cout << "C" << endl;

                  chi_sq += pow(((0.2 - weight_C_core) / 0.05), 2.0);
                  num_cons += 1;
                  if (write_file == 1) cout << "and C in core" << endl;*/
        }

        chi_sq = chi_sq * (1 + chi_sq_multiplier);
        int dof;
        if (num_cons - NP > 0) dof = num_cons - NP;
        else dof = 1;
        chi_sq = chi_sq / dof; // divide by dof if you want
        if (write_file == 1 || fit == 1) cout << "chi sq = " << fixed << setprecision(2) << chi_sq << endl;
        if (write_file == 1) cout << "red chi sq = " << fixed << setprecision(2) << chi_sq << tab << "#cons = " <<  num_cons << tab << "#param = " << NP << tab << "#dof = " << dof << endl;
      }

      // CHI SQUARED CALCULATION INCLUDES BOTH EARTH AND MARS
      // if ((write_file == 1) && (kk == earth_planet))
      //   cout << endl << endl << "Earth ELEMENTS USED IN CHI SQUARE CALCULATION:" << endl;
      // if (kk == earth_planet)
      // {
      //   chi_sq = 0.0;
      //   for (int i = 0; i < n; i++)
      //   {
      //     if (fit_elem[i] == 'y')
      //     {
      //       chi_sq += pow(((p_mantle[i] - weight[i]) / (p_mantle[i] * sig_p_mantle[i])), 2.0);
      //       num_cons += 1;
      //       if (write_file == 1) for (int j = 0; j < 2; j++) cout << el_sym[i][j];
      //       if (write_file == 1) cout << endl;
      //     }
      //   }

      //   double nb_ta(0); // Contribution of Nb/Ta ratio to chi squared
      //   int nb(0);
      //   for (int i = 0; i < n; i++)
      //   {
      //     if ((el_sym[i][0] == 'N') && (el_sym[i][1] == 'b'))
      //       nb = i;
      //   }
      //   if (fit_elem[nb] == 'y' && fit_elem[nb + 1] == 'y')
      //   {
      //     nb_ta = (weight[nb] / weight[nb + 1]);
      //     chi_sq += pow(((14.0 - nb_ta) / 0.3), 2.0); // Mantle value is 14.0 +/-0.3 (Mnker et al)
      //     num_cons += 1;
      //   }

      //   /*
      //   chi_sq += pow(((0.0715 - weight_H2O) / 0.01), 2.0);
      //   num_cons += 1;
      //   if (write_file == 1) cout << "H2O" << endl;
      //   */

      //   if (fit_cond == 2)
      //   {
      //     //chi_sq += pow(((0.0225 - weight_S_mantle) / 0.0050), 2.0);
      //     num_cons += 1;
      //     if (write_file == 1) cout << "S" << endl;

      //     //chi_sq += pow(((1.7 - weight_S_core) / 0.3), 2.0);
      //     num_cons += 1;
      //     if (write_file == 1) cout << "S in core" << endl;

      //     /*        chi_sq += pow(((0.0075 - weight_C_mantle) / 0.0020), 2.0);
      //             num_cons += 1;
      //             if (write_file == 1) cout << "C" << endl;

      //             chi_sq += pow(((0.2 - weight_C_core) / 0.05), 2.0);
      //             num_cons += 1;
      //             if (write_file == 1) cout << "and C in core" << endl;*/
      //   }
      // }
      // if ((write_file == 1) && (kk == mars_planet))
      //   cout << endl << endl << "Mars ELEMENTS USED IN CHI SQUARE CALCULATION:" << endl;
      // if (kk == mars_planet)
      // {
      //   for (int i = 0; i < n; i++)
      //   {
      //     if (fit_elem[i] == 'y')
      //     {
      //       chi_sq += pow(((m_mantle[i] - weight[i]) / (m_mantle[i] * sig_m_mantle[i])), 2.0);
      //       num_cons += 1;
      //       if (write_file == 1) for (int j = 0; j < 2; j++) cout << el_sym[i][j];
      //       if (write_file == 1) cout << endl;
      //     }
      //   }

      //   chi_sq = chi_sq * (1 + chi_sq_multiplier);
      //   int dof;
      //   if (num_cons - NP > 0) dof = num_cons - NP;
      //   else dof = 1;
      //   chi_sq = chi_sq / dof; // divide by dof if you want
      //   if (write_file == 1 || fit == 1) cout << "chi sq = " << fixed << setprecision(5) << chi_sq << endl;
      //   if (write_file == 1) cout << "red chi sq = " << fixed << setprecision(5) << chi_sq << tab << "#cons = " <<  num_cons << tab << "#param = " << NP << tab << "#dof = " << dof << endl;

      // }
    }
    // FIT ONLY MOON (Moon v1)
    // Here, earth number is referring to the
    else if (fit_planet == 3)
      // CHI SQUARED CALCULATION INCLUDES EARTH PLANET
    {
      if ((write_file == 1) && (kk == moon_body)) // moon_body is actually the moon in this case
        cout << endl << endl << "ELEMENTS USED IN CHI SQUARE CALCULATION:" << endl;
      if (kk == moon_body)
      {
        for (int i = 0; i < n; i++) { // Calculate the value of num_cons
            if (fit_elem[i] == 'y') 
            {
              num_cons += 1;
            }  
          }
        if (fit_cond == 3) num_cons += 1; // Add additional num_con for lunar core size fitting, Moon ver. 
        
        int dof; // initialize dof value
        
        chi_sq = 0.0;
        for (int i = 0; i < n; i++)
        {
          double old_chisq(0); //added to see individual chi contribution
          if (fit_elem[i] == 'y')
          {
            if (num_cons - NP > 0) dof = num_cons - NP; // calculate dof value (originally was in loop below at chi_sq = chi_sq * (1+mult))
            else dof = 1;
            old_chisq = chi_sq; //added to see individual chi contribution
            chi_sq += pow(((moon_mantle[i] - weight[i]) / (moon_mantle[i] * sig_moon_mantle[i])), 2.0);
            if (write_file == 1) for (int j = 0; j < 2; j++) cout << el_sym[i][j];
            if (write_file == 1) cout << " contributes " << (chi_sq - old_chisq) / dof;
            if (write_file == 1) cout << endl;
          }
        }

        double nb_ta(0); // Contribution of Nb/Ta ratio to chi squared
        int nb(0);
        double old_chisq(0); //added to see individual chi contribution
        for (int i = 0; i < n; i++)
        {
          if ((el_sym[i][0] == 'N') && (el_sym[i][1] == 'b'))
            nb = i;
        }
        if (fit_elem[nb] == 'y' && fit_elem[nb + 1] == 'y')
        {
          nb_ta = (weight[nb] / weight[nb + 1]);
          old_chisq = chi_sq; //added to see individual chi contribution
          chi_sq += pow(((14.0 - nb_ta) / 0.3), 2.0); // Mantle value is 14.0 +/-0.3 (Mnker et al)
          num_cons += 1;
          if (write_file == 1) cout << "Nb_Ta contributes " << (chi_sq - old_chisq) / dof; //added to see individual chi contribution
          if (write_file == 1) cout << endl; //added to see individual chi contribution
        }
        // New chi sq contribution condition for the mass of the core
        // Moon ver.
        // Core mass if 2.0% +- 1.0% from Gaffney et al. (2023)
        if (fit_cond == 3) 
        {
          double metal_mass = 0.0; 
          double lunar_core_mass = 2.0; // Lunar masses
          double lunar_core_mass_err = 1.0; // Lunar masses
          metal_mass = sum_met * 100 / (sum_met + sum_sil); // calculate metal mass using same formula as metal mass calculation above
          old_chisq = chi_sq;
          chi_sq += pow((lunar_core_mass - metal_mass)/ lunar_core_mass_err, 2.0);
          if (write_file == 1) cout << "Lunar core mass contributes " << (chi_sq - old_chisq)/ dof << endl;
        }
        
        /*
        chi_sq += pow(((0.0715 - weight_H2O) / 0.01), 2.0);
        num_cons += 1;
        if (write_file == 1) cout << "H2O" << endl;
        */

        if (fit_cond == 2)
        {
          //chi_sq += pow(((0.0225 - weight_S_mantle) / 0.0050), 2.0);
          num_cons += 1;
          if (write_file == 1) cout << "S" << endl;

          //chi_sq += pow(((1.7 - weight_S_core) / 0.3), 2.0);
          num_cons += 1;
          if (write_file == 1) cout << "S in core" << endl;

          /*        chi_sq += pow(((0.0075 - weight_C_mantle) / 0.0020), 2.0);
                  num_cons += 1;
                  if (write_file == 1) cout << "C" << endl;

                  chi_sq += pow(((0.2 - weight_C_core) / 0.05), 2.0);
                  num_cons += 1;
                  if (write_file == 1) cout << "and C in core" << endl;*/
        }

        chi_sq = chi_sq * (1 + chi_sq_multiplier);

        chi_sq = chi_sq / dof; // divide by dof if you want
        if (write_file == 1 || fit == 1) {
          cout << "chi sq = " << fixed << setprecision(2) << chi_sq << endl;
        }
        if (write_file == 1) cout << "red chi sq = " << fixed << setprecision(2) << chi_sq << tab << "#cons = " <<  num_cons << tab << "#param = " << NP << tab << "#dof = " << dof << endl;
      }
    }
  }
//_________________________________________________________________________

  if (write_file == 1)
  {
    cout << endl << "FINAL PLANETS" << endl;
    for (int i = 1; i < n_emb; i++) // OUTPUT FINAL PLANET DETAILS TO SCREEN
      if (b_status[i] == 2)
      {
        cout << "planet# = " << i << tab << "No of collisions = " << e_coll[i] << tab << "Mass = "
             << e_mass[i][e_coll[i]];
        cout << endl;
      }

// READ FILE "final_planets.txt" to obtain final distances and to order planets by increasing distance in outfiles
    int no_planets = 0;
    string sim_title;
    int p_index[n_emb] = {0};
    int p_no[n_emb] = {0};
    double p_dist[n_emb] = {0};
    float mass;
    int jj = 0;

    ifstream infile4("final_planets.txt"); // declaration of 4th input file
    if (infile4.is_open())
    {
      getline (infile4, sim_title); // Reads simulation # as string
      infile4 >> no_planets; // Final number of planets in the simulation
      for (jj = 0; jj < no_planets; jj++)
      {
        infile4 >> p_index[jj] >> p_no[jj] >> mass >> p_dist[jj];
        getline (infile4, line);
        planet_no[jj] = p_no[jj];
      }
    }

    infile4.close();

    // OUTPUT FINAL PLANET COMPOSITION DETAILS TO FILES
    outfile1 << "XSIL FOR PLANETESIMAL COLLISIONS (min, av., max):" << tab << xsil_pl_min << ", " << xsil_pl / n_pl << ", " << xsil_pl_max << endl;
    outfile1 << "XSIL FOR EMBRYO COLLISIONS (min, av., max):" << tab << tab << xsil_em_min << ", " << xsil_em / n_em << ", " << xsil_em_max << endl;
    outfile1 << "Numbers of (a) planetesimal-embryo and (b) embryo-embryo collisions: (a) " << n_pl << tab << "(b) " << n_em << endl;
    outfile1 << "Number of collisions either without equilibration or without metal being added: " << n0_pl << endl;
    outfile2 << "XSIL FOR PLANETESIMAL COLLISIONS (min, av., max):" << tab << xsil_pl_min << ", " << xsil_pl / n_pl << ", " << xsil_pl_max << endl;
    outfile2 << "XSIL FOR EMBRYO COLLISIONS (min, av., max):" << tab << tab << xsil_em_min << ", " << xsil_em / n_em << ", " << xsil_em_max << endl;
    outfile2 << "Numbers of (a) planetesimal-embryo and (b) embryo-embryo collisions: (a) " << n_pl << tab << "(b) " << n_em << endl;
    outfile2 << "Number of collisions either without equilibration or without metal being added: " << n0_pl << endl;

    for (jj = 0; jj < no_planets; jj++)   // Runs through each successive planet
    {
      outfile1 << endl;
      outfile2 << endl;
      outfile3 << endl;
      outfile6 << endl;
      outfile7 << endl;
      outfile8 << endl;
      outfile9 << endl;
      outfile12 << endl;
      // outfile13 << endl;
      outfile1 << "PLANET #" << p_no[jj] << endl;
      outfile1 << "Final Distance = " << p_dist[jj] << " Earth radii" << endl;
      outfile2 << "PLANET #" << p_no[jj] << endl;
      outfile2 << "Final Distance  = " << p_dist[jj] << " Earth radii" << endl;
      outfile3 << "PLANET #" << p_no[jj] << endl;
      outfile3 << "Final Distance = " << p_dist[jj] << " Earth radii" << endl;
      outfile6 << "PLANET #" << p_no[jj] << endl;
      outfile6 << "Final Distance = " << p_dist[jj] << " Earth radii" << endl;
      outfile7 << "PLANET #" << p_no[jj] << endl;
      outfile7 << "Final Distance = " << p_dist[jj] << " Earth radii" << endl;
      outfile8 << "PLANET #" << p_no[jj] << endl;
      outfile8 << "Final Distance = " << p_dist[jj] << " Earth radii" << endl;
      outfile9 << "PLANET #" << p_no[jj] << endl;
      outfile9 << "Final Distance = " << p_dist[jj] << " Earth radii" << endl;
      outfile12 << "PLANET #" << p_no[jj] << endl;
      outfile12 << "Final Distance = " << p_dist[jj] << " Earth radii" << endl;
      // if (jj == 27)
      // {
      //   outfile13 << "PLANET #" << p_no[jj] << endl;
      //   outfile13 << "Final Distance = " << p_dist[jj] << " AU" << endl;
      // }
      // outfile13 << "PLANET #" << p_no[jj] << endl;
      // outfile13 << "Final Distance = " << p_dist[jj] << " AU" << endl;

      tm_col = 1;
      write_headings();
      int imp(0);
      j = p_no[jj];    // Corresponds to the surviving embryo numbers
      // kk is simulation and jj is planet

      double core_moles[n + m] = {0.0, 0.0};
      double core_moles_oxygen = 0.0;
      double core_moles_H = {0.0};
      double core_moles_S = {0.0};
      double core_moles_C = {0.0};

      for (imp = 0; imp < (e_coll[j] + 1); imp++) // Runs through each successive collision for planet j
      { // summing up core-forming components - in order to later determine incremental mass fractions
        for (int i = n; i < (n + m); i++) core_moles[i] += incr_moles[j][i][imp];
        core_moles_oxygen += incr_moles_oxygen[j][imp];
        core_moles_S += incr_moles_S[j][imp];
        core_moles_C += incr_moles_C[j][imp];
        core_moles_H += incr_moles_H[j][imp];
      }

      double sum_core = 0.0;  // Determine final mass for each core
      for (int i = n; i < (n + m); i++) sum_core += core_moles[i] * AW[i - m];
      sum_core += core_moles_oxygen * AW_O; // Adds oxygen
      sum_core += core_moles_H * AW_H;
      sum_core += core_moles_S * AW_S; // Adds S
      sum_core += core_moles_C * AW_C; // Adds S

      double total_moles[n + m] = {0.0};
      double total_moles_H2O = 0.0;
      double total_moles_oxygen = 0.0;
      double total_moles_H = 0.0;
      double total_moles_S_mantle = 0.0;
      double total_moles_S_core = 0.0;
      double total_moles_C_mantle = 0.0;
      double total_moles_C_core = 0.0;
      for (imp = 0; imp < (e_coll[j] + 1); imp++)
      {
        for (int i = 0; i < (n + m); i++) total_moles[i] += s_moles[j][i][imp];
        total_moles_H2O += s_moles_H2O[j][imp];
        total_moles_oxygen += s_moles_oxygen[j][imp];
        total_moles_H += s_moles_H[j][imp];
        total_moles_S_mantle += s_moles_S_mantle[j][imp];
        total_moles_S_core += s_moles_S_core[j][imp];
        total_moles_C_mantle += s_moles_C_mantle[j][imp];
        total_moles_C_core += s_moles_C_core[j][imp];
      }

      for (imp = 0; imp < (e_coll[j] + 1); imp++) // RUNS THROUGH EACH SUCCESSIVE COLLISION FOR PLANET j, and  outputs composition
      {

        outfile1 << setw(10) << left << format(e_time[j][imp] / 1.e6, 8) << tab; // Moon_ver ; origin was set(8) format ( , 6)
        outfile1 << setw(8) << left << format(e_mass[j][imp], 6) << tab;

        outfile2 << setw(8) << left << format(e_time[j][imp] / 1.e6, 6) << tab;
        outfile2 << setw(8) << left << format(e_mass[j][imp], 6) << tab;
        outfile2 << setw(8) << left << format(e_P[j][imp], 6) << tab;
        outfile2 << setw(8) << left << format(e_T[j][imp], 6) << tab;
        outfile2 << setw(8) << left << format(e_fO2[j][imp], 6) << tab;

        for (int i = 0; i < (n + m); i++) moles[i] = e_moles[j][i][imp];
        moles_H2O = e_moles_H2O[j][imp];
        moles_oxygen = e_moles_oxygen[j][imp];
        moles_H = e_moles_H[j][imp];
        moles_S_mantle = e_moles_S_mantle[j][imp];
        moles_S_core = e_moles_S_core[j][imp];
        moles_C_mantle = e_moles_C_mantle[j][imp];
        moles_C_core = e_moles_C_core[j][imp];

        weight_comp(moles, moles_oxygen, moles_H2O, moles_H, moles_S_mantle, moles_S_core, moles_C_mantle, moles_C_core, AW, MW);

        // Now output compositions of core-formation increments
        if (e_fO2[j][imp] != 0.0) // Zero fO2 value indicates no core forming event
        {
          outfile3 << setw(8) << left << format(e_time[j][imp] / 1.e6, 6) << tab;
          outfile3 << setw(8) << left << format(e_mass[j][imp], 6) << tab;
          outfile3 << setw(8) << left << format(e_P[j][imp], 6) << tab;
          outfile3 << setw(8) << left << format(e_T[j][imp], 6) << tab;
          outfile3 << setw(8) << left << format(e_fO2[j][imp], 6) << tab;

          for (int i = n; i < (n + m); i++) moles[i] = incr_moles[j][i][imp];
          moles_oxygen = incr_moles_oxygen[j][imp];
          moles_H = incr_moles_H[j][imp];
          moles_S_core = incr_moles_S[j][imp];
          moles_C_core = incr_moles_C[j][imp];
          incr_comp(moles, moles_oxygen, moles_H, moles_S_core, moles_C_core, sum_core, AW);
        }

        if (imp > 0 )
        {
          outfile6 << setw(8) << left << format(e_time[j][imp] / 1.e6, 6) << tab;
          outfile6 << setw(8) << left << format(1.e2 * i_mass[j][imp], 6) << tab;
          outfile6 << setw(8) << left << format(1.e6 * e_equil[j][imp], 6) << tab;

          outfile7 << setw(8) << left << format(e_time[j][imp] / 1.e6, 6) << tab;
          outfile7 << setw(8) << left << format(1.e2 * i_mass[j][imp], 6) << tab;
          outfile7 << setw(8) << left << format(1.e2 * i_equil[j][imp], 6) << tab;

          for (int i = 0; i < (n + m); i++) moles[i] = i_moles[j][i][imp];
          moles_H2O = i_moles_H2O[j][imp];
          moles_oxygen = i_moles_oxygen[j][imp];
          moles_H = i_moles_H[j][imp];
          moles_S_mantle = i_moles_S_mantle[j][imp];
          moles_S_core = i_moles_S_core[j][imp];
          moles_C_mantle = i_moles_C_mantle[j][imp];
          moles_C_core = i_moles_C_core[j][imp];

          weight_comp_imp(moles, moles_oxygen, moles_H2O, moles_H, moles_S_mantle, moles_S_core, moles_C_mantle, moles_C_core, AW, MW);
        }

        outfile8 << setw(8) << left << format(e_time[j][imp] / 1.e6, 6) << tab;
        outfile8 << setw(8) << left << format(e_mass[j][imp], 6) << tab;

        outfile9 << setw(8) << left << format(e_time[j][imp] / 1.e6, 6) << tab;
        outfile9 << setw(8) << left << format(e_mass[j][imp], 6) << tab;

        for (int i = 0; i < (n + m); i++) moles[i] = s_moles[j][i][imp];
        moles_H2O = s_moles_H2O[j][imp];
        moles_oxygen = s_moles_oxygen[j][imp];
        moles_H = s_moles_H[j][imp];
        moles_S_mantle = s_moles_S_mantle[j][imp];
        moles_S_core = s_moles_S_core[j][imp];
        moles_C_mantle = s_moles_C_mantle[j][imp];
        moles_C_core = s_moles_C_core[j][imp];

        moles_comp_source(moles, moles_oxygen, moles_H2O, moles_H, moles_S_mantle, moles_S_core, moles_C_mantle, moles_C_core,
                          total_moles, total_moles_oxygen, total_moles_H2O, total_moles_H, total_moles_S_mantle, total_moles_S_core, total_moles_C_mantle, total_moles_C_core);

        // write to COMPOSITION_FE.txt:
        outfile12 << setw(8) << left << format(e_time[j][imp] / 1.e6, 6) << tab; // time
        outfile12 << setw(8) << left << format(e_mass[j][imp], 6) << tab; // mass

         // if (no_1 == 12 && write_file == 1 && e_time[no_1][e_coll[no_1]]/1e6 < 0.37 && e_time[no_1][e_coll[no_1]]/1e6 > 0.337  ) cout << endl << "the time is: "<< e_time[no_1][e_coll[no_1]]/1e6 <<  " and this is e_fe_sil SPOT 3a = " << no_1_Fe_sil_isotope[1] <<  " and e_Fe_sil_isotope[no_1][1][e_coll[no_1]]:" << e_Fe_sil_isotope[no_1][1][e_coll[no_1]] << endl;

        outfile12 << setw(24) << left << format(e_Fe_sil_isotope[j][0][imp], 18) << tab; // Fe54_sil
        outfile12 << setw(24) << left << format(e_Fe_sil_isotope[j][1][imp], 18) << tab;// Fe56_sil
        outfile12 << setw(24) << left << format(e_Fe_sil_isotope[j][2][imp], 18) << tab;// Fe57_sil

        outfile12 << setw(24) << left << format(e_Fe_met_isotope[j][0][imp], 18) << tab; // Fe54_met
        outfile12 << setw(24) << left << format(e_Fe_met_isotope[j][1][imp], 18) << tab; // Fe56_met
        outfile12 << setw(24) << left << format(e_Fe_met_isotope[j][2][imp], 18) << tab; // Fe57_met

        
        double mantle_d56 = ((e_Fe_sil_isotope[j][1][imp] / e_Fe_sil_isotope[j][0][imp]) / (Fe56_IRMM / Fe54_IRMM) - 1.0) * 1000.0 ;
        outfile12 << setw(13) << left << format(mantle_d56, 6) << tab; // mantle_d56
        
        double core_d56 = ((e_Fe_met_isotope[j][1][imp] / e_Fe_met_isotope[j][0][imp]) / (Fe56_IRMM / Fe54_IRMM) - 1.0) * 1000.0;
        outfile12 << setw(13) << left << format(core_d56, 6) << tab;// core_d56
        
        double DELTA_56 = mantle_d56 - core_d56;
        outfile12 << setw(13) << left << format(DELTA_56, 6) << tab; // DELTA_56
        
        
        double mantle_d57 = ((e_Fe_sil_isotope[j][2][imp] / e_Fe_sil_isotope[j][0][imp]) / (Fe57_IRMM / Fe54_IRMM) - 1.0) * 1000.0 ;
        outfile12 << setw(13) << left << format(mantle_d57, 6) << tab; // mantle_d57
        
        double core_d57 = ((e_Fe_met_isotope[j][2][imp] / e_Fe_met_isotope[j][0][imp]) / (Fe57_IRMM / Fe54_IRMM) - 1.0) * 1000.0;
        outfile12 << setw(13) << left << format(core_d57, 6) << tab; // core_d57
        
        double DELTA_57 = mantle_d57 - core_d57;
        outfile12 << setw(13) << left << format(DELTA_57, 6) << tab; // DELTA_57

        outfile12 << endl;

        // if (no_1 == 12 && write_file == 1 && e_time[no_1][e_coll[no_1]]/1e6 < 0.37 && e_time[no_1][e_coll[no_1]]/1e6 > 0.337  ) cout << endl << "the time is: "<< e_time[no_1][e_coll[no_1]]/1e6 <<  " and this is e_fe_sil SPOT 3b = " << no_1_Fe_sil_isotope[1] <<  " and e_Fe_sil_isotope[no_1][1][e_coll[no_1]]:" << e_Fe_sil_isotope[no_1][1][e_coll[no_1]] << endl;

        if (j == 28)
        {
          // outfile13 << "Time" << "Mass" << "Sil54" << "Sil57" << "Met54" << "Met57" << " mantle_d57" << "core_d57" << "DELTA57";
          outfile13 << setw(8) << left << format(e_time[j][imp] / 1.e6, 6) << tab;
          outfile13 << setw(8) << left << format(e_mass[j][imp], 6) << tab;

          // print Fe54_sil, Fe57_sil
          outfile13 << setw(24) << left << format(e_Fe_sil_isotope[j][0][imp], 18) << tab;
          outfile13 << setw(24) << left << format(e_Fe_sil_isotope[j][1][imp], 18) << tab;

          // print Fe54_met, Fe57_met
          outfile13 << setw(24) << left << format(e_Fe_met_isotope[j][0][imp], 18) << tab;
          outfile13 << setw(24) << left << format(e_Fe_met_isotope[j][1][imp], 18) << tab;

          //mantle_d57
          double mantle_d57 = ((e_Fe_sil_isotope[j][1][imp] / e_Fe_sil_isotope[j][0][imp]) / (Fe57_IRMM / Fe54_IRMM) - 1.0) * 1000.0 ;
          outfile13 << setw(13) << left << format(mantle_d57, 6) << tab;

          double core_d57 = ((e_Fe_met_isotope[j][1][imp] / e_Fe_met_isotope[j][0][imp]) / (Fe57_IRMM / Fe54_IRMM) - 1.0) * 1000.0;
          outfile13 << setw(13) << left << format(core_d57, 6) << tab;

          double DELTA_57 = mantle_d57 - core_d57;
          outfile13 << setw(13) << left << format(DELTA_57, 6) << tab;
          outfile13 << endl;
        }
        if (j == 28)
        {
          // outfile14 << "Time" << "Mass" << "Sil54" << "Sil57" << "Met54" << "Met57" << " mantle_d57" << "core_d57" << "DELTA57";
          outfile14 << setw(8) << left << format(e_time[j][imp] / 1.e6, 6) << tab;
          outfile14 << setw(8) << left << format(e_mass[j][imp], 6) << tab;

          // print Fe54_sil, Fe57_sil
          outfile14 << setw(24) << left << format(e_Fe_sil_isotope[j][0][imp], 18) << tab;
          outfile14 << setw(24) << left << format(e_Fe_sil_isotope[j][1][imp], 18) << tab;

          // print Fe54_met, Fe57_met
          outfile14 << setw(24) << left << format(e_Fe_met_isotope[j][0][imp], 18) << tab;
          outfile14 << setw(24) << left << format(e_Fe_met_isotope[j][1][imp], 18) << tab;

          //mantle_d57
          double mantle_d57 = ((e_Fe_sil_isotope[j][1][imp] / e_Fe_sil_isotope[j][0][imp]) / (Fe57_IRMM / Fe54_IRMM) - 1.0) * 1000.0 ;
          outfile14 << setw(13) << left << format(mantle_d57, 6) << tab;

          double core_d57 = ((e_Fe_met_isotope[j][1][imp] / e_Fe_met_isotope[j][0][imp]) / (Fe57_IRMM / Fe54_IRMM) - 1.0) * 1000.0;
          outfile14 << setw(13) << left << format(core_d57, 6) << tab;

          double DELTA_57 = mantle_d57 - core_d57;
          outfile14 << setw(13) << left << format(DELTA_57, 6) << tab;
          outfile14 << endl;
        }

      }

      outfile1 << "____________________________________________________________________________________________________________________________________________________________________________________________________" << endl;
      outfile2 << "____________________________________________________________________________________________________________________________________________________________________________________________________" << endl;
      outfile3 << "____________________________________________________________________________________________________________________________________________________________________________________________________" << endl;
      outfile6 << "____________________________________________________________________________________________________________________________________________________________________________________________________" << endl;
      outfile7 << "____________________________________________________________________________________________________________________________________________________________________________________________________" << endl;
      outfile8 << "____________________________________________________________________________________________________________________________________________________________________________________________________" << endl;
      outfile9 << "____________________________________________________________________________________________________________________________________________________________________________________________________" << endl;
      outfile12 << "____________________________________________________________________________________________________________________________________________________________________________________________________" << endl;

      outfile1 << setw(16) << left << "PROP ERRORS:" << tab;
      for (int i = 0; i < n; i++)
      {
        if (sil_c[i] == 0) // output in wt%
          outfile1 << setw(8) << left << format(err[p_no[jj]][i], 6) << tab;

        else if (sil_c[i] == 1) // output in ppm
        {
          outfile1 << setw(8) << left << format(1.e4 * err[p_no[jj]][i], 6) << tab;
        }
        else if (sil_c[i] == 2) // output in ppb
        {
          outfile1 << setw(8) << left << format(1.e7 * err[p_no[jj]][i], 6) << tab;
        }
      }
      outfile1 << endl;

      if (fit_planet == 3)
      {
        outfile1 << endl << setw(16) << left << "MOON MANTLE:" << tab;
        for (int i = 0; i < n; i++)
        {
          if (i < 6) outfile1 << setw(8) << left << format(MoonM[i], 6) << tab;
          else outfile1 << setw(8) << left << format(MoonM[i], 6) << tab;
        }
        outfile1 << setw(8) << left << "14 (0.3)" << tab << setw(8) << left << "1000" << tab << setw(8) << left << "200" << endl;

      }
      else 
      {
        outfile1 << endl << setw(16) << left << "MOON MANTLE:" << tab;
        for (int i = 0; i < n; i++)
        {
          if (i < 6) outfile1 << setw(8) << left << format(EM[i], 6) << tab;
          else outfile1 << setw(8) << left << format(EM[i], 6) << tab;
        }
        outfile1 << setw(8) << left << "14 (0.3)" << tab << setw(8) << left << "1000" << tab << setw(8) << left << "200" << endl;

      }
      
      outfile1 << setprecision(4) << endl << "Start of LATE VENEER at " << t0_LV[p_no[jj]] / 1e6 << " My *** ";
      if (matte_HSE == 0) outfile1 << "No FeS segregation" << endl;
      else if (matte_HSE == 1) outfile1 << "Earliest time of FeS segregation = " << t_matte[p_no[jj]] / 1e6 << " My" << endl;

      //Determine final mass based on mass balance:
      double object_mass(0);
      for (int i = 0; i < (n + m); i++) moles[i] = e_moles[j][i][e_coll[j]];
      moles_H2O = e_moles_H2O[j][e_coll[j]];
      moles_oxygen = e_moles_oxygen[j][e_coll[j]];
      moles_H = e_moles_H[j][e_coll[j]];
      moles_S_mantle = e_moles_S_mantle[j][e_coll[j]];
      moles_S_core = e_moles_S_core[j][e_coll[j]];
      moles_C_mantle = e_moles_C_mantle[j][e_coll[j]];
      moles_C_core = e_moles_C_core[j][e_coll[j]];

      object_mass = mass_f(moles, moles_H2O, moles_oxygen, moles_H, moles_S_mantle, moles_S_core, moles_C_mantle, moles_C_core, AW, MW);

      outfile1 << "Final mass based on accumulated mass balance results = " << object_mass << endl;
      outfile1 << "Mass accreted from embryos = " << e_em_mass[j] << tab << "Mass accreted from planetesimals = " << e_pl_mass[j] << endl;
      outfile1 << "FRACTION OF MASS ACCRETED FROM WATER-BEARING COMPOSITIONS = " << 100.0 * e_SN_mass[j] / object_mass << " %" << endl;
      outfile2 << "Mass accreted from embryos = " << e_em_mass[j] << tab << "Mass accreted from planetesimals = " << e_pl_mass[j] << endl;
      // cout << "Mass beyond SN = " << e_SN_mass[j] << endl;
      outfile1 << endl;

      // IN THE CASE OF THE EARTH-LIKE PLANET, OUTPUT FINAL DETAILS AND CHI SQ:
      if (p_no[jj] == moon_body)
      {
        outfile10 << endl;
        outfile10 << "FINAL WT% COMPOSITIONS FOR PLANET #" << moon_body << " (with observed mantle values from Palme & O'Neill 2013):" << endl;
        outfile10 << "Components included in chi square calculation indicated by *" << endl;
        outfile10 << setw(8) << left << "MANTLE" << tab << setw(25) << left << "CALCULATED (UNCERTAINTY)" << tab << setw(8) << " " << tab << setw(25) << left << "OBSERVED (UNCERTAINTY)" << tab << setw(8) << " " << tab << setw(25) << left << "CALC % ERROR (OBS % ERROR)" << endl;

        string elem_name;
        for (int i = 0; i < n; i++)
        {
          elem_name.clear();
          for (int j = 0; j < 2; j++)
          {
            if ((j == 1) && ((el_sym[i][0] == 'V') || (el_sym[i][0] == 'W'))) continue;
            elem_name += el_sym[i][j];
          }
          if (el_sym[i][0] == 'A') elem_name += "2O3";
          else if ((el_sym[i][0] == 'S') && (el_sym[i][1] == 'i')) elem_name += "O2";
          else if ((i < (n - m + 2)) && (el_sym[i][0] != 'N')) elem_name += "O";
          if (fit_elem[i] == 'y') elem_name += "*";
          outfile10 << setw(8) << left << elem_name;

          outfile10 << tab << scientific << setw(11) << setprecision(5) << left << wtbest[moon_body][i] << " (" << err[moon_body][i] << ")" << tab << setw(8) << " " << tab << setw(11) << setprecision(5) << left << p_mantle[i] << " (" << sig_p_mantle[i]*p_mantle[i] << ")" << tab << setw(8) << " " << tab << setw(8) << format(100.0 * (wtbest[moon_body][i] - p_mantle[i]) / p_mantle[i], 6) << " (" << format(100.0 * sig_p_mantle[i], 6) << ")" << endl;
        }
        elem_name.clear();
        elem_name += el_sym_H2O;
        outfile10 << setw(8) << left << elem_name;
        outfile10 << tab << scientific << setw(11) << setprecision(5) << left << wtbest_H2O[moon_body] << setw(14) << " " << tab << setw(8) << " " << tab << setw(11) << setprecision(5) << left << 1.00000e-01 << " (" << 3.00000e-2 << ")" << tab << setw(8) << " " << tab << setw(8) << format(100.0 * (wtbest_H2O[moon_body] - 0.1) / 0.1, 6) << endl;
        elem_name.clear();
        elem_name += el_sym_S;
        elem_name += "*";
        outfile10 << setw(8) << left << elem_name;
        outfile10 << tab << scientific << setw(11) << setprecision(5) << left << wtbest_S_mantle[moon_body] << setw(14) << " " << tab << setw(8) << " " << tab << setw(11) << setprecision(5) << left << 2.00000e-02 << " (" << 5.00000e-3 << ")" << tab << setw(8) << " " << tab << setw(8) << format(100.0 * (wtbest_S_mantle[moon_body] - 0.02) / 0.02, 6) << endl;
        elem_name.clear();
        elem_name += el_sym_C;
        elem_name += "*";
        outfile10 << setw(8) << left << elem_name;
        outfile10 << tab << scientific << setw(11) << setprecision(5) << left << wtbest_C_mantle[moon_body] << setw(14) << " " << tab << setw(8) << " " << tab << setw(11) << setprecision(5) << left << 0.0075 << " (" << 0.0020 << ")" << tab << setw(8) << " " << tab << setw(8) << format(100.0 * (wtbest_C_mantle[moon_body] - 0.0075) / 0.0075, 6) << endl;
        outfile10 << setw(8) << left << "Chi Sq = " << tab << setw(8) << format(chi_sq, 6) << endl << endl;
        outfile2 << setw(8) << left << "Chi Sq = " << tab << setw(8) << format(chi_sq, 6) << endl << endl;

        outfile10 << "CORE" << endl;
        for (int i = (n - m); i < n; i++)
        {
          elem_name.clear();
          for (int j = 0; j < 2; j++)
          {
            if ((j == 1) && ((el_sym[i][0] == 'V') || (el_sym[i][0] == 'W'))) continue;
            elem_name += el_sym[i][j];
          }
          outfile10 << setw(8) << left << elem_name;


          outfile10 << tab << setprecision(5) << setw(11) << wtbest[moon_body][i + m] << " (" << err[moon_body][i] << ")" << endl;
        }
        outfile10 << setw(8) << left << "O" << tab << setw(11) << setprecision(5) << wtbest_O[moon_body] << endl;
        outfile10 << setw(8) << left << "H" << tab << setw(11) << setprecision(5) << wtbest_H[moon_body] << endl;
        outfile10 << setw(8) << left << "S" << tab << setw(11) << setprecision(5) << wtbest_S_core[moon_body] << setw(14) << " " << tab << setw(8) << " " << tab << setw(11) << setprecision(5) << left << 1.7 << " (" << 0.3 << ")" << tab << setw(8) << " " << tab << setw(8) << format(100.0 * (wtbest_S_core[moon_body] - 1.7) / 1.7, 6) << endl;
        outfile10 << setw(8) << left << "C" << tab << setw(11) << setprecision(5) << wtbest_C_core[moon_body] << setw(14) << " " << tab << setw(8) << " " << tab << setw(11) << setprecision(5) << left << 0.2 << " (" << 0.05 << ")" << tab << setw(8) << " " << tab << setw(8) << format(100.0 * (wtbest_C_core[moon_body] - 0.2) / 0.2, 6) << endl;
      }

      // commented out mars_planet related information in this noMarsNumber
      // IN THE CASE OF THE MARS-LIKE PLANET, OUTPUT FINAL DETAILS AND CHI SQ:
      // if (p_no[jj] == mars_planet)
      // {
      //   outfile11 << endl;
      //   outfile11 << "FINAL WT% COMPOSITIONS FOR PLANET #" << mars_planet << " (with observed mantle values from Taylor):" << endl;
      //   outfile11 << "Components included in chi square calculation indicated by *" << endl;
      //   outfile11 << setw(8) << left << "MANTLE" << tab << setw(25) << left << "CALCULATED (UNCERTAINTY)" << tab << setw(8) << " " << tab << setw(25) << left << "OBSERVED (UNCERTAINTY)" << tab << setw(8) << " " << tab << setw(25) << left << "CALC % ERROR (OBS % ERROR)" << endl;

      //   string elem_name;
      //   for (int i = 0; i < n; i++)
      //   {
      //     elem_name.clear();
      //     for (int j = 0; j < 2; j++)
      //     {
      //       if ((j == 1) && ((el_sym[i][0] == 'V') || (el_sym[i][0] == 'W'))) continue;
      //       elem_name += el_sym[i][j];
      //     }
      //     if (el_sym[i][0] == 'A') elem_name += "2O3";
      //     else if ((el_sym[i][0] == 'S') && (el_sym[i][1] == 'i')) elem_name += "O2";
      //     else if ((i < (n - m + 2)) && (el_sym[i][0] != 'N')) elem_name += "O";
      //     if (fit_elem[i] == 'y') elem_name += "*";
      //     outfile11 << setw(8) << left << elem_name;

      //     outfile11 << tab << scientific << setw(11) << setprecision(5) << left << wtbest[mars_planet][i] << " (" << err[mars_planet][i] << ")" << tab << setw(8) << " " << tab << setw(11) << setprecision(5) << left << m_mantle[i] << " (" << sig_m_mantle[i]*m_mantle[i] << ")" << tab << setw(8) << " " << tab << setw(8) << format(100.0 * (wtbest[mars_planet][i] - m_mantle[i]) / m_mantle[i], 6) << " (" << format(100.0 * sig_m_mantle[i], 6) << ")" << endl;
      //   }
      //   elem_name.clear();
      //   elem_name += el_sym_H2O;
      //   outfile11 << setw(8) << left << elem_name;
      //   outfile11 << tab << scientific << setw(11) << setprecision(5) << left << wtbest_H2O[mars_planet] << setw(14) << " " << tab << setw(8) << " " << tab << setw(11) << setprecision(5) << left << 1.00000e-01 << " (" << 3.00000e-2 << ")" << tab << setw(8) << " " << tab << setw(8) << format(100.0 * (wtbest_H2O[mars_planet] - 0.1) / 0.1, 6) << endl;
      //   // elem_name.clear();
      //   // elem_name += el_sym_S;
      //   // elem_name += "";
      //   // outfile11 << setw(8) << left << elem_name;
      //   // outfile11 << tab << scientific << setw(11) << setprecision(5) << left << wtbest_S_mantle[mars_planet] << setw(14) << " " << tab << setw(8) << " " << tab << setw(11) << setprecision(5) << left << 2.00000e-02 << " (" << 5.00000e-3 << ")" << tab << setw(8) << " " << tab << setw(8) << format(100.0 * (wtbest_S_mantle[mars_planet] - 0.02) / 0.02, 6) << endl;
      //   // elem_name.clear();
      //   // elem_name += el_sym_C;
      //   // elem_name += "";
      //   // outfile11 << setw(8) << left << elem_name;
      //   // outfile11 << tab << scientific << setw(11) << setprecision(5) << left << wtbest_C_mantle[mars_planet] << setw(14) << " " << tab << setw(8) << " " << tab << setw(11) << setprecision(5) << left << 0.0075 << " (" << 0.0020 << ")" << tab << setw(8) << " " << tab << setw(8) << format(100.0 * (wtbest_C_mantle[mars_planet] - 0.0075) / 0.0075, 6) << endl;
      //   outfile11 << setw(8) << left << "Chi Sq = " << tab << setw(8) << format(chi_sq, 6) << endl << endl;
      //   outfile2 << setw(8) << left << "Chi Sq = " << tab << setw(8) << format(chi_sq, 6) << endl << endl;

      //   outfile11 << "CORE" << endl;
      //   for (int i = (n - m); i < n; i++)
      //   {
      //     elem_name.clear();
      //     for (int j = 0; j < 2; j++)
      //     {
      //       if ((j == 1) && ((el_sym[i][0] == 'V') || (el_sym[i][0] == 'W'))) continue;
      //       elem_name += el_sym[i][j];
      //     }
      //     outfile11 << setw(8) << left << elem_name;


      //     outfile11 << tab << setprecision(5) << setw(11) << wtbest[mars_planet][i + m] << " (" << err[mars_planet][i] << ")" << endl;
      //   }
      //   outfile11 << setw(8) << left << "O" << tab << setw(11) << setprecision(5) << wtbest_O[mars_planet] << endl;
      //   outfile11 << setw(8) << left << "H" << tab << setw(11) << setprecision(5) << wtbest_H[mars_planet] << endl;
      //   //outfile11 << setw(8) << left << "S" << tab << setw(11) << setprecision(5) << wtbest_S_core[mars_planet] << setw(14) << " " << tab << setw(8) << " " << tab << setw(11) << setprecision(5) << left << 1.7 << " (" << 0.3 << ")" << tab << setw(8) << " " << tab << setw(8) << format(100.0 * (wtbest_S_core[mars_planet] - 1.7) / 1.7, 6) << endl;
      //   //outfile11 << setw(8) << left << "C" << tab << setw(11) << setprecision(5) << wtbest_C_core[mars_planet] << setw(14) << " " << tab << setw(8) << " " << tab << setw(11) << setprecision(5) << left << 0.2 << " (" << 0.05 << ")" << tab << setw(8) << " " << tab << setw(8) << format(100.0 * (wtbest_C_core[mars_planet] - 0.2) / 0.2, 6) << endl;
      // }
    }

    // OUTPUT EMBRYO COMPOSITION DETAILS TO MANTLE AND CORE FILES
    outfile4 << endl;
    outfile5 << endl;
    outfile8 << endl;
    outfile9 << endl;

    // cout << "em_max = " << em_max << endl; // EMBRYO MANTLE DEBUGGING COMMENT
    // cout << "write file value" << write_file << endl; // EMBRYO MANTLE DEBUGGING COMMENT
    // THIS IS WHERE EMBRYO_MANTLE OUTFILE IS WRITTEN
    for (jj = 0; jj < em_max; jj++)  // Runs through each successive embryo
    {
      if (em_merge_time[jj] != 0.0) // Excludes the final planets (for which zero is the default value)
      {
        outfile4 << endl << "EMBRYO #" << em_no[jj] << endl;
        outfile4 << "Merges with #" << em_target[jj] << " at " << em_merge_time[jj] << "years" << endl;
        outfile5 << endl << "EMBRYO #" << em_no[jj] << endl;
        outfile5 << "Merges with #" << em_target[jj] << " at " << em_merge_time[jj] << "years" << endl;
        outfile8 << endl << "EMBRYO #" << em_no[jj] << endl;
        outfile8 << "Merges with #" << em_target[jj] << " at " << em_merge_time[jj] << "years" << endl;
        outfile9 << endl << "EMBRYO #" << em_no[jj] << endl;
        outfile9 << "Merges with #" << em_target[jj] << " at " << em_merge_time[jj] << "years" << endl;
        write_headings_em();

        int imp(0);
        j = em_no[jj]; // Corresponds to the surviving embryo numbers
        // cout << "Surviving Embryo Numbers: " << j << endl; // EMBRYO MANTLE DEBUGGING COMMENT

        double total_moles[n + m] = {0.0};
        double total_moles_H2O = 0.0;
        double total_moles_oxygen = 0.0;
        double total_moles_H = 0.0;
        double total_moles_S_mantle = 0.0;
        double total_moles_S_core = 0.0;
        double total_moles_C_mantle = 0.0;
        double total_moles_C_core = 0.0;
        for (imp = 0; imp < (e_coll[j] + 1); imp++)
        {
          // cout << "e_coll[j]: " << e_coll[j] << endl; // EMBRYO MANTLE DEBUGGING COMMENT
          for (int i = 0; i < (n + m); i++) total_moles[i] += s_moles[j][i][imp];
          total_moles_H2O += s_moles_H2O[j][imp];
          total_moles_oxygen += s_moles_oxygen[j][imp];
          total_moles_H += s_moles_H[j][imp];
          total_moles_S_mantle += s_moles_S_mantle[j][imp];
          total_moles_S_core += s_moles_S_core[j][imp];
          total_moles_C_mantle += s_moles_C_mantle[j][imp];
          total_moles_C_core += s_moles_C_core[j][imp];
        }

        for (imp = 0; imp < (e_coll[j]); imp++) // RUNS THROUGH EACH SUCCESSIVE COLLISION FOR EMBRYO j, and  outputs composition
        {
          for (int i = 0; i < (n + m); i++) moles[i] = e_moles[j][i][imp];
          moles_H2O = e_moles_H2O[j][imp];
          moles_oxygen = e_moles_oxygen[j][imp];
          moles_H = e_moles_H[j][imp];
          moles_S_mantle = e_moles_S_mantle[j][imp];
          moles_S_core = e_moles_S_core[j][imp];
          moles_C_mantle = e_moles_C_mantle[j][imp];
          moles_C_core = e_moles_C_core[j][imp];
          // EMBRYO MANTLE DEBUGGING COMMENT
          outfile4 << setw(8) << left << format(e_time[j][imp] / 1.e6, 6) << tab;
          outfile4 << setw(8) << left << format(e_mass[j][imp], 6) << tab;

          outfile5 << setw(8) << left << format(e_time[j][imp] / 1.e6, 6) << tab;
          outfile5 << setw(8) << left << format(e_mass[j][imp], 6) << tab;

          outfile5 << setw(8) << left << format(e_P[j][imp], 6) << tab;
          outfile5 << setw(8) << left << format(e_T[j][imp], 6) << tab;
          outfile5 << setw(8) << left << format(e_fO2[j][imp], 6) << tab;

          if (imp == e_coll[j] - 1) norm_conc = 1; // Outputs final errors and moon mantle concentrations at end for each planet
          else norm_conc = 0;

          weight_comp_em(moles, moles_oxygen, moles_H2O, moles_H, moles_S_mantle, moles_S_core, moles_C_mantle, moles_C_core, AW, MW);

          if (norm_conc == 1)
          {
            outfile4 << "_______________________________________________________________________________________________________________________________________________________________________" << endl;
            outfile4 << setw(16) << "PROP ERRORS:" << tab;

            for (int i = 0; i < n; i++)
            {
              if (sil_c[i] == 0) // output in wt%
                outfile4 << setw(8) << left << format(err[j][i], 6) << tab;

              else if (sil_c[i] == 1) // output in ppm
              {
                outfile4 << setw(8) << left << format(1e4 * err[j][i], 6) << tab;
              }
              else if (sil_c[i] == 2) // output in ppb
              {
                outfile4 << setw(8) << left << format(1e7 * err[j][i], 6) << tab;
              }
            }

            outfile4 << setprecision(4) << endl << "Start of LATE VENEER at " << t0_LV[p_no[jj]] / 1e6 << " My *** ";
            if (matte_HSE == 0) outfile4 << "No FeS segregation" << endl;
            else if (matte_HSE == 1) outfile4 << "Earliest time of FeS segregation = " << t_matte[p_no[jj]] / 1e6 << " My" << endl;

          }

          outfile8 << setw(8) << left << format(e_time[j][imp] / 1.e6, 6) << tab;
          outfile8 << setw(8) << left << format(e_mass[j][imp], 6) << tab;

          outfile9 << setw(8) << left << format(e_time[j][imp] / 1.e6, 6) << tab;
          outfile9 << setw(8) << left << format(e_mass[j][imp], 6) << tab;

          for (int i = 0; i < (n + m); i++) moles[i] = s_moles[j][i][imp];
          moles_H2O = s_moles_H2O[j][imp];
          moles_oxygen = s_moles_oxygen[j][imp];
          moles_H = s_moles_H[j][imp];
          moles_S_mantle = s_moles_S_mantle[j][imp];
          moles_S_core = s_moles_S_core[j][imp];
          moles_C_mantle = s_moles_C_mantle[j][imp];
          moles_C_core = s_moles_C_core[j][imp];

          moles_comp_source(moles, moles_oxygen, moles_H2O, moles_H, moles_S_mantle, moles_S_core, moles_C_mantle, moles_C_core,
                            total_moles, total_moles_oxygen, total_moles_H2O, total_moles_H, total_moles_S_mantle, total_moles_S_core, total_moles_C_mantle, total_moles_C_core);


        }

        outfile8 << "_______________________________________________________________________________________________________________________________________________________________________" << endl;
        outfile9 << "_______________________________________________________________________________________________________________________________________________________________________" << endl;

        //Determine final mass based on mass balance:
        double object_mass(0);
        for (int i = 0; i < (n + m); i++)
          moles[i] = e_moles[j][i][e_coll[j] - 1];
        moles_H2O = e_moles_H2O[j][e_coll[j] - 1];
        moles_oxygen = e_moles_oxygen[j][e_coll[j] - 1];
        moles_H = e_moles_H[j][e_coll[j] - 1];
        moles_S_mantle = e_moles_S_mantle[j][e_coll[j] - 1];
        moles_S_core = e_moles_S_core[j][e_coll[j] - 1];
        moles_C_mantle = e_moles_C_mantle[j][e_coll[j] - 1];
        moles_C_core = e_moles_C_core[j][e_coll[j] - 1];

        object_mass = mass_f(moles, moles_H2O, moles_oxygen, moles_H, moles_S_mantle, moles_S_core, moles_C_mantle, moles_C_core, AW, MW);

        outfile4 << "Final mass based on accumulated mass balance results = " << object_mass << endl;
        outfile4 << "Mass accreted from embryos = " << e_em_mass[j] << tab << "Mass accreted from planetesimals = " << e_pl_mass[j] << endl;
        outfile4 << "FRACTION OF MASS ACCRETED FROM WATER-BEARING COMPOSITIONS = " << 100.0 * e_SN_mass[j] / object_mass << " %" << endl;
        outfile5 << "Mass accreted from embryos = " << e_em_mass[j] << tab << "Mass accreted from planetesimals = " << e_pl_mass[j] << endl;
        // cout << "Mass beyond SN = " << e_SN_mass[j] << endl;
        outfile4 << endl;
      }
    }
    cout << "N-BODY SIMULATION No: " << simulation_no << endl;
  }                                         // END OF FINAL OUTPUT (when write_file==1)

  return chi_sq;
}                       // END OF func

//____________________________________________________________________________________
// WEIGHT_COMP FUNCTION that converts moles of elements in silicate and metal to wt% values for final planets
// and writes results to output files
void weight_comp(double moles[], double moles_O, double moles_H2O, double moles_H, double moles_S_mantle, double moles_S_core, double moles_C_mantle, double moles_C_core, const double AW[], const double MW[])
{
  double weight[n + m + 1] = {0}; // All elements including oxygen
  double weight_H2O = 0.0;
  double weight_H = 0.0;
  double weight_S_mantle = 0.0;
  double weight_S_core = 0.0;
  double weight_C_mantle = 0.0;
  double weight_C_core = 0.0;

  // Calculate silicate components
  double sum_sil = 0;
  for (int i = 0; i < n; i++)
  {
    weight[i] = moles[i] * MW[i];
    sum_sil += weight[i];
  }
  weight_H2O = moles_H2O * MW_H2O;
  sum_sil += weight_H2O;
  weight_S_mantle = moles_S_mantle * AW_S;
  sum_sil += weight_S_mantle;
  weight_C_mantle = moles_C_mantle * AW_C;
  sum_sil += weight_C_mantle;
  for (int i = 0; i < n; i++) weight[i] = weight[i] * 100.0 / sum_sil;
  weight_H2O = weight_H2O * 100.0 / sum_sil;
  weight_S_mantle = weight_S_mantle * 100.0 / sum_sil;
  weight_C_mantle = weight_C_mantle * 100.0 / sum_sil;

  // Calculate metal components
  double sum_met = 0;
  for (int i = n; i < (n + m); i++)
  {
    weight[i] = moles[i] * AW[i - m];
    sum_met += weight[i];
  }

  weight[n + m] = moles_O * AW_O;     // Adds oxygen
  sum_met += weight[n + m];
  weight_H = moles_H * AW_H;    // Adds H
  sum_met += weight_H;
  weight_S_core = moles_S_core * AW_S; // Adds S
  sum_met += weight_S_core;
  weight_C_core = moles_C_core * AW_C; // Adds S
  sum_met += weight_C_core;

  if (sum_met > 0.)
  {
    for (int i = n; i < n + m + 1; i++)
      weight[i] = weight[i] * 100.0 / sum_met;
    weight_H = weight_H * 100.0 / sum_met; //
    weight_S_core = weight_S_core * 100.0 / sum_met;
    weight_C_core = weight_C_core * 100.0 / sum_met;
  }
  else
  {
    for (int i = n; i < n + m + 1; i++)
      weight[i] = 0.0;
    weight_H = 0.0; //
    weight_S_core = 0.0;
    weight_C_core = 0.0;
  }

  if (write_file == 1)
  {
    for (int i = 0; i < n; i++)
    {
      if (sil_c[i] == 0)   // output in wt%
      {
        if ((el_sym[i][0] == 'N') && (el_sym[i][1] == 'a')) // Elemental concentration of sodium
          outfile1 << setw(8) << left << format(weight[i]*AW[i] / MW[i], 6) << tab;
        else outfile1 << setw(8) << left << format(weight[i], 6) << tab;
      }
      else if (sil_c[i] == 1) // output in ppm
      {
        outfile1 << setw(8) << left << format(1.e4 * weight[i]*AW[i] / MW[i], 6) << tab;
      }
      else if (sil_c[i] == 2) // output in ppb
      {
        outfile1 << setw(8) << left << format(1e7 * weight[i]*AW[i] / MW[i], 6) << tab;
      }
    }

    double nb_ta(0);
    for (int i = 0; i < n; i++)                     // determines Nb/Ta ratio
    {
      if ((el_sym[i][0] == 'N') && (el_sym[i][1] == 'b'))
        nb_ta = (weight[i] * AW[i] / MW[i]) / (weight[i + 1] * AW[i + 1] / MW[i + 1]);
    }
    outfile1 << setw(8) << left << format(nb_ta, 6) << tab;
    outfile1 << setw(8) << left << format(1.e4 * weight_H2O, 6) << tab;
    outfile1 << setw(8) << left << format(1.e4 * weight_S_mantle, 6) << tab;
    outfile1 << setw(8) << left << format(1.e4 * weight_C_mantle, 6) << endl;

    //Output METAL composition to FILE

    for (int i = n; i < n + m + 1; i++)
    {
      if ((met_c[i - m] == 0) || (i == n + m))     // output in wt%
        outfile2 << setw(8) << left << format(weight[i], 6) << tab;
      else if (met_c[i - m] == 1)               // output in ppm
      {
        outfile2 << setw(8) << left << format(1.e4 * weight[i], 6) << tab;
      }
      else if (met_c[i - m] == 2)               // output in ppb
      {
        outfile2 << setw(8) << left << format(1.e7 * weight[i], 6) << tab;
      }
    }
    outfile2 << setw(8) << left << format(1.e4 * weight_H, 6) << tab;
    outfile2 << setw(8) << left << format(weight_S_core, 6) << tab;
    outfile2 << setw(8) << left << format(weight_C_core, 6) << tab;
    outfile2 << setw(8) << left << format(sum_met * 100.0 / (sum_met + sum_sil), 6) << endl; // Mass fraction of core
  }

  return;
}                 // END of weight_comp

// WEIGHT_COMP_EM FUNCTION that converts moles of elements in silicate and metal to wt% values for embryos
// and writes results to output files
void weight_comp_em(double moles[], double moles_O, double moles_H2O, double moles_H, double moles_S_mantle, double moles_S_core, double moles_C_mantle, double moles_C_core, const double AW[], const double MW[])
{
  double weight[n + m + 1] = {0};
  double weight_H2O = 0.0;
  double weight_H = 0.0;
  double weight_S_mantle = 0.0;
  double weight_S_core = 0.0;
  double weight_C_mantle = 0.0;
  double weight_C_core = 0.0;

  // Calculate silicate components
  double sum_sil = 0;
  for (int i = 0; i < n; i++)
  {
    weight[i] = moles[i] * MW[i];
    // cout << "moles " << i << " " << moles[i] <<endl; // EMBRYO MANTLE DEBUGGING COMMENT
    // cout << "MW: " << i << " " << MW[i] << endl; // EMBRYO MANTLE DEBUGGING COMMENT
    sum_sil += weight[i];
  }
  weight_H2O = moles_H2O * MW_H2O;
  sum_sil += weight_H2O;
  weight_S_mantle = moles_S_mantle * AW_S;
  sum_sil += weight_S_mantle;
  weight_C_mantle = moles_C_mantle * AW_C;
  sum_sil += weight_C_mantle;
  for (int i = 0; i < n; i++) weight[i] = weight[i] * 100.0 / sum_sil;
  weight_H2O = weight_H2O * 100.0 / sum_sil;
  weight_S_mantle = weight_S_mantle * 100.0 / sum_sil;
  weight_C_mantle = weight_C_mantle * 100.0 / sum_sil;

  // Calculate metal components
  double sum_met = 0;
  for (int i = n; i < (n + m); i++)
  {
    weight[i] = moles[i] * AW[i - m];
    sum_met += weight[i];
  }

  weight[n + m] = moles_O * AW_O;     // Adds oxygen
  sum_met += weight[n + m];
  weight_H = moles_H * AW_H;    // Adds H
  sum_met += weight_H;
  weight_S_core = moles_S_core * AW_S; // Adds S
  sum_met += weight_S_core;
  weight_C_core = moles_C_core * AW_C; // Adds C
  sum_met += weight_C_core;

  if (sum_met > 0.)
  {
    for (int i = n; i < n + m + 1; i++)
      weight[i] = weight[i] * 100.0 / sum_met;
    weight_H = weight_H * 100.0 / sum_met; //
    weight_S_core = weight_S_core * 100.0 / sum_met;
    weight_C_core = weight_C_core * 100.0 / sum_met;
  }
  else
  {
    for (int i = n; i < n + m + 1; i++)
      weight[i] = 0.0;
    weight_H = 0.0; //
    weight_S_core = 0.0;
    weight_C_core = 0.0;
  }

  if (write_file == 1)
  {
    //Output SILICATE composition to FILE

    for (int i = 0; i < n; i++)
    {
      if (sil_c[i] == 0)   // output in wt%
      {
        if ((el_sym[i][0] == 'N') && (el_sym[i][1] == 'a')) outfile4 << setw(8) << left << format(weight[i]*AW[i] / MW[i], 6) << tab;
        else outfile4 << setw(8) << left << format(weight[i], 6) << tab;
        // cout << "outfile4 weight[i]: " << weight[i] << endl; // EMBRYO MANTLE DEBUGGING COMMENT
      }
      else if (sil_c[i] == 1) // output in ppm
      {
        outfile4 << setw(8) << left << format(1.e4 * weight[i]*AW[i] / MW[i], 6) << tab;
      }
      else if (sil_c[i] == 2) // output in ppb
      {
        outfile4 << setw(8) << left << format(1e7 * weight[i]*AW[i] / MW[i], 6) << tab;
      }
    }

    double nb_ta(0);
    for (int i = 0; i < n; i++)                     // determines Nb/Ta ratio
    {
      if ((el_sym[i][0] == 'N') && (el_sym[i][1] == 'b'))
        nb_ta = (weight[i] * AW[i] / MW[i]) / (weight[i + 1] * AW[i + 1] / MW[i + 1]);
    }
    outfile4 << setw(8) << left << format(nb_ta, 6) << tab;
    outfile4 << setw(8) << left << format(1.e4 * weight_H2O, 6) << tab;
    outfile4 << setw(8) << left << format(1.e4 * weight_S_mantle, 6) << tab;
    outfile4 << setw(8) << left << format(1.e4 * weight_C_mantle, 6) << endl;

    //Output METAL composition to FILE

    for (int i = n; i < n + m + 1; i++)
    {
      if ((met_c[i - m] == 0) || (i == n + m))     // output in wt%
        outfile5 << setw(8) << left << format(weight[i], 6) << tab;
      else if (met_c[i - m] == 1)               // output in ppm
      {
        outfile5 << setw(8) << left << format(1.e4 * weight[i], 6) << tab;
      }
      else if (met_c[i - m] == 2)               // output in ppb
      {
        outfile5 << setw(8) << left << format(1.e7 * weight[i], 6) << tab;
      }
    }
    outfile5 << setw(8) << left << format(1.e4 * weight_H, 6) << tab;
    outfile5 << setw(8) << left << format(weight_S_core, 6) << tab;
    outfile5 << setw(8) << left << format(weight_C_core, 6) << tab;
    outfile5 << setw(8) << left << format(sum_met * 100.0 / (sum_met + sum_sil), 6) << endl; // Mass fraction of core
  }

  return;
}                 // END of weight_comp_em
//_________________________________________________________________________________


// WEIGHT_COMP_IMP FUNCTION that converts moles of elements in silicate and metal to wt% values for each impactor on the final planets
// and writes results to output files
void weight_comp_imp(double moles[], double moles_O, double moles_H2O, double moles_H, double moles_S_mantle, double moles_S_core, double moles_C_mantle, double moles_C_core, const double AW[], const double MW[])
{
  double weight[n + m + 1] = {0}; // All elements including oxygen
  double weight_H2O = 0.0;
  double weight_H = 0.0;
  double weight_S_mantle = 0.0;
  double weight_S_core = 0.0;
  double weight_C_mantle = 0.0;
  double weight_C_core = 0.0;

  // Calculate silicate components
  double sum_sil = 0;
  for (int i = 0; i < n; i++)
  {
    weight[i] = moles[i] * MW[i];
    sum_sil += weight[i];
  }
  weight_H2O = moles_H2O * MW_H2O;
  sum_sil += weight_H2O;
  weight_S_mantle = moles_S_mantle * AW_S;
  sum_sil += weight_S_mantle;
  weight_C_mantle = moles_C_mantle * AW_C;
  sum_sil += weight_C_mantle;
  for (int i = 0; i < n; i++) weight[i] = weight[i] * 100.0 / sum_sil;
  weight_H2O = weight_H2O * 100.0 / sum_sil;
  weight_S_mantle = weight_S_mantle * 100.0 / sum_sil;
  weight_C_mantle = weight_C_mantle * 100.0 / sum_sil;

  // Calculate metal components
  double sum_met = 0;
  for (int i = n; i < (n + m); i++)
  {
    weight[i] = moles[i] * AW[i - m];
    sum_met += weight[i];
  }

  weight[n + m] = moles_O * AW_O;     // Adds oxygen
  sum_met += weight[n + m];
  weight_H = moles_H * AW_H;    // Adds H
  sum_met += weight_H;
  weight_S_core = moles_S_core * AW_S; // Adds S
  sum_met += weight_S_core;
  weight_C_core = moles_C_core * AW_C; // Adds S
  sum_met += weight_C_core;

  if (sum_met > 0.)
  {
    for (int i = n; i < n + m + 1; i++) weight[i] = weight[i] * 100.0 / sum_met;
    weight_H = weight_H * 100.0 / sum_met; //
    weight_S_core = weight_S_core * 100.0 / sum_met;
    weight_C_core = weight_C_core * 100.0 / sum_met;
  }
  else
  {
    for (int i = n; i < n + m + 1; i++) weight[i] = 0.0;
    weight_H = 0.0; //
    weight_S_core = 0.0;
    weight_C_core = 0.0;
  }

  if (write_file == 1)
  {
    for (int i = 0; i < n; i++)
    {
      if (sil_c[i] == 0)   // output in wt%
      {
        if ((el_sym[i][0] == 'N') && (el_sym[i][1] == 'a')) // Elemental concentration of sodium
          outfile6 << setw(8) << left << format(weight[i]*AW[i] / MW[i], 6) << tab;
        else outfile6 << setw(8) << left << format(weight[i], 6) << tab;
      }
      else if (sil_c[i] == 1) // output in ppm
      {
        outfile6 << setw(8) << left << format(1.e4 * weight[i]*AW[i] / MW[i], 6) << tab;
      }
      else if (sil_c[i] == 2) // output in ppb
      {
        outfile6 << setw(8) << left << format(1e7 * weight[i]*AW[i] / MW[i], 6) << tab;
      }
    }

    double nb_ta(0);
    for (int i = 0; i < n; i++)                     // determines Nb/Ta ratio
    {
      if ((el_sym[i][0] == 'N') && (el_sym[i][1] == 'b'))
        nb_ta = (weight[i] * AW[i] / MW[i]) / (weight[i + 1] * AW[i + 1] / MW[i + 1]);
    }
    outfile6 << setw(8) << left << format(nb_ta, 6) << tab;
    outfile6 << setw(8) << left << format(1.e4 * weight_H2O, 6) << tab;
    outfile6 << setw(8) << left << format(1.e4 * weight_S_mantle, 6) << tab;
    outfile6 << setw(8) << left << format(1.e4 * weight_C_mantle, 6) << endl;

    //Output METAL composition to FILE

    for (int i = n; i < n + m + 1; i++)
    {
      if ((met_c[i - m] == 0) || (i == n + m))     // output in wt%
        outfile7 << setw(8) << left << format(weight[i], 6) << tab;
      else if (met_c[i - m] == 1)               // output in ppm
      {
        outfile7 << setw(8) << left << format(1.e4 * weight[i], 6) << tab;
      }
      else if (met_c[i - m] == 2)               // output in ppb
      {
        outfile7 << setw(8) << left << format(1.e7 * weight[i], 6) << tab;
      }
    }
    outfile7 << setw(8) << left << format(1.e4 * weight_H, 6) << tab;
    outfile7 << setw(8) << left << format(weight_S_core, 6) << tab;
    outfile7 << setw(8) << left << format(weight_C_core, 6) << tab;
    outfile7 << setw(8) << left << format(sum_met * 100.0 / (sum_met + sum_sil), 6) << endl; // Mass fraction of core
  }

  return;
}                 // END of weight_comp_imp

// MOLES_COMP_SOURCE FUNCTION that converts moles of elements in silicate and metal to wt% values for each impactor on the final planets
// and writes results to output files
void moles_comp_source(double moles[], double moles_O, double moles_H2O, double moles_H, double moles_S_mantle, double moles_S_core, double moles_C_mantle, double moles_C_core,
                       double total_moles[], double total_moles_O, double total_moles_H2O, double total_moles_H, double total_moles_S_mantle, double total_moles_S_core, double total_moles_C_mantle, double total_moles_C_core)
{
  for (int i = 0; i < n + m; i++) if (total_moles[i] == 0.0) total_moles[i] = 1.;
  if (total_moles_H2O == 0.0) total_moles_H2O = 1.;
  if (total_moles_S_mantle == 0.0) total_moles_S_mantle = 1.;
  if (total_moles_C_mantle == 0.0) total_moles_C_mantle = 1.;
  if (total_moles_O == 0.0) total_moles_O = 1.;
  if (total_moles_H == 0.0) total_moles_H = 1.;
  if (total_moles_S_core == 0.0) total_moles_S_core = 1.;
  if (total_moles_C_core == 0.0) total_moles_C_core = 1.;

  if (write_file == 1)
  {
    for (int i = 0; i < n; i++) outfile8 << setw(8) << left << format(100.0 * moles[i] / total_moles[i], 6) << tab;
    outfile8 << setw(8) << left << format(100.0 * moles_H2O / total_moles_H2O, 6) << tab;
    outfile8 << setw(8) << left << format(100.0 * moles_S_mantle / total_moles_S_mantle, 6) << tab;
    outfile8 << setw(8) << left << format(100.0 * moles_C_mantle / total_moles_C_mantle, 6) << endl;

    for (int i = n; i < n + m; i++) outfile9 << setw(8) << left << format(100.0 * moles[i] / total_moles[i], 6) << tab;
    outfile9 << setw(8) << left << format(100.0 * moles_O / total_moles_O, 6) << tab;
    outfile9 << setw(8) << left << format(100.0 * moles_H / total_moles_H, 6) << tab;
    outfile9 << setw(8) << left << format(100.0 * moles_S_core / total_moles_S_core, 6) << tab;
    outfile9 << setw(8) << left << format(100.0 * moles_C_core / total_moles_C_core, 6) << endl;
  }

  return;
}                 // END of moles_comp_source

// FUNCTION incr_comp. This outputs composition of each successive batch of metal that is added to cores
void incr_comp(double moles[], double moles_oxygen, double moles_H, double moles_S_core, double moles_C_core, double sum_core, const double AW[])
{
  double weight[n + m + 1] = {0.0};
  double weight_H = 0.0;
  double weight_S = 0.0;
  double weight_C = 0.0;

  double sum_met = 0.0;
  for (int i = n; i < (n + m); i++)
  {
    weight[i] = moles[i] * AW[i - m];
    sum_met += weight[i];
  }

  weight[n + m] = moles_oxygen * AW_O; // Adds oxygen
  sum_met += weight[n + m];
  weight_H = moles_H * AW_H;    // Adds hydrogen
  sum_met += weight_H;
  weight_S = moles_S_core * AW_S;     // Adds S
  sum_met += weight_S;
  weight_C = moles_C_core * AW_C;     // Adds C
  sum_met += weight_C;

  if (sum_met > 0.0)
  {
    for (int i = n; i < n + m + 1; i++) weight[i] = weight[i] * 100.0 / sum_met;
    weight_H = weight_H * 100.0 / sum_met; //
    weight_S = weight_S * 100.0 / sum_met;
    weight_C = weight_C * 100.0 / sum_met;
  }
  else
  {
    for (int i = n; i < n + m + 1; i++) weight[i] = 0.0;
    weight_H = 0.0; //
    weight_S = 0.0;
    weight_C = 0.0;
  }

  for (int i = n; i < n + m + 1; i++)
  {
    if ((i > (n + 3)) && (i < (m + n)))          // Output wt%, ppm or ppb values for elements from Co to O
    {
      if ((el_sym[i - m][0] == 'C') && (el_sym[i - m][1] == 'r'))
        outfile3 << setw(8) << left << format(weight[i], 6) << tab; // wt% for Cr

      else if ((el_sym[i - m][0] == 'V') || ((el_sym[i - m][0] == 'M') && (el_sym[i - m][1] == 'n'))
               || ((el_sym[i - m][0] == 'G') && (el_sym[i - m][1] == 'a')) || ((el_sym[i - m][0] == 'Z') && (el_sym[i - m][1] == 'n'))
               || ((el_sym[i - m][0] == 'C') && (el_sym[i - m][1] == 'u')))
      {
        outfile3 << setw(8) << left << format(1.e4 * weight[i], 6) << tab;
      }
      else
      {
        outfile3 << setw(8) << left << format(1.e7 * weight[i], 6) << tab;
      }
    }
    else
      outfile3 << setw(8) << left << format(weight[i], 6) << tab; // Output wt% values elements
  }

  outfile3 << setw(8) << left << format(1.e4 * weight_H, 6) << tab;

  outfile3 << setw(8) << left << format(weight_S, 6) << tab;
  outfile3 << setw(8) << left << format(weight_C, 6) << tab;
  outfile3 << setw(8) << left << format(sum_met * 100.0 / sum_core, 6) << endl; // Mass fraction of core

  return;
}                 // END of incr_comp


//_________________________________________________________________________________

void write_headings()
{
  if (tm_col == 1)
  {
    outfile1 << setw(8) << left << "Time" << tab << setw(8) << left << "Mass" << tab;
    outfile2 << setw(8) << left << "Time" << tab << setw(8) << left << "Mass" << tab << setw(8) << left << "P" << tab << setw(8) << left << "T" << tab << setw(8) << left << "logfO2" << tab;
    outfile3 << setw(8) << left << "Time" << tab << setw(8) << left << "Mass" << tab << setw(8) << left << "P" << tab << setw(8) << left << "T" << tab << setw(8) << left << "logfO2" << tab;
    outfile6 << setw(8) << left << "Time" << tab << setw(8) << left << "Mass" << tab << setw(8) << left << "Xsil" << tab;
    outfile7 << setw(8) << left << "Time" << tab << setw(8) << left << "Mass" << tab << setw(8) << left << "Xmet" << tab;
    outfile8 << setw(8) << left << "Time" << tab << setw(8) << left << "Mass" << tab;
    outfile9 << setw(8) << left << "Time" << tab << setw(8) << left << "Mass" << tab;
    // outfile12 << setw(8) << left << "Time" << tab << setw(8) << left << "Mass" << tab << setw(13) << left << "Mantle Fe54" << tab << setw(9) << tab << left << "Mantle Fe57" << tab << setw(13) << tab << left << "Core Fe54" << tab << setw(11) << tab << left << "  Core Fe57" << tab << setw(13) << tab << left << "Mantle d57/54" << tab << setw(13) << left << "Core d57/54" << tab << setw(13) << left << "DELTA57 Mantle-Core" << endl;
    // outfile13 << setw(8) << left << "Time" << tab << setw(8) << left << "Mass" << tab << endl;
    // outfile12 << setw(8) << left  << "time" << tab << setw(9) << "mass" << tab << setw(9) << "Fe54_sil" << tab << setw(9) << "Fe56_sil" << tab << setw(9) << "Fe57_sil" << tab << setw(9) << "Fe54_met" << tab << setw(9) << "Fe56_met" << tab << setw(9) << "Fe57_met" << tab << setw(9) << "mantle_d56" << tab << setw(9) << "core_d56" << tab << setw(9) << "DELTA_56" << tab << setw(9) << "mantle_d57" << tab << setw(9) << "core_d57" << tab << setw(9) << "DELTA_57" << endl;
    outfile12 << setw(8) << left  << "time" << tab << setw(9) << "mass" << tab << setw(24) << "Fe54_sil" << tab << setw(24) << "Fe56_sil" << tab << setw(24) << "Fe57_sil" << tab << setw(24) << "Fe54_met" << tab << setw(24) << "Fe56_met" << tab << setw(24) << "Fe57_met" << tab << setw(13) << "mantle_d56" << tab << setw(13) << "core_d56" << tab << setw(13) << "DELTA_56" << tab << setw(13) << "mantle_d57" << tab << setw(13) << "core_d57" << tab << setw(13) << "DELTA_57"   << endl; 

  }

  string elem_name;
  // Mantle composition.
  for (int i = 0; i < n; i++)
  {
    elem_name.clear();
    for (int j = 0; j < 2; j++)
    {
      if ((j == 1) && ((el_sym[i][0] == 'V') || (el_sym[i][0] == 'W'))) continue;
      elem_name += el_sym[i][j];
    }

    if ((el_sym[i][0] == 'A') && (el_sym[i][1] == 'l')) elem_name += "2O3";
    else if ((el_sym[i][0] == 'S') && (el_sym[i][1] == 'i')) elem_name += "O2";
    else if ((i < (n - m + 2)) && (el_sym[i][0] != 'N')) elem_name += "O"; // Na is excluded as oxide
    outfile1 << setw(8) << left << elem_name << tab;
  }
  outfile1 << setw(8) << left << "Nb/Ta" << tab;
  outfile1 << setw(8) << left << "H2O" << tab;
  outfile1 << setw(8) << left << el_sym_S << tab;
  outfile1 << setw(8) << left << el_sym_C << endl;

  if (tm_col == 1)
    outfile1 << setw(8) << left << "Ma" << tab << setw(8) << left << "Me" << tab;
  for (int i = 0; i < n; i++)
    if (sil_c[i] == 0) outfile1 << setw(8) << left << "wt%" << tab;
    else if (sil_c[i] == 1) outfile1 << setw(8) << left << "ppm" << tab;
    else if (sil_c[i] == 2) outfile1 << setw(8) << left << "ppb" << tab;
  outfile1 << setw(8) << left << " " << tab << setw(8) << left << "ppm" << tab << setw(8) << left << "ppm" << tab << setw(8) << left << "ppm" << endl;

//
  for (int i = 0; i < n; i++)
  {
    elem_name.clear();
    for (int j = 0; j < 2; j++)
    {
      if ((j == 1) && ((el_sym[i][0] == 'V') || (el_sym[i][0] == 'W'))) continue;
      elem_name += el_sym[i][j];
    }

    if ((el_sym[i][0] == 'A') && (el_sym[i][1] == 'l')) elem_name += "2O3";
    else if ((el_sym[i][0] == 'S') && (el_sym[i][1] == 'i')) elem_name += "O2";
    else if ((i < (n - m + 2)) && (el_sym[i][0] != 'N')) elem_name += "O"; // Na is excluded as oxide
    outfile6 << setw(8) << left << elem_name << tab;
  }
  outfile6 << setw(8) << left << "Nb/Ta" << tab;
  outfile6 << setw(8) << left << "H2O" << tab;
  outfile6 << setw(8) << left << el_sym_S << tab;
  outfile6 << setw(8) << left << el_sym_C << endl;

  if (tm_col == 1)
    outfile6 << setw(8) << left << "Ma" << tab << setw(8) << left << "cMe" << tab << setw(8) << left << "ppm" << tab;
  for (int i = 0; i < n; i++)
    if (sil_c[i] == 0) outfile6 << setw(8) << left << "wt%" << tab;
    else if (sil_c[i] == 1) outfile6 << setw(8) << left << "ppm" << tab;
    else if (sil_c[i] == 2) outfile6 << setw(8) << left << "ppb" << tab;
  outfile6 << setw(8) << left << " " << tab << setw(8) << left << "ppm" << tab << setw(8) << left << "ppm" << tab << setw(8) << left << "ppm" << endl;

//
  for (int i = 0; i < n; i++)
  {
    elem_name.clear();
    for (int j = 0; j < 2; j++)
    {
      if ((j == 1) && ((el_sym[i][0] == 'V') || (el_sym[i][0] == 'W'))) continue;
      elem_name += el_sym[i][j];
    }

    if ((el_sym[i][0] == 'A') && (el_sym[i][1] == 'l')) elem_name += "2O3";
    else if ((el_sym[i][0] == 'S') && (el_sym[i][1] == 'i')) elem_name += "O2";
    else if ((i < (n - m + 2)) && (el_sym[i][0] != 'N')) elem_name += "O"; // Na is excluded as oxide
    outfile8 << setw(8) << left << elem_name << tab;
  }
  outfile8 << setw(8) << left << "H2O" << tab;
  outfile8 << setw(8) << left << el_sym_S << tab;
  outfile8 << setw(8) << left << el_sym_C << endl;

  if (tm_col == 1) outfile8 << setw(8) << left << "Ma" << tab << setw(8) << left << "Me" << tab;
  for (int i = 0; i < n; i++) outfile8 << setw(8) << left << "m%" << tab;
  outfile8 << setw(8) << left << "m%" << tab << setw(8) << left << "m%" << tab << setw(8) << left << "m%" << endl;

  // Core composition:

  for (int i = (n - m); i < n; i++)
  {
    elem_name.clear();
    for (int j = 0; j < 2; j++)
    {
      if ((j == 1) && ((el_sym[i][0] == 'V') || (el_sym[i][0] == 'W'))) continue;
      elem_name += el_sym[i][j];
    }
    outfile2 << setw(8) << left << elem_name << tab;
  }
  outfile2 << setw(8) << left << "O" << tab;
  outfile2 << setw(8) << left << "H" << tab;
  outfile2 << setw(8) << left << el_sym_S << tab;
  outfile2 << setw(8) << left << el_sym_C << tab;
  outfile2 << setw(8) << left << "Metal mass" << endl;
  if (tm_col == 1)
    outfile2 << setw(8) << left << "Ma" << tab << setw(8) << left << "Me" << tab << setw(8) << left << "GPa" << tab << setw(8) << left << "K" << tab << setw(8) << left << "D-IW" << tab;
  for (int i = (n - m); i < n; i++)
    if (met_c[i] == 0) outfile2 << setw(8) << left << "wt%" << tab;
    else if (met_c[i] == 1) outfile2 << setw(8) << left << "ppm" << tab;
    else if (met_c[i] == 2) outfile2 << setw(8) << left << "ppb" << tab;
  outfile2 << setw(8) << left << "wt%" << tab << setw(8) << left << "ppm" << tab << setw(8) << left << "wt%" << tab << setw(8) << left << "wt%" << tab << setw(8) << left << "wt%" << endl;

//
  for (int i = (n - m); i < n; i++)
  {
    elem_name.clear();
    for (int j = 0; j < 2; j++)
    {
      if ((j == 1) && ((el_sym[i][0] == 'V') || (el_sym[i][0] == 'W'))) continue;
      elem_name += el_sym[i][j];
    }
    outfile7 << setw(8) << left << elem_name << tab;
  }
  outfile7 << setw(8) << left << "O" << tab;
  outfile7 << setw(8) << left << "H" << tab;
  outfile7 << setw(8) << left << el_sym_S << tab;
  outfile7 << setw(8) << left << el_sym_C << tab;
  outfile7 << setw(8) << left << "Metal mass" << endl;
  if (tm_col == 1)
    outfile7 << setw(8) << left << "Ma" << tab << setw(8) << left << "cMe" << tab << setw(8) << left << "wt%" << tab;
  for (int i = (n - m); i < n; i++)
    if (met_c[i] == 0) outfile7 << setw(8) << left << "wt%" << tab;
    else if (met_c[i] == 1) outfile7 << setw(8) << left << "ppm" << tab;
    else if (met_c[i] == 2) outfile7 << setw(8) << left << "ppb" << tab;
  outfile7 << setw(8) << left << "wt%" << tab << setw(8) << left << "ppm" << tab << setw(8) << left << "wt%" << tab << setw(8) << left << "wt%" << tab << setw(8) << left << "wt%" << endl;

//
  for (int i = (n - m); i < n; i++)
  {
    elem_name.clear();
    for (int j = 0; j < 2; j++)
    {
      if ((j == 1) && ((el_sym[i][0] == 'V') || (el_sym[i][0] == 'W'))) continue;
      elem_name += el_sym[i][j];
    }
    outfile9 << setw(8) << left << elem_name << tab;
  }
  outfile9 << setw(8) << left << "O" << tab;
  outfile9 << setw(8) << left << "H" << tab;
  outfile9 << setw(8) << left << el_sym_S << tab;
  outfile9 << setw(8) << left << el_sym_C << endl;

  if (tm_col == 1) outfile9 << setw(8) << left << "Ma" << tab << setw(8) << left << "Me" << tab;
  for (int i = (n - m); i < n; i++) outfile9 << setw(8) << left << "m%" << tab;
  outfile9 << setw(8) << left << "m%" << tab << setw(8) << left << "m%" << tab << setw(8) << left << "m%" <<  tab << setw(8) << left << "m%" << endl;

  // Compositions of core-forming increments.

  if (tm_col == 1)
  {
    for (int i = (n - m); i < n; i++)
    {
      elem_name.clear();
      for (int j = 0; j < 2; j++)
      {
        if ((j == 1) && ((el_sym[i][0] == 'V') || (el_sym[i][0] == 'W'))) continue;
        elem_name += el_sym[i][j];
      }
      outfile3 << setw(8) << left << elem_name << tab;
    }
    outfile3 << setw(8) << left << "O" << tab;
    outfile3 << setw(8) << left << "H" << tab;
    outfile3 << setw(8) << left << el_sym_S << tab;
    outfile3 << setw(8) << left << el_sym_C << tab;
    outfile3 << setw(8) << left << "Core frac" << endl;
    outfile3 << setw(8) << left << "Ma" << tab << setw(8) << left << "Me" << tab << setw(8) << left << "GPa" << tab << setw(8) << left << "K" << tab << setw(8) << left << "D-IW" << tab;
    for (int i = (n - m); i < n; i++)
      if (met_c[i] == 0) outfile3 << setw(8) << left << "wt%" << tab;
      else if (met_c[i] == 1) outfile3 << setw(8) << left << "ppm" << tab;
      else if (met_c[i] == 2) outfile3 << setw(8) << left << "ppb" << tab;
    outfile3 << setw(8) << left << "wt%" << tab << setw(8) << left << "ppm" << tab << setw(8) << left << "wt%" << tab << setw(8) << left << "wt%" << tab << setw(8) << left << "wt%" << endl;
  }

  return;
}

//____________________________________________________________________________________
//_________________________________________________________________________________

void write_headings_em()
{
  if (tm_col == 1)
  {
    outfile4 << setw(8) << left << "Time" << tab << setw(8) << left << "Mass" << tab;
    outfile5 << setw(8) << left << "Time" << tab << setw(8) << left << "Mass" << tab << setw(8) << left << "P" << tab << setw(8) << left << "T" << tab << setw(8) << left << "logfO2" << tab;
  }

  string elem_name;
  // Mantle composition.
  for (int i = 0; i < n; i++)
  {
    elem_name.clear();
    for (int j = 0; j < 2; j++)
    {
      if ((j == 1) && ((el_sym[i][0] == 'V') || (el_sym[i][0] == 'W'))) continue;
      elem_name += el_sym[i][j];
    }

    if ((el_sym[i][0] == 'A') && (el_sym[i][1] == 'l')) elem_name += "2O3";
    else if ((el_sym[i][0] == 'S') && (el_sym[i][1] == 'i')) elem_name += "O2";
    else if ((i < (n - m + 2)) && (el_sym[i][0] != 'N')) elem_name += "O"; // Na is excluded as oxide
    outfile4 << setw(8) << left << elem_name << tab;
  }
  outfile4 << setw(8) << left << "Nb/Ta" << tab;
  outfile4 << setw(8) << left << "H2O" << tab;
  outfile4 << setw(8) << left << el_sym_S << tab;
  outfile4 << setw(8) << left << el_sym_C << endl;

  if (tm_col == 1)
    outfile4 << setw(8) << left << "Ma" << tab << setw(8) << left << "Me" << tab;
  for (int i = 0; i < n; i++)
    if (sil_c[i] == 0) outfile4 << setw(8) << left << "wt%" << tab;
    else if (sil_c[i] == 1) outfile4 << setw(8) << left << "ppm" << tab;
    else if (sil_c[i] == 2) outfile4 << setw(8) << left << "ppb" << tab;
  outfile4 << setw(8) << " " << tab << setw(8) << left << "ppm" << tab << setw(8) << left << "ppm" << tab << setw(8) << left << "ppm" << endl;

  // Core composition:

  for (int i = (n - m); i < n; i++)
  {
    elem_name.clear();
    for (int j = 0; j < 2; j++)
    {
      if ((j == 1) && ((el_sym[i][0] == 'V') || (el_sym[i][0] == 'W'))) continue;
      elem_name += el_sym[i][j];
    }
    outfile5 << setw(8) << left << elem_name << tab;
  }
  outfile5 << setw(8) << left << "O" << tab;
  outfile5 << setw(8) << left << "H" << tab;
  outfile5 << setw(8) << left << el_sym_S << tab;
  outfile5 << setw(8) << left << el_sym_C << tab;
  outfile5 << setw(8) << left << "Metal mass" << endl;
  if (tm_col == 1)
    outfile5 << setw(8) << left << "Ma" << tab << setw(8) << left << "Me" << tab << setw(8) << left << "GPa" << tab << setw(8) << left << "K" << tab << setw(8) << left << "D-IW" << tab;
  for (int i = (n - m); i < n; i++)
    if (met_c[i] == 0) outfile5 << setw(8) << left << "wt%" << tab;
    else if (met_c[i] == 1) outfile5 << setw(8) << left << "ppm" << tab;
    else if (met_c[i] == 2) outfile5 << setw(8) << left << "ppb" << tab;
  outfile5 << setw(8) << left << "wt%" << tab << setw(8) << left << "ppm" << tab << setw(8) << left << "wt%" << tab << setw(8) << left << "wt%" << tab << setw(8) << left << "wt%" << endl;

  return;
}

//_____________________________________________________________________________________________________________

double pressure(double object_mass)
{
  const double earth_mass = 5.97219E24; // Mass of Earth in kg
  const double moon_mass = 7.34767E22;      // Mass of Moon in kg
  object_mass = object_mass * moon_mass;
  const double pi = 3.141593;
  const double g_const = 6.673e-11;   // Gravitational constant
  double r_ratio = 0.5;
  double density_ratio = 2.5;
  double bulk_density = 3739.3 + 2.9803e-22 * object_mass;
  double radius = pow((3.0 * object_mass / bulk_density / 4.0 / pi), 0.33333);
  double mantle_density = bulk_density / (1.0 + (density_ratio - 1) * pow(r_ratio, 3.0));
  double p_cmb = 2.0 * pi * g_const * radius * radius * mantle_density * (1.0 - r_ratio) / 3.0;
  p_cmb = p_cmb * (mantle_density * (1.0 + r_ratio) + 2.0 * (density_ratio * mantle_density - mantle_density)
                   * r_ratio * r_ratio);
  p_cmb = p_cmb / 1.0e9;

  return p_cmb;
}

//____________________________________________________________________________________

double temperature(double P)
{
  // This new formulation is based on Frost & Liebske (2012), Tronnes & Frost (2002) and Zhang & Herzberg
  double temp(0);
  if (P < 24.0)
    temp = 1874 + 55.43 * P - 1.739 * P * P + 0.01932 * P * P * P;
  else if (P >= 24.0)
    temp = 1249 + 58.28 * P - 0.395 * P * P + 0.001 * P * P * P;
  return temp;
}

//____________________________________________________________________________________



//____________________________________________________________________________________

// void distrib_coeff(double (&kd)[n], double P, double T, double Xs)
// {
//   double sig_logKd[n] = {0.0};

//   for (int i = (n - m + 1); i < n; i++)
//   {
//     if (i == (n - m + 1)) // this means Si; Value of sig_si specified in file input_2*.txt - fraction of sigma added/subtracted
//     {
//       sig_logKd[i] = sqrt(pow(sig_a[i], 2.0) + pow (sig_b[i] / T, 2.0)
//                           + pow(25.0 * P / T, 2.0) + pow(0.84 * P * P / T, 2.0)  + pow(0.007 * P * P * P / T, 2.0));
//       kd[i] = exp((a[i] + b[i] / T + (-155.0 * P + 2.26 * P * P - 0.011 * P * P * P) / T  + sig_si * sig_logKd[i]) * 2.3026);// Kd(Si)
//     }

//     else if  (i != (n - m + 1) && el_type[i] != 'h')                                    // Ni and trace elements but not HSEs
//     {
//       if (P < 5.0)                                                // low pressure <5GPa, use low pressure coeff for Ni, Co
//       {
//           sig_logKd[i] = sqrt(pow(lowp_sig_a[i], 2.0) + pow (lowp_sig_b[i] / T, 2.0) + pow(lowp_sig_c[i] * P / T, 2.0)); // Does not yet include the S term

//       if (error == 0) //not using HSEs currently
//         kd[i] = exp(   (  lowp_a[i] + (lowp_b[i] / T) + (lowp_c[i] * P / T) + d[i] * log10(1) ) * 2.3026); // no d array for low pressure
//       else if (error == 1)
//         kd[i] = exp((lowp_a[i] + (lowp_b[i] / T) + (lowp_c[i] * P / T) + sig_logKd[i]) * 2.3026); // Currently no S effect because "d" values are all zero for non-HSEs
//       else if (error == -1)
//         kd[i] = exp((lowp_a[i] + (lowp_b[i] / T) + (lowp_c[i] * P / T) - sig_logKd[i]) * 2.3026);
//       }

//       else if (P >= 5.0)                                           // high pressure >5GPa, use high pressure coeff for Ni,Co
//       {
//         sig_logKd[i] = sqrt(pow(sig_a[i], 2.0) + pow (sig_b[i] / T, 2.0) + pow(sig_c[i] * P / T, 2.0)); // Does not yet include the S term

//       if (error == 0)
//         kd[i] = exp(   (  a[i] + (b[i] / T) + (c[i] * P / T) + d[i] * log10(1 - Xs) ) * 2.3026);
//       else if (error == 1)
//         kd[i] = exp((a[i] + (b[i] / T) + (c[i] * P / T) + sig_logKd[i]) * 2.3026); // Currently no S effect because "d" values are all zero for non-HSEs
//       else if (error == -1)
//         kd[i] = exp((a[i] + (b[i] / T) + (c[i] * P / T) - sig_logKd[i]) * 2.3026);
//     }

//     else if  (i != (n - m + 1) && el_type[i] == 'h' && Xs <= 0.35)       // HSEs; metal composition below Fe-FeS eutectic (at 11 GPa)
//     {
//       sig_logKd[i] = sqrt(pow(sig_a[i], 2.0) + pow (sig_b[i] / T, 2.0) + pow(sig_c[i] * P / T, 2.0)); // Does not yet include the S term

//       if (error == 0)
//         kd[i] = exp(   (  a[i] + (b[i] / T) + (c[i] * P / T) + d[i] * log10(1 - Xs) ) * 2.3026);
//       else if (error == 1)
//         kd[i] = exp((a[i] + (b[i] / T) + (c[i] * P / T) + sig_logKd[i]) * 2.3026);
//       else if (error == -1)
//         kd[i] = exp((a[i] + (b[i] / T) + (c[i] * P / T) - sig_logKd[i]) * 2.3026);
//     }

//     else if  (i != (n - m + 1) && el_type[i] == 'h' && Xs > 0.35)        // HSEs; metal/sulfide composition above Fe-FeS eutectic (at 11 GPa)
//     {
//       sig_logKd[i] = sqrt(pow(sulf_sig_a[i], 2.0) + pow (sulf_sig_b[i] / T, 2.0) + pow(sulf_sig_c[i] * P / T, 2.0));

//       if (error == 0)
//         kd[i] = sulf_a[i] + sulf_b[i] / T + sulf_c[i] * P / T;
//       else if (error == 1)
//         kd[i] = (sulf_a[i] + sulf_b[i] / T + sulf_c[i] * P / T) + sig_logKd[i];
//       else if (error == -1)
//         kd[i] = (sulf_a[i] + sulf_b[i] / T + sulf_c[i] * P / T) - sig_logKd[i];
//     }
//   }

//   return;
// }

void distrib_coeff(double (&kd)[n], double P, double T, double Xs)
{
  double sig_logKd[n] = {0.0};

  for (int i = (n - m + 1); i < n; i++)
  {
    if (i == (n - m + 1)) // this means Si; Value of sig_si specified in file input_2*.txt - fraction of sigma added/subtracted
    {
      sig_logKd[i] = sqrt(pow(sig_a[i], 2.0) + pow (sig_b[i] / T, 2.0)
                          + pow(25.0 * P / T, 2.0) + pow(0.84 * P * P / T, 2.0)  + pow(0.007 * P * P * P / T, 2.0));
      kd[i] = exp((a[i] + b[i] / T + (-155.0 * P + 2.26 * P * P - 0.011 * P * P * P) / T  + sig_si * sig_logKd[i]) * 2.3026);// Kd(Si)
    }

    else if  (i != (n - m + 1) && el_type[i] != 'h')                                    // Ni and trace elements but not HSEs
    {
      //   if (P < 5.0)                                                // low pressure <5GPa, use low pressure coeff for Ni, Co
      //   {
      //       sig_logKd[i] = sqrt(pow(lowp_sig_a[i], 2.0) + pow (lowp_sig_b[i] / T, 2.0) + pow(lowp_sig_c[i] * P / T, 2.0)); // Does not yet include the S term

      //   if (error == 0)
      //     kd[i] = exp((lowp_a[i] + (lowp_b[i] / T) + (lowp_c[i] * P / T)) * 2.3026); // no d array for low pressure
      //   else if (error == 1)
      //     kd[i] = exp((lowp_a[i] + (lowp_b[i] / T) + (lowp_c[i] * P / T) + sig_logKd[i]) * 2.3026); // Currently no S effect because "d" values are all zero for non-HSEs
      //   else if (error == -1)
      //     kd[i] = exp((lowp_a[i] + (lowp_b[i] / T) + (lowp_c[i] * P / T) - sig_logKd[i]) * 2.3026);
      //   }

      //else if (P >= 5.0)                                           // high pressure >5GPa, use high pressure coeff for Ni,Co
      //{
      sig_logKd[i] = sqrt(pow(sig_a[i], 2.0) + pow (sig_b[i] / T, 2.0) + pow(sig_c[i] * P / T, 2.0)); // Does not yet include the S term

      if (error == 0)
        kd[i] = exp((a[i] + (b[i] / T) + (c[i] * P / T) + d[i] * log10(1 - Xs) ) * 2.3026);
      else if (error == 1)
        kd[i] = exp((a[i] + (b[i] / T) + (c[i] * P / T) + sig_logKd[i]) * 2.3026); // Currently no S effect because "d" values are all zero for non-HSEs
      else if (error == -1)
        kd[i] = exp((a[i] + (b[i] / T) + (c[i] * P / T) - sig_logKd[i]) * 2.3026);
      //}
    }

    else if  (i != (n - m + 1) && el_type[i] == 'h' && Xs <= 0.35)       // HSEs; metal composition below Fe-FeS eutectic (at 11 GPa)
    {
      sig_logKd[i] = sqrt(pow(sig_a[i], 2.0) + pow (sig_b[i] / T, 2.0) + pow(sig_c[i] * P / T, 2.0)); // Does not yet include the S term

      if (error == 0)
        kd[i] = exp((a[i] + (b[i] / T) + (c[i] * P / T) + d[i] * log10(1 - Xs) ) * 2.3026);
      else if (error == 1)
        kd[i] = exp((a[i] + (b[i] / T) + (c[i] * P / T) + sig_logKd[i]) * 2.3026);
      else if (error == -1)
        kd[i] = exp((a[i] + (b[i] / T) + (c[i] * P / T) - sig_logKd[i]) * 2.3026);
    }

    else if  (i != (n - m + 1) && el_type[i] == 'h' && Xs > 0.35)        // HSEs; metal/sulfide composition above Fe-FeS eutectic (at 11 GPa)
    {
      sig_logKd[i] = sqrt(pow(sulf_sig_a[i], 2.0) + pow (sulf_sig_b[i] / T, 2.0) + pow(sulf_sig_c[i] * P / T, 2.0));

      if (error == 0)
        kd[i] = sulf_a[i] + sulf_b[i] / T + sulf_c[i] * P / T;
      else if (error == 1)
        kd[i] = (sulf_a[i] + sulf_b[i] / T + sulf_c[i] * P / T) + sig_logKd[i];
      else if (error == -1)
        kd[i] = (sulf_a[i] + sulf_b[i] / T + sulf_c[i] * P / T) - sig_logKd[i];
    }
  }
  return;
}
//____________________________________________________________________________________


// Oxygen partitioning model (Frost et al 2010, JGR)
double oxygen_kd_Frost(double xfeo_mw, double P, double T)
{
  const double gas_const = 8.314;
  const double kprime_liq = 4.397;
  const double kprime_mw = 4.3;
  const double v0_liq = 13.244;
  const double v0_mw = 12.25;
  const double a0_liq = 4.923E-5;
  const double a1_liq = 2.968E-9;
  const double a2_liq = -0.0806;
  const double a3_liq = -0.0014437;
  const double a0_mw = 3.481E-5;
  const double a1_mw = 2.968E-9;
  const double a2_mw = -0.0806;
  const double a3_mw = -0.0014437;

  double g0_feo_mw(0), g0_feo_liq(0), vt_liq(0), k0_mw(0), vt_mw(0), k0_liq(0);
  double int_vdp_liq(0), int_vdp_mw(0), wg_mw(0), rtlngamma_mw(0), rtlngamma_liq(0);
  double gfeo_mw(0), wfe_feo(0), w_feo_fe;
  double o_kd(0);

  // pressure values are multiplied by 10^4 to convert from GPa to bars

  g0_feo_mw = -279318 + 252.848 * T - 46.12826 * T * log(T) - 0.0057402984 * pow(T, 2.0);
  g0_feo_liq = -279318 + 252.848 * T - 46.12826 * T * log(T) - 0.0057402984 * T * T + 34008 - 20.969 * T;
  vt_liq = (v0_liq / 10.0) * exp(a0_liq * (T - 298.15) + (a1_liq / 2.0) * (T * T - 298.15 * 298.15)
                                 - a2_liq * ((1 / T) - (1 / 298.15)) + a3_liq * (log(T) - log(298.15)));
  vt_mw = (v0_mw / 10.0) * exp(a0_mw * (T - 298.15) + a1_mw / 2.0 * (T * T - 298.15 * 298.15) - a2_mw
                               * ((1 / T) - (1 / 298.15)) + a3_mw * (log(T) - log(298.15)));
  k0_mw = 1.5E6 - 264.0 * (T - 298.0) + 0.01906 * (T - 298.0) * (T - 298.15); // should be 298.15!
  k0_liq = 802655 - 100.0 * (T - 298.15);
  int_vdp_liq = vt_liq * (k0_liq / (kprime_liq - 1.0) * (exp(((kprime_liq -
                          1.0) / kprime_liq) * log(1 + kprime_liq / k0_liq * P * 1E4)) - 1.0));
  int_vdp_mw = vt_mw * (k0_mw / (kprime_mw - 1) * (exp(((kprime_mw - 1) / kprime_mw)
                        * log(1 + kprime_mw / k0_mw * P * 1E4)) - 1.0));
  wg_mw = 11000 - 0.17 * (T - 1400.0) + 0.011 * P * 1E4;
  rtlngamma_mw = wg_mw * (1.0 - xfeo_mw) * (1.0 - xfeo_mw);
  gfeo_mw = g0_feo_mw + gas_const * T * log(xfeo_mw) + rtlngamma_mw + int_vdp_mw;
  wfe_feo = 83307.0 - 8.978 * T - 0.0959422008665705 * P * 1E4;
  w_feo_fe = 135943.0 - 31.122 * T - 0.0596202529287227 * P * 1E4;

  double xfeo_liq_guess = 0.01; // arbitrary starting guess for iterations
  double xfeo_liq_guess_new(0);
  double increm = 0.01;
  double deltag_feo = 1E6; // arbitrary initial free energy difference
  double old_deltag_feo(0);
  double gfeo_liq(0);
  double xfeo_liq = 0;

  int count = 0;
  while (fabs(deltag_feo) > 0.01)
  {
    count++;
    if (count > 10000)
    {
      //cout << endl << "Oxygen Kd iterations have reached 10000 at time = " << t << endl;
      //cout << "Fraction of total Fe in metal = " << setprecision(4) << xfe_met << endl;
      o_kd = 1e20;
      return o_kd;
    }
    xfeo_liq_guess_new = xfeo_liq_guess - increm;

    if (xfeo_liq_guess_new < 0.000000000) do //ensures that xFeO does not become negative
      {
        increm /= 2.0;
        xfeo_liq_guess_new = xfeo_liq_guess - increm;
      } while (xfeo_liq_guess_new < 0.000000000);

    xfeo_liq_guess = xfeo_liq_guess_new;
    old_deltag_feo = deltag_feo;

    rtlngamma_liq = (1 - xfeo_liq_guess) * (1 - xfeo_liq_guess) *
                    (w_feo_fe + (2 * (wfe_feo - w_feo_fe) * xfeo_liq_guess));

    gfeo_liq = g0_feo_liq + (gas_const * T * log(xfeo_liq_guess)) + rtlngamma_liq
               + int_vdp_liq;
    deltag_feo = gfeo_mw - gfeo_liq;

    if (fabs(deltag_feo) > fabs(old_deltag_feo))
      increm /= -2.0;
  }

  xfeo_liq = xfeo_liq_guess;
  o_kd = xfeo_liq / ((1 + xfeo_liq) * (1 + xfeo_liq) * xfeo_mw);

  return o_kd;
}
//_________________________________________________________________________________

double oxygen_kd_Rubie(double moles[], double moles_oxygen, double moles_S_met,
                       double moles_S_sil, double moles_C_met, double moles_C_sil, double pkd[],
                       double feo_new, double & nio_new, double & sio2_new,
                       double & fe_new, double & ni_new, double & si_new, double & o_new,
                       double & xfeo_mw, double & xsi_liq, double & xfe_liq, double & xo_liq, double & xni_liq)
{
// The partition coefficient of O is the return variable
  double o_kd_new;

  double moles_TOL = 1.0e-15;
// Total moles of lithophile elements
  double moles_lith = 0.0;
  for (int i = 0; i < (n - m); i++) moles_lith += moles[i];
  moles_lith += moles_S_sil;
  moles_lith += moles_C_sil;
// For reference, below are the indexes for elements in elements.h
// Element index: 0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22
// Element name: Al Mg Ca Na Fe Si Ni Co Nb Ta  V Cr Pt Pd Ru Ir  W Mo Mn Ga Zn In Cu
// Element index:            18 19 20 21 22 23 24 25 26 27 28 29 30 31
// Metals start at Fe at index n then continue to n+m where m = n - 4
// If n = 17, then only up to W is included

// For convenience, S and C don't change
  double fe_old = moles[n];
  double si_old = moles[n + 1];
  double ni_old = moles[n + 2];
  double o_old = moles_oxygen;
  double feo_old = moles[n - m];
  double sio2_old = moles[n - m + 1];
  double nio_old = moles[n - m + 2];
  double s_old = moles_S_met;
  double c_old = moles_C_met;
  double so2_old = moles_S_sil;
  double co2_old = moles_C_sil;

  double s_new = s_old;
  double c_new = c_old;
  double so2_new = so2_old;
  double co2_new = co2_old;

// feo_new is established by a guess input to oxygen_kd_Rubie
// Check to see if it is a bad guess
  if (feo_new < 0.0 || feo_new > feo_old + fe_old)
  {
    nio_new = -1.0; sio2_new = -1.0; fe_new = -1.0; ni_new = -1.0; si_new = -1.0; o_new = -1.0; xfeo_mw = -1.0; xsi_liq = -1.0; xfe_liq = -1.0; xo_liq = -1.0; xni_liq = -1.0; o_kd_new = -1.0e10;
    return o_kd_new;
  }
  if (feo_new < moles_TOL) feo_new = moles_TOL;

// fe_new is established by Eqs S-14a (Rubie et al. 2011)
  fe_new = feo_old + fe_old - feo_new;
// Check to see if there are not enough Fe atoms
  if (fe_new < 0.0)
  {
    nio_new = -1.0; sio2_new = -1.0; fe_new = -1.0; ni_new = -1.0; si_new = -1.0; o_new = -1.0; xfeo_mw = -1.0; xsi_liq = -1.0; xfe_liq = -1.0; xo_liq = -1.0; xni_liq = -1.0; o_kd_new = -1.0e10;
    return o_kd_new;
  }
  else if (fe_new < moles_TOL)
  {
    fe_new = moles_TOL;
    feo_new = feo_old + fe_old - fe_new;
  }

// nio_new is established by Eq S-15 (Rubie et al. 2011):
  nio_new = feo_new * (nio_old + ni_old) / ((feo_old + fe_old - feo_new) * pkd[n - m + 2] + feo_new);
// Check to see if there are not enough Ni atoms
  if (nio_new < 0.0 || nio_new > nio_old + ni_old)
  {
    nio_new = -1.0; sio2_new = -1.0; fe_new = -1.0; ni_new = -1.0; si_new = -1.0; o_new = -1.0; xfeo_mw = -1.0; xsi_liq = -1.0; xfe_liq = -1.0; xo_liq = -1.0; xni_liq = -1.0; o_kd_new = -1.0e10;
    return o_kd_new;
  }

// ni_new is established by Eqs S-14b (Rubie et al. 2011):
  ni_new = nio_old + ni_old - nio_new;
// Check to see if there are not enough Ni atoms
  if (ni_new < 0.0 || ni_new > nio_old + ni_old)
  {
    nio_new = -1.0; sio2_new = -1.0; fe_new = -1.0; ni_new = -1.0; si_new = -1.0; o_new = -1.0; xfeo_mw = -1.0; xsi_liq = -1.0; xfe_liq = -1.0; xo_liq = -1.0; xni_liq = -1.0; o_kd_new = -1.0e10;
    return o_kd_new;
  }

// sio2_new is established by solving Eq 17 (Rubie et al. 2011):
  double alpha = sio2_old + si_old;
  double gamma = fe_new + ni_new + s_new + c_new + feo_old + nio_old + 3.0 * sio2_old + o_old - feo_new - nio_new + si_old;
  double sigma = feo_new + nio_new + so2_new + co2_new + moles_lith;
// Calculate quadratic coefficients:
  double ax1 = 3.0 * feo_new * feo_new - fe_new * fe_new * pkd[n - m + 1];
  double ax2 = -(((gamma + 3.0 * alpha) * feo_new * feo_new) + sigma * pkd[n - m + 1] * fe_new * fe_new);
  double ax3 = alpha * gamma * feo_new * feo_new;
// quadratic solution follows:
  double ax4 = -ax2 - sqrt(ax2 * ax2 - 4.0 * ax1 * ax3); // minus before 2nd term is OK
  sio2_new = ax4 / (2 * ax1);
// Check to see if there are not enough Si atoms
  if (sio2_new < 0.0 || sio2_new > sio2_old + si_old )
  {
    nio_new = -1.0; sio2_new = -1.0; fe_new = -1.0; ni_new = -1.0; si_new = -1.0; o_new = -1.0; xfeo_mw = -1.0; xsi_liq = -1.0; xfe_liq = -1.0; xo_liq = -1.0; xni_liq = -1.0; o_kd_new = -1.0e10;
    return o_kd_new;
  }

// Next are mass balance equations 9c and 9d:
  si_new = si_old + sio2_old - sio2_new;
  if (si_new < 0.0 || si_new > sio2_old + si_old )
  {
    nio_new = -1.0; sio2_new = -1.0; fe_new = -1.0; ni_new = -1.0; si_new = -1.0; o_new = -1.0; xfeo_mw = -1.0; xsi_liq = -1.0; xfe_liq = -1.0; xo_liq = -1.0; xni_liq = -1.0; o_kd_new = -1.0e10;
    return o_kd_new;
  }

//
  o_new = feo_old + nio_old + 2 * sio2_old + o_old - feo_new - nio_new - 2 * sio2_new;
  if (o_new < 0.0)
  {
    nio_new = -1.0; sio2_new = -1.0; fe_new = -1.0; ni_new = -1.0; si_new = -1.0; o_new = -1.0; xfeo_mw = -1.0; xsi_liq = -1.0; xfe_liq = -1.0; xo_liq = -1.0; xni_liq = -1.0; o_kd_new = -1.0e10;
    return o_kd_new;
  }

  // Next calculate XFeO(MW) using Eq 5:
  double xfeo_new = feo_new / (feo_new + nio_new + sio2_new + moles_lith + so2_new + co2_new);
  xfeo_mw = 1.1478 * xfeo_new + 1.3189 * pow(xfeo_new, 2.0);
  // SAJ: rough fix
  if (xfeo_mw > 1.0) xfeo_mw = 1.0;

  // Calculate oxygen Kd using Frost et al (2010) model:
  xsi_liq = si_new / (si_new + fe_new + ni_new + o_new + s_new + c_new);
  xfe_liq = fe_new / (si_new + fe_new + ni_new  + o_new + s_new + c_new);
  xo_liq = o_new / (si_new + fe_new + ni_new + o_new + s_new + c_new);
  xni_liq = ni_new / (si_new + fe_new + ni_new + o_new + s_new + c_new);
// Next 3 expressions calculate lefthand side of Eq 6:
  double a1 = fe_new * o_new;
  double a2 = xfeo_mw * pow((fe_new + ni_new + o_new + si_new + s_new + c_new), 2.0);
  o_kd_new = a1 / a2;

  if (o_kd_new != o_kd_new || o_kd_new < 0 || nio_new < 0 || sio2_new < 0 || fe_new < 0 || ni_new < 0 || si_new < 0 || o_new < 0 || xfeo_mw > 1 || xfe_liq > 1 || xo_liq > 1 || xni_liq > 1 || xsi_liq > 1 || xfeo_mw < 0 || xfe_liq < 0 || xo_liq < 0 || xni_liq < 0  || xsi_liq < 0)
  {
    nio_new = -1.0; sio2_new = -1.0; fe_new = -1.0; ni_new = -1.0; si_new = -1.0; o_new = -1.0; xfeo_mw = -1.0; xsi_liq = -1.0; xfe_liq = -1.0; xo_liq = -1.0; xni_liq = -1.0; o_kd_new = -1.0e10;
    return o_kd_new;
  }

  return o_kd_new;
}


// MASS BALANCE CALCULATION "REACT"
void react(double (&new_moles)[nm3], double moles[], double moles_oxygen, double moles_S_met, double moles_S_sil, double moles_C_met, double moles_C_sil, double P, double T, double pkd[])
{
  double moles_TOL = 1.0e-15;
// react counter; counts ALL reacts
  react_call_no++;

// Triggers print and exit statements at a specific react
  int react_call_no_test = -1;
  if (react_call_no == react_call_no_test) cout << endl;
  if (react_call_no == react_call_no_test) cout << t << endl;
  ofstream reactfile, reactfile2, reactfile3, reactfile4;
  if (react_call_no == react_call_no_test) reactfile.open ("reactfile.txt");
  if (react_call_no == react_call_no_test) reactfile2.open ("reactfile2.txt");
  if (react_call_no == react_call_no_test) reactfile3.open ("reactfile3.txt");
  if (react_call_no == react_call_no_test) reactfile4.open ("reactfile4.txt");


// Declarations of new mole values:
  double feo_new(0.0), nio_new(0.0), sio2_new(0.0);  // silicate values
  double fe_new(0.0), ni_new(0.0), si_new(0.0), o_new(0.0); //metal values
  double te_new[m - 3] = {0.0}, teo_new[m - 3] = {0.0}; // "3" - Fe, Si and Ni
  double xfeo_mw(0.0);
  double xsi_liq(-1.0), xfe_liq(-1.0), xo_liq(-1.0), xni_liq(-1.0);
  double zero_xfe_liq = 1.0;

  double o_kd, o_kd_new;
  double delta, delta_neg, delta_pos;
  double old_delta, delta_max, delta_min;
  double step;

// Total moles of lithophile elements
  double moles_lith = 0.0;
  for (int i = 0; i < (n - m); i++) moles_lith += moles[i];

// moles[n-m] gives value for FeO in original equilibrating liquid
  double moles_feo_old = moles[n - m];
// moles[n] gives value for Fe in original equilibrating liquid
  double moles_fe_old = moles[n];
// guess for equilibrated FeO abundance
  double moles_feo_guess;

// Map delta if write statements are set
  if (react_call_no == react_call_no_test)
  {
    moles_feo_guess = 1.0e-8;
    while (moles_feo_guess < 1.0e-5)
    {
      o_kd_new = oxygen_kd_Rubie(moles, moles_oxygen, moles_S_met, moles_S_sil, moles_C_met, moles_C_sil, pkd, moles_feo_guess, nio_new, sio2_new, fe_new, ni_new, si_new, o_new, xfeo_mw, xsi_liq, xfe_liq, xo_liq, xni_liq);
      o_kd = oxygen_kd_Frost(xfeo_mw, P, T);
      delta = o_kd_new - o_kd;
      if (o_kd_new >= 0 && nio_new >= 0.0 && sio2_new >= 0.0 && fe_new >= 0.0 && ni_new >= 0.0 && si_new >= 0.0 && o_new >= 0.0 && xfeo_mw >= 0.0 && xsi_liq >= 0.0 && xfe_liq >= 0.0 && xo_liq >= 0.0 && xni_liq >= 0.0 && xfeo_mw <= 1.0 && xfe_liq <= 1.0 && xo_liq <= 1.0 && xni_liq <= 1.0 && xsi_liq <= 1.0 )
      {
        reactfile << 1 << tab << scientific << setprecision(20) << moles_feo_guess << tab << o_kd_new << tab << o_kd << tab << delta << endl;
        //reactfile << 1 << tab << scientific << setprecision(20) << moles_feo_guess << tab << o_kd_new << tab << o_kd << tab << delta << tab << nio_new << tab << sio2_new << tab << fe_new << tab << ni_new << tab << si_new << tab << o_new << tab << xfeo_mw << tab << xsi_liq << tab << xfe_liq << tab << xo_liq << tab << xni_liq << endl;
      }
      else
      {
        reactfile << 0 << tab << scientific << setprecision(20) << moles_feo_guess << tab << o_kd_new << tab << o_kd << tab << delta << endl;
        //reactfile << 0 << tab << scientific << setprecision(20) << moles_feo_guess << tab << o_kd_new << tab << o_kd << tab << delta << tab << nio_new << tab << sio2_new << tab << fe_new << tab << ni_new << tab << si_new << tab << o_new << tab << xfeo_mw << tab << xsi_liq << tab << xfe_liq << tab << xo_liq << tab << xni_liq << endl;
      }
      moles_feo_guess = moles_feo_guess * pow(10, 0.0001);
    }
  }

// This is the highest possible value of FeO in the equilibrating liquid
  double max_moles_feo_guess = moles_feo_old + moles_fe_old;
// This is the lowest possible value of FeO in the equilibrating liquid
  double min_moles_feo_guess = 0.0;

// In the next section the goal is just to get ONE good guess between the min and max possible
// For later convenience, we start guessing high and keep trying lower and lower
  moles_feo_guess = max_moles_feo_guess;
  o_kd_new = oxygen_kd_Rubie(moles, moles_oxygen, moles_S_met, moles_S_sil, moles_C_met, moles_C_sil, pkd, moles_feo_guess, nio_new, sio2_new, fe_new, ni_new, si_new, o_new, xfeo_mw, xsi_liq, xfe_liq, xo_liq, xni_liq);
  int count = 0;
  while (o_kd_new < 0 || nio_new < 0 || sio2_new < 0 || fe_new < 0 || ni_new < 0 || si_new < 0 || o_new < 0 || xfeo_mw > 1 || xfe_liq > 1 || xo_liq > 1 || xni_liq > 1 || xsi_liq > 1 || xfeo_mw < 0 || xfe_liq < 0 || xo_liq < 0 || xni_liq < 0  || xsi_liq < 0)
  {
    count++;
// This steps through a number of different step sizes allowing the code to run faster
    if (count < 100)
    {
      moles_feo_guess = moles_feo_guess / pow(10.0, 0.1);
    }
    else if (count == 100) moles_feo_guess = max_moles_feo_guess / pow(10.0, 0.01);
    else if (count < 1100)
    {
      moles_feo_guess = moles_feo_guess / pow(10.0, 0.01);
    }
    else if (count == 1100) moles_feo_guess = max_moles_feo_guess / pow(10.0, 0.001);
    else if (count < 11100)
    {
      moles_feo_guess = moles_feo_guess / pow(10.0, 0.001);
    }
    else
    {
      cout << "First count in react() exceeded 11100 at " << react_call_no << endl;
      exit(1);
    }
    o_kd_new = oxygen_kd_Rubie(moles, moles_oxygen, moles_S_met, moles_S_sil, moles_C_met, moles_C_sil, pkd, moles_feo_guess, nio_new, sio2_new, fe_new, ni_new, si_new, o_new, xfeo_mw, xsi_liq, xfe_liq, xo_liq, xni_liq);
  }
  if (react_call_no == react_call_no_test) cout << "One good guess of moles FeO: " << scientific << setprecision(20) << moles_feo_guess << tab << "Okd: " << o_kd_new << endl;

// Search for the maximum guess of FeO for which oxygen_kd_Rubie returns a physically plausible result
  double max_moles_feo_good_guess = moles_feo_guess;
// Placemarker for the minimum no good guess greater than the max_moles_feo_good_guess
  double max_moles_feo_nogood_guess = max_moles_feo_guess;
  count = 0;
  while (abs(max_moles_feo_nogood_guess - max_moles_feo_good_guess) > moles_TOL)
  {
    count++;
    if (count > 10000)
    {
      cout << "Second count in react() exceeded 10000 at " << react_call_no << endl;
      exit(1);
    }
    moles_feo_guess = (max_moles_feo_nogood_guess + max_moles_feo_good_guess) / 2.0;
    o_kd_new = oxygen_kd_Rubie(moles, moles_oxygen, moles_S_met, moles_S_sil, moles_C_met, moles_C_sil, pkd, moles_feo_guess, nio_new, sio2_new, fe_new, ni_new, si_new, o_new, xfeo_mw, xsi_liq, xfe_liq, xo_liq, xni_liq);
    if (o_kd_new >= 0.0 && nio_new >= 0.0 && sio2_new >= 0.0 && fe_new >= 0.0 && ni_new >= 0.0 && si_new >= 0.0 && o_new >= 0.0 && xfeo_mw >= 0.0 && xsi_liq >= 0.0 && xfe_liq >= 0.0 && xo_liq >= 0.0 && xni_liq >= 0.0 && xfeo_mw <= 1.0 && xfe_liq <= 1.0 && xo_liq <= 1.0 && xni_liq <= 1.0 && xsi_liq <= 1.0 )
    {
// The guess is good
      max_moles_feo_good_guess = moles_feo_guess;
    }
    else
    {
// The guess is no good
      max_moles_feo_nogood_guess = moles_feo_guess;
    }
  }
// delta_max is the difference at the maximum good guess FeO
// if negative, then a zero delta exists
// if positive, then the zero delta is likely unresolved, and the answer is likely that all Fe is in the silicate
  o_kd = oxygen_kd_Frost(xfeo_mw, P, T);
  delta_max = o_kd_new - o_kd;
  if (react_call_no == react_call_no_test) cout << "Maximum good guess of moles FeO: " << scientific << setprecision(20) << max_moles_feo_good_guess << tab << "delta_max: " << delta_max << endl;

// Search for the minimum guess of FeO for which oxygen_kd_Rubie returns a physically plausible result
  double min_moles_feo_good_guess = moles_feo_guess;
// Placemarker for the maximum no good guess less than the min_moles_feo_good_guess
  double min_moles_feo_nogood_guess = min_moles_feo_guess;
  count = 0;
  while (abs(min_moles_feo_good_guess - min_moles_feo_nogood_guess) > moles_TOL)
  {
    count++;
    if (count > 10000)
    {
      cout << "Third count in react() exceeded 10000 at " << react_call_no << endl;
      exit(1);
    }
    moles_feo_guess = (min_moles_feo_nogood_guess + min_moles_feo_good_guess) / 2.0;
    o_kd_new = oxygen_kd_Rubie(moles, moles_oxygen, moles_S_met, moles_S_sil, moles_C_met, moles_C_sil, pkd, moles_feo_guess, nio_new, sio2_new, fe_new, ni_new, si_new, o_new, xfeo_mw, xsi_liq, xfe_liq, xo_liq, xni_liq);
    if (o_kd_new >= 0.0 && nio_new >= 0.0 && sio2_new >= 0.0 && fe_new >= 0.0 && ni_new >= 0.0 && si_new >= 0.0 && o_new >= 0.0 && xfeo_mw >= 0.0 && xsi_liq >= 0.0 && xfe_liq >= 0.0 && xo_liq >= 0.0 && xni_liq >= 0.0 && xfeo_mw <= 1.0 && xfe_liq <= 1.0 && xo_liq <= 1.0 && xni_liq <= 1.0 && xsi_liq <= 1.0)
    {
// The guess is good
      min_moles_feo_good_guess = moles_feo_guess;
    }
    else
    {
// The guess is no good
      min_moles_feo_nogood_guess = moles_feo_guess;
    }
  }
// delta_min is the difference at the minimum good guess FeO
// if positive, then a zero delta exists
// if negative, then the zero delta is likely unresolved and the answer is likely that all Fe is in the metal
  o_kd = oxygen_kd_Frost(xfeo_mw, P, T);
  delta_min = o_kd_new - o_kd;
  if (react_call_no == react_call_no_test) cout << "Minimum good guess of moles FeO: " << scientific << setprecision(20) << min_moles_feo_good_guess << tab << "delta_min: " << delta_min << endl;

// If this is true, then the two estimates of Okd are likely do not cross
  if (delta_min * delta_max >= 0.0)
  {
// If this is true then they are closest when the FeO content is very high
    if (abs(delta_max) <= abs(delta_min))
    {
// Return the answer that the equilibrating liquid is 100% silicate
      zero_xfe_liq = 0.0;          // Indicates no metal
      for (int i = 0; i < n; i++)
      {
        if (i < n - m) new_moles[i] = moles[i]; // Lithophile elements
        else new_moles[i] = moles[i] + moles[i + m]; // Siderophile elements
      }
      for (int i = n; i < (n + m); i++) new_moles[i] = 0.0; // No metal components
      new_moles[n + m] = 0.0; // No oxygen in metal
      new_moles[n + m + 1] = 0.0; // LOGfO2 - arbitarily assigned to be zero when no metal-silicate equilibration
      new_moles[n + m + 2] = zero_xfe_liq;
      no_react_call_no++;
      if (react_call_no == react_call_no_test) exit(1);
      return;
    }
// Conversely, then they are closest when the FeO content is very low
    else
    {
      cout << "here" << tab << react_call_no << endl;
      moles_feo_guess = min_moles_feo_good_guess;
      o_kd_new = oxygen_kd_Rubie(moles, moles_oxygen, moles_S_met, moles_S_sil, moles_C_met, moles_C_sil, pkd, moles_feo_guess, nio_new, sio2_new, fe_new, ni_new, si_new, o_new, xfeo_mw, xsi_liq, xfe_liq, xo_liq, xni_liq);
    }
  }
// In this case, the delta goes through zero and it's just a matter of converging to the best estimate of moles_feo_guess
  else
  {
//
    double moles_feo_neg_guess, moles_feo_pos_guess;
    if (react_call_no == react_call_no_test) cout << "An answer exists between the two extremes" << endl;

    if (delta_min >= 0.0 && delta_max <= 0.0)
    {
      delta_neg = delta_max;
      moles_feo_neg_guess = max_moles_feo_good_guess;
      delta_pos = delta_min;
      moles_feo_pos_guess = min_moles_feo_good_guess;
    }
    else if (delta_min <= 0.0 && delta_max >= 0.0)
    {
      delta_neg = delta_min;
      moles_feo_neg_guess = min_moles_feo_good_guess;
      delta_pos = delta_max;
      moles_feo_pos_guess = max_moles_feo_good_guess;
    }
    count = 0;
    while (abs(moles_feo_pos_guess - moles_feo_neg_guess) > moles_TOL)
    {
      count++;
      if (count > 10000)
      {
        cout << "Fourth count in react() exceeded 10000 at " << react_call_no << endl;
        exit(1);
      }
      moles_feo_guess = (moles_feo_pos_guess + moles_feo_neg_guess) / 2.0;
      o_kd_new = oxygen_kd_Rubie(moles, moles_oxygen, moles_S_met, moles_S_sil, moles_C_met, moles_C_sil, pkd, moles_feo_guess, nio_new, sio2_new, fe_new, ni_new, si_new, o_new, xfeo_mw, xsi_liq, xfe_liq, xo_liq, xni_liq);
      o_kd = oxygen_kd_Frost(xfeo_mw, P, T);
      delta = o_kd_new - o_kd;
      if (react_call_no == react_call_no_test) reactfile2 << scientific << setprecision(20) << moles_feo_guess << tab << delta << tab << o_kd_new << tab << o_kd << endl;
      if (o_kd_new >= 0.0 && nio_new >= 0.0 && sio2_new >= 0.0 && fe_new >= 0.0 && ni_new >= 0.0 && si_new >= 0.0 && o_new >= 0.0 && xfeo_mw >= 0.0 && xsi_liq >= 0.0 && xfe_liq >= 0.0 && xo_liq >= 0.0 && xni_liq >= 0.0 && xfeo_mw <= 1.0 && xfe_liq <= 1.0 && xo_liq <= 1.0 && xni_liq <= 1.0 && xsi_liq <= 1.0)
      {
        if (delta > 0.0)
        {
          delta_pos = delta;
          moles_feo_pos_guess = moles_feo_guess;
        }
        else if (delta == 0.0) break;
        else
        {
          delta_neg = delta;
          moles_feo_neg_guess = moles_feo_guess;
        }
      }
      else
      {
        cout << "Bad FeO guess of " << scientific << setw(20) << moles_feo_guess << " of " << max_moles_feo_good_guess << " at " << react_call_no << endl;
        cout << endl << scientific << setw(20) << min_moles_feo_good_guess << tab << max_moles_feo_good_guess << endl;
        cout << endl << scientific << setw(20) << moles_feo_pos_guess << tab << moles_feo_neg_guess << tab << moles_TOL << endl;
        exit(1);
      }
    }
  }
  if (react_call_no == react_call_no_test) cout << "The answer is: " << scientific << setprecision(20) << moles_feo_guess << tab << "delta: " << delta << endl;
// result of mass balance - moles of FeO
  feo_new = moles_feo_guess;

//  if (react_call_no == react_call_no_test) cout << react_call_no << tab << feo_new << tab << fixed << fe_new * 55.85 / (fe_new * 55.85 + ni_new * 58.71 + si_new * 28.09 + o_new * 16.00) * 100.0 <<  endl;

// END OF ITERATION SECTION - MAJOR ELEMENT CONCENTRATIONS ARE NOW DEFINED
// _________________________________________________________________________________________________________
// If there is too little Fe and Ni in the final metal composition, then its all put in the mantle and assumed to be oxidized there
  if ( (fe_new * 55.85 + ni_new * 58.71 ) / (fe_new * 55.85 + ni_new * 58.71 + si_new * 28.09 + o_new * 16.00) < 0.5)
  {
    if (react_call_no == react_call_no_test) cout << "ended up here" << fixed << setprecision(4) << tab <<  (fe_new * 55.85 + ni_new * 58.71 ) / (fe_new * 55.85 + ni_new * 58.71 + si_new * 28.09 + o_new * 16.00) << endl;
    zero_xfe_liq = 0.0;          // Indicates no metal
    for (int i = 0; i < n; i++)
    {
      if (i < n - m) new_moles[i] = moles[i]; // Lithophile elements
      else new_moles[i] = moles[i] + moles[i + m]; // Siderophile elements
    }
    for (int i = n; i < (n + m); i++) new_moles[i] = 0.0; // No metal components
    new_moles[n + m] = 0.0; // No oxygen in metal
    new_moles[n + m + 1] = 0.0; // LOGfO2 - arbitarily assigned to be zero when no metal-silicate equilibration
    new_moles[n + m + 2] = zero_xfe_liq;
    no_react_call_no++;
    if (react_call_no == react_call_no_test) exit(1);
    return;
  }
// Otherwise the final composition is determined
  else
  {
    double total_moles_sil = 2.0 * moles[0] + moles[1] + moles[2] + moles[3] + feo_new + sio2_new + nio_new + moles_S_sil + moles_C_sil;
    double xalo = moles[0] * 2.0 / total_moles_sil;
    double xsio2 = sio2_new / total_moles_sil;
    double xmgo = moles[1] / total_moles_sil;
    double xcao = moles[2] / total_moles_sil;
    double xfeo = feo_new / total_moles_sil;

// Determine trace element concentrations:
    double tot_sil = moles_lith + feo_new + sio2_new + nio_new + moles_S_sil + moles_C_sil; //total moles of silicate
    double tot_met = fe_new + ni_new + o_new + si_new + moles_S_met + moles_C_met;    //total moles of metal

// Following lines calculate trace element concentrations in metal and silicate
    for (int i = 0; i < (m - 3); i++) // "3" limits the elements to Fe, Si and Ni
    {
// Eq.3.1:
      double te_k = (pkd[n - m + 3 + i]) * pow(((fe_new / tot_met) / (feo_new / tot_sil)), (val[n - m + 3 + i] / 2.0) );
// Eq. 3.3:
      double cc = te_k * (moles[n - m + 3 + i] + moles[n + 3 + i]);
      double dd = te_k + (tot_sil / tot_met);
      te_new[i] = cc / dd;
      teo_new[i] = moles[n - m + 3 + i] + moles[n + 3 + i] - te_new[i];
// SAJ: quick fix to prevent negative values
      if (teo_new[i] < 0.0)
      {
        teo_new[i] = 0.0;
        te_new[i] = moles[n - m + 3 + i] + moles[n + 3 + i];
      }
      else if (te_new[i] < 0.0)
      {
        te_new[i] = 0.0;
        teo_new[i] = moles[n - m + 3 + i] + moles[n + 3 + i];
      }

    }

// calculate fO2:
    double xfe = fe_new / (fe_new + ni_new + si_new + o_new + moles_S_met + moles_C_met);
    double logfo2 = (2 * log(xfeo_mw / xfe)) / 2.30259;
    if (react_call_no == react_call_no_test) cout << endl << t << tab << xfeo_mw << tab << xfe << tab << logfo2 << tab << 5 << endl;
// Write equilibrated silicate and metal molar compositions to array:
    for (int i = 0; i < (n - m); i++) new_moles[i] = moles[i]; // lithophile elements
    new_moles[n - m] = moles_feo_guess; // silicate components
    new_moles[n - m + 1] = sio2_new;
    new_moles[n - m + 2] = nio_new;
    for (int i = 0; i < (m - 3); i++) new_moles[n - m + 3 + i] = teo_new[i];
    new_moles[n] = fe_new;        // metal components
    new_moles[n + 1] = si_new;
    new_moles[n + 2] = ni_new;
    for (int i = 0; i < (m - 3); i++) new_moles[n + 3 + i] = te_new[i];
    new_moles[n + m] = o_new;
    new_moles[n + m + 1] = logfo2; // LOGfO2
    new_moles[n + m + 2] = zero_xfe_liq; // indicator of presence (=1) or absence (=0) of metal

    int react_failed = 0;
    for (int i = 0; i < n + m + 1; i++) if (new_moles[i] < 0.0) react_failed++;
    if (react_failed > 0)
    {
      cout << endl << "REACT HAS FAILED at react call no: " << react_call_no <<  endl;
      for (int i = 0; i < n + m + 1; i++) cout << i << tab << scientific << moles[i] << tab << new_moles[i] << endl;
      exit(1);
    }

    if (react_call_no == react_call_no_test) exit(1);
    return;          // EQUILIBRATED COMPOSITIONS + 2 other parameters
  }
}
//_____________________________________________________________________________________________________

// Function that calculates Ds and the distribution of W and Mo between metal and silicate using Jennings notes
void WMo_partitioning(double (& moles)[nm3], double moles_oxygen, double moles_S_met, double moles_S_sil, double moles_C_met, double moles_C_sil, double P, double T)
{
// For reference, below are the indexes for elements in elements.h
// Element index: 0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22
// Element name: Al Mg Ca Na Fe Si Ni Co Nb Ta  V Cr Pt Pd Ru Ir  W Mo Mo Ga Zn In Cu
// Element index:             n +1 +2 +3 +4 +5 +6 +7 +8 +9+10+11+12+13+14+15+16+17+18
// Metals start at Fe at index n then continue to n+m where m = n - 4
// If n = 17, then only up to W is included
//  MOLAR COMPOSITIONS NEEDED TO BE DEFINED ARE:
//            moles[x]          Number of moles of the xth oxide and metal species where x varies from 0 to n+m-1
//            moles_H2O         Number of moles of water
//            moles_S_sil       Number of moles of sulfur as oxide
//            moles_C_sil       Number of moles of carbon as oxide
//            moles_oxygen      Number of moles of oxygen as metal
//            moles_H           Number of moles of hydrogen as metal
//            moles_S_met       Number of moles of sulfur as metal
//            moles_C_met       Number of moles of carbon as metal
//__________________________________________________________________________________________________________
  int n_comp = 8;     // NUMBER OF COMPONENTS
  double XFeO;
  double XFe, XSi, XS, XO, XNi, XC, XW, XMo;
  double elem[8];

  double moles_sil = 0.0;
  for (int i = 0; i < n; i++)  moles_sil += moles[i];
  moles_sil += moles_S_sil;
  moles_sil += moles_C_sil;

  double moles_met = 0.0;
  for (int i = n; i < (n + m); i++)  moles_met += moles[i];
  moles_met += moles_oxygen;        // oxygen
  moles_met += moles_S_met;
  moles_met += moles_C_met;

  XFeO = moles[n - m] / moles_sil;
  XFe = moles[n] / moles_met;
  XSi = moles[n + 1] / moles_met;
  XNi = moles[n + 2] / moles_met;
  XS = moles_S_met / moles_met;
  XO = moles_oxygen / moles_met;
  XC = moles_C_met / moles_met;
  XW = moles[n + 12] / moles_met;
  XMo = moles[n + 13] / moles_met;

  elem[0] = XFe; // Fe
  elem[1] = XSi; // Si
  elem[2] = XS; // S
  elem[3] = XO; // O
  elem[4] = XNi; // Ni
  elem[5] = XC; // C
  elem[6] = XW; // W
  elem[7] = XMo; // Mo

// declare epsilon and ln values:
  double epsilon[8][8] = {
//    Fe       Si        S        O       Ni        C        W       Mo
    { 0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0},  // Fe
    { 0.0,    8.60,   8.955,   -7.15,    7.50,   9.767,    18.3,   22.00},  // Si
    { 0.0,   8.955,  -5.659, -17.168,    2.17,   6.091,    6.47,    2.27},  // S
    { 0.0,   -7.15, -17.168, -10.507,   1.374, -20.408,    4.16,   1.261},  // O
    { 0.0,    7.50,    2.17,   1.374,   0.119,   2.306,     0.0,     0.0},  // Ni
    { 0.0,   9.767,   6.091, -20.408,   2.036,  12.826,   -7.38,   -7.03},  // C
    { 0.0,    18.3,    6.47,    4.16,     0.0,   -7.38,     0.0,     0.0},  // W
    { 0.0,   22.00,    2.27,   1.261,     0.0,   -7.03,     0.0,    3.94},  // Mo
  };

// Adjust epsilon values for temperature:
  double ep[8][8];           // temperature-corrected epsilon value
  for (int i = 0; i < n_comp; i++) {
    for (int j = 0; j < n_comp; j++) {
      ep[i][j] = epsilon[i][j] * 1873.0 / T;
    }
  }

// CALCULATE ln gamma(Fe) (Eq(23) in Ma 2001)
// (t1...t5 are the 5 successive summation terms in this equation)
  double t1 = 0.0;
  for (int i = 2; i < n_comp; i++) {
    t1 += ep[i][i] * (elem[i] + log(1 - elem[i]));
  }

  double t2 = 0.0;
  for (int j = 2; j < (n_comp - 1); j++) {
    for (int k = (j + 1); k < n_comp; k++) {
      if ((elem[j] > 0) && (elem[k] > 0)) t2 += -ep[j][k] * elem[j] * elem[k] * (1 + (log(1 - elem[j])) / elem[j] + (log(1 - elem[k])) / elem[k]);
    }
  }

  double t3 = 0.0;
  for (int i = 2; i < n_comp; i++) {
    for (int k = 2; k < n_comp; k++) {
      if ((k != i) && (elem[k] > 0)) t3 += ep[i][k] * elem[i] * elem[k] * (1 + (log(1 - elem[k])) / elem[k] - 1 / (1 - elem[i]));
    }
  }

  double t4 = 0.0;
  for (int j = 2; j < (n_comp - 1); j++) {
    for (int k = (j + 1); k < n_comp; k++) {
      t4 += 0.5 * ep[j][k] * elem[j] * elem[j] * elem[k] * elem[k] * (1 / (1 - elem[j]) + 1 / (1 - elem[k]) - 1);
    }
  }

  double t5 = 0.0;
  for (int i = 2; i < n_comp; i++) {
    for (int k = 2; k < n_comp; k++) {
      if (k != i) t5 += -ep[i][k] * elem[i] * elem[i] * elem[k] * elem[k] * (1 / (1 - elem[i]) + 1 / (1 - elem[k]) + elem[i] / (2 * (1 - elem[i]) * (1 - elem[i])) - 1);
    }
  }

  double ln_gamma_fe = t1 + t2 + t3 + t4 + t5;
  //cout << "ln gamma(Fe) = " << setprecision(8) << ln_gamma_fe << tab << "gamma(Fe) = " << exp(ln_gamma_fe) << endl;

// CALCULATE ln gamma of solute components (Eq(24) in Ma 2001)
// (t1...t3 are the 3 successive summation terms in this equation)

// W
  int iW = 6;
// gamma_W_0 is unknown and assigned to 1 in Jennings note
  double gamma_W_0 = 3.0;

  double t1_W = 0.0;
  t1_W = ln_gamma_fe + (1873.0 / T) * log(gamma_W_0) - ep[iW][iW] * log(1 - elem[iW]);

  double t2_W = 0.0;
  for (int k = 2; k < n_comp; k++) {
    if ((k != iW) && (elem[k] > 0)) t2_W += - ep[iW][k] * elem[k] * (1 + (log(1 - elem[k])) / elem[k] - 1 / (1 - elem[iW]));
  }

  double t3_W = 0.0;
  for (int k = 2; k < n_comp; k++) {
    if (k != iW) t3_W += ep[iW][k] * elem[k] * elem[k] * elem[iW] * (1 / (1 - elem[iW]) + 1 / (1 - elem[k]) + elem[iW] / (2 * (1 - elem[iW]) * (1 - elem[iW])) - 1);
  }

  double ln_gamma_W = t1_W + t2_W + t3_W;

  //cout << "ln gamma(W) = " << setprecision(8) << ln_gamma_W << tab << "gamma(W) = " << exp(ln_gamma_W) << endl;

// Fit for K an Kd of W using fits from Jennings note
  double a_W = 0.61; //0.592;
  double b_W = -4091; //-6487;
  double d_W = 0.0; //-9.28;
  double K_W = pow(10.0, a_W + b_W / T + d_W * log10(1 - XFeO));
  double Kd_W = exp(log(K_W) - ln_gamma_W + 3.0 * ln_gamma_fe);
  double D_W = Kd_W * pow(XFe, 3.0) / pow(XFeO, 3.0);
  double total_moles_W = moles[16] + moles[n + 12];
  double moles_W_met = (total_moles_W * D_W * moles_met / moles_sil) / (1.0 + D_W * moles_met / moles_sil);
  double moles_W_sil = total_moles_W - moles_W_met;
  moles[16] = moles_W_sil;
  moles[n + 12] = moles_W_met;

  if (moles[16] != moles[16]) cout << "Kd_W = " << Kd_W << tab << "ln_gamma_W = " << ln_gamma_W << endl;

// Mo
  int iMo = 7;
// gamma_Mo_0 is unknown and assigned to 1 in Jennings note
  double gamma_Mo_0 = 1.00;

  double t1_Mo = 0.0;
  t1_Mo = ln_gamma_fe + (1873.0 / T) * log(gamma_Mo_0) - ep[iMo][iMo] * log(1 - elem[iMo]);

  double t2_Mo = 0.0;
  for (int k = 2; k < n_comp; k++) {
    if ((k != iMo) && (elem[k] > 0)) t2_Mo += - ep[iMo][k] * elem[k] * (1 + (log(1 - elem[k])) / elem[k] - 1 / (1 - elem[iMo]));
  }


  double t3_Mo = 0.0;
  for (int k = 2; k < n_comp; k++) {
    if (k != iMo) t3_Mo += ep[iMo][k] * elem[k] * elem[k] * elem[iMo] * (1 / (1 - elem[iMo]) + 1 / (1 - elem[k]) + elem[iMo] / (2 * (1 - elem[iMo]) * (1 - elem[iMo])) - 1);
  }

  double ln_gamma_Mo = t1_Mo + t2_Mo + t3_Mo;
  //cout << "ln gamma(Mo) = " << setprecision(8) << ln_gamma_Mo << tab << "gamma(Mo) = " << exp(ln_gamma_Mo) << endl;


// Fit for K an Kd of Mo using fits from Jennings note
  double a_Mo = 1.47;
  double b_Mo = -1448;
  double c_Mo = -67.1;
  double K_Mo = pow(10.0, a_Mo + b_Mo / T + c_Mo * P / T);
  double Kd_Mo = exp(log(K_Mo) - ln_gamma_Mo + 2.0 * ln_gamma_fe);
  double D_Mo = Kd_Mo * pow(XFe, 2.0) / pow(XFeO, 2.0);
  double total_moles_Mo = moles[17] + moles[n + 13];
  double moles_Mo_met = (total_moles_Mo * D_Mo * moles_met / moles_sil) / (1.0 + D_Mo * moles_met / moles_sil);
  double moles_Mo_sil = total_moles_Mo - moles_Mo_met;
  moles[17] = moles_Mo_sil;
  moles[n + 13] = moles_Mo_met;

  //testfile << D_Mo << tab << log(XC)

  return;
}

//_______________________________________________________________________________________________________

// Function that calculates Ds and the distribution of S between metal and silicate using Boujibar et al (2014) (which has some problems, especially with effect of oxygen)
double S_partitioning_blanchard(double moles[], double moles_oxygen, double moles_S_met, double moles_S_sil, double moles_C_met, double moles_C_sil, double P, double T)
{
// For reference, below are the indexes for elements in elements.h
// Element index: 0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22
// Element name: Al Mg Ca Na Fe Si Ni Co Nb Ta  V Cr Pt Pd Ru Ir  W Mo Mn Ga Zn In Cu
// Element index:             n +1 +2 +3 +4 +5 +6 +7 +8 +9+10+11+12+13+14+15+16+17+18
// Metals start at Fe at index n then continue to n+m where m = n - 4
// If n = 17, then only up to W is included
//  MOLAR COMPOSITIONS NEEDED TO BE DEFINED ARE:
//            moles[x]          Number of moles of the xth oxide and metal species where x varies from 0 to n+m-1
//            moles_H2O         Number of moles of water
//            moles_S_sil       Number of moles of sulfur as oxide
//            moles_C_sil       Number of moles of carbon as oxide
//            moles_oxygen      Number of moles of oxygen as metal
//            moles_H           Number of moles of hydrogen as metal
//            moles_S_met       Number of moles of sulfur as metal
//            moles_C_met       Number of moles of carbon as metal
//__________________________________________________________________________________________________________

  int i = 0;
  double X[n + m] = {0.0};      // mole fractions
  double total_moles_sil = 0.0;
  double total_moles_met = 0.0;
  double total_moles_S = moles_S_sil + moles_S_met;

// First calculate mole fractions of Al, Mg, Ca, Fe, Na, H, and Si in silicate melt, as required for calculating sulfide capacity
// (According to Smythe et al., 2017, Am Min, 102,795-803)

  for (i = 0; i < n; i++)
    total_moles_sil += moles[i];
  for (i = 0; i < n; i++)
    X[i] = moles[i] / total_moles_sil;

  for (i = n; i < n + m; i++)
    total_moles_met += moles[i];
  total_moles_met = total_moles_met + moles_oxygen;
  for (i = n; i < n + m; i++)
    X[i] = moles[i] / total_moles_met;
  double X_oxygen = moles_oxygen / total_moles_met;

// Constants from Smythe et al (Table 2) for defining Cs:
  double A[n] = {0.0};
  A[0] = 9295.0;
  A[1] = 13767.0;
  A[2] = 19893.0;
  A[3] = 14197.0;
  A[4] = -7080.0;
  double A_SiFe = 117827.0;
  double A0 = 8.03;


  /* The following are test values from Ingrid that overwrite the actual values
    X[0] = 0.03;    // Al in silicate
    X[1] = 0.34;    // Mg in silicate
    X[2] = 0.03;    // Ca in silicate Al in silicate
    x[3] = 0.02;    // Na in silicate
    X[4] = 0.27;    // Fe in silicate
    X[5] = 0.32;    // Si in silicate
    X_oxygen = 0.15;  // O in metal
    X[n-1] = 0.08;
    P = 30.0; T = 2962.0;

    cout << endl << "P = " << P << "   T = " << T << "  XFeO = "<< X[4] << endl;
  */

  double lnCs = 0.0;    // Natural log of sulfide calacity
  for (i = 0; i < 5; i++)
    lnCs += X[i] * A[i] / T;
  lnCs = A0 + lnCs + A_SiFe * X[4] * X[5] / T;

  double Cs = exp(lnCs);  // Sulfide capacity according to Smyth et al

// Metal-silicate partition coefficient Ds of sulfur based on Ingrid's fit:

  const double a = 6.0, b = 6037.0, c = -87.0, d = -17.0;
  double log_Ds = 0.0, Ds = 0.0;
  log_Ds = a + b / T + c * P / T + log10(X[4]) - log10(Cs) + d * (1873.0 / T) * log10(1.0 - X_oxygen) + 9.3 * (1873.0 / T) * log10(1.0 - X[n + 1]);
  Ds = pow(10.0, log_Ds);

//  cout << "Results of Ingrid's model:" << endl;
//  cout << "   ln Cs = " << lnCs << "   log Ds = " << log_Ds << endl;

//  Following is old code based on Boujibar et al from version v-29
  /*
     double moles_sil = 0.0;
     double moles_met = 0.0;


  // Calculate wt% FeO in silicate melt:
     double sum_sil = 0;
     double weight[n] = { 0.0 };
     for (i = 0; i < n; i++)
     {
        weight[i] = moles[i] * MW[i];
        sum_sil += weight[i];
     }

     double weight_S_sil = moles_S_sil * AW_S;
     sum_sil += weight_S_sil;
     double weight_C_sil = moles_C_sil * AW_C;
     sum_sil += weight_C_sil;
     for (i = 0; i < n; i++)  weight[i] = weight[i] * 100.0 / sum_sil;
     weight_S_sil = weight_S_sil * 100.0 / sum_sil;
     double wtFeO = weight[4];
    cout << "wt % FeO = " << wtFeO << endl;

   /*
     for (i = 0; i < n; i++)  moles_sil += moles[i];
     moles_sil += moles_S_sil;
     moles_sil += moles_C_sil;

     for (i = n; i < (n + m); i++)  moles_met += moles[i];
     moles_met += moles_oxygen;        // oxygen
     moles_met += moles_S_met;
     moles_met += moles_C_met;

     double Xfeo = 0.0;        // silicate mole fractions...
     double Xcao = 0.0;
     double Xmgo = 0.0;
     double Xfe = 0.0;       // metal mole fractions...
     double Xsi = 0.0;
     double Xni = 0.0;

  // Calculate sulfide capacity according to Boujibar et al
     double log_Cs = 0.0;
     double Xfeo = X[4];
     double Xcao = X[2];
     double Xmgo = X[1];
     log_Cs = -5.704 + 3.15 * Xfeo + 2.65 * Xcao + 0.12 * Xmgo;

  // Calculate log Ds(met-sil)
     const double B = 405.0, C = 136.0, D = 32.0, E = 181.0, F = 305.0, H = 1.13, I = 10.7, K = -3.72;
  //   double log_Ds = 0.0;
     double Xfe = moles[n];
     double Xsi = X[n + 1];
     double Xni = X[n + 2];

     log_Ds = log10(wtFeO) - log_Cs + B / T + C * P / T + D * log10(1.0 - Xsi) + E * pow(log10(1.0 - Xsi), 2.0) + F * pow(log10(1.0 - Xsi), 3.0)
              + H * log10(1.0 - Xfe) + I * log10(1.0 - Xni) + K;

     Ds = pow(10.0, log_Ds);

     cout << "Results of Boujibar's model:" << endl;
    cout << "   log Cs = " << log_Cs << "   log Ds = " << log_Ds << endl;
      int ggg; cin >> ggg;
  */

  /*  if (t==3.57000E+03 || t== 2.11842E+07 || t==1.12998E+08)
    {
      cout << endl << "P = " << P << "  T = " << T << "  Ds = " << Ds << endl;
      int hhhh; cin >> hhhh;
    }
  */
//  if(no_1==12 && no_2==27) cout << "Ds = " << Ds << endl;
//   if (P>5.0) cout << "P=" << P << tab << "Ds=" << Ds << endl;



  moles_S_met = (total_moles_S * Ds * total_moles_met / total_moles_sil) / (1.0 + Ds * total_moles_met / total_moles_sil); // This seems OK

//   double moles_S_met_2 = (Ds * total_moles_S / total_moles_sil) / (1.0 / total_moles_met + Ds / total_moles_sil);

//  if (t==1.12998E+08) cout << "1: " << moles_S_met << " 2: " << moles_S_met_2 << endl;
//  if (t==1.12998E+08) cout << "Direct D in S partitioning = " << (moles_S_met/moles_met)/((total_moles_S-moles_S_met)/moles_sil) << endl;

  return moles_S_met;
}

// Function that calculates Ds and the distribution of S between metal and silicate using Boujibar et al (2014) (which has some problems, especially with effect of oxygen)
double S_partitioning_Boujibar(double moles[], double moles_oxygen, double moles_S_met, double moles_S_sil, double moles_C_met, double moles_C_sil, double P, double T)
{
// For reference, below are the indexes for elements in elements.h
// Element index: 0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22
// Element name: Al Mg Ca Na Fe Si Ni Co Nb Ta  V Cr Pt Pd Ru Ir  W Mo Mn Ga Zn In Cu
// Element index:             n +1 +2 +3 +4 +5 +6 +7 +8 +9+10+11+12+13+14+15+16+17+18
// Metals start at Fe at index n then continue to n+m where m = n - 4
// If n = 17, then only up to W is included
//  MOLAR COMPOSITIONS NEEDED TO BE DEFINED ARE:
//            moles[x]          Number of moles of the xth oxide and metal species where x varies from 0 to n+m-1
//            moles_H2O         Number of moles of water
//            moles_S_sil       Number of moles of sulfur as oxide
//            moles_C_sil       Number of moles of carbon as oxide
//            moles_oxygen      Number of moles of oxygen as metal
//            moles_H           Number of moles of hydrogen as metal
//            moles_S_met       Number of moles of sulfur as metal
//            moles_C_met       Number of moles of carbon as metal
//__________________________________________________________________________________________________________
  int i = 0;
  double moles_sil = 0.0;
  double moles_met = 0.0;
  double old_moles_S_met = 0.0;
  double total_moles_S = moles_S_sil + moles_S_met;
  int itt = 0;

  do {
    itt++;
// Calculate wt% FeO in silicate melt:
    double sum_sil = 0;
    double weight[n] = { 0.0 };
    for (i = 0; i < n; i++)
    {
      weight[i] = moles[i] * MW[i];
      sum_sil += weight[i];
    }

    double weight_S_sil = moles_S_sil * AW_S;
    sum_sil += weight_S_sil;
    double weight_C_sil = moles_C_sil * AW_C;
    sum_sil += weight_C_sil;
    for (i = 0; i < n; i++)  weight[i] = weight[i] * 100.0 / sum_sil;
    weight_S_sil = weight_S_sil * 100.0 / sum_sil;
    double wtFeO = weight[4];
// cout << "wt % FeO = " << wtFeO << endl;

    for (i = 0; i < n; i++)  moles_sil += moles[i];
    moles_sil += moles_S_sil;
    moles_sil += moles_C_sil;

    for (i = n; i < (n + m); i++)  moles_met += moles[i];
    moles_met += moles_oxygen;        // oxygen
    moles_met += moles_S_met;
    moles_met += moles_C_met;

    double Xfeo = 0.0;        // silicate mole fractions...
    double Xcao = 0.0;
    double Xmgo = 0.0;
    double Xfe = 0.0;       // metal mole fractions...
    double Xsi = 0.0;
    double Xni = 0.0;

// Calculate sulfide capacity
    double log_Cs = 0.0;
    Xfeo = moles[4] / moles_sil;
    Xcao = moles[2] / moles_sil;
    Xmgo = moles[1] / moles_sil;
    log_Cs = -5.704 + 3.15 * Xfeo + 2.65 * Xcao + 0.12 * Xmgo;

// Calculate log Ds(met-sil)
    const double B = 405.0, C = 136.0, D = 32.0, E = 181.0, F = 305.0, H = 1.13, I = 10.7, K = -3.72;
    double log_Ds = 0.0;
    Xfe = moles[n] / moles_met;
    Xsi = moles[n + 1] / moles_met;
    Xni = moles[n + 2] / moles_met;

    log_Ds = log10(wtFeO) - log_Cs + B / T + C * P / T + D * log10(1.0 - Xsi) + E * pow(log10(1.0 - Xsi), 2.0) + F * pow(log10(1.0 - Xsi), 3.0)
             + H * log10(1.0 - Xfe) + I * log10(1.0 - Xni) + K;

    double Ds = pow(10.0, log_Ds);

    old_moles_S_met = moles_S_met;
    moles_S_met = (total_moles_S * Ds * moles_met / moles_sil) / (1.0 + Ds * moles_met / moles_sil); // This seems OK

//    double moles_S_met_2 = (Ds * total_moles_S / moles_sil) / (1.0 / moles_met + Ds / moles_sil);

    // testfile << itt << scientific << tab << moles_S_met << tab << (old_moles_S_met - moles_S_met) / total_moles_S  << endl;
  }
  while (fabs(old_moles_S_met - moles_S_met) > 1.0e-6);

  return moles_S_met;
}
//_______________________________________________________________________________________________________

// Function that calculates Dc and the distribution of C between metal and silicate using Hirschmann et al (2016)
double C_partitioning(double moles[], double moles_oxygen, double moles_S_met, double moles_S_sil, double moles_C_met, double moles_C_sil, double P, double T)
{
  int i = 0;
  double moles_sil = 0.0;
  double moles_met = 0.0;
  double total_moles_C = moles_C_sil + moles_C_met;

// Calculate wt% FeO in silicate melt:
  for (i = 0; i < n; i++) moles_sil += moles[i];
  moles_sil += moles_S_sil;
  moles_sil += moles_C_sil;

  for (i = n; i < (n + m); i++)  moles_met += moles[i];
  moles_met += moles_oxygen;        // oxygen
  moles_met += moles_S_met;
  moles_met += moles_C_met;

  moles_C_met = (total_moles_C * Dc * moles_met / moles_sil) / (1.0 + Dc * moles_met / moles_sil); // This seems OK

  return moles_C_met;
}
//_______________________________________________________________________________________________________

// Function that calculates HSE partitioning between FeS and silicate liquid
void FeS_HSE(double h_moles[], double h_moles_H2O, double h_moles_S, double h_moles_FeS, double (& HSE_moles)[2 * n_hse], double P)
{
  int i = 0;
  double P_FeS = P * matte_P_factor; // effective pressure of FeS-silicate equilibration
  double HSE_met[n_hse] = {0.0};
  double HSE_sil[n_hse] = {0.0};

  double sulf_logKd = 0.0;
  double sig_sulf_logKd = 0.0;

  // Determine equilibrated HSE concentrations (immediately prior to sulphide segregation)
  double tot_sil = 0.0;
  double tot_met = 0.0;

  for (int i = 0; i < n; i++) tot_sil += h_moles[i]; //total moles of silicate
  tot_sil += h_moles_H2O;
  tot_sil += h_moles_S;
  tot_met = h_moles_FeS;     // total moles of FeS in mantle

  int j = 0;              // j counts the HSEs and runs from 0 to <n_hse (currently 4)
  for (int i = 0; i < n; i++)
  {
    if (el_type[i] == 'h' || el_type[i] == 'H')
    {
      // log Kd from Vera's FeS-silicate liquid fits

      sig_sulf_logKd = sqrt(pow(sulf_sig_a[i], 2.0) + pow (sulf_sig_b[i] / temperature(P_FeS), 2.0) + pow(sulf_sig_c[i] * P_FeS / temperature(P_FeS), 2.0));

      if (error == 0)
        sulf_logKd = sulf_a[i] + sulf_b[i] / temperature(P_FeS) + sulf_c[i] * P_FeS / temperature(P_FeS);

      else if (error == 1)
        sulf_logKd = (sulf_a[i] + sulf_b[i] / temperature(P_FeS) + sulf_c[i] * P_FeS / temperature(P_FeS)) + sig_sulf_logKd;


      else if (error == -1)
        sulf_logKd = (sulf_a[i] + sulf_b[i] / temperature(P_FeS) + sulf_c[i] * P_FeS / temperature(P_FeS)) - sig_sulf_logKd;

      // Eq.3.1 of trace element calculation:
      double te_D = exp(sulf_logKd * 2.3026) * pow((0.5 / (h_moles[4] / tot_sil)), (val[i] / 2.0) );
      // Eq. 3.3:
      double cc = te_D * (h_moles[i]);
      double dd = te_D + (tot_sil / tot_met);
      HSE_met[j] = cc / dd;
      HSE_sil[j] = h_moles[i] - HSE_met[j];
      j++;
    }
  }

  for (int i = 0; i < n_hse; i++) HSE_moles[i] = HSE_met[i];
  for (int i = n_hse; i < 2 * n_hse; i++) HSE_moles[i] = HSE_sil[i - n_hse];

  return;          // EQUILIBRATED HSE CONCENTRATIONS - first FeS and second silicate concentrations
}
//_____________________________________________________________________________________________________

// This function determines ppm S concentration in the mantle and the degree of mantle supersaturation of S
double S_conc(double moles[], double moles_H2O, double moles_S_mantle, double moles_C_mantle, double P)
{
  double S_ppm(0);
  double S_frac_matte(0);
  // cout << "moles_S_mantle = " << moles_S_mantle << endl;
  // cout << "moles_Fe = " << moles[4] << endl;
  int i = 0;
  double sum_sil = 0;
  for (int i = 0; i < n; i++) sum_sil += moles[i] * MW[i];
  sum_sil += moles_H2O * MW_H2O;
  sum_sil += moles_S_mantle * AW_S;
  sum_sil += moles_C_mantle * AW_C;
// S concentration prior to Hadean matte segregation
  S_ppm = moles_S_mantle * AW_S * 1e6 / sum_sil;
  double P_FeS = P * matte_P_factor;

  // S saturation from Laurenz et al (2016)
  double S_sat = exp(14.2 - 11032.0 / temperature(P_FeS) - 379.0 * P_FeS / temperature(P_FeS));

  // S saturation from Mavrogenes & O'Neill (1999):
  // double S_sat = exp(8.29 - 0.019*P*1e4/temperature(P));
  // Old model without Antje's data:
  //double S_sat = exp(15.6455 - 13616.14/temperature(P_FeS) - 441.413 *P_FeS/temperature(P_FeS));

  return S_frac_matte = S_ppm / S_sat;
}

//____________________________________________________________________________________

// Function that calculates object mass based on number of accumulated moles (from the mass balance)
double mass_f(double moles[], double moles_H2O, double moles_O, double moles_H, double moles_S_mantle, double moles_S_core, double moles_C_mantle, double moles_C_core, const double AW[], const double MW[])
{
  double object_mass(0);
  int i = 0;

  double sum_sil = 0;
  for (int i = 0; i < n; i++)
    sum_sil += moles[i] * MW[i];
  sum_sil += moles_H2O * MW_H2O;
  sum_sil += moles_S_mantle * AW_S;
  sum_sil += moles_C_mantle * AW_C;

  double sum_met = 0;
  for (int i = n; i < (n + m); i++)
    sum_met += moles[i] * AW[i - m];
  sum_met += moles_O * AW_O;
  sum_met += moles_S_core * AW_S;
  sum_met += moles_C_core * AW_C;
  sum_met += moles_H * AW_H;
  object_mass = (sum_sil + sum_met);
  return object_mass;
}

//______________________________________________________________________________________________________________________________________________
void write_outfile_header(ofstream & outfile, string simulation_no, int dist_comp, double p_factor, double em_equil, double pf_em, double dist_1,
                          double dist_2, double dist_3, double erf_width_1, double erf_width_2, double xfe_met_1, double xfe_met_2, double xfe_met_3, double xsi_met_0, double wtci_H2O,
                          double mass_adjust, double xmet_pl, double xmet_em, double x_mantle, double em_xsil, double t_overstep, const char el_type[],
                          char el_sym[][2], double dist_0[], double dist_0_S, double dist_0_C, double matte_HSE, double matte_P_factor, double pf,
                          double moon_xfe_met_1, double moon_xfe_met_2, double moon_dist_1, double moon_dist_2, double moon_disk_start, double moon_disk_end) // Moon_version 
{
// WRITE HEADER INFORMATION FOR THE VARIOUS OUPUT FILES
  outfile << endl << "__________________________________________________________________" << endl;
  outfile << "N-BODY SIMULATION No: " << simulation_no << endl;
  if (dist_comp == 1)
    outfile << "<CONSTANT> DISTANCE-COMPOSITION MODEL" << endl;
  else if (dist_comp == 2)
    outfile << "<STEP> DISTANCE-COMPOSITION MODEL" << endl;
  else if (dist_comp == 3)
    outfile << "<GRADIENT> DISTANCE-COMPOSITION MODEL" << endl;
  else if (dist_comp == 4)
    outfile << "<ERF> DISTANCE-COMPOSITION MODEL" << endl;
  else if (dist_comp == 5)
    outfile << "<MOON> DISTANCE-COMPOSITION MODEL" << endl;
  if (em_equil == 0)
  {
    outfile << "PRESSURE FACTOR IS THE SAME FOR EMBRYOS AND PLANETESIMALS" << endl;
    outfile << setprecision(3);
    outfile << "Met-Sil equilibration pressure = " << p_factor << " x P(CMB)" << endl;
  }
  else if (em_equil == 1)
  {
    outfile << "PRESSURE FACTOR IS DIFFERENT FOR EMBRYOS AND PLANETESIMALS" << endl;
    outfile << setprecision(3);
    outfile << "For embryo collisions, Met-Sil equilibration pressure = " << pf_em << " x P(CMB) (fixed parameter)" << endl;
    outfile << "For planetesimal impacts, Met-Sil equilibration pressure = " << p_factor << " x P(CMB) (refined parametr)" << endl;
  }
  if (dist_comp == 1)  // CONSTANT model
    outfile << "Value of  xfe_met in H2O-free compositions: " << xfe_met_2 << endl;
  else if (dist_comp == 2) // STEP Model
  {
    outfile << "First distance value (dist_1): " << dist_1 << " AU" << endl;
    outfile << "Second distance value (dist_2): " << dist_2 << " AU" << endl;
    outfile << "Third distance value (dist_3): " << dist_3 << " AU" << endl;
    outfile << "Value of xfe_met from 0 AU to dist_1: " << xfe_met_1 << endl;
    outfile << "Fraction of Si initially in metal from 0 AU to dist_1: " << xsi_met_0 << endl;
    outfile << "Value of  xfe_met from dist_1 to dist_2: " << xfe_met_2 << endl;
    outfile << "Value of  xfe_met from dist_2 to dist_3: " << xfe_met_3 << endl;
  }
  else if (dist_comp == 3) // GRADIENT model
  {
    outfile << "Distance beyond which XSi(met)=zero (D1): " << dist_1 << " AU" << endl;
    outfile << "Distance beyond which XFe(met)is close to zero (D2): " << dist_2 << " AU" << endl;
    outfile << "Fraction of Si initially in metal in reduced composition (X0): " << xsi_met_0 << endl;
    outfile << "Value of xfe_met in reduced composition (X1) " << xfe_met_1 << endl;
    outfile << "Value of xfe_met at D2 (X2): " << xfe_met_2 << endl;
    outfile << "Maximum value of xfe_met in oxidised, H2O-free compositions (X3): " << xfe_met_3 << endl;
  }
  else if (dist_comp == 4) // ERROR FUNCTION model
  {
    outfile << "Distance at half of maximum value for Si metal: " << dist_1 << " AU" << endl;
    outfile << "Width of transition region between max and min value: " << erf_width_1 << " AU" << endl;
    outfile << "Distance at half of maximum value for Fe metal:: " << dist_2 << " AU" << endl;
    outfile << "Width of transition region between max and min value: " << erf_width_2 << " AU" << endl;
    outfile << "Fraction of Si initially in metal in reduced composition (X0): " << xsi_met_0 << endl;
    outfile << "Value of xfe_met in reduced composition (X1) " << xfe_met_1 << endl;
    outfile << "Maximum value of xfe_met in oxidised, H2O-free compositions (X3): " << xfe_met_3 << endl;
  }
  else if (dist_comp == 5) // MOON GRADIENT model
  {
    outfile << "Value of moon_xfe_met_1: " << moon_xfe_met_1 << endl;
    outfile << "Value of moon_xfe_met_2: " << moon_xfe_met_2 << endl;
    outfile << "moon_dist_1: " << moon_dist_1 << " (Earth radii)" << endl;
    outfile << "moon_dist_2: " << moon_dist_2 << " (Earth radii)" << endl;
    outfile << "moon_disk_start: " << moon_disk_start << " (Earth radii)" << endl;
    outfile << "moon_disk_end: " << moon_disk_end << " (Earth radii)" << endl;

  }
  outfile << "Distance beyond which compositions are H2O-bearing (D3): " << dist_3 << " AU" << endl;
  outfile << "H2O concentration in H2O-bearing composition: " << wtci_H2O << " wt%" << endl;
  outfile << "Mass adjustment factor: " << mass_adjust << endl;
  outfile << "Fraction of impacting planetesimal core that equilibrates: " << xmet_pl << endl;
  outfile << "Fraction of impacting embryo core that equilibrates: " << xmet_em << endl;
  if (x_mantle == 0.0)
    outfile << "Factor by which xsil is changed relative to Deguen model for EMBRYO impacts: " << em_xsil << endl;
  else if (x_mantle > 0.0)
    outfile << "Fraction of mantle that equilibrates with metal for all impacts: " << x_mantle << endl;
  outfile << "T overstep (relative to peridotite liquidus) = " << setprecision(3) << t_overstep << " K" << endl;
  outfile << "VOLATILE ELEMENT PARAMETERS - inner distances that define linear concentration gradients from D3: " << endl;
  for (int i = 0; i < n; i++)
  {
    if (el_type[i] == 'v' || el_type[i] == 'H')
    {
      for (int j = 0; j < 2; j++)
        outfile << el_sym[i][j];
      outfile << ": " << dist_0[i] << " AU;  ";
    }
  }
  outfile << "S: " << dist_0_S << " AU" << tab;
  outfile << "C: " << dist_0_C << " AU" << endl;
  outfile << "Dc: " << Dc << endl;
  if (matte_HSE == 0) outfile << "NO FeS(+ HSE) SEGREGATION" << endl;
  else if (matte_HSE == 1)
  {
    outfile << "INCLUDES FeS(+ HSE) SEGREGATION" << endl;
    outfile << "matte_P_factor = " << matte_P_factor << "   ks = matte_P_factor*pressure factor(pf) = " << matte_P_factor*pf << endl;
  }
  return;
}

//______________________________________________________________________________________________________________________________________________
double Fe_metalfraction(int dist_comp, double dist, double dist_1, double dist_2, double dist_3, double xfe_met_1,  double xfe_met_2, double xfe_met_3, double erf_width_1, double erf_width_2,
                        double moon_xfe_met_1, double moon_xfe_met_2, double moon_dist_1, double moon_dist_2, double moon_disk_start, double moon_disk_end) // Moon_version 
{
//DETERMINE STARTING FE METAL FRACTION BASED ON HELIOCENTRIC DISTANCE (AU)
  double xfe_met = 0.0;
//__________________________________
// Constant composition-distance model
  if (dist_comp == 1)
  {
// Within dist_3, the fraction of Fe in the metal is constant at xfe_met_2
    if (dist < dist_3) xfe_met = xfe_met_2;
// Exterior to dist_3, there is no Fe in the metal
    else if (dist >= dist_3) xfe_met = 0.0;
  }
//__________________________________
// Step composition-distance model
  else if (dist_comp == 2)
  {
// Interior to dist_1, the fraction of Fe in the metal is constant at xfe_met_1
    if (dist < dist_1) xfe_met = xfe_met_1;
// Between dist_1 and dist_2, the fraction of Fe in the metal is constant at xfe_met_2
    else if ((dist >= dist_1) && (dist < dist_2)) xfe_met = xfe_met_2;
// Between dist_2 and dist_3, the fraction of Fe in the metal is constant at xfe_met_3
    else if ((dist >= dist_2) && (dist < dist_3)) xfe_met = xfe_met_3;
// Exterior to dist_3, there is no Fe in the metal
    else if (dist >= dist_3) xfe_met = 0.0;
  }
//__________________________________
// Gradient composition-distance model
  else if (dist_comp == 3)
  {
// Interior to dist_1, the fraction of Fe in the metal is constant at xfe_met_1
    if (dist < dist_1) xfe_met = xfe_met_1;
// Between dist_1 and dist_2, the fraction of Fe in the metal is a slope between xfe_met_1 and xfe_met_2
    else if ((dist >= dist_1) && (dist < dist_2))
    {
      xfe_met = ((xfe_met_1 - xfe_met_2) / (dist_1 - dist_2)) * (dist - dist_1) + xfe_met_1;
    }
// Between dist_2 and dist_3, the fraction of Fe in the metal is a slope between xfe_met_2 and 0.0
    else if ((dist >= dist_2) && (dist < dist_3))
    {
      xfe_met = ((xfe_met_2) / (dist_2 - dist_3)) * (dist - dist_2) + xfe_met_2;
    }
// Exterior to dist_3, there is no Fe in the metal
    else if (dist >= dist_3) xfe_met = 0.0;
  }
//__________________________________
// Error function model
  else if (dist_comp == 4)
  {
    if (dist < dist_3)
    {
      xfe_met = xfe_met_1 * (0.5 - 0.5 * erf(2.0 * (dist - dist_2) / erf_width_2));
    }
    else
    {
      xfe_met = 0.0;
    }
  }
//__________________________________
// Moon ver
// MOON Gradient composition-distance model
  else if (dist_comp == 5)
  {
// Interior to moon_dist_1, the fraction of Fe in the metal is constant at moon_xfe_met_1
    if (dist < moon_dist_1) xfe_met = moon_xfe_met_1;
// Between moondist_1 and moondist_2, the fraction of Fe in the metal is a slope between moon_xfe_met_1 and moon_xfe_met_2
    else if ((dist >= moon_dist_1) && (dist < moon_dist_2))
    {
      xfe_met = ((moon_xfe_met_1 - moon_xfe_met_2) / (moon_dist_1 - moon_dist_2)) * (dist - moon_dist_1) + moon_xfe_met_1;
    }
// Exterior to moon_dist_2, fraction of Fe in metl is constant at moon_xfe_met_2
    else if (dist >= moon_dist_2) xfe_met = moon_xfe_met_2;
  }
  return xfe_met;
}

//______________________________________________________________________________________________________________________________________________
double Si_metalfraction(int dist_comp, double dist, double dist_1, double dist_2, double dist_3, double xsi_met_0, double erf_width_1, double erf_width_2)
{
  double xsi_met = 0.0;
//DETERMINE STARTING SI METAL FRACTION BASED ON HELIOCENTRIC DISTANCE (AU)
//__________________________________
// Constant composition-distance model
  if (dist_comp == 1)
  {
// No Si ever starts as a metal
    xsi_met = 0.0;
  }
//__________________________________
// Step composition-distance model
  else if (dist_comp == 2)
  {
// Interior to dist_1, the fraction of Si in the metal is xsi_met_0
    if (dist < dist_1) xsi_met = xsi_met_0;
// Exterior to dist_1, there is no Si in the metal
    else if (dist >= dist_1) xsi_met = 0.0;
  }
//__________________________________
// Gradient composition-distance model
  else if (dist_comp == 3)
  {
// Interior to dist_1, the fraction of Si in the metal is xsi_met_0
    if (dist < dist_1) xsi_met = xsi_met_0;
// Exterior to dist_1, there is no Si in the metal
    else if (dist >= dist_1) xsi_met = 0.0;
  }
//__________________________________
// Error function model
  else if (dist_comp == 4)
  {
    if (dist < dist_3)
    {
      xsi_met = xsi_met_0 * (0.5 - 0.5 * erf(2.0 * (dist - dist_1) / erf_width_1));
    }
    else
    {
      xsi_met = 0.0;
    }
  }
//_________________________________
// Moon ver.
// MOON Gradient model
  else if (dist_comp == 5)
  {
    xsi_met = 0.05;
  }
  return xsi_met;
}

//______________________________________________________________________________________________________________________________________________
void starting_molar_proportions(double dist,
                                double (&moles)[nm], double & moles_H2O, double & moles_S_sil, double & moles_C_sil,
                                double & moles_oxygen, double & moles_H, double & moles_S_met, double & moles_C_met)
{
//DETERMINE STARTING MOLAR COMPOSITIONS FROM CI composition defined in INPUT_*.txt, elements.h and heliocentric distance
// For reference, below are the indexes for elements in elements.h
// Element index: 0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22
// Element name: Al Mg Ca Na Fe Si Ni Co Nb Ta  V Cr Pt Pd Ru Ir  W Mo Mn Ga Zn In Cu
// Element index:            17 18 19 20 21 22 23 24 25 26 27 28 29
// Metals start at Fe at index n then continue to n+m where m = n - 4
// If n = 17, then only up to W is included
//  MOLAR COMPOSITIONS NEEDED TO BE DEFINED ARE:
//            moles[x]          Number of moles of the xth oxide and metal species where x varies from 0 to n+m-1
//            moles_H2O         Number of moles of water
//            moles_S_sil       Number of moles of sulfur as oxide
//            moles_C_sil       Number of moles of carbon as oxide
//            moles_oxygen      Number of moles of oxygen as metal
//            moles_H           Number of moles of hydrogen as metal
//            moles_S_met       Number of moles of sulfur as metal
//            moles_C_met       Number of moles of carbon as metal
//__________________________________________________________________________________________________________
// CALCULATE STARTING MOLES OF ELEMENTS AS OXIDES IN SILICATE (formulated as XO(valence/2), except for H2O)
  for (int i = 0; i < n; i++)
  {
//______________________________________________________________________
// Volatile elements as oxides, e.g. Na, but not S
// THESE FOLLOW A LINEAR GRADIENT BETWEEN 0 AND CI COMPOSITION ACCORDING TO HELIOCENTRIC DISTANCE
    if (el_type[i] == 'v' || el_type[i] == 'H')
    {
      if (dist <= dist_0[i]) moles[i] = 0.0;
      else if (dist >= dist_3) moles[i] = wtci[i] / AW[i];
      else moles[i] = (wtci[i] * (dist - dist_0[i]) / (dist_3 - dist_0[i])) / AW[i];
    }
// Si as an oxide
    else if (i == 5) moles[i] = wtci[i] * en_dep[i] * (1.0 - xsi_met) / AW[i];
// Lithophile elements as oxides; these are entirely oxides
    else if (i < 4) moles[i] = wtci[i] * en_dep[i] / AW[i];
// Fe as an oxide
    else if (i == 4)
    {
      moles[i] = wtci[i] * en_dep[i] * (1.0 - xfe_met) / AW[i];
    } 
// Siderophile elements as oxides; assume that 1% of their total abundance is initially as an oxide
    else moles[i] = wtci[i] * en_dep[i] * (1.0 - 0.99) / AW[i];
  }

//______________________________________________________________________
// H2O only beyond dist_3
  if (dist >= dist_3) moles_H2O = wtci_H2O / MW_H2O;
  else moles_H2O = 0.0;

//______________________________________________________________________
// NO SULFUR STARTS AS OXIDE
  moles_S_sil = 0.0;

//______________________________________________________________________
// NO CARBON STARTS AS OXIDE
  moles_C_sil = 0.0;

//__________________________________________________________________________________________________________
// CALCULATE STARTING MOLES OF ELEMENTS IN METAL
  for (int i = n; i < (n + m); i++)
  {
// Fe as a metal
    if (i == n) 
    {
      moles[i] = wtci[i - m] * en_dep[i - m] * xfe_met / AW[i - m];
    }
// Si as a metal
    else if (i == n + 1) moles[i] = wtci[i - m] * en_dep[i - m] * xsi_met / AW[i - m];
// All other siderophile elements; assumes that 1% of their total abundance is initially as a metal
    else moles[i] = wtci[i - m] * en_dep[i - m] * 0.99 / AW[i - m];
  }

//______________________________________________________________________
// NO O STARTS AS A METAL
  moles_oxygen = 0.0;

//______________________________________________________________________
// NO H STARTS AS A METAL
  moles_H = 0.0;

//______________________________________________________________________
// ALL SULFUR IS A METAL AND FOLLOWS A LINEAR GRADIENT BETWEEN 0 AND CI COMPOSITION ACCORDING TO HELIOCENTRIC DISTANCE
  if (dist <= dist_0_S) moles_S_met = 0.0;
  else if (dist >= dist_3) moles_S_met = wtci_S / AW_S;
  else moles_S_met = (wtci_S * (dist - dist_0_S) / (dist_3 - dist_0_S)) / AW_S;

//______________________________________________________________________
// ALL CARBON IS A METAL AND FOLLOWS A LINEAR GRADIENT BETWEEN 0 AND CI COMPOSITION ACCORDING TO HELIOCENTRIC DISTANCE
  if (dist <= dist_0_C) moles_C_met = 0.0;
  else if (dist >= dist_3) moles_C_met = wtci_C / AW_C;
  else moles_C_met = (wtci_C * (dist - dist_0_C) / (dist_3 - dist_0_C)) / AW_C;

  return;
}

//______________________________________________________________________________________________________________________________________________
double get_mass_factor(double moles[], double moles_H2O, double moles_S_sil, double moles_C_sil,
                       double moles_oxygen, double moles_H, double moles_S_met, double moles_C_met)
{
// Calculate "mass factor", which is the relation between the total moles and a unit mass
// The code assumes that this unit mass is an Earth mass
  double mass_factor;
  double sum_sil = 0.0;
  double sum_met = 0.0;

// Sum the total moles of oxides
  for (int i = 0; i < n; i++) sum_sil += moles[i] * MW[i];
  sum_sil += moles_H2O * MW_H2O;
  sum_sil += moles_S_sil * MW_S2O;
  sum_sil += moles_C_sil * MW_CO;

// Sum the total moles of metals
  for (int i = n; i < (n + m); i++) sum_met += moles[i] * AW[i - m];
  sum_met += moles_oxygen * AW_O;
  sum_met += moles_H * AW_H;
  sum_met += moles_S_met * AW_S;
  sum_met += moles_C_met * AW_C;
// Calculate mass factor
  mass_factor = 1.0 / (sum_sil + sum_met);

  return mass_factor;
}

//______________________________________________________________________________________________________________________________________________
void undifferentiated_composition(double mass, double mass_factor,
                                  double (&moles)[nm], double & moles_H2O, double & moles_S_sil, double & moles_C_sil,
                                  double & moles_oxygen, double & moles_H, double & moles_S_met, double & moles_C_met)
{
// GET THE COMPOSITION OF AN UNDIFFERENTIATED, METAL-FREE BODY
// For reference, below are the indexes for elements in elements.h
// Element index: 0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22
// Element name: Al Mg Ca Na Fe Si Ni Co Nb Ta  V Cr Pt Pd Ru Ir  W Mo Mn Ga Zn In Cu
// Element index:            18 19 20 21 22 23 24 25 26 27 28 29 30 31
// Metals start at Fe at index n then continue to n+m where m = n - 4
// If n = 18, then only up to Mo is included
//______________________________________________________________________
// UNDIFFERENTIATED BODIES HAVE ALL ATOMS AS OXIDES
  for (int i = 0; i < n - m; i++) moles[i] = moles[i] * mass_factor * mass;
  for (int i = n - m; i < n; i++) moles[i] = moles[i] * mass_factor * mass + moles[i + m] * mass_factor * mass;
  moles_H2O = (moles_H2O + moles_H) * mass_factor * mass;
  moles_S_sil = (moles_S_sil + moles_S_met) * mass_factor * mass;
  moles_C_sil = (moles_C_sil + moles_C_met) * mass_factor * mass;
//______________________________________________________________________
// UNDIFFERENTIATED BODIES HAVE NO ATOMS AS METALS
  for (int i = n; i < (n + m); i++) moles[i] = 0.0;
  moles_oxygen = 0.0;
  moles_H = 0.0;
  moles_S_met = 0.0;
  moles_C_met = 0.0;

  return;
}

//______________________________________________________________________________________________________________________________________________
void differentiated_composition(double mass, double mass_factor,
                                double (&moles)[nm], double & moles_H2O, double & moles_S_sil, double & moles_C_sil,
                                double & moles_oxygen, double & moles_H, double & moles_S_met, double & moles_C_met)
{
// GET THE COMPOSITION OF AN DIFFERENTIATED, METAL-FREE BODY
// For reference, below are the indexes for elements in elements.h
// Element index: 0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22
// Element name: Al Mg Ca Na Fe Si Ni Co Nb Ta  V Cr Pt Pd Ru Ir  W Mo Mn Ga Zn In Cu
// Element index:            17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35
// Metals start at Fe at index n then continue to n+m where m = n - 4
// If n = 17, then only up to W is included
//______________________________________________________________________
// OXIDES
  for (int i = 0; i < n; i++) moles[i] = moles[i] * mass_factor * mass;
  moles_H2O = moles_H2O * mass_factor * mass;
  moles_S_sil = moles_S_sil * mass_factor * mass;
  moles_C_sil = moles_C_sil * mass_factor * mass;
//______________________________________________________________________
// METALS
  for (int i = n; i < (n + m); i++) moles[i] = moles[i] * mass_factor * mass; ;
  moles_oxygen = moles_oxygen * mass_factor * mass;
  moles_H = moles_H * mass_factor * mass;
  moles_S_met = moles_S_met * mass_factor * mass;
  moles_C_met = moles_C_met * mass_factor * mass;

  return;
}

//______________________________________________________________________________________________________________________________________________
void metal_silicate_equilibration(double mass,
                                  double (&mseq_moles)[nm], double & mseq_moles_H2O, double & mseq_moles_S_mantle, double & mseq_moles_C_mantle,
                                  double & mseq_moles_O, double & mseq_moles_H, double & mseq_moles_S_core, double & mseq_moles_C_core,
                                  double & P, double & T, double & fo2, double Fe_isotope, int & zero_xfe_liq) //GN fractionate
{
// GET THE COMPOSITION OF A DIFFERENTIATED, METAL-BEARING BODY
// For reference, below are the indexes for elements in elements.h
// Element index: 0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22
// Element name: Al Mg Ca Na Fe Si Ni Co Nb Ta  V Cr Pt Pd Ru Ir  W Mo Mn Ga Zn In Cu
// Element index:            17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35
// Metals start at Fe at index n then continue to n+m where m = n - 4
// If n = 17, then only up to W is included
//__________________________________________________________________________________________________________
// Input for react()
  double react_moles[n + m + 3] = {0.0};
  double react_moles_H2O = 0.0;
  double react_moles_S_mantle = 0.0;
  double react_moles_C_mantle = 0.0;
  double react_moles_O = 0.0;
  double react_moles_H = 0.0;
  double react_moles_S_core = 0.0;
  double react_moles_C_core = 0.0;
  fo2 = 0.0;

  mseq_call_no++;

  /*  if (mseq_moles_S_core < 0) cout << endl << "S1: " << mseq_call_no << tab << react_call_no << endl;
    if (mseq_moles_O < 0) cout << endl << "O1: " << mseq_call_no << tab << react_call_no << endl;
    if (mseq_moles_H2O < 0) cout << endl << "H2O1: " << mseq_call_no << tab << react_call_no << endl;
    if (mseq_moles_H < 0) cout << endl << "H1: " << mseq_call_no << tab << react_call_no << endl;*/

// Set the react variables for react()
  for (int i = 0; i < (n + m); i++) react_moles[i] = mseq_moles[i];
  react_moles[n + m] = react_moles_O = mseq_moles_O + mseq_moles_H2O;
  if (react_moles[n] > 1.0e-20) react_moles[n + m + 2] = 1.0;
  react_moles_H2O = 0.0;
  react_moles_S_mantle = mseq_moles_S_mantle;
  react_moles_C_mantle = mseq_moles_C_mantle;
  react_moles_H = mseq_moles_H + 2.0 * mseq_moles_H2O;
  react_moles_S_core = mseq_moles_S_core;
  react_moles_C_core = mseq_moles_C_core;

// Total S moles
  double total_moles_S = react_moles_S_mantle + react_moles_S_core;
  double total_moles_C = react_moles_C_mantle + react_moles_C_core;
// Sum of metal moles
  double total_moles_met = 0.0;
  for (int i = n; i < (n + m); i++) total_moles_met += react_moles[i];
  total_moles_met += react_moles_O;
  total_moles_met += react_moles_S_core;
  total_moles_met += react_moles_C_core;

// Sum of all non-oxygen moles
  double total_moles = 0.0;
  for (int i = 0; i < (n + m); i++) total_moles += react_moles[i];
  total_moles += react_moles_S_mantle + react_moles_S_core;
  total_moles += react_moles_C_mantle + react_moles_C_core;

// Ratio of metal S to all metal atoms
  double mole_frac_S_met = react_moles_S_core / total_moles_met;
  double old_mole_frac_S_met = 1.0;
  double older_mole_frac_S_met = 1.0;
// Ratio of metal C to all metal atoms
  double mole_frac_C_met = react_moles_C_core / total_moles_met;
  double old_mole_frac_C_met = 1.0;
// Get the pressure and temperature for react from the embryo mass
  P = p_factor * pressure(mass);
  T = temperature(P) + t_overstep;
// DIFFERENTIATION WITH S in equilibrating composition:
  if (total_moles_S > 0.0)
  {
    int it = 0;
    while (fabs(old_mole_frac_S_met - mole_frac_S_met) * total_moles_met / total_moles > 1e-6 || fabs(old_mole_frac_C_met - mole_frac_C_met) * total_moles_met / total_moles > 1e-6 ) //&& it<6
    {
      it++;
      // if (react_call_no < -1) testfile << it << tab << scientific << fabs(old_mole_frac_S_met - mole_frac_S_met) << tab << old_mole_frac_S_met << tab << mole_frac_S_met << tab << react_moles_S_core << tab << total_moles_met << tab << total_moles << endl;
      // Determine distribution coefficients for siderophile elements (except oxygen):
      double pkd[n];                 // Declation of pointer to kd array
      distrib_coeff(pkd, P, T, mole_frac_S_met);
      older_mole_frac_S_met = old_mole_frac_S_met;
      old_mole_frac_S_met = mole_frac_S_met;
      old_mole_frac_C_met = mole_frac_C_met;

      // ALTERNATIVE THIRD REACT CALL:
      double reacted_moles[nm3];
      react(reacted_moles, react_moles, react_moles_O, react_moles_S_core, react_moles_S_mantle, react_moles_C_core, react_moles_C_mantle, P, T, pkd);
      for (int i = 0; i < n + m; i++) react_moles[i] = reacted_moles[i];
      react_moles_O = react_moles[n + m] = reacted_moles[n + m];
      fo2 = react_moles[n + m + 1] = reacted_moles[n + m + 1];
      react_moles[n + m + 2] = reacted_moles[n + m + 2];

      zero_xfe_liq = int(react_moles[n + m + 2] + 0.5);

      if (zero_xfe_liq == 0) // Indicates no metal
      {
        react_moles_S_mantle = total_moles_S;
        react_moles_S_core = 0.0;
        react_moles_C_mantle = total_moles_C;
        react_moles_C_core = 0.0;
        mole_frac_S_met = old_mole_frac_S_met; // This terminates the iterations when there is no metal left after "react"
      }
      else
      {
        // Now calculate distribution of S between metal and silicate
        react_moles_S_core = S_partitioning_blanchard(react_moles, react_moles_O, react_moles_S_core, react_moles_S_mantle, react_moles_C_core, react_moles_C_mantle, P, T);
        react_moles_S_mantle = total_moles_S - react_moles_S_core;

        react_moles_C_core = C_partitioning(react_moles, react_moles_O, react_moles_S_core, react_moles_S_mantle, react_moles_C_core, react_moles_C_mantle, P, T);
        react_moles_C_mantle = total_moles_C - react_moles_C_core;

        WMo_partitioning(react_moles, react_moles_O, react_moles_S_core, react_moles_S_mantle, react_moles_C_core, react_moles_C_mantle, P, T);

        total_moles_met = 0.0;
        for (int i = n; i < (n + m); i++) total_moles_met += react_moles[i];
        total_moles_met += react_moles_O;
        total_moles_met += react_moles_S_core;
        total_moles_met += react_moles_C_core;
        mole_frac_S_met = react_moles_S_core / total_moles_met;
        mole_frac_C_met = react_moles_C_core / total_moles_met;

        if (fabs(older_mole_frac_S_met - mole_frac_S_met) * total_moles_met / total_moles < 1e-6) mole_frac_S_met = old_mole_frac_S_met;
      }

      if (it > 10000) //Generally 12 iterations or less
      {
        mole_frac_S_met = old_mole_frac_S_met;
        mole_frac_C_met = old_mole_frac_C_met;
      }
    }
  }
// DIFFERENTIATION WITH NO S in equilibrating composition:
  else
  {
    int it = 0;
    while ( fabs(old_mole_frac_C_met - mole_frac_C_met) * total_moles_met / total_moles > 1e-6 ) //&& it<6
    {
      it++;
      // if (react_call_no < -1 ) testfile << it << tab << scientific << old_mole_frac_C_met - mole_frac_C_met << endl;
      if (it > 10000)
      {
        cout << endl << "second it = " << it << " and react call no = " << react_call_no << endl; //Generally 12 iterations or less
        cout << "C differed = " << fabs(old_mole_frac_C_met - mole_frac_C_met) << endl;
        exit(1);
      }

      // Determine distribution coefficients for siderophile elements (except oxygen):
      double pkd[n];                       // Declation of pointer to kd array
      distrib_coeff(pkd, P, T, mole_frac_S_met);
      old_mole_frac_C_met = mole_frac_C_met;

      // Mass balance/partitioning calculation
      // THIRD REACT CALL:
      double reacted_moles[nm3];
      react(reacted_moles, react_moles, react_moles_O, react_moles_S_core, react_moles_S_mantle, react_moles_C_core, react_moles_C_mantle, P, T, pkd);
      for (int i = 0; i < n + m; i++) react_moles[i] = reacted_moles[i];
      react_moles_O = react_moles[n + m] = reacted_moles[n + m];
      fo2 = react_moles[n + m + 1] = reacted_moles[n + m + 1];
      react_moles[n + m + 2] = reacted_moles[n + m + 2];
      react_moles_S_core = 0.0;
      react_moles_S_mantle = 0.0;
      react_moles_C_core = C_partitioning(react_moles, react_moles_O, react_moles_S_core, react_moles_S_mantle, react_moles_C_core, react_moles_C_mantle, P, T);
      react_moles_C_mantle = total_moles_C - react_moles_C_core;
      WMo_partitioning(react_moles, react_moles_O, react_moles_S_core, react_moles_S_mantle, react_moles_C_core, react_moles_C_mantle, P, T);

      total_moles_met = 0.0;
      for (int i = n; i < (n + m); i++) total_moles_met += react_moles[i];
      total_moles_met += react_moles_O;
      total_moles_met += react_moles_S_core;
      total_moles_met += react_moles_C_core;
      mole_frac_C_met = react_moles_C_core / total_moles_met;

    }
  }

  zero_xfe_liq = int(react_moles[n + m + 2] + 0.5);
  if (zero_xfe_liq == 1)
  {
    react_moles_H2O = 0.0;
  }
  else if (zero_xfe_liq == 0) // Indicates no metal
  {
    react_moles_H2O = 0.5 * react_moles_H;
    react_moles_H = 0.0;
    react_moles_S_mantle = total_moles_S;
    react_moles_S_core = 0.0;
    react_moles_C_mantle = total_moles_C;
    react_moles_C_core = 0.0;
  }

// Set final compositions
  for (int i = 0; i < (n + m); i++) mseq_moles[i] = react_moles[i];
  mseq_moles_H2O = react_moles_H2O;
  mseq_moles_S_mantle = react_moles_S_mantle;
  mseq_moles_C_mantle = react_moles_C_mantle;
  mseq_moles_O = react_moles_O;
  mseq_moles_H = react_moles_H;
  mseq_moles_S_core = react_moles_S_core;
  mseq_moles_C_core = react_moles_C_core;

//  if (mseq_call_no == 10322) cout << 10322 << tab << fo2 << endl;

  /*
    if (mseq_moles_S_core < 0) cout << endl << "S2: " << mseq_call_no << tab << react_call_no << endl;
    if (mseq_moles_O < 0) cout << endl << "O2: " << mseq_call_no << tab << react_call_no << endl;
    if (mseq_moles_H2O < 0) cout << endl << "H2O2: " << mseq_call_no << tab << react_call_no << endl;
    if (mseq_moles_H < 0) cout << endl << "H2: " << mseq_call_no << tab << react_call_no << endl;
  */

  return;
}

//______________________________________________________________________________________________________________________________________________
void initiate_bulk_composition(double mass, int e_no, int coll_no,
                               double moles[], double moles_H2O, double moles_S_sil, double moles_C_sil,
                               double moles_oxygen, double moles_H, double moles_S_met, double moles_C_met,
                               double P, double T, double fo2)
{

// Define mantle and core concentrations of differentiated body

  for (int i = 0; i < (n + m); i++) e_moles[e_no][i][coll_no] = moles[i];
  e_moles_H2O[e_no][coll_no] = moles_H2O;
  e_moles_S_mantle[e_no][coll_no] = moles_S_sil;
  e_moles_C_mantle[e_no][coll_no] = moles_C_sil;
  e_moles_oxygen[e_no][coll_no] = moles_oxygen;
  e_moles_H[e_no][coll_no] = moles_H;
  e_moles_S_core[e_no][coll_no] = moles_S_met;
  e_moles_C_core[e_no][coll_no] = moles_C_met;
  e_P[e_no][coll_no] = P;
  e_T[e_no][coll_no] = T;
  e_fO2[e_no][coll_no] = fo2;

  // // Print out e_moles number to check what the initial moles/ where the moles for each element went
  // for (int i = 0; i < (n+m); i++) cout << "emoles " << i << ", " << moles[i] << endl; 

  for (int i = 0; i < (n + m); i++) s_moles[e_no][i][coll_no] = moles[i];
  s_moles_H2O[e_no][coll_no] = moles_H2O;
  s_moles_S_mantle[e_no][coll_no] = moles_S_sil;
  s_moles_C_mantle[e_no][coll_no] = moles_C_sil;
  s_moles_oxygen[e_no][coll_no] = moles_oxygen;
  s_moles_H[e_no][coll_no] = moles_H;
  s_moles_S_core[e_no][coll_no] = moles_S_met;
  s_moles_C_core[e_no][coll_no] = moles_C_met;

  // Store compositions of metal that are being added to proto-core (composition of core-forming increments)
  for (int i = n; i < (n + m); i++) incr_moles[e_no][i][coll_no] = e_moles[e_no][i][coll_no];
  incr_moles_oxygen[e_no][coll_no] = e_moles_oxygen[e_no][coll_no];
  incr_moles_H[e_no][coll_no] = e_moles_H[e_no][coll_no];
  incr_moles_S[e_no][coll_no] = e_moles_S_core[e_no][coll_no];
  incr_moles_C[e_no][coll_no] = e_moles_C_core[e_no][coll_no];

  return;
}


//______________________________________________________________________________________________________________________________________________
void initiate_isotopic_composition(int e_no, int coll_no, double Fe_sil_isotope[], double Fe_met_isotope[])
{
   

  // [0] Fe54
  e_Fe_sil_isotope[e_no][0][coll_no] = Fe_sil_isotope[0];
  e_Fe_met_isotope[e_no][0][coll_no] = Fe_met_isotope[0];
  
  // [1] Fe57
  e_Fe_sil_isotope[e_no][1][coll_no] = Fe_sil_isotope[1]; //3isotope
  e_Fe_met_isotope[e_no][1][coll_no] = Fe_met_isotope[1]; //3isotope

  // [2] Fe57
  e_Fe_sil_isotope[e_no][2][coll_no] = Fe_sil_isotope[2]; //3isotope
  e_Fe_met_isotope[e_no][2][coll_no] = Fe_met_isotope[2]; //3isotope

  
  return;
}

//______________________________________________________________________________________________________________________________________________
void update_bulk_composition(int e_no, int coll_no,
                             double equil_moles[], double equil_moles_H2O, double equil_moles_S_sil, double equil_moles_C_sil,
                             double equil_moles_oxygen, double equil_moles_H, double equil_moles_S_met, double equil_moles_C_met,
                             double noneq_moles[], double noneq_moles_H2O, double noneq_moles_S_sil, double noneq_moles_C_sil,
                             double noneq_moles_oxygen, double noneq_moles_H, double noneq_moles_S_met, double noneq_moles_C_met,
                             double P, double T, double fo2)
{
// Define mantle and core concentrations of differentiated body
  for (int i = 0; i < (n + m); i++) e_moles[e_no][i][coll_no] = equil_moles[i] + noneq_moles[i];
  e_moles_H2O[e_no][coll_no] = equil_moles_H2O + noneq_moles_H2O;
  e_moles_S_mantle[e_no][coll_no] = equil_moles_S_sil + noneq_moles_S_sil;
  e_moles_C_mantle[e_no][coll_no] = equil_moles_C_sil + noneq_moles_C_sil;
  e_moles_oxygen[e_no][coll_no] = equil_moles_oxygen + noneq_moles_oxygen;
  e_moles_H[e_no][coll_no] = equil_moles_H + noneq_moles_H;
  e_moles_S_core[e_no][coll_no] = equil_moles_S_met + noneq_moles_S_met;
  e_moles_C_core[e_no][coll_no] = equil_moles_C_met + noneq_moles_C_met;
  e_P[e_no][coll_no] = P;
  e_T[e_no][coll_no] = T;
  e_fO2[e_no][coll_no] = fo2;

  return;
}


//______________________________________________________________________________________________________________________________________________
void update_Fe_composition(int e_no, int coll_no,
                           double equil_moles[], double equil_Fe_sil_isotope[], double equil_Fe_met_isotope[],
                           double noneq_moles[], double noneq_Fe_sil_isotope[], double noneq_Fe_met_isotope[])
{
   
  // [0] Fe54
  e_Fe_sil_isotope[e_no][0][coll_no] = (equil_Fe_sil_isotope[0] * equil_moles[n - m] + noneq_Fe_sil_isotope[0] * noneq_moles[n - m]) / (equil_moles[n - m] + noneq_moles[n - m]);
  e_Fe_met_isotope[e_no][0][coll_no] = (equil_Fe_met_isotope[0] * equil_moles[n] + noneq_Fe_met_isotope[0] * noneq_moles[n]) / (equil_moles[n] + noneq_moles[n]);

  // [1] Fe56 //3isotope
  e_Fe_sil_isotope[e_no][1][coll_no] = (equil_Fe_sil_isotope[1] * equil_moles[n - m] + noneq_Fe_sil_isotope[1] * noneq_moles[n - m]) / (equil_moles[n - m] + noneq_moles[n - m]);
  e_Fe_met_isotope[e_no][1][coll_no] = (equil_Fe_met_isotope[1] * equil_moles[n] + noneq_Fe_met_isotope[1] * noneq_moles[n]) / (equil_moles[n] + noneq_moles[n]);
 
  // [2] Fe57 //3isotope
  e_Fe_sil_isotope[e_no][2][coll_no] = (equil_Fe_sil_isotope[2] * equil_moles[n - m] + noneq_Fe_sil_isotope[2] * noneq_moles[n - m]) / (equil_moles[n - m] + noneq_moles[n - m]);
  e_Fe_met_isotope[e_no][2][coll_no] = (equil_Fe_met_isotope[2] * equil_moles[n] + noneq_Fe_met_isotope[2] * noneq_moles[n]) / (equil_moles[n] + noneq_moles[n]);

 
  return;
}

//______________________________________________________________________________________________________________________________________________
void update_contributions(int no_1, int no_2, double mass_1, double mass_2)
{
// First occurence of this body
  if (b_status[no_1] < 2)
  {
    e_em_mass[no_1] = mass_1;
// if no_2 is a planetesimal
    if (mass_2 <= max_pl * mass_adjust) e_pl_mass[no_1] = mass_2;
// if no_2 is an embryo
    else if (mass_2 > max_pl * mass_adjust)
    {
      if (b_status[no_2] < 2) e_em_mass[no_1] += mass_2;
      else if (b_status[no_1] > 1)
      {
        e_em_mass[no_1] += e_em_mass[no_2];
        e_pl_mass[no_1] += e_pl_mass[no_2];
      }
    }
  }
// #1 already has a collisional history
  else if (b_status[no_1] > 1)
  {
// #2 is a planetesimal
    if (mass_2 <= max_pl * mass_adjust) e_pl_mass[no_1] += mass_2;
// #2 is an embryo
    else if (mass_2 > max_pl * mass_adjust)
    {
      if (b_status[no_2] < 2) e_em_mass[no_1] += mass_2;
      else if (b_status[no_2] > 1)
      {
        e_em_mass[no_1] += e_em_mass[no_2];
        e_pl_mass[no_1] += e_pl_mass[no_2];
      }
    }
  }

// Sum masses of H2O-bearing bodies accreted from beyond <dist_3>
  e_SN_mass[no_1] = e_SN_mass[no_1] + e_SN_mass[no_2];

  return;
}

//______________________________________________________________________________________________________________________________________________
void accretion_no_equilibration(double (& moles_target)[nm], double & moles_H2O_target, double & moles_S_mantle_target, double & moles_C_mantle_target,
                                double & moles_O_target, double & moles_H_target, double & moles_S_target, double & moles_C_target,
                                double moles_projectile[], double moles_H2O_projectile, double moles_S_mantle_projectile, double moles_C_mantle_projectile,
                                double moles_O_projectile, double moles_H_projectile, double moles_S_projectile, double moles_C_projectile)
{
// Mantle is combined with mantle and core with core
  for (int i = 0; i < (n + m); i++) moles_target[i] = moles_target[i] + moles_projectile[i];
  moles_H2O_target = moles_H2O_target + moles_H2O_projectile;
  moles_S_mantle_target = moles_S_mantle_target + moles_S_mantle_projectile;
  moles_C_mantle_target = moles_C_mantle_target + moles_C_mantle_projectile;
  moles_O_target = moles_O_target + moles_O_projectile;
  moles_H_target = moles_H_target + moles_H_projectile;
  moles_S_target = moles_S_target + moles_S_projectile;
  moles_C_target = moles_C_target + moles_C_projectile;
  return;
}

//______________________________________________________________________________________________________________________________________________
double silicate_fraction(double mass_1, double mass_2,
                         double no_1_moles[], double no_1_moles_H2O, double no_1_moles_S_mantle, double no_1_moles_C_mantle,
                         double no_1_moles_oxygen, double no_1_moles_H, double no_1_moles_S_core, double no_1_moles_C_core,
                         double no_2_moles[], double no_2_moles_H2O, double no_2_moles_S_mantle, double no_2_moles_C_mantle,
                         double no_2_moles_oxygen, double no_2_moles_H, double no_2_moles_S_core, double no_2_moles_C_core)
{
// Calculates fraction of target's silicate mantle that equilibrates with impactor's core using Eq.(6) of Deguen et al (2011, EPSL)
// Assumed ratio metal-density/silicate-density (2.29 WORKS WELL FOR EARTH)
  const double theta = 2.29;
// constant in Eq.(6) of Deguen et al (2011) (0.2 gives slightly better results than 0.25)
  const double alpha = 0.25;
// Mass of Earth (kg)
  const double earth_mass = 5.97219E24;
  const double moon_mass = 7.34767E22;      // Mass of Moon in kg
  double mass[2] = {0.0};
  double met_frac[2] = {0.0};

  mass[0] = mass_1 * moon_mass;
  mass[1] = mass_2 * moon_mass;

// CALCULATE METAL FRACTION FOR EACH OF THE TWO BODIES (i = 0, 1):
  for (int i = 0; i < 2; i++)
  {
    double weight[n + m + 1] = {0.0};
    double moles[n + m] = {0.0};
    double moles_H2O = 0.0;
    double moles_S_mantle = 0.0;
    double moles_C_mantle = 0.0;
    double moles_O = 0.0;
    double moles_H = 0.0;
    double moles_S_core = 0.0;
    double moles_C_core = 0.0;

    if (i == 0)
    {
      for (int j = 0; j < n + m; j++) moles[j] = no_1_moles[j];
      moles_H2O = no_1_moles_H2O;
      moles_S_mantle = no_1_moles_S_mantle;
      moles_C_mantle = no_1_moles_C_mantle;
      moles_O = no_1_moles_oxygen;
      moles_H = no_1_moles_H;
      moles_S_core = no_1_moles_S_core;
      moles_C_core = no_1_moles_C_core;
    }
    else if (i == 1)
    {
      for (int j = 0; j < n + m; j++) moles[j] = no_2_moles[j];
      moles_H2O = no_2_moles_H2O;
      moles_O = no_2_moles_oxygen;
      moles_H = no_2_moles_H;
      moles_S_mantle = no_2_moles_S_mantle;
      moles_S_core = no_2_moles_S_core;
      moles_C_mantle = no_2_moles_C_mantle;
      moles_C_core = no_2_moles_C_core;
    }

// Calculate silicate chemical components
    double sum_sil = 0;
    for (int j = 0; j < n; j++)
    {
      weight[j] = moles[j] * MW[j];
      sum_sil += weight[j];
    }
    sum_sil += moles_H2O * MW_H2O;
    sum_sil += moles_S_mantle * AW_S;
    sum_sil += moles_C_mantle * AW_C;

// Calculate metal chemical components
    double sum_met = 0;
    for (int j = n; j < (n + m); j++)
    {
      weight[j] = moles[j] * AW[j - m];
      sum_met += weight[j];
    }
    weight[n + m] = moles_O * AW_O;
    sum_met += weight[n + m];
    sum_met += moles_H * AW_H;
    sum_met += moles_S_core * AW_S;
    sum_met += moles_C_core * AW_C;


    met_frac[i] = sum_met / (sum_met + sum_sil);
  }



  double eq_depth = 0.0;
  double imp_core_rad = 0.0;
  double imp_core_vol = 0.0;
  double target_mantle_vol = 0.0;

  for (int i = 0; i < 2; i++)
  {
    double bulk_density = 3739.3 + 2.9803e-22 * mass[i]; // kg m-3
    double body_rad = pow((3.0 * mass[i] / (4.0 * pi * bulk_density)), (1.0 / 3.0));
    double total_vol = 4.0 * pi * pow(body_rad, 3.0) / 3.0;
    double mass_met = met_frac[i] * mass[i];

// Calculate radius of core:
    double rad_cubed = 3.0 * mass_met / (4.0 * pi * theta * (bulk_density - mass_met * (1.0 - 1.0 / theta) / total_vol));
    double core_rad = pow(rad_cubed, (1.0 / 3.0));
    double core_vol = 4.0 * pi * pow(core_rad, 3.0) / 3.0;
    double mantle_vol = total_vol - core_vol;
// Target body
    if (i == 0)
    {
      eq_depth = (body_rad - core_rad) * p_factor;
      target_mantle_vol = mantle_vol;
    }
// Impactor
    if (i == 1)
    {
      imp_core_rad = core_rad;
      imp_core_vol = core_vol;
    }
  }
// Vol fraction of metal (Eq.(6) of Deguen et al 2011)
  double phi = pow((1 + alpha * eq_depth / imp_core_rad), -3.0);
  double vol_eq_mantle = imp_core_vol * (1.0 - phi) / phi;

// fraction of target's mantle that equilibrates
  return (vol_eq_mantle / target_mantle_vol);
}

//______________________________________________________________________________________________________________________________________________
void record_silicate_fraction(double mass_2, double xsil)
{
// Record the min, avg, and max xsil for planetesimals and embryos
// Also record the number of planetesimal and embryo accretion events
// No_2 is a planetesimal
  if (mass_2 <= max_pl * mass_adjust)
  {
    n_pl++;
    xsil_pl += xsil;
    if (xsil > xsil_pl_max) xsil_pl_max = xsil;
    if (xsil < xsil_pl_min) xsil_pl_min = xsil;
  }
// No_2 is an embryo
  else if (mass_2 > max_pl * mass_adjust)
  {
    n_em++;
    xsil_em += xsil;
    if (xsil > xsil_em_max) xsil_em_max = xsil;
    if (xsil < xsil_em_min) xsil_em_min = xsil;
  }
  return;
}

void identify_molar_reservoirs(double xsil, double xmet,
                               double no_1_moles[], double no_1_moles_H2O, double no_1_moles_S_mantle, double no_1_moles_C_mantle,
                               double no_1_moles_oxygen, double no_1_moles_H, double no_1_moles_S_core, double no_1_moles_C_core,
                               double no_1_Fe_sil_isotope[3], double no_1_Fe_met_isotope[3],
                               double no_2_moles[], double no_2_moles_H2O, double no_2_moles_S_mantle, double no_2_moles_C_mantle,
                               double no_2_moles_oxygen, double no_2_moles_H, double no_2_moles_S_core, double no_2_moles_C_core,
                               double no_2_Fe_sil_isotope[3], double no_2_Fe_met_isotope[3],
                               double (&equil_moles)[nm], double & equil_moles_H2O, double & equil_moles_S_mantle, double & equil_moles_C_mantle,
                               double & equil_moles_oxygen, double & equil_moles_H, double & equil_moles_S_core, double & equil_moles_C_core,
                               double (&equil_Fe_sil_isotope)[3], double (&equil_Fe_met_isotope)[3],
                               double (&noneq_moles)[nm], double & noneq_moles_H2O, double & noneq_moles_S_mantle, double & noneq_moles_C_mantle,
                               double & noneq_moles_oxygen, double & noneq_moles_H, double & noneq_moles_S_core, double & noneq_moles_C_core,
                               double (&noneq_Fe_sil_isotope)[3], double (&noneq_Fe_met_isotope)[3])
{

// Equilibrating moles
  for (int i = 0; i < n; i++) equil_moles[i] = xsil * no_1_moles[i];
  equil_moles_H2O = 0.0;
  equil_moles_S_mantle = xsil * no_1_moles_S_mantle;
  equil_moles_C_mantle = xsil * no_1_moles_C_mantle;
  for (int i = n; i < (n + m); i++) equil_moles[i] = xmet * no_2_moles[i];
  equil_moles_oxygen = xmet * no_2_moles_oxygen + xsil * no_1_moles_H2O;
  equil_moles_H = xmet * no_2_moles_H + 2.0 * xsil * no_1_moles_H2O;
  equil_moles_S_core = xmet * no_2_moles_S_core;
  equil_moles_C_core = xmet * no_2_moles_C_core;
  
 
// equilibrating moles isotopic: //3isotope ToDo -- is this correct????
  equil_Fe_sil_isotope[0] = no_1_Fe_sil_isotope[0];
  equil_Fe_sil_isotope[1] = no_1_Fe_sil_isotope[1];
  equil_Fe_sil_isotope[2] = no_1_Fe_sil_isotope[2];
  equil_Fe_met_isotope[0] = no_2_Fe_met_isotope[0];
  equil_Fe_met_isotope[1] = no_2_Fe_met_isotope[1];
  equil_Fe_met_isotope[2] = no_2_Fe_met_isotope[2];

// Nonequilibrating moles
  for (int i = 0; i < n; i++) noneq_moles[i] = (1.0 - xsil) * no_1_moles[i] + no_2_moles[i];
  noneq_moles_H2O = (1.0 - xsil) * no_1_moles_H2O + no_2_moles_H2O;
  noneq_moles_S_mantle = (1.0 - xsil) * no_1_moles_S_mantle + no_2_moles_S_mantle;
  noneq_moles_C_mantle = (1.0 - xsil) * no_1_moles_C_mantle + no_2_moles_C_mantle;
  for (int i = n; i < (n + m); i++) noneq_moles[i] = no_1_moles[i] + (1.0 - xmet) * no_2_moles[i];
  noneq_moles_oxygen = no_1_moles_oxygen + (1.0 - xmet) * no_2_moles_oxygen;
  noneq_moles_H = no_1_moles_H + (1.0 - xmet) * no_2_moles_H;
  noneq_moles_S_core = no_1_moles_S_core + (1.0 - xmet) * no_2_moles_S_core;
  noneq_moles_C_core = no_1_moles_C_core + (1.0 - xmet) * no_2_moles_C_core;

// non_equilibrating moles isotopic //3isotope ToDo
  noneq_Fe_sil_isotope[0] = (no_1_Fe_sil_isotope[0] * no_1_moles[n - m] * (1.0 - xsil) + no_2_Fe_sil_isotope[0] * no_2_moles[n - m]) / (no_1_moles[n - m] * (1.0 - xsil) + no_2_moles[n - m]);
  noneq_Fe_sil_isotope[1] = (no_1_Fe_sil_isotope[1] * no_1_moles[n - m] * (1.0 - xsil) + no_2_Fe_sil_isotope[1] * no_2_moles[n - m]) / (no_1_moles[n - m] * (1.0 - xsil) + no_2_moles[n - m]);
  noneq_Fe_sil_isotope[2] = (no_1_Fe_sil_isotope[2] * no_1_moles[n - m] * (1.0 - xsil) + no_2_Fe_sil_isotope[2] * no_2_moles[n - m]) / (no_1_moles[n - m] * (1.0 - xsil) + no_2_moles[n - m]); //3isotope
  noneq_Fe_met_isotope[0] = (no_1_Fe_met_isotope[0] * no_1_moles[n] + no_2_Fe_met_isotope[0] * no_2_moles[n] * (1.0 - xmet)) / (no_1_moles[n] + no_2_moles[n] * (1.0 - xmet));
  noneq_Fe_met_isotope[1] = (no_1_Fe_met_isotope[1] * no_1_moles[n] + no_2_Fe_met_isotope[1] * no_2_moles[n] * (1.0 - xmet)) / (no_1_moles[n] + no_2_moles[n] * (1.0 - xmet));
  noneq_Fe_met_isotope[2] = (no_1_Fe_met_isotope[2] * no_1_moles[n] + no_2_Fe_met_isotope[2] * no_2_moles[n] * (1.0 - xmet)) / (no_1_moles[n] + no_2_moles[n] * (1.0 - xmet)); //3isotope

      
  return;
}

void identify_Fe_isotopic_reservoirs()
{
  return;
}
void update_sources(int e_no, int coll_no, double xsil, double xmet,
                    double no_1_moles[], double no_1_moles_H2O, double no_1_moles_S_mantle, double no_1_moles_C_mantle,
                    double no_1_moles_oxygen, double no_1_moles_H, double no_1_moles_S_core, double no_1_moles_C_core,
                    double no_2_moles[], double no_2_moles_H2O, double no_2_moles_S_mantle, double no_2_moles_C_mantle,
                    double no_2_moles_oxygen, double no_2_moles_H, double no_2_moles_S_core, double no_2_moles_C_core,
                    double equil_moles[], double equil_moles_H2O, double equil_moles_S_mantle, double equil_moles_C_mantle,
                    double equil_moles_oxygen, double equil_moles_H, double equil_moles_S_core, double equil_moles_C_core,
                    double noneq_moles[], double noneq_moles_H2O, double noneq_moles_S_mantle, double noneq_moles_C_mantle,
                    double noneq_moles_oxygen, double noneq_moles_H, double noneq_moles_S_core, double noneq_moles_C_core)
{
// GET THE COMPOSITION OF AN UNDIFFERENTIATED, METAL-FREE BODY
// For reference, below are the indexes for elements in elements.h
// Element index: 0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22
// Element name: Al Mg Ca Na Fe Si Ni Co Nb Ta  V Cr Pt Pd Ru Ir  W Mo Mn Ga Zn In Cu
// Element index:            17 18 19 20 21 22 23 24 25 26 27 28 29
// Metals start at Fe at index n then continue to n+m where m = n - 4
// If n = 17, then only up to W is included
//______________________________________________________________________
// EVERYTHING IS ALWAYS DONE PROPORTIONALLY, I.E. ASSUMING PERFECT MIXING

// Add nonequilibrating moles from mantle and core of impactor
  for (int i = 0; i < n; i++) s_moles[e_no][i][coll_no] = no_2_moles[i];
  s_moles_H2O[e_no][coll_no] = no_2_moles_H2O;
  s_moles_S_mantle[e_no][coll_no] = no_2_moles_S_mantle;
  s_moles_C_mantle[e_no][coll_no] = no_2_moles_C_mantle;
  for (int i = n; i < (n + m); i++) s_moles[e_no][i][coll_no] = (1.0 - xmet) * no_2_moles[i];
  s_moles_oxygen[e_no][coll_no] = (1.0 - xmet) * no_2_moles_oxygen;
  s_moles_H[e_no][coll_no] = (1.0 - xmet) * no_2_moles_H;
  s_moles_S_core[e_no][coll_no] = (1.0 - xmet) * no_2_moles_S_core;
  s_moles_C_core[e_no][coll_no] = (1.0 - xmet) * no_2_moles_C_core;

  double d_sid[m] = {0.0};
  double sili_O = 0.0;
  double d_O = 0.0;
  double d_H = 0.0;
  double d_S = 0.0;
  double d_C = 0.0;
  for (int i = 0; i < m; i++) if (equil_moles[i + n] > 0.0) d_sid[i] = equil_moles[i + n] / (equil_moles[i + n] + equil_moles[i + n - m]);
  if (equil_moles_H > 0.0) d_H = equil_moles_H / (equil_moles_H + equil_moles_H2O);
  if (equil_moles_S_core > 0.0) d_S = equil_moles_S_core / (equil_moles_S_core + equil_moles_S_mantle);
  if (equil_moles_C_core > 0.0) d_C = equil_moles_C_core / (equil_moles_C_core + equil_moles_C_mantle);
  sili_O = 1.5 * equil_moles[0] + equil_moles[1] + equil_moles[2] + equil_moles[3] + equil_moles[4] + 2.0 * equil_moles[5] + 2.0 * equil_moles_H2O;
  if (equil_moles_oxygen > 0.0) d_O = equil_moles_oxygen / (equil_moles_oxygen + sili_O);

// Add equilibrating moles from core of impactor to mantle and core
  for (int i = (n - m); i < n; i++) s_moles[e_no][i][coll_no] += (1.0 - d_sid[i - (n - m)]) * xmet * no_2_moles[i + m];
  s_moles_H2O[e_no][coll_no] += (1.0 - d_H) * xmet * no_2_moles_H;
  s_moles_S_mantle[e_no][coll_no] += (1.0 - d_S) * xmet * no_2_moles_S_core;
  s_moles_C_mantle[e_no][coll_no] += (1.0 - d_C) * xmet * no_2_moles_C_core;
  for (int i = n; i < (n + m); i++) s_moles[e_no][i][coll_no] += d_sid[i - n] * xmet * no_2_moles[i];
  s_moles_oxygen[e_no][coll_no] += d_O * xmet * no_2_moles_oxygen;
  s_moles_H[e_no][coll_no] += d_H * xmet * no_2_moles_H;
  s_moles_S_core[e_no][coll_no] += d_S * xmet * no_2_moles_S_core;
  s_moles_C_core[e_no][coll_no] += d_C * xmet * no_2_moles_C_core;

// Move moles from mantle sources to the core during this core formation event
  for (int source = 0; source < coll_no; source++)
  {
// Add to core
    for (int i = n; i < n + m; i++) s_moles[e_no][i][source] += d_sid[i - n] * xsil * s_moles[e_no][i - m][source];
    s_moles_oxygen[e_no][source] += d_O * xsil * (1.5 * s_moles[e_no][0][source] + s_moles[e_no][1][source] + s_moles[e_no][2][source] + s_moles[e_no][3][source] + s_moles[e_no][4][source] + 2.0 * s_moles[e_no][5][source] + 2.0 * s_moles_H2O[e_no][source]);
    s_moles_H[e_no][source] += d_H * xsil * s_moles_H2O[e_no][source];
    s_moles_S_core[e_no][source] += d_S * xsil * s_moles_S_mantle[e_no][source];
    s_moles_C_core[e_no][source] += d_C * xsil * s_moles_C_mantle[e_no][source];
    // Take from mantle
    for (int i = (n - m); i < n; i++) s_moles[e_no][i][source] -= d_sid[i - (n - m)] * xsil * s_moles[e_no][i][source];
    s_moles_H2O[e_no][source] -= d_H * xsil * s_moles_H2O[e_no][source];
    s_moles_S_mantle[e_no][source] -= d_S * xsil * s_moles_S_mantle[e_no][source];
    s_moles_C_mantle[e_no][source] -= d_C * xsil * s_moles_C_mantle[e_no][source];
  }

  return;
}

void update_sources_hse(int e_no, int coll_no, double h_moles[], double h_moles_S, double h_moles_FeS, double pmoles_HSE[])
{
// GET THE COMPOSITION OF AN UNDIFFERENTIATED, METAL-FREE BODY
// For reference, below are the indexes for elements in elements.h
// Element index: 0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22
// Element name: Al Mg Ca Na Fe Si Ni Co Nb Ta  V Cr Pt Pd Ru Ir  W Mo Mn Ga Zn In Cu
// Element index:            17 18 19 20 21 22 23 24 25 26 27 28 29
// Metals start at Fe at index n then continue to n+m where m = n - 4
// If n = 17, then only up to W is included
//______________________________________________________________________
// EVERYTHING IS ALWAYS DONE PROPORTIONALLY, I.E. ASSUMING PERFECT MIXING

// New concentrations of S and Fe after extraction of FeS to the core
  double d_S = 0.0;
  double d_Fe = 0.0;
  double d_hse[n_hse] = {0.0};
  if (h_moles_FeS > 0.0) d_S = h_moles_FeS / (h_moles_FeS + h_moles_S);
  if (h_moles_FeS > 0.0) d_Fe = h_moles_FeS / (h_moles_FeS + h_moles[n - m]);
  for (int hse_no = 0; hse_no < n_hse; hse_no++) if (pmoles_HSE[hse_no] > 0.0) d_hse[hse_no] = pmoles_HSE[hse_no] / (pmoles_HSE[hse_no] + pmoles_HSE[hse_no + n_hse]);

// Move moles from mantle sources to the core during this core formation event
  for (int source = 0; source < coll_no + 1; source++)
  {
// Add to core
    s_moles[e_no][n][source] += d_Fe * s_moles[e_no][n - m][source];
    s_moles_S_core[e_no][source] += d_S * s_moles_S_mantle[e_no][source];
    int hse_no = 0;
    for (int i = n; i < (n + m); i++)
    {
      if (el_type[i] == 'h' || el_type[i] == 'H')
      {
        s_moles[e_no][i][source] += d_hse[hse_no] * s_moles[e_no][i - (n - m)][source]; // Subtract the moles of HSEs contained in FeS from the mantle:
        hse_no++;
      }
    }
// Take from mantle
    s_moles[e_no][n - m][source] -= d_Fe * s_moles[e_no][n - m][source];
    s_moles_S_mantle[e_no][source] -= d_S * s_moles_S_mantle[e_no][source];
    hse_no = 0;
    for (int i = n - m; i < n; i++)
    {
      if (el_type[i] == 'h' || el_type[i] == 'H')
      {
        s_moles[e_no][i][source] -= d_hse[hse_no] * s_moles[e_no][i][source]; // Subtract the moles of HSEs contained in FeS from the mantle:
        hse_no++;
      }
    }
  }

  return;
}


//______________________________________________________________________________________________________________________________________________
void starting_isotopic_proportions(double (&Fe_sil_isotope)[3], double (&Fe_met_isotope)[3])//GN fractionate
{
  // GIVE initial constant STARTING ISOTOPIC COMPOSITIONS
  double Fe54, Fe56, Fe57;
  double Fe_Total;
  // INITIAL CHONDRITIC Fe isotopic composition == IRMM-014
  // from Isotope Catalogue 2008
  // This applies to non-differentiated bodies
  // Fe54 = 5.845;
  // Fe56 = 91.754; //3isotope
  // Fe57 = 2.119; //3isotope
  Fe54 = Fe54_IRMM;
  Fe56 = Fe56_IRMM;
  Fe57 = Fe57_IRMM;
  // 58Fe=0.281
  Fe_Total = Fe54 + Fe56 + Fe57; // to normalize //3isotope

  //silicate
  Fe_sil_isotope[0] = Fe54 / Fe_Total; // initiate 54Fe value in silicate
  Fe_sil_isotope[1] = Fe56 / Fe_Total; // initiate 56Fe value in silicate //3isotope
  Fe_sil_isotope[2] = Fe57 / Fe_Total; // initiate 57Fe value in silicate //3isotope
  
  //metal 
  Fe_met_isotope[0] = Fe54 / Fe_Total; // initiate 54Fe value in metal
  Fe_met_isotope[1] = Fe56 / Fe_Total; // initiate 56Fe value in metal //3isotope
  Fe_met_isotope[2] = Fe57 / Fe_Total; // initiate 57Fe value in metal //3isotope
  
  // cout << "------------ starting_isotopic_proportions() ------------" << endl;
  // cout << "Fe_sil_isotope[0]: " << Fe_sil_isotope[0] << endl;
  // cout << "Fe_sil_isotope[1]: " << Fe_sil_isotope[1] << endl;
  // cout << "Fe_sil_isotope[2] : " << Fe_sil_isotope[2]  << endl;
  
  
  // cout << "Fe54_IRMM: " << Fe54_IRMM << endl;
  // cout << "Fe56_IRMM : " << Fe56_IRMM  << endl;
  // cout << "Fe57_IRMM : " << Fe57_IRMM  << endl;
  // double mantle_d56_here = ((Fe_sil_isotope[1] / Fe_sil_isotope[0]) / (Fe56_IRMM / Fe54_IRMM) - 1.0) * 1000.0 ;
  // double core_d56_here = ((Fe_met_isotope[1] / Fe_met_isotope[0]) / (Fe56_IRMM / Fe54_IRMM) - 1.0) * 1000.0 ;
  // double DELTA_56_here = mantle_d56_here - core_d56_here;
  // cout << "mantle_d56_here: " << mantle_d56_here << endl;
  // cout << "core_d56_here: " << core_d56_here << endl;
  // cout << "DELTA_56_here: " << DELTA_56_here << endl;
  // cout << endl;

  // double mantle_d57_here = ((Fe_sil_isotope[2] / Fe_sil_isotope[0]) / (Fe57_IRMM / Fe54_IRMM) - 1.0) * 1000.0 ;
  // double core_d57_here = ((Fe_met_isotope[2] / Fe_met_isotope[0]) / (Fe57_IRMM / Fe54_IRMM) - 1.0) * 1000.0 ;
  // double DELTA_57_here = mantle_d57_here - core_d57_here;

  // cout << "mantle_d57_here: " << mantle_d57_here << endl;
  // cout << "core_d57_here: " << core_d57_here << endl;
  // cout << "DELTA_57_here: " << DELTA_57_here << endl;
  // cout << "------------ ------------------------------ ------------" << endl;
  // cout << endl;

  return;
}

double Liu_Partition(double T, double P)
{
  double F_a;
  double F_b;
  double log_beta_a;
  double log_beta_b;
  double delta;

  if (P < 26.0)
  {
    F_b = 1.29 * P + 184.1;

  }
  else if (P > 35.0)
  {
    F_b = 1.81 * P + 235.8;
  }
  else
  {
    F_b = 8.819 * P - 9.01155; // interpolation between high and low pressure
  }
  F_a = 2.03 * P + 174.2; //Fe86.8
  log_beta_b = (2904 / 1000) * F_b / pow(T, 2);
  log_beta_a = (2904 / 1000) * F_a / pow(T, 2);
  delta = 1000 * (log_beta_b - log_beta_a);
  return delta;
}




void Fe_fractionate(double (&Fe_sil_isotope)[3], double (&Fe_met_isotope)[3],
                    double equil_moles[], double equil_moles_H2O, double equil_moles_S_sil, double equil_moles_C_sil,
                    double equil_moles_oxygen, double equil_moles_H, double equil_moles_S_met, double equil_moles_C_met,
                    double P, double T) //GN fractionate pass through double Fe_isotope[]... moles[] as placeholder for now
{
 
  
  // For reference, below are the indexes for elements in elements.h
  // Element index: 0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22
  // Element name: Al Mg Ca Na Fe Si Ni Co Nb Ta  V Cr Pt Pd Ru Ir  W Mo Mn Ga Zn In Cu
  // Element index:            18 19 20 21 22 23 24 25 26 27 28 29 30 31
  // Metals start at Fe at index n then continue to n+m where m = n - 4
  // If n = 17, then only up to W is included

  //double Fe_isotope[n_emb][1][n_imps] = {0.0, 0.0, 0.0};
  //cout << "called Fe_fractionate" << endl;

  //From react():
  // moles[n] gives value for Fe in original equilibrating liquid
  // here, we want to use mseq_moles--but how does mseq_moles work? where is the array?
  // double moles_fe_old = moles[n];

  // calculate amount of Fe54 to Fe56 from total Fe

  // Inputs
  // mseq_Fe,
  // mseq_FeO
  // double mseqFe_isotope[],
  // // double T --> like w/distrib_coeff; // Get the pressure and temperature for react from the embryo mass
  // P = p_factor * pressure(mass);
  // T = temperature(P) + t_overstep;
  // double P

  // Input/Outputs &pass_by_reference
  // double (&Fe_sil_isotope)[1],
  // double (&Fe_met_isotope)[1]

  // Outputs

  // Internals
  // delta variables
  int delta_version;
  
  // 6 knowns:
  double Fe54_tot; //t54 in mathematica
  double Fe56_tot; //t56 in mathematica
  double Fe57_tot; //t57 in mathematica
  double Dfe; // DD in mathematica   --> phase partition coeff of Fe for equilibration event, renaming variable 
  double delta; // BB in mathematica --> Delta56-54Fe_{silicate - metal} // should be for Fe57 ToDo
  double theta56 = 0.678; // CC in mathematica // ratio of delta56Fe to delta57Fe
  
  //double alpha; // dependent on delta, but also equal to Fe_sil_56/Fe_sil_54 - Fe_met_56/Fe_met_54 in system of equations // still necessary? ToDo
  // double delta_prime; //check difference bwtn delta and deltaprime
  
  // Other internals:
  double Fe_tot; // total Fe
  // double AA, BB, CC; //DD ToDo
  double avar, bvar, cvar, dvar;
  double mseq_FeO = equil_moles[n - m];
  double mseq_Fe = equil_moles[n];


  fractionate_call_no++;
  no_fractionate_call_no = -1;

  // Choosing Delta: adjust which version of delta is used by changing variable delta_version
  delta_version = delta_version_Fe; // possible values: 1,2,3
  
  if (delta_version == 1) //constant value for delta
  {
    delta = 0.01;
  }
  
  else if (delta_version == 2) // delta = linear function of bridgemanite
  {
    double F_a, F_b, log_beta_a, log_beta_b;
    F_b = 208.051 + 2.51551 * P; // bridgmanite bestfit line taken from Liu et al 2017--from Rustad & Yin

    F_a = 2.03 * P + 174.2; //Fe_86.8 Ni_8.6 Si_4.6
    log_beta_b = (2904.0 / 1000.0) * F_b / pow(T, 2);
    log_beta_a = (2904.0 / 1000.0) * F_a / pow(T, 2);
    delta = (1000.0 * (log_beta_b - log_beta_a) / theta56); //note-- THE GIVEN DELTA FROM THIS LIU PAPER IS capDelta56Fe. To convert that to capDelta57, need to divide by theta56. 
  }
  
  else if (delta_version == 3) // use Liu partition, with discontinuity filled in with
  {
    delta = (Liu_Partition(T, P) / theta56); //note-- THE GIVEN DELTA FROM THIS LIU PAPER IS capDelta56Fe. To convert that to capDelta57, need to divide by theta56. 
  }

  else if (delta_version == 4) // use constant value for delta, except it's now dictated by the INPUT_36.txt parameter delta_amount
  {
    delta = delta_amount; 
  }

  else if (delta_version == 5) // use Ni and Shahar <F> for pyrolite vs <F> for hcp Fe
  {
    if (nishahar_version == 0) // Force constant
    {
      double F_a, F_b, log_beta_a, log_beta_b;
    // Equation (1) from Ni and Shahar 2022:
    // for silicate: they use the following for pyrolite
    // F_b = <F_pyrolite> = (Fe2+/Fe_net)(a1 + a2*P) + (Fe3+/Fe_net)(b1 + b2*P)
    // they calculate a1,a2,b1,b2 in the paper and I use their values below
    // using FeO and Fe as proxies for Fe2+ and Fe3+, respectively, without considering any more complex disproportionation
    // This is not correct. Fe is a proxy for Fe+0 oxidation state so I need to figure out how to approximate the Fe2+ and Fe3+ proportion
    // F_b = (mseq_FeO / (mseq_FeO + mseq_Fe)) * ((150.4) + 4.64 * P) + (mseq_Fe / (mseq_FeO + mseq_Fe)) * (284.5 - 1.13 * P);
    
    // instead (as an approximation) we can consider their lowest Fe3+ state, and that is their Pyrolite_3 with an Fe3+/sum(Fe)=0.087, which is at IW+1.4
    // this is an extreme upper bound on the IW for our model 
    F_b = (1.0 - 0.087) * ((150.4) + 4.64 * P) + (0.087) * (284.5 - 1.13 * P);

    // For metal: they use Murphy et al (2013) force constant <F> for hcp Fe as the metal reference point
    // to get the equation for that, I pulled the equation from the figure 1 in Ni and Shahar
    // <F> as a function of Pressure at 3500 K equation of line: <F> = 2.16*P + 182.38
    F_a = 2.16 * P + 182.38; //hcp pure Fe

    log_beta_b = (2904.0 / 1000.0) * F_b / pow(T, 2);
    log_beta_a = (2904.0 / 1000.0) * F_a / pow(T, 2);
    delta = 1000.0 * (log_beta_b - log_beta_a); 
    }

    else if (nishahar_version == 1) // cap-Delta57_sil-met from figure Ni and Shahar 2022, Figure 2
    {
      // delta = function of pressure from figure 2, keeping temperature constant at 3500K 
      delta =  0.0009 * P - 0.0118; // Fe3+/Total_Fe = 0.0
    }
    else if (nishahar_version == 2)
    {
      delta =  0.0005 * P - 0.0015; // Fe3+/Total_Fe = 0.2
    }
    else if (nishahar_version == 3) 
    {
      delta =  0.0001 * P + 0.0089; // Fe3+/Total_Fe = 0.4
    }
    else if (nishahar_version == 4) 
    {
      delta = -0.0004 * P + 0.0191; // Fe3+/Total_Fe = 0.6
    }
    else if (nishahar_version == 5) 
    {
      delta = -0.0012 * P + 0.0392; // Fe3+/Total_Fe = 1
    }
  }
  

  double old_delta_silicate = 1.0E3 * ((Fe_sil_isotope[2] / Fe_sil_isotope[0]) / (R57) - 1.0); //delta ratio defined by Fe57/54 //3isotope
  double old_delta_metal = 1.0E3 * ((Fe_met_isotope[2] / Fe_met_isotope[0]) / (R57) - 1.0); //delta ratio defined by Fe57/54 //3isotope

  double moles_sil = 0;
  for (int i = 0; i < n; i++)
  {
    moles_sil += equil_moles[i];
  }
  moles_sil += equil_moles_H2O;
  moles_sil += equil_moles_S_sil;
  moles_sil += equil_moles_C_sil;

  double moles_met = 0;
  for (int i = n; i < (n + m); i++)
  {
    moles_met += equil_moles[i];
  }
  moles_met += equil_moles_oxygen;
  moles_met += equil_moles_H;
  moles_met += equil_moles_S_met;
  moles_met += equil_moles_C_met;

  // pass values to 4 known values: Fe54_tot, Fe556_tot, Fe57_tot, DFe
  Fe54_tot = Fe_sil_isotope[0] * mseq_FeO + Fe_met_isotope[0] * mseq_Fe; // total 54Fe in moles //3isotope 
  Fe56_tot = Fe_sil_isotope[1] * mseq_FeO + Fe_met_isotope[1] * mseq_Fe; // total 56Fe in moles //3isotope
  Fe57_tot = Fe_sil_isotope[2] * mseq_FeO + Fe_met_isotope[2] * mseq_Fe; // total 57Fe in moles //3isotope
  

  Dfe = (mseq_Fe / moles_met) / (mseq_FeO / moles_sil); // partition coefficient as defined by mols of equilibrating Fe metal divided by equilibrating FeO silicate //3isotope

  // pass values to Other internals:
  Fe_tot = Fe54_tot + Fe56_tot + Fe57_tot; // total fe in moles //3isotope
  // -------------------------------------------------------------------------------------------------------------------------------- //
  // ------------------------------------------- Solution set  ---------------------------------------------------------------------- //
  // -------------------------------------------------------------------------------------------------------------------------------- /

  // Constant avar, bvar, cvar, dvar:
  // mathematica variable --> cpp variable:
  //   t54 given by--> Fe54_tot
  //   t56 given by--> Fe56_tot
  //   t57 given by--> Fe57_tot
  //    DD given by--> Dfe 
  //    BB given by--> delta
  //    CC given by--> theta56 
  //   R56 defined as global variable
  //   R57 defined as global variable
  // TotEC stands for "Total Equilibrating Core", given by--> moles_met (total equilibrating core metal moles, for clarity, note that:  mseq_Fe are just for Fe)
  // TotEM stands for "Total Equilibrating Mantle", given by--> moles_sil (total equilibrating mantle silicate moles, for clarity, note that: mseq_FeO are just for Fe)
  
  //variable conversions:
  const double t54 = Fe54_tot; 
  const double t56 = Fe56_tot; 
  const double t57 = Fe57_tot; 
  const double DD = Dfe; // for DD here */;  
  const double BB = delta;/* Set your value for BB here */;   
  const double CC = theta56; /* Set your value for CC here */;   
  const double TotEM = moles_met; 
  const double TotEC = moles_sil; 

  double avarNumer = sqrt(((DD  * TotEC + TotEM) * (-4000. * BB * pow(R57,2) * (CC * R56 + R57) * t54 * (t54 + t56 + t57) * TotEM + pow((R57 * (1000. + (1000. + (-1000. + BB) * CC) * R56 +  BB * R57) * t54 + 1000. * (CC * R56 + R57) * t57),2) * (DD * TotEC +  TotEM))));
  double avarDenom = 2. * BB * pow((CC * R56 + R57),2) * (DD * TotEC + TotEM);
  avar = avarNumer/avarDenom;
  
  bvar = (500. / (BB * (CC * R56 + R57))) * ((R57 * CC * R56 * t54 - R57 * t54 - R57 * R56 * t54 - t57 * (CC * R56 + R57)) / R57);
    
  cvar = (-pow(R57, 2) * (1000. * pow(R56, 2) * t54 + R56 * (1000. + BB * R57) * t54 - 2. * BB * R57 * t56 + 1000. * R56 * t57) * (DD * TotEC + TotEM) + CC * R56 * R57 * (DD * (R56 * (1000. - (-2000. + BB) * R57) * t54 + (1. + R57) * (1000. + BB * R57) * t54 + 4. * BB * R57 * t56 + 1000. * (1. - R56 + R57) * t57) * TotEC + ((1000. * (1. + R56) - (-1000. + BB + (-2000. + BB) * R56) * R57 + BB * pow(R57, 2)) * t54 + 2. * BB * R57 * (t56 - t57) - 1000. * (-1. + R56 - R57) * t57) * TotEM) + pow(CC, 2) * pow(R56, 2) * (DD * ((-1000. + BB) * R57 * (1. + R57) * t54 + 2. * BB * R57 * t56 + 1000. * (1. + R57) * t57) * TotEC + (R57 * (BB * (-1. + R57) - 1000. * (1. + R57)) * t54 + 2. * (500. - (-500. + BB) * R57) * t57) * TotEM)) / (2. * BB * R57 * pow(CC * R56 + R57, 2) * (DD * TotEC + TotEM));
    
  dvar = 1. / (2. * BB * pow(CC * R56 + R57, 2) * (DD * TotEC + TotEM)) * (DD * (-(-1. + (-1. + CC) * R56) * R57 * (1000. + (1000. + (-1000. + BB) * CC) * R56 + BB * R57) * t54 + 2. * (CC * R56 + R57) * (500. + (500. + (-500. + BB) * CC) * R56 + BB * R57) * t57) * TotEC + (-R57 * (-1000. + (-1. + CC) * (1000. + (-1000. + BB) * CC) * pow(R56, 2) + BB * R57 + R56 * (-2000. + (2000. + BB) * CC + BB * (-1. + CC) * R57)) * t54 - 2. * BB * R57 * (CC * R56 + R57) * t56 + 2. * (500. + (500. + (-500. + BB) * CC) * R56) * (CC * R56 + R57) * t57) * TotEM);
    
  // Use solution set 2 since it produces physical answers: 
  // -------------------------------------------------------------------------------------------------------------------------------- //
  // Calculations, from simplified-3isotope.nb, and from solving-6-equations.nb mathematica Solve[]
  double m54 = t54 / 2. + bvar + avar * ((CC * R56 + R57) / R57);
  double s54 = t54 / 2. - bvar - avar * ((CC * R56 + R57) / R57);

  double m56 = cvar + avar * (R56 - CC * R56 - CC * R56 / R57);
  double s56 = t56 - cvar + avar * (-R56 + CC * R56 + (CC * R56) / R57);

  double m57 = dvar + (CC * R56 - 1. - R56) * avar;
  double s57 = t57 - dvar + (R56 + 1. - CC * R56) * avar;

  // -------------------------------------------------------------------------------------------------------------------------------- //
  // NOW HAVE TO NORMALIZE SO THE total of Fe_met[0,1,2] = 1 and Fe_sil[0,1,2]=1
  Fe_met_isotope[0] = m54 / (m54 + m56 + m57); // update 54Fe value, in fractional abundance, not moles -- send n_emb and n_imp to Fe_isotope
  Fe_met_isotope[1] = m56 / (m54 + m56 + m57); // update 56Fe value -- send n_emb and n_imp to Fe_isotope //3isotope
  Fe_met_isotope[2] = m57 / (m54 + m56 + m57); // update 57Fe value -- send n_emb and n_imp to Fe_isotope //3isotope
  
  Fe_sil_isotope[0] = s54 / (s54 + s56 + s57);
  Fe_sil_isotope[1] = s56 / (s54 + s56 + s57);//3isotope
  Fe_sil_isotope[2] = s57 / (s54 + s56 + s57);//3isotope
  

  double new_delta56_silicate = 1.0E3 * ((Fe_sil_isotope[1] / Fe_sil_isotope[0]) / (R56) - 1.0); //3isotope
  double new_delta56_metal = 1.0E3 * ((Fe_met_isotope[1] / Fe_met_isotope[0]) / (R56) - 1.0); //3isotope

  double new_delta57_silicate = 1.0E3 * ((Fe_sil_isotope[2] / Fe_sil_isotope[0]) / (R57) - 1.0); //3isotope
  double new_delta57_metal = 1.0E3 * ((Fe_met_isotope[2] / Fe_met_isotope[0]) / (R57) - 1.0); //3isotope

  // if fe_sil_isotope54 and fesilisotope56 dont add to 1, print out: 54sil 54met 56sil 56met
  if (abs(Fe_sil_isotope[0] + Fe_sil_isotope[1] + Fe_sil_isotope[2] - 1.0 ) > 0.0001)
  {
    cout << "darn" << endl;
    cout << Fe_sil_isotope[0] << " " << Fe_met_isotope[0] << " " << Fe_sil_isotope[1] << " " << Fe_met_isotope[1] << " " << fractionate_call_no << endl;
    }
  if (no_fractionate_call_no == fractionate_call_no) exit(0);

     
  return;

}



std::string format(double f, int n_digits)
{
  std::stringstream ss;
  if (f == 0) {
    ss << std::fixed << std::setprecision(n_digits - 1) << 0.0;
    return ss.str();
  }
  int n_decimal = int(::ceil(::log10(f < 0 ? -f : f))); /*digits before decimal point*/
  double order = ::pow(10., 1.0 * (n_digits - n_decimal));
  if (n_decimal > 0) {
    ss << std::fixed << std::setprecision(std::max(n_digits - n_decimal, 0)) << ceil(f * order + 0.5) / order;
  }
  else {
    order = ::pow(10., 1.0 * n_digits);
    ss << std::fixed << std::setprecision(n_digits - 1) << ceil(f * order + 0.5) / order;
  }
  return ss.str();
}
