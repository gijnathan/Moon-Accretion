// FUNCTION input_parameters - Defines an array of imput parameters within core-formation function "func"
#include <iostream>
#include <iomanip>
#include <fstream>
#include "nr.h"
using namespace std;


void fit_mse_parameters(double (&param)[8], int dist_comp, int fit, double p_factor, double dist_1, double dist_2, double xsi_met_0, double xfe_met_2, double xfe_met_3, double erf_width_1, double erf_width_2, Vec_I_DP &x);
void fit_hse_parameters(double (&param)[4], double dist_0_S, double dist_0_C, double matte_P_factor, double time_magma, Vec_I_DP &x);

void fit_mse_parameters(double (&param)[8], int dist_comp, int fit, double p_factor, double dist_1, double dist_2, double xsi_met_0, double xfe_met_2, double xfe_met_3, double erf_width_1, double erf_width_2, Vec_I_DP &x)
{
   param[0] = p_factor;
   param[1] = dist_1;
   param[2] = dist_2;
   param[3] = xsi_met_0;
   param[4] = xfe_met_2;
   param[5] = xfe_met_3;
   param[6] = erf_width_1;
   param[7] = erf_width_2;

   if (fit == 1)
// List the required fit parameters (above-listed values are overwritten).
//Must correspond to values listed in "input.txt" file
   {
      if (dist_comp == 1)     // CONSTANT
      {
         param[0] = x[0];  // p_factor
         param[4] = x[1];  // xfe_met_2
      }
      else if (dist_comp == 2)   // STEP
      {
         param[0] = x[0];  // p_factor
         param[1] = x[1];  // dist_1
         param[2] = x[2];  // dist_2
         param[3] = x[3];  // xsi_met_0
         param[4] = x[4];  // xfe_met_2
         param[5] = x[5];  // xfe_met_3
      }
      else if (dist_comp == 3)   // GRADIENT
      {
         param[0] = x[0];  // p_factor
         param[1] = x[1];  // dist_1
         param[2] = x[2];  // dist_2
         param[3] = x[3];  // xsi_met_0
         param[4] = x[4];  // xfe_met_2
      }
      else if (dist_comp == 4)   // ERF
      {
         param[0] = x[0];  // p_factor
         param[1] = x[1];  // dist_1
         param[2] = x[2];  // dist_2
         param[3] = x[3];  // xsi_met_0
         param[6] = x[4];  // erf_width_1
         param[7] = x[5];  // erf_width_2
      }
   }
   return;
}

void fit_hse_parameters(double (&param)[4], double dist_0_S, double dist_0_C, double matte_P_factor, double time_magma, Vec_I_DP &x)
{
   param[0] = dist_0_S;
   param[1] = dist_0_C;
   param[2] = matte_P_factor;
   param[3] = time_magma;

   param[0] = x[0];
   param[1] = x[1];
   param[2] = x[2];
   param[3] = x[3];

   return;
}


/*
double* dist_comp_parameters(int vol_grad, int fit_vol_grad, double dist_0_S, double dist_0_C, double matte_P_factor, double time_magma, Vec_I_DP &x)
{
   double* param = new double[13];
   if (!param)
   {
      cout << "Memory allocation failure in INPUT_PARAMETERS.";
      exit(1);
   }

   param[0] = matte_P_factor;
   param[1] = dist_0_S;
   param[2] = dist_0_C;
   param[3] = time_magma;

   if (fit_vol_grad == 1)  // List the required fit parameters (above-listed values are overwritten).
      //Must correspond to values listed in "input.txt" file
   {
      if (vol_grad == 1)     // CONSTANT
      {
         param[0] = x[0];  // p_factor
         param[3] = x[3];  // xfe_met_2
      }
      else if (vol_grad == 2)   // STEP
      {
         param[0] = x[0];  // p_factor
         param[1] = x[1];  // dist_1
         param[6] = x[2];  // xfe_met_1
         param[3] = x[3];  // xfe_met_2
         param[4] = x[4];  // xsi_met_0
      }
      else if (vol_grad == 3)   // GRADIENT
      {
         param[0] = x[0];  // p_factor
         param[1] = x[1];  // dist_1
         param[2] = x[2];  // dist_2
         param[4] = x[3];  // xsi_met_0
         param[3] = x[4];  // xfe_met_2
      }
      else if (vol_grad == 4)   // GRADIENT with S and HSEs
      {
         param[1] = x[0];  // dist_1
         param[2] = x[1];  // dist_2
         param[5] = x[2];  // dist_3
         param[11] = x[3];   // dist_0_S
         param[4] = x[4];  // xsi_met_0
         param[3] = x[5];  // xfe_met_2
         param[0] = x[6];  // p_factor
         param[9] = x[7];    // matte_P_factor
         param[10] = x[8];   // time_magma
      }
      else if (vol_grad == 5)   // GRADIENT with S, C, HSEs
      {
         param[1] = x[0];  // dist_1
         param[2] = x[1];  // dist_2
         param[5] = x[2];  // dist_3
         param[11] = x[3];   // dist_0_S
         param[12] = x[4];   // dist_0_C
         param[4] = x[5];  // xsi_met_0
         param[3] = x[6];  // xfe_met_2
         param[0] = x[7];  // p_factor
         param[9] = x[8];    // matte_P_factor
         param[10] = x[9];   // time_magma
      }
   }
   return param;
}
*/