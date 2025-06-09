//List of element details used in calculations. The constant "n" declared as a global constant in main program must not exceed the
// total number of elements listed here.

using namespace std;

// Ni-Co partitioning from Kegler et al.
// HSE parrameters from Vera's corrected fits of 10th September 2015

// LIST LITHOPHILES FIRST; THEN FE, SI, NI (in that order), THEN TRACE SIDEROPHILES; Ta should be listed immediately after Nb in order to obtain correct Nb/Ta ratio
// atomic weights and oxide molecular weights defined for the N elements:
// Formula unit XO(val/2) is used for all oxides (eg NbO(2.5))

//Element index #:                             0         1         2         3         4         5         6         7         8         9        10        11        12        13        14        15        16        17        18        19        20        21        22
//                                            Al        Mg        Ca        Na        Fe        Si        Ni        Co        Nb        Ta         V        Cr        Pt        Pd        Ru        Ir         W        Mo        Mn        Ga        Zn        In        Cu
extern const double AW[] =            {    26.98,    24.31,    40.08,    22.99,    55.85,    28.09,    58.71,    58.93,    92.91,   180.95,    50.94,     52.0,   195.08,   106.42,   101.07,    192.2,   183.85,    95.94,    54.94,    69.72,    65.37,   114.82,    63.55};
extern const double MW[] =            {    50.98,    40.32,    56.08,    30.99,    71.85,    60.09,    74.71,    74.93,   132.91,   220.95,    74.94,     68.0,   211.08,   122.42,   117.07,    208.2,   231.85,   127.94,    70.94,    85.72,    81.37,   122.82,    71.55};

// Partitioning (log Kd=a+b/T+cP/T+log(1-Xs)) parameters and valence (incomplete W values from Wade et al 2012, Eq.22):
// Added low pressure a, b, c, sig_a, sig_b, sig_c values for Ni, Co for P < 5 GPa
// testing estimate low pressure vals for Ta for lowp_a, lowp_b, lowp_c
extern const double a[] =             {      0.0,      0.0,      0.0,      0.0,      0.0,     2.98,     1.06,   0.1261,     2.66,     0.84,    -0.48,    -0.17, -3.80381,  1.31537,  0.83869, 0.231867,     0.61,    0.611,     0.54,     0.52,     -0.3,     1.03,     0.68};
extern const double sig_a[] =         {      0.0,      0.0,      0.0,      0.0,      0.0,     0.02,      0.4,     0.38,     0.11,     0.09,     0.38,     0.13, 0.057894, 0.043893, 0.215401, 0.242633,     0.28,    0.303,     0.05,     0.45,     0.36,      0.3,      0.0};
extern const double lowp_a[] =        {      0.0,      0.0,      0.0,      0.0,      0.0,     2.98,    0.414,   0.3798,     2.66,     0.84,    -0.48,    -0.17, -3.80381,  1.31537,  0.83869, 0.231867,     0.61,    0.611,     0.54,     0.52,     -0.3,     1.03,     0.68};
extern const double lowp_sig_a[] =    {      0.0,      0.0,      0.0,      0.0,      0.0,     0.02,     0.28,     0.24,     0.11,     0.09,     0.38,     0.13, 0.057894, 0.043893, 0.215401, 0.242633,     0.28,    0.303,     0.05,     0.45,     0.36,      0.3,      0.0};

extern const double b[] =             {      0.0,      0.0,      0.0,      0.0,      0.0, -15934.0,   1553.1,   2057.8, -14032.0, -13806.0,  -5063.0,  -2730.0,  22589.3,  7843.92,  12396.7,  17253.4,  -0.4091,      0.0,  -7142.0,   -514.0,  -1574.0,    302.0,   2300.0};
extern const double sig_b[] =         {      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,   921.05,   823.25,      0.0,      0.0,    894.0,    312.0,  54.8674,  36.9069,  31.9987,  54.9997,    670.0,      0.0,      0.0,    914.0,    739.0,    621.0,      0.0};
extern const double lowp_b[] =        {      0.0,      0.0,      0.0,      0.0,      0.0, -15934.0,   3744.7,   1501.0, -14032.0, -13806.0,  -5063.0,  -2730.0,  22589.3,  7843.92,  12396.7,  17253.4,  -0.4091,      0.0,  -7142.0,   -514.0,  -1574.0,    302.0,   2300.0};
extern const double lowp_sig_b[] =    {      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,   470.63,   400.94,      0.0,      0.0,    894.0,    312.0,  54.8674,  36.9069,  31.9987,  54.9997,    670.0,      0.0,      0.0,    914.0,    739.0,    621.0,      0.0};

extern const double c[] =             {      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,    -97.9,   -57.48,   -199.0,   -101.0,      0.0,      0.0,  138.477,  -124.42,  20.2824, -44.4699,      0.0,      0.0,     52.0,    -78.0,    -51.0,   -106.0,    -37.0};
extern const double sig_c[] =         {      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,     7.66,     9.53,     16.0,     15.0,      0.0,      0.0,  10.9322,  7.91161,    33.74,  50.9154,      0.0,      0.0,      9.0,     13.0,     12.0,     10.0,      0.0};
extern const double lowp_c[] =        {      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,   -271.7,   -16.75,   -199.0,   -101.0,      0.0,      0.0,  138.477,  -124.42,  20.2824, -44.4699,      0.0,      0.0,     52.0,    -78.0,    -51.0,   -106.0,    -37.0};
extern const double lowp_sig_c[] =    {      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,    31.88,    30.93,     16.0,     15.0,      0.0,      0.0,  10.9322,  7.91161,    33.74,  50.9154,      0.0,      0.0,      9.0,     13.0,     12.0,     10.0,      0.0};

extern const double d[] =             {      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,  14.4292,  9.42348,  16.3602,  24.6429,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0};
extern const double sig_d[] =         {      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,  0.25259, 0.182349,  1.02227,  3.17365,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0};

// Partitioning (log Kd=a+b/T+cP/T) parameters for FeS - silicate liquid partitioning of HSEs
extern const double sulf_a[] =        {      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,   2.4318,  2.74257, -1.33538, -1.28861,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0};
extern const double sulf_sig_a[] =    {      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0, 0.189945, 0.227032, 0.578663, 0.645279,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0};

extern const double sulf_b[] =        {      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,  3305.93,  1686.99,  13315.5,  14058.1,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0};
extern const double sulf_sig_b[] =    {      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,  520.551,   585.43,  1430.39,   1797.9,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0};

extern const double sulf_c[] =        {      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,  71.1881, -18.8536,  -48.935, -10.2483,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0};
extern const double sulf_sig_c[] =    {      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,  5.57156,  7.35806,  6.60265,  39.5069,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0};

//VALENCE:5
extern const int val[] =              {        3,        2,        2,        1,        2,        4,        2,        2,        5,        5,        3,        2,        1,        1,        2,        2,        6,        4,        2,         2,       2,        1,        1};

// Following two arrays specify output concentrations (0=wt%, 1=ppm, 2=ppb) in mantle and core files respectively:
extern const int sil_c[] =            {        0,        0,        0,        0,        0,        0,        1,        1,        2,        2,        1,        1,        2,        2,        2,        2,        2,        2,        1,         1,       1,        2,        1};
extern const int met_c[] =            {       -1,       -1,       -1,       -1,        0,        0,        0,        0,        2,        2,        1,        0,        2,        2,        2,        2,        2,        2,        1,         1,       1,        2,        1};

//Element index #:                             0         1         2         3         4         5         6         7         8         9        10        11        12        13        14        15        16        17        18        19        20        21        22
//                                            Al        Mg        Ca        Na        Fe        Si        Ni        Co        Nb        Ta         V        Cr        Pt        Pd        Ru        Ir         W        Mo        Mn        Ga        Zn        In        Cu
// Following array specifies the elements to be included in least square calculation (y=included, n=not included)
// fit_mse[] used in Rubie 2015 =
extern const char fit_mse[] =         {      'n',      'n',      'n',      'n',      'y',      'y',      'y',      'y',      'y',      'y',      'y',      'y',      'n',      'n',      'n',      'n',      'n',      'n',      'n',       'n',     'n',      'n',      'n'};
extern const char fit_na[] =          {      'n',      'n',      'n',      'y',      'n',      'n',      'n',      'n',      'n',      'n',      'n',      'n',      'n',      'n',      'n',      'n',      'n',      'n',      'n',       'n',     'n',      'n',      'n'};
extern const char fit_hse[] =         {      'n',      'y',      'n',      'n',      'y',      'y',      'y',      'y',      'y',      'y',      'y',      'y',      'y',      'y',      'y',      'y',      'n',      'n',      'n',       'n',     'n',      'n',      'n'};

// fit_mse[] used for my fitmars runs:
// extern const char fit_mse[] =         {      'n',      'y',      'n',      'n',      'y',      'y',      'y',      'y',      'y',      'n',      'n',      'y',      'n',      'n',      'n',      'n',      'n',      'n',      'n',       'n',     'n',      'n',      'n'};

///major elements of mars el 0-5 (no Na)

// Following array specifies type of element
// (lithophile: 'l'; siderophile: 's'; volatile: 'v'
// highly siderophile: 'h'; highly siderophile and volatile: 'H'
extern const char el_type[] =         {      'l',      'l',      'l',      'v',      's',      's',      's',      's',      's',      's',      's',      's',      'h',      'h',      'h',       'h',     'w',      'w',      'v',       'v',     'v',      'v',      'v'};

// Oxygen atomic weight:
extern const double AW_O = 16.0;
// Water molecular weight
extern const double MW_H2O = 18.0;
extern const double MW_S2O = 64.06;
extern const double MW_CO = 28.01;
extern const double AW_S = 32.065;      //Sulfur
extern const double AW_C = 12.0107;    // carbon
extern const double AW_H = 1.0;         //Hydrogen

//Element index #:                             0         1         2         3         4         5         6         7         8         9        10        11        12        13        14        15        16        17        18        19        20        21        22
// Earth mantle oxide/element concentrations  Al        Mg        Ca        Na        Fe        Si        Ni        Co        Nb        Ta         V        Cr        Pt        Pd        Ru        Ir         W        Mo        Mn        Ga        Zn        In        Cu
// WT%, PPM OR PPB (Palme & O'Neill 2014):    NOTE: need to label, but these values don't quite match the labels given in Palme ONeill 2014, need to find out conversion btwn ppm and wt% and ppb
extern const double EM[] =            {     4.49,    36.77,     3.65,    0.259,      8.1,     45.4,   1860.0,    102.0,    595.0,     43.0,     86.0,   2520.0,      7.6,      7.1,      7.4,      3.5,     12.0,     47.0,   1050.0,      4.4,     53.5,     18.0,     20.0};

// Earth mantle oxide/element concentrations, ALL as WT%, with FRACTIONAL uncertainties(from Palme & O'Neill 2013):
extern const double p_mantle[] =      {     4.49,    36.77,     3.65,    0.259,      8.1,     45.4,    0.186,   0.0102, 595.0e-7,  43.0e-7,  86.0e-4,    0.252,   7.6e-7,   7.1e-7,   7.4e-7,   3.5e-7,  12.0e-7,  47.0e-7,    0.105,  44.0e-5, 535.0e-5,  18.0e-7,   2.0e-3};
// Fractional uncertainties on mantle concentrations:
extern const double sig_p_mantle[] =  {     0.08,     0.01,     0.08,     0.05,  0.00617,    0.007,     0.05,     0.05,      0.2,     0.05,     0.05,      0.1,      0.2,      0.2,      0.2,      0.1,      0.3,      0.4,      0.1,     0.05,     0.05,      0.2,      0.5};

//Element index #:                             0         1         2         3         4         5         6         7         8         9        10        11        12        13        14        15        16        17        18        19        20        21        22
//                                            Al        Mg        Ca        Na        Fe        Si        Ni        Co        Nb        Ta         V        Cr        Pt        Pd        Ru        Ir         W        Mo        Mn        Ga        Zn        In        Cu
//Mars mantle concentrations from Taylor 2013, table 4 -- V is calculated from 1.9x CI valule from Lodders (2003), Mo is from Yoshizaki (2020)
extern const double MM[] =      	  {     2.61,     30.68,    2.43,      0.4,    18.14,    43.85,     330.,      71.,     501.,      27.2,      105.,     4990.,    3.1,      2.4,      2.6,      2.0,       74.,     0.5,     0.34,     6.6,     18.9,      6.9,      2.};

// Mars mantle oxide/element concentrations, ALL as WT%, with FRACTIONAL uncertainties(from Table 4 (for oxides 0 thru 5) and Table 5 for elements 6 thru 22 from Taylor 2013):
extern const double m_mantle[] =      {     2.61,     30.68,    2.43,      0.4,    18.14,    43.85,  330.e-4,   71.e-4,  501.e-7,   27.2e-7, 105.83e-4,  4990.e-4, 3.1e-7,   2.4e-7,   2.6e-7,   2.0e-7,    74.e-7,  0.5e-4,     0.34,  6.6e-4,   18.9e-4,  6.9e-7,   2.e-4};

// Fractional uncertainties on mantle concentrations, calculated values from Table 4 Taylor 2013: 
//extern const double sig_m_mantle[] =  {    0.021,     0.016,   0.011,      0.1,    0.028,    0.024,     0.17,     0.18,     0.02,      0.02,      0.02,      0.04,   0.13,     0.17,     0.17,     0.25,      0.21,     0.8,     0.04,    0.06,      0.08,    0.16,    0.17};
// to test using same values as Earth uncertainty use this uncertainty array:
//extern const double sig_m_mantle[] =  {     0.08,     0.01,     0.08,     0.05,  0.00617,    0.007,     0.05,     0.05,      0.2,     0.05,     0.05,      0.1,      0.2,      0.2,      0.2,      0.1,      0.3,      0.4,      0.1,     0.05,     0.05,      0.2,      0.5};
//highest uncertainty between taylor 2013 values and earth values:
extern const double sig_m_mantle[] = {      0.08,    0.016,     0.08,      0.1,    0.028,    0.024,     0.17,     0.18,      0.2,      0.05,     0.05,      0.1,     0.2,       0.2,     0.2,       0.25,      0.3,     0.8,     0.1,     0.06,     0.08,      0.2,      0.5};


//Element index #:                             0         1         2         3         4         5         6         7         8         9        10        11        12        13        14        15        16        17        18        19        20        21        22
//                                            Al        Mg        Ca        Na        Fe        Si        Ni        Co        Nb        Ta         V        Cr        Pt        Pd        Ru        Ir         W        Mo        Mn        Ga        Zn        In        Cu
//Mars mantle concentrations from Yoshizaki 2020 (In is from Taylor 2013)
//extern const double MM[] =      	  {     3.03,     31.01,    2.88,     1.18,    14.67,    45.56,     360.,      96.,     640.,      33.9,     123.,     6000.,    9.6,      7.4,      6.2,      3.7,       69.,     0.5,     0.288,      8.7,     45.,      6.9,      2.6};
// Yoshizaki 2020 Mars mantle oxide/element concentrations, ALL as WT%:
//extern const double m_mantle[] =      {     3.03,     31.01,    2.88,     1.18,    14.67,    45.56,  360.e-4,    96.e-4,   640.e-7,   33.9e-7,  123.e-4, 6000.e-4, 9.6e-7, 7.40e-7, 6.20e-7,   3.70e-7,   69.e-7,   0.5e-4,     0.288,   8.7e-4,   45.e-4,   6.9e-7,  2.6e-4};
//highest uncertainty between Yoshizaki 2020 and earth values:
//extern const double sig_m_mantle[] = {      0.1,       0.06,     0.1,     0.23,     0.07,     0.04,     0.26,      0.46,     0.2,      0.14,     0.1,       0.17,     0.8,     0.7,      0.6,      0.8,       0.54,     0.8,     0.18,     0.22,     0.33,      0.2,    0.5};



// NOTE Kd(W) & Kd(Mo) parameters listed above are from Wade et al (GCA, 2012, Eqs.22 & 16) for a "MORB-pyrolite" composition. Kd(W) & Kd(Mo) are recalculated
// in function "REACT" based on actual silicate liquid compositions



//Wade & Wood:
//Ni   0.75   3097.0   -123.0
//Co   0.17   2511.0   -45.0
//Kegler
//Ni   1.06   1553.0   -98.0
//Co   0.13   2057.0   -57.0
//(   a      b      c   )
//Siebert (2011)
//Ni   0.35   2934.0   -83.0
//Co   0.20   1892.0   -66.0

//Rose et al
//Ni
