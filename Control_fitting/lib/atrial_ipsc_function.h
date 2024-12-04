/****
******* Ni-Grandi atrial-like iPSC-CM model
*****/

// 17:07:35, Mon, 02-March-2020, By Haibo
// adding IKur from my new atrial cell model

// 17:45:17, Thu, 12-December-2019, By Haibo
// adding cvode



// To use the original Kernik et al. Model, 
// undefine USE_V1
// set CaL_v_shift=0;




// Kernik-Clancy iPSC-CM model
//**********************************************
//Kernik DC, Morotti S, Wu H, Garg P, Duff HJ, Kurokawa J, Jalife J, Wu JC, Grandi E, Clancy CE.
//A computational model of induced pluripotent stem-cell derived cardiomyocytes
//incorporating experimental variability from multiple data sources"
//J Physiol. 2019 Jul 6. doi: 10.1113/JP277724
//**********************************************
//
//Converted to C-code by Mao-Tsuen Jeng
//
//Colleen Clancy Lab @ UC davis
//
//May-21-2019



#ifndef ATRILA_IPSC_FUNCTION__H
#define ATRILA_IPSC_FUNCTION__H
#include <cmath>
#include <iostream>
// void atrial_ipsc_function( double t, double *Y, double *model_parameter_inputs, double *dY, double *currents );
void atrial_ipsc_function( double t, double *Y, double *model_parameter_inputs, double *dY, double *currents );


#endif