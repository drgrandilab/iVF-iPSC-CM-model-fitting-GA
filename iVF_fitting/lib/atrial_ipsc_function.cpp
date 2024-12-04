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

#include "atrial_ipsc_function.h"

// #define USE_V1

// void atrial_ipsc_function( double t, double *Y, double *model_parameter_inputs, double *dY, double *currents );
void atrial_ipsc_function( double t, double *Y, double *model_parameter_inputs, double *dY, double *currents ) {


	//% differential equations for Kernik-Clancy iPSC-CM model
	//
	//%% State variable definitions:
	//% 1: Vm (millivolt)
	//
	//% Ionic Flux: -------------------------------------------------------------
	//% 2: Ca_SR (millimolar)
	//% 3: Cai (millimolar)
	//% 4: Nai (millimolar)
	//% 5: Ki (millimolar)
	//% 6: Ca_ligand (millimolar)
	//
	//% Current Gating (dimensionless):------------------------------------------
	//% 7: d     (activation in i_CaL)
	//% 8: f1    (inactivation in i_CaL)
	//% 9: fCa   (calcium-dependent inactivation in i_CaL)
	//% 10: Xr1   (activation in i_Kr)
	//% 11: Xr2  (inactivation in i_Kr
	//% 12: Xs   (activation in i_Ks)
	//% 13: h    (inactivation in i_Na)
	//% 14: j    (slow inactivation in i_Na)
	//% 15: m    (activation in i_Na)
	//% 16: Xf   (inactivation in i_f)
	//% 17: s    (inactivation in i_to)
	//% 18: r    (activation in i_to)
	//% 19: dCaT (activation in i_CaT)
	//% 20: fCaT (inactivation in i_CaT)
	//% 21: R (in Irel)
	//% 22: O (in Irel)
	//% 23: I (in Irel)
	// 17:10:03, Mon, 02-March-2020, By Haibo
	//% 24: IKur - a
	//% 25: IKur - i1 fast
	//% 26: IKur - i2 slow inactivation
	// % 27: o_k2p - IK2P, opening
	// % 28: o_SK - ISK, opening
	//% -------------------------------------------------------------------------------


	//%% Parameter inputs:

	//%Current parameter values:

	double Gto_atrial = 1;
	double GK1_atrial = 1;//0.75/1.1;
	double GNaca_atrial = 1;//0.5;
	double GNaK_atrial = 1;//0.8;
	double CaL_v_shift = 0;//3.114;
	double GCat_atrial = 1;
	double Gserca = 1;
	double Kmf_Scale = 1;
	double G_CaL_atrial = 1;

	double Gkur_atrial = 0;
	double GSK_atrial = 0;
	double GK2P_atrial = 0;
	double GKAch_atrial = 0;
	double Scale_volume = 1;
	double GKr_atrial = 1;
	double GNa_atrial = 1;
	double Gf_atrial = 1.0;

	double GbNa_atrial = 1.0;
	double GPCa_atrial = 1.0;
	double GbCa_atrial = 1.0;

	double Buffer_Scale = 1.0;
	double ec50SR_Scale = 1.0;
	// G_CaL_atrial = 0*10.5;
	// GCat_atrial = 0*1;

	const bool atrial_like = 0;

	// if atrial-like
	// from Paci 2013


	double * x_scale_conductance = &model_parameter_inputs[0]; // (1:16);
	double * x_K1 = &model_parameter_inputs[16]; // (17:22);
	double * x_KR = &model_parameter_inputs[22]; // (23:33);
	double * x_IKS = &model_parameter_inputs[33]; // (34:39);
	double * xTO = &model_parameter_inputs[39]; // (40:50);
	double * x_cal = &model_parameter_inputs[50]; // (51:61);
	double * x_cat = &model_parameter_inputs[61]; // (62);
	double * x_NA = &model_parameter_inputs[62]; // (63:76);
	double * x_F = &model_parameter_inputs[76]; // (77:82);


	//    Flags:
	double stim_flag = model_parameter_inputs[82]; // (83);  // % dimensionless (in stim_mode)
	// double stim_flag = 0;
	double voltageclamp = model_parameter_inputs[84]; // (85); // %square pulses if =1


	double * scale_para = &model_parameter_inputs[87];

	// debug
#ifdef DEBUG
	for (int i = 0; i < 16; ++i)
	{
		std::cout << scale_para[i] << std::endl;
	}
#endif

	double Ina_shift = 0;
	double RyR_scale = 1.;



	// #ifdef USE_V1

	// std::cout << "test";
	if (atrial_like) {
		Scale_volume = /*scale_para[12] * 1;//*/1;

		Gto_atrial   = scale_para[0] * 1;//1.2*2;
		GK1_atrial   = scale_para[1] * 1;//0.8* 0.9;// / 1.1;
		GNaca_atrial = scale_para[2] * 1;//0.8;
		Gserca       = scale_para[3] * 1;//1;
		GNaK_atrial  = scale_para[4] * 1;//1;
		CaL_v_shift  = 5*(scale_para[5]-1);//3.114;
		G_CaL_atrial = scale_para[6] * 1;//0.5;
		GCat_atrial  = scale_para[7] * 1;//1.5;
		Gkur_atrial  = 0* scale_para[8] * 1;
		GSK_atrial   = 0* scale_para[9] * 0.5;
		GK2P_atrial  = 0* scale_para[10] * 1;
		GKAch_atrial = 0*scale_para[11] * 0.5;
		GbNa_atrial = scale_para[12];
		GNa_atrial   = scale_para[13] * 1;//1.1;
		Gf_atrial    = scale_para[14] * 1.0;
		GKr_atrial   = scale_para[15] * 1;
		GPCa_atrial = scale_para[16];
		GbCa_atrial = scale_para[17];

		Ina_shift = 5*(scale_para[11]-1);

		// RyR_scale = scale_para[2];

	} else {


		Gto_atrial   = scale_para[0] * 1;//1.2*2;
		GK1_atrial   = scale_para[1] * 1;//0.8* 0.9;// / 1.1;
		GNaca_atrial = scale_para[2] * 1;//0.8;
		Gserca       = scale_para[3] * 1;//1;
		GNaK_atrial  = scale_para[4] * 1;//1;
		CaL_v_shift  = 5*(scale_para[5]-1);//3.114;  // shift is zero at default (para=1)
		G_CaL_atrial = scale_para[6] * 1;//0.5;
		GCat_atrial  = scale_para[7] * 1;//1.5;

		Buffer_Scale = scale_para[8];  // SR BUFFER
		ec50SR_Scale = scale_para[9];
		Kmf_Scale = scale_para[10];   // not used if USE_V1 not defined
		Ina_shift =  5*(scale_para[11]-1) ;   // shift is zero at default (para=1) // not used if USE_V1 not defined
		// Gkur_atrial  = 0* scale_para[8] * 1;
		// GSK_atrial   = 0* scale_para[9] * 0.5;
		// GK2P_atrial  = 0* scale_para[10] * 1;
		// GKAch_atrial = 0*scale_para[11] * 0.5;
		GbNa_atrial = scale_para[12];
		GNa_atrial   = scale_para[13] * 1;//1.1;
		Gf_atrial    = scale_para[14] * 1;
		GKr_atrial   = scale_para[15] * 1;
		GPCa_atrial = scale_para[16];
		GbCa_atrial = scale_para[17];


	}

	// #endif
	// to simulate orignal model, set CaL_v_shift to zero
	// CaL_v_shift = 0;

	// -------------------------------------------------------------------------------
	//  Constants for flag protocols:

	// for stim:
	double cyclelength = 1852;          // 1000ms = 1hz beating
	double i_stim_Amplitude = 3.;       // pA/pF (in stim_mode)
	double i_stim_End = 1000e3;           // milisecond (in stim_mode)
	double i_stim_PulseDuration = 5.;   // milisecond (in stim_mode)
	double i_stim_Start = 0;            // milisecond (in stim_mode)

	// for square pulse voltage clamp:
	double v_clamp_step = 0;
	double v_clamp_rest = -65;
	double steplength = 100;
	double R_clamp = 0.02;
	// -------------------------------------------------------------------------------
	// Cell geometry
	double Cm = Scale_volume * 60; // pF
	double V_tot = Scale_volume * 3960; // um^3
	double Vc_tenT = Scale_volume * 16404;
	double VSR_tenT = Scale_volume * 1094;
	double V_tot_tenT = Vc_tenT + VSR_tenT; // V_total data from Hwang et al., V_c and V_SR  proportionally scaled from Ten Tusscher 2004 values
	double Vc = V_tot * ( Vc_tenT / V_tot_tenT ); // =3712.4 um^3 (93.7% total volume)
	double V_SR = V_tot * ( VSR_tenT / V_tot_tenT ); // =247.6 um^3 (6.3% total volume)

	// -------------------------------------------------------------------------------
	// Constants
	double T = 310.0;   // kelvin (in model_parameters)'
	double R = 8.314472;   // joule_per_mole_kelvin (in model_parameters)
	double F = 96.4853415;   // coulomb_per_mmole (in model_parameters)

	double Ko = 5.4;   // millimolar (in model_parameters)
	double Cao = 1.8;  // millimolar (in model_parameters
	double Nao = 140.0;   // millimolar (in model_parameters)

	// -------------------------------------------------------------------------------
	// Reversal Potentials:
	double E_Ca = 0.5 * R * T / F * log( Cao / Y[2] ); // millivolt
	double E_Na = R * T / F * log( Nao / Y[3] );  // millivolt
	double E_K = R * T / F * log( Ko / Y[4] );  // millivolt

	// -------------------------------------------------------------------------------
	// Inward Rectifier K+ current (Ik1):
	// define parameters from x_K1
	double xK11 = x_K1[1]; // (2);
	double xK12 = x_K1[2]; // (3);
	double xK13 = x_K1[3]; // (4);
	double xK14 = x_K1[4]; // (5);
	double xK15 = x_K1[5]; // (6);

	double alpha_xK1 = xK11 * exp( ( Y[0] + xK13 ) / xK12 );
	double beta_xK1 = exp( ( Y[0] + xK15 ) / xK14 );
	double XK1_inf = alpha_xK1 / ( alpha_xK1 + beta_xK1 );

	// Current:
	double g_K1 = GK1_atrial * x_K1[0] * x_scale_conductance[0];
	double i_K1 = g_K1 * XK1_inf * ( Y[0] - E_K ) * sqrt( Ko / 5.4 );

	// -------------------------------------------------------------------------------
	// Rapid Delayed Rectifier Current (Ikr):
	// define parameters from x_KR
	double Xr1_1 = x_KR[1]; // (2);
	double Xr1_2 = x_KR[2]; // (3);
	double Xr1_5 = x_KR[3]; // (4);
	double Xr1_6 = x_KR[4]; // (5);
	double Xr2_1 = x_KR[5]; // (6);
	double Xr2_2 = x_KR[6]; // (7);
	double Xr2_5 = x_KR[7]; // (8);
	double Xr2_6 = x_KR[8]; // (9);
	// parameter-dependent values:
	double Xr1_3 = Xr1_5 * Xr1_1;
	double Xr2_3 = Xr2_5 * Xr2_1;
	double Xr1_4 = 1. / ( ( 1. / Xr1_2 ) + ( 1. / Xr1_6 ) );
	double Xr2_4 = 1. / ( ( 1. / Xr2_2 ) + ( 1. / Xr2_6 ) );

	// 10: Xr1 (dimensionless) (activation in i_Kr_Xr1)
	double alpha_Xr1 = Xr1_1 * exp( ( Y[0] ) / Xr1_2 );
	double beta_Xr1 = Xr1_3 * exp( ( Y[0] ) / Xr1_4 );
	double Xr1_inf = alpha_Xr1 / ( alpha_Xr1 + beta_Xr1 );
	double tau_Xr1 = ( ( 1. / ( alpha_Xr1 + beta_Xr1 ) ) + x_KR[9] ); // (10));
	dY[9] = ( Xr1_inf - Y[9] ) / tau_Xr1;

	// 11: Xr2 (dimensionless) (inactivation in i_Kr_Xr2)
	double alpha_Xr2 = Xr2_1 * exp( ( Y[0] ) / Xr2_2 );
	double beta_Xr2 = Xr2_3 * exp( ( Y[0] ) / Xr2_4 );
	double Xr2_inf = alpha_Xr2 / ( alpha_Xr2 + beta_Xr2 );
	double tau_Xr2 = ( ( 1. / ( alpha_Xr2 + beta_Xr2 ) ) + x_KR[10] ) ; // (11));
	// dY(11) = (Xr2_inf-Y(11))./tau_Xr2;
	dY[10] = ( Xr2_inf - Y[10] ) / tau_Xr2;

	// Current:
	double g_Kr = GKr_atrial * x_KR[0] * x_scale_conductance[1] ; // nS_per_pF (in i_Kr)
	double i_Kr = g_Kr * ( Y[0] - E_K ) * Y[9] * Y[10] * sqrt( Ko / 5.4 );

	// ----------------------------------------------------------------------------
	// Slow delayed rectifier current (IKs):
	// define parameters from x_IKS:
	double ks1 = x_IKS[1]; // (2);
	double ks2 = x_IKS[2]; // (3);
	double ks5 = x_IKS[3]; // (4);
	double ks6 = x_IKS[4]; // (5);
	double tauks_const = x_IKS[5]; // (6);

	// parameter-dependent values:
	double ks3 = ks5 * ks1;
	double ks4 = 1. / ( ( 1. / ks2 ) + ( 1. / ks6 ) );

	// 12: Xs (dimensionless) (activation in i_Ks)
	double alpha_Xs = ks1 * exp( ( Y[0] ) / ks2 );
	double beta_Xs = ks3 * exp( ( Y[0] ) / ks4 );
	double Xs_inf = alpha_Xs / ( alpha_Xs + beta_Xs );
	double tau_Xs = ( 1. / ( alpha_Xs + beta_Xs ) ) + tauks_const;
	dY[11] = ( Xs_inf - Y[11] ) / tau_Xs;

	// Current:
	double g_Ks = x_IKS[0] * x_scale_conductance[2] ;   // nS_per_pF (in i_Ks)
	double i_Ks = g_Ks * ( Y[0] - E_K ) * Y[11] * Y[11] ;





	double i_Kur, i_SK, i_KAch, i_K2P;
	i_Kur= i_SK= i_KAch= i_K2P= 0;

	if(atrial_like) {
	// NEw IKur model here; // 16:37:09, Mon, 25-February-2019, By Haibo
	double K_Q10 = 3.0;//3.5308257834747638;
	double inf = 1.0 / (1 + exp(-(Y[0] + 5.52) / (8.6))); // change to simple version, haibo.  removed IKur mutations here.
	double tau = (45.67 / (1.0 + exp((Y[0] + 9.23) / 10.53)) + 4.27)
	             * ( 0.55 / (1.0 + exp((Y[0] + 60.87) / -22.88)) + 0.03); // incorporate deactivation Thu 20 Jul 2017 23:24:01 BST Haibo
	tau = tau / K_Q10;
	// y[42] = inf + (y[42] - inf) * exp(-(dt) / tau);
	dY[23] = (inf - Y[23]) / tau;
	// y[43]
	inf = 1.0 / (1.0 + exp( (Y[0] + 7.5) / 6.67)); // Thu 20 Jul 2017 23:24:21 BST
	tau = 1300.0 / (exp((Y[0] - 133.0) / 97.0) * 5.04 + 0.2 * exp((Y[0] + 10.9233071) / (-12.0))) + 370.0; // at 37deg  Thu 20 Jul 2017 23:24:35 BST haibo Feng et al. 1998
	// y[43] = inf + (y[43] - inf) * exp(-(dt) / tau);
	dY[24] = (inf - Y[24]) / tau;
	tau = 1400.0 / (exp((Y[0] + 0.1) / 8.86) + 1.0) + 478.0 / (exp((Y[0] - 41.0) / -3.0) + 1.0) + 6500.0; // // at 37deg  Thu 20 Jul 2017 23:24:35 BST haibo Feng et al. 1998
	// y[44] = inf + (y[44] - inf) * exp(-(dt) / tau);
	dY[25] = (inf - Y[25]) / tau;
	double CNZ_gkur = /*(1.0 + 0.2 * para.kPKA_Ikur) * (1 - 0.5 * AF) **/ 0.006398 * 0.9;
	i_Kur = Gkur_atrial * CNZ_gkur * (4.5128 + 1.899769 / (1.0 + exp((Y[0] - 20.5232) / (-8.26597))))
	               * Y[23] * (0.65 * Y[24] + 0.35 * Y[25]) * (Y[0] - E_K);



	// IK2p model from CZ model // 11:19:03, Mon, 31-August-2020, By Haibo
	double i_K2P;
	double V = Y[0];  // vm
	double GK2p = 0.005;
	double o_k2p = Y[26];

	i_K2P = GK2P_atrial * GK2p * o_k2p * (V - E_K);

	inf = 0.2 + 0.8 / ( 1 + exp(-(V - (10/*+het.IK2p_Shift*/)) / 14.0));

	tau = 2.0 + 40 / (1.0 + exp((V + 25) * (V + 25) / 80.0));

	// o_k2p = inf + (o_k2p - inf) * exp(-(dt) / tau);  //
	dY[26] = (inf - Y[26]) / tau;

	// ISK model from CZ model // 11:19:03, Mon, 31-August-2020, By Haibo

	double KCa_on = 47e6 / 1000.0; // change from s to ms.
	double KCa_off = 13 / 1000.0;   // change from s to ms.
	double Cass = Y[2];
	double O_SK = Y[27];
	double d_O = (1 - O_SK) * KCa_on * Cass * Cass - KCa_off * O_SK;

	// double EC_50 = 0.43 * 1.0e-3;
	double GSK = GSK_atrial * 0.06 * 3.6;
	// O_SK = O_SK + dt * d_O;
	dY[27] = d_O;

	// double tau = 3.0;
	// double inf = 1.0 / (1.0 + pow(EC_50 / Cai, 3.6));


	i_SK =  GSK * O_SK * (1.0 / (1.0 + exp( (V - E_K + 120) / 45.0))) * (V - E_K);


	// %% I_k,ach: Muscarinic-Receptor-Activated K Current  // 11:57:36, Mon, 31-August-2020, By Haibo
	double Ach = 1e-4;

	// IKach from https://github.com/drgrandilab/Morotti-et-al-2016-Human-Atrial-Model-with-Updated-INa/blob/master/morotti_et_al_ham_ina_ran_model.m;
	// note the difference between 0.4 and 0.04 from Grandi et al. 2011
	// 0.4 seems to better reproduce the original plot Fig VI of Grandi et al., 2011
	i_KAch = GKAch_atrial * 1.0 / (1 + pow(0.03 / Ach, 2.1)) * (0.08 + 0.4 / (1 + exp((Y[0] + 91) / 12.0))) * (Y[0] - E_K);

	} else {
		dY[23] = dY[24] = dY[25] = dY[26] = dY[27] = 0;
	}


	// -------------------------------------------------------------------------------
	// Transient outward current (Ito):
	// define parameters from xTO
	double r1 = xTO[1]; // (2);
	double r2 = xTO[2]; // (3);
	double r5 = xTO[3]; // (4);
	double r6 = xTO[4]; // (5);
	double s1 = xTO[5]; // (6);
	double s2 = xTO[6]; // (7);
	double s5 = xTO[7]; // (8);
	double s6 = xTO[8]; // (9);
	double tau_r_const = xTO[9]; // (10);
	double tau_s_const = xTO[10]; // (11);

	// parameter-dependent values:
	double r3 = r5 * r1;
	double r4 = 1. / ( ( 1. / r2 ) + ( 1. / r6 ) );
	double s3 = s5 * s1;
	double s4 = 1. / ( ( 1. / s2 ) + ( 1. / s6 ) );

	// 17: s (dimensionless) (inactivation in i_to)
	double alpha_s = s1 * exp( ( Y[0] ) / s2 );
	double beta_s = s3 * exp( ( Y[0] ) / s4);
	double s_inf = alpha_s / ( alpha_s + beta_s );
	double tau_s = ( ( 1. / ( alpha_s + beta_s ) ) + tau_s_const );

	if (atrial_like)
	{
		s_inf = 1.0 / ( 1 + exp( -(Y[0] + 1.0) / 11.0 ) ) ;
		tau_s = 3.5 * exp(-((Y[0] / 30.0) * (Y[0] / 30.0))) + 1.5;
	}

	dY[16]  = ( s_inf - Y[16] ) / tau_s;

	// 18: r (dimensionless) (activation in i_to)
	double alpha_r = r1 * exp( ( Y[0] ) / r2 );
	double beta_r = r3 * exp( ( Y[0] ) / r4 );
	double r_inf = alpha_r / ( alpha_r + beta_r );
	double tau_r = ( 1. / ( alpha_r + beta_r ) ) + tau_r_const + atrial_like * 10;

	if (atrial_like) {

		r_inf = (1.0) / ( 1 + exp( (Y[0] + 40.5) / 11.5) ) ;
		tau_r = 25.635 * exp(-(((Y[0] + 52.45) / 15.8827) * ((Y[0] + 52.45) / 15.8827))) + 24.14; //% 14.14
	}
	dY[17] = ( r_inf - Y[17] ) / tau_r;

	// Current:
	double g_to = Gto_atrial * xTO[0] * x_scale_conductance[3] ; // nS_per_pF (in i_to)
	double i_to = g_to * ( Y[0] - E_K ) * Y[16] * Y[17] ;

	// -------------------------------------------------------------------------------
	// L-type Ca2+ current (ICaL):
	// define parameters from x_cal
	double d1 = x_cal[1]; // (2);
	double d2 = x_cal[2]; // (3);
	double d5 = x_cal[3]; // (4);
	double d6 = x_cal[4]; // (5);
	double f1 = x_cal[5]; // (6);
	double f2 = x_cal[6]; // (7);
	double f5 = x_cal[7]; // (8);
	double f6 = x_cal[8]; // (9);
	double taud_const = x_cal[9]; // (10);
	double tauf_const = x_cal[10]; // (11);

	// parameter-dependent values:
	double d3 = d5 * d1;
	double d4 = 1. / ( ( 1. / d2 ) + ( 1. / d6 ) );
	double f3 = f5 * f1;
	double f4 = 1. / ( ( 1. / f2 ) + ( 1. / f6 ) );


	// 7: d (dimensionless) (activation in i_CaL)
	double alpha_d = d1 * exp( ( ( Y[0] - CaL_v_shift) ) / d2 );
	double beta_d = d3 * exp( ( ( Y[0] - CaL_v_shift ) ) / d4 );
	double d_inf = alpha_d / ( alpha_d + beta_d );
	double tau_d = ( ( 1. / ( alpha_d + beta_d ) ) + taud_const );
	dY[6] = ( d_inf - Y[6] ) / tau_d;

	// 8: f (dimensionless) (inactivation  i_CaL)
	double alpha_f = f1 * exp( ( ( Y[0] ) ) / f2 );
	double beta_f = f3 * exp( ( ( Y[0] ) ) / f4 );
	double f_inf = alpha_f / ( alpha_f + beta_f );
	double tau_f = ( ( 1. / ( alpha_f + beta_f ) ) + tauf_const );
	dY[7] = ( f_inf - Y[7] ) / tau_f;

	// 9: fCa (dimensionless) (calcium-dependent inactivation in i_CaL)
	// from Ten tusscher 2004
	double scale_Ical_Fca_Cadep = 1.2;
	double alpha_fCa = 1.0 / ( 1.0 + pow( ( ( scale_Ical_Fca_Cadep * Y[2] ) / 0.000325 ), 8. ) ) ; //^8.0);
	double beta_fCa = 0.1 / ( 1.0 + exp( ( scale_Ical_Fca_Cadep * Y[2] - 0.0005 ) / 0.0001 ) );
	double gamma_fCa = 0.2 / ( 1.0 + exp( ( scale_Ical_Fca_Cadep * Y[2] - 0.00075 ) / 0.0008 ) );

	double fCa_inf = ( ( alpha_fCa + beta_fCa + gamma_fCa + 0.23 ) / ( 1.46 ) );
	double tau_fCa = 2. ; // ms

	double k_fca;
	if ( fCa_inf > Y[8] && Y[0] > -60. ) {
		k_fca = 0;
	} else {
		k_fca = 1;
	} //end
	dY[8] = k_fca * ( fCa_inf - Y[8] ) / tau_fCa;



	// Current:
	double p_CaL =  x_cal[0] * x_scale_conductance[4]; // nS_per_pF (in i_CaL)
	double p_CaL_shannonCa = 5.4e-4;
	double p_CaL_shannonNa = 1.5e-8;
	double p_CaL_shannonK = 2.7e-7;
	double p_CaL_shannonTot = p_CaL_shannonCa + p_CaL_shannonNa + p_CaL_shannonK;
	double p_CaL_shannonCap = p_CaL_shannonCa / p_CaL_shannonTot;
	double p_CaL_shannonNap = p_CaL_shannonNa / p_CaL_shannonTot;
	double p_CaL_shannonKp = p_CaL_shannonK / p_CaL_shannonTot;

	double p_CaL_Ca = G_CaL_atrial * p_CaL_shannonCap * p_CaL;
	double p_CaL_Na = G_CaL_atrial * p_CaL_shannonNap * p_CaL;
	double p_CaL_K = G_CaL_atrial * p_CaL_shannonKp * p_CaL;

	double ibarca = ( p_CaL_Ca * 4.0 * Y[0] * F * F
	                  / ( R * T )
	                  * ( 0.341 * Y[2] * exp( 2.0 * Y[0] * F
	                          / ( R * T ) )
	                      - 0.341 * Cao )
	                  / ( exp( 2.0 * Y[0] * F
	                           / ( R * T ) ) - 1.0 ) );
	double i_CaL_Ca =  ibarca * Y[6] * Y[7] * Y[8];


	double ibarna = ( p_CaL_Na * Y[0] * F * F
	                  / ( R * T )
	                  * ( 0.75 * Y[3]
	                      * exp( Y[0] * F / ( R * T ) )
	                      - 0.75 * Nao )
	                  / ( exp( Y[0] * F / ( R * T ) )
	                      - 1.0 ) );
	double i_CaL_Na =  ibarna * Y[6] * Y[7] * Y[8];

	double ibark = ( p_CaL_K * Y[0] * F * F
	                 / ( R * T )
	                 * ( 0.75 * Y[4]
	                     * exp( Y[0] * F / ( R * T ) )
	                     - 0.75 * Ko )
	                 / ( exp( Y[0] * F / ( R * T ) )
	                     - 1.0 ) );
	double i_CaL_K = ibark * Y[6] * Y[7] * Y[8] ;

	/*if (atrial_like) {

		// grandi et al. 2011 model
		double ISO = 0;
		double dss = 1 / (1 + exp(-(Y[0] + 3 * ISO + 9) / 6)); //in Maleckar v1/2=-9 S=6 (mV); Courtemanche v1/2=-9 S=5.8 (mV)
		double taud = 1 * dss * (1 - exp(-(Y[0] + 3 * ISO + 9) / 6)) / (0.035 * (Y[0] + 3 * ISO + 9));
		double fss = 1 / (1 + exp((Y[0] + 3 * ISO + 30) / 7)) + 0.2 / (1 + exp((50 - Y[0] - 3 * ISO) / 20)); // in Maleckar v1/2=-27.4 S=7.1 (mV); Courtemanche v1/2=-28 S=6.9 (mV)
		double tauf = 1 / (0.0197 * exp( -pow(0.0337 * (Y[0] + 3 * ISO + 25), 2) ) + 0.02);
		dY[6] = (dss - Y[6]) / taud;
		dY[7] = (fss - Y[7]) / tauf;
		dY[8] = 1.7 * Y[2] * (1 - Y[8]) - 1 * 11.9e-3 * Y[8]; // fCa_junc   koff!!!!!!!!
		// ydot(7) = 1.7*y(37)*(1-y(7))-1*11.9e-3*y(7); // fCa_sl
		i_CaL_Na =  ibarna * Y[6] * Y[7] * (1-Y[8]);
		i_CaL_Ca =  ibarca * Y[6] * Y[7] * (1-Y[8]);

		i_CaL_K = ibark * Y[6] * Y[7] * (1-Y[8]) ;

	}*/

	if (atrial_like) {

		double tau_d;//, tau_f, tau_f_Ca;

		if (fabs(Y[0] + 10.0) < 1.0e-10)
			tau_d = 4.579 / (1.0 + exp((Y[0] + 10.0) / -6.24));
		else
			tau_d = (1.0 - exp((Y[0] + 10.0) / -6.24)) / (0.035 * (Y[0] + 10.0) * (1.0 + exp((Y[0] + 10.0) / -6.24)));
		double d_infinity = pow(1.0 + exp((Y[0] + 10.0) / -8.0), -1.0);

		dY[6] = (d_infinity - Y[6]) / tau_d;

		double f_infinity = exp(-(Y[0] + 28.0) / 6.9) / (1.0 + exp(-(Y[0] + 28.0) / 6.9));
		double tau_f = 9.0 * pow(0.0197 * exp(-pow(0.0337, 2.0) * pow(Y[0] + 10.0, 2.0)) + 0.02, -1.0);
		dY[7] = (f_infinity - Y[7]) / tau_f;
		// dstate[6] = (f_infinity-state[6])/tau_f;
		// dstate[4] = (d_infinity-state[4])/tau_d;
		// state[4] = Euler_inf(dt, state[4], d_infinity, tau_d);
		// state[6] = Euler_inf(dt, state[6], f_infinity, tau_f);
		double f_Ca_infinity = pow(1.0 + Y[2] / 0.00035, -1.0);
		double tau_f_Ca = 2.0;
		// dstate[5] = (f_Ca_infinity-state[5])/tau_f_Ca;
		// state[5] = Euler_inf(dt, state[5], f_Ca_infinity, tau_f_Ca);
		dY[8] = (f_Ca_infinity - Y[8]) / tau_f_Ca;

		double g_Ca_L = 0.12375;

		i_CaL_Ca = G_CaL_atrial* g_Ca_L * Y[6] * Y[7] * Y[8] * (Y[0] - 65.0);

		i_CaL_K = ibark * Y[6] * Y[7] * Y[8] ;
		i_CaL_Na = ibarna * Y[6] * Y[7] * Y[8] ;
	}

	double i_CaL = i_CaL_Ca + i_CaL_Na + i_CaL_K;
	// -------------------------------------------------------------------------------
	// T-type Calcium Current (ICaT):
	// SAN T-TYPE CA2+ model (Demir et al., Maltsev-Lakatta ), G_CaT determined by fit to Kurokawa IV:

	// 19: dCaT (activation in i_CaT)
	double dcat_inf = 1. / ( 1. + exp( -( ( Y[0] ) + 26.3 ) / 6. ) );
	double tau_dcat = 1. / ( 1.068 * exp( ( ( Y[0] ) + 26.3 ) / 30. )
	                         + 1.068 * exp( -( ( Y[0] ) + 26.3 ) / 30. ) );
	dY[18] = ( dcat_inf - Y[18] ) / tau_dcat;

	// 20: fCaT (inactivation in i_CaT)
	double fcat_inf = 1. / ( 1. + exp( ( ( Y[0] ) + 61.7 ) / 5.6 ) ) ;
	double tau_fcat = 1. / ( 0.0153 * exp( -( ( Y[0] ) + 61.7 ) / 83.3 )
	                         + 0.015 * exp( ( ( Y[0] ) + 61.7 ) / 15.38 ) );
	dY[19] = ( fcat_inf - Y[19] ) / tau_fcat;

	double g_CaT = GCat_atrial * x_cat[0] * x_scale_conductance[5] ; // nS_per_pF (in i_CaT)
	double i_CaT = g_CaT * ( Y[0] - E_Ca ) * Y[18] * Y[19];
	// double i_CaT = g_CaT * ( Y[0] - (45) ) * Y[18] * Y[19];

	// -------------------------------------------------------------------------------
	// Sodium Current (INa):
	// define parameters from x_Na
	double m1 = x_NA[1]; // (2);
	double m2 = x_NA[2]; // (3);
	double m5 = x_NA[3]; // (4);
	double m6 = x_NA[4]; // (5);
	double h1 = x_NA[5]; // (6);
	double h2 = x_NA[6]; // (7);
	double h5 = x_NA[7]; // (8);
	double h6 = x_NA[8]; // (9);
	double j1 = x_NA[9]; // (10);
	double j2 = x_NA[10]; // (11);
	double tau_m_const = x_NA[11]; // (12);
	double tau_h_const = x_NA[12]; // (13);
	double tau_j_const = x_NA[13]; // (14);

	// parameter-dependent values:
	double m3 = m5 * m1;
	double m4 = 1. / ( ( 1. / m2 ) + ( 1. / m6 ) );
	double h3 = h5 * h1;
	double h4 = 1. / ( ( 1. / h2 ) + ( 1. / h6 ) );
	double j5 = h5;
	double j6 = h6;
	double j3 = j5 * j1;
	double j4 = 1. / ( ( 1. / j2 ) + ( 1. / j6 ) );

	// 13: h (dimensionless) (inactivation in i_Na)
	double alpha_h = h1 * exp( ( Y[0] ) / h2 );
	double beta_h = h3 * exp( ( Y[0] ) / h4 );
	double h_inf = ( alpha_h / ( alpha_h + beta_h ) );
	double tau_h = ( ( 1. / ( alpha_h + beta_h ) ) + tau_h_const );
	dY[12] = ( h_inf - Y[12] ) / tau_h;

	// 14: j (dimensionless) (slow inactivation in i_Na)
	double alpha_j = j1 * exp( ( Y[0] ) / j2 );
	double beta_j = j3 * exp( ( Y[0] ) / j4 );
	double j_inf = ( alpha_j / ( alpha_j + beta_j ) );
	double tau_j = ( ( 1. / ( alpha_j + beta_j ) ) + tau_j_const );
	dY[13] = ( j_inf - Y[13] ) / tau_j;

	// 15: m (dimensionless) (activation in i_Na)
	double alpha_m = m1 * exp( ( Y[0] ) / m2 );
	double beta_m = m3 * exp( ( Y[0] ) / m4 );
	double m_inf = alpha_m / ( alpha_m + beta_m );
	double tau_m = ( ( 1. / ( alpha_m + beta_m ) ) + tau_m_const );
	dY[14] = ( m_inf - Y[14] ) / tau_m;




	#ifdef USE_V1
	// double Ina_h_shift = /*0.5 * */(5 - 2 * 5 * 1 / (1 + exp(-(para.NaV_CKp - 0.12) / 0.1))); // h gate dependent on INa camKII
	double Vm = Y[0] + Ina_shift;
	// double GNa = 9/*9*/ ; // [mS/uF] Courtemanche (adjusted for dV/dt max) MOD1

	double am = (Vm == -47.13) ? 3.2 : (0.32 * (Vm + 47.13) / (1 - exp( -0.1 * (Vm + 47.13))));

	double bm = 0.08 * exp(-Vm / 11.0);

	// double Vm = y[38] - Ina_h_shift;

	double ah = (Vm >= -40) ? 0 : 0.135 * exp( -(Vm + 80) / 6.8 );

	double bh = (Vm >= -40) ? (1.0 / (0.13 * (1 + exp( -(Vm + 10.66) / 11.1 )))) : (3.56 * exp( 0.079 * Vm) + 3.1e5 * exp(0.35 * Vm));

	double aj = (Vm >= -40) ? 0 : ( (-127140 * exp(0.2444 * Vm) - 3.474e-5 * exp(-0.04391 * Vm)) * (Vm + 37.78)) / (1 + exp( 0.311 * (Vm + 79.23) ) );
	double bj = (Vm >= -40) ? ((0.3 * exp(-2.535e-7 * Vm)) / (1 + exp( -0.1 * (Vm + 32) ))) : (0.1212 * exp( -0.01052 * Vm)) / (1 + exp( -0.1378 * (Vm + 40.14) ));

	dY[14] = am * (1 - Y[14]) - bm * Y[14];
	dY[13]  = ah * (1 - Y[13]) - bh * Y[13];
	dY[12]  = aj * (1 - Y[12]) - bj * Y[12];
	Vm = Y[0];
	#endif

	// Current:
	double g_Na = GNa_atrial * x_NA[0] * x_scale_conductance[6] ; // nS_per_pF (in i_Na)
	double i_Na = g_Na * Y[14] * Y[14] * Y[14] * Y[12] * Y[13] * ( Y[0] - E_Na );

	// -------------------------------------------------------------------------------
	// -------------------------------------------------------------------------------
	// Funny/HCN current (If):
	// define parameters from x_F
	double xF1 = x_F[1]; // (2);
	double xF2 = x_F[2]; // (3);
	double xF5 = x_F[3]; // (4);
	double xF6 = x_F[4]; // (5);
	double xF_const = x_F[5]; // (6);

	// parameter-dependent values:
	double xF3 = xF5 * xF1;
	double xF4 = 1. / ( ( 1. / xF2 ) + ( 1. / xF6 ) );

	// 16: Xf (dimensionless) (inactivation in i_f)
	double alpha_Xf = xF1 * exp( ( Y[0] ) / xF2 );
	double beta_Xf = xF3 * exp( ( Y[0] ) / xF4 );
	double Xf_inf = alpha_Xf / ( alpha_Xf + beta_Xf );
	double tau_Xf = ( ( 1. / ( alpha_Xf + beta_Xf ) ) + xF_const );
	dY[15] = ( Xf_inf - Y[15] ) / tau_Xf;

	// Current:

	// fprintf(stderr, "%f\n", Gf_atrial);
	double g_f = Gf_atrial * x_F[0] * x_scale_conductance[7] ; // nS_per_pF (in i_f)
	double NatoK_ratio = 0.491; // Verkerk et al. 2013
	double Na_frac = NatoK_ratio / ( NatoK_ratio + 1. );
	double i_fNa = Na_frac * g_f * Y[15] * ( Y[0] - E_Na );
	double i_fK = ( 1. - Na_frac ) * g_f * Y[15] * ( Y[0] - E_K );
	double i_f = i_fNa + i_fK;


	#ifdef USE_V1

	double ih_pK = 0.6167; // ih_PNa = 0.3833;
	double slope_gh = 13.3;
	double vhalf_gh = 114.4;
	//y_inf = 1.0/(1.0 + exp((v+vhalf_gh)/16.3));
	double y_inf = 1.0 / (1.0 + exp((Vm + vhalf_gh) / slope_gh));
	double tau_y_1_2 = 1 * (1.5049 / (exp(-(Vm + 590.3) * 0.01094) + exp((Vm - 85.1) / 17.2)));
	dY[15] = ( y_inf - Y[15] ) / tau_y_1_2;

	// double y_1_2_dot = (y_inf - y_1_2) / tau_y_1_2;
	i_fK  = ih_pK *  g_f * Y[15]  * (Vm - E_K);
	i_fNa = (1 - ih_pK) * g_f * Y[15]  * (Vm - E_Na);
	i_f =  i_fNa + i_fK;
	#endif
	// -------------------------------------------------------------------------------
	// Na+/Ca2+ Exchanger current (INaCa):
	// Ten Tusscher formulation
	double KmCa = 1.38;         // Cai half-saturation constant millimolar (in i_NaCa)
	double KmNai = 87.5;        // Nai half-saturation constnat millimolar (in i_NaCa)
	double Ksat = 0.1;          // saturation factor dimensionless (in i_NaCa)
	double gamma = 0.35 * 2;    // voltage dependence parameter dimensionless (in i_NaCa)
	double alpha = 2.5 * 1.1;   // factor to enhance outward nature of inaca dimensionless (in i_NaCa)
	double kNaCa = GNaca_atrial * 1000. * 1.1 * x_scale_conductance[8] ;  // maximal inaca pA_per_pF (in i_NaCa)

	double i_NaCa = ( kNaCa * ( ( exp( gamma * Y[0] * F / ( R * T ) )
	                              * ( Y[3] * Y[3] * Y[3] ) * Cao )
	                            - ( exp( ( gamma - 1.0 ) * Y[0] * F / ( R * T ) )
	                                * ( Nao * Nao * Nao ) * Y[2] * alpha ) )
	                  / ( ( ( KmNai * KmNai * KmNai ) + ( Nao * Nao * Nao ) )
	                      * ( KmCa + Cao )
	                      * ( 1.0 + Ksat * exp( ( gamma - 1.0 ) * Y[0] * F / ( R * T ) ) ) ) );

	// -------------------------------------------------------------------------------
	// Na+/K+ pump current (INaK):
	// Ten Tusscher formulation
	double Km_K = 1.0;   // Ko half-saturation constant millimolar (in i_NaK)
	double Km_Na = 40.0;   //  Nai half-saturation constant millimolar (in i_NaK)
	double PNaK = GNaK_atrial * 1.362 * 1.818 * x_scale_conductance[12] ;   // maxiaml nak pA_per_pF (in i_NaK)
	double i_NaK = ( PNaK * ( ( Ko * Y[3] )
	                          / ( ( Ko + Km_K ) * ( Y[3] + Km_Na )
	                              * ( 1.0
	                                  + 0.1245 * exp( -0.1 * Y[0] * F / ( R * T ) )
	                                  + 0.0353 * exp( -Y[0] * F / ( R * T ) ) ) ) ) );

	// -------------------------------------------------------------------------------
	// SR Uptake/SERCA (J_up):
	// Ten Tusscher formulation
	double Kup = 0.00025 * 0.702;   // millimolar (in calcium_dynamics)
	double VmaxUp = Gserca * 0.000425 * 0.26 * x_scale_conductance[9] ; // millimolar_per_milisecond (in calcium_dynamics)
	double i_up = VmaxUp / ( 1.0 + Kup * Kup / ( Y[2] * Y[2] ) );

	#ifdef USE_V1
	// from grandi et al. 2011
	double hillSRCaP = 1.787;       // [mM]
	double Vmax_SRCaP = Gserca * 0.2 * 0.25 * 0.5 * 5.3114e-3; //% [mM/msec] (286 umol/L cytosol/sec)
	double Kmf = Kmf_Scale * Kup;//(2.5 - 1.25 * 0) * 0.246e-3;    //% [mM] default
	double Kmr = 1.7;               //% [mM]L cytosol
	i_up = Vmax_SRCaP * (pow(Y[2] / Kmf, hillSRCaP) - pow(Y[1] / Kmr, hillSRCaP) )
	       / (1 + pow(Y[2] / Kmf, hillSRCaP) + pow(Y[1] / Kmr, hillSRCaP) );
	#endif
	// -------------------------------------------------------------------------------
	// SR Leak (J_leak):
	// Ten Tusscher formulation
	double V_leak = x_scale_conductance[11] * 0.00008 * 0.02;   // per_millisecond (in calcium_dynamics)
	double i_leak = ( Y[1] - Y[2] ) * V_leak;

	// -------------------------------------------------------------------------------
	// SR Release/RYR (J_rel):
	// re-fit parameters. scaled to account for differences in calcium concentration in
	// cleft (cleft is used in shannon-bers model geometry, not in this model geometry)
	double ks = 12.5 * x_scale_conductance[10];     // [1/ms]
	double koCa = 56320. * 11.43025;                // [mM^-2 1/ms]
	double kiCa = 54. * 0.3425;                     // [1/mM/ms]
	double kom = 1.5 * 0.1429;                      // [1/ms]
	double kim = 0.001 * 0.5571;                    // [1/ms]
	double ec50SR = ec50SR_Scale * 0.45;
	double MaxSR = 15.;
	double MinSR = 1.;

	double kCaSR = MaxSR - ( MaxSR - MinSR ) / ( 1. + pow( ( ec50SR / Y[1] ), 2.5 ) );
	double koSRCa = koCa / kCaSR;
	double kiSRCa = kiCa * kCaSR;
	double RI = 1. - Y[20] - Y[21] - Y[22] ;

	dY[20] = ( ( kim * RI - kiSRCa * Y[2] * Y[20] )
	           - ( koSRCa * Y[2] * Y[2] * Y[20]  - kom * Y[21] ) );   // R
	dY[21] = ( ( koSRCa * Y[2] * Y[2] * Y[20] - kom * Y[21] )
	           - ( kiSRCa * Y[2] * Y[21] - kim * Y[22] ) ); // O
	dY[22] = ( ( kiSRCa * Y[2] * Y[21] - kim * Y[22] )
	           - ( kom * Y[22] - koSRCa * Y[2] * Y[2] * RI ) );   // I

	double i_rel = RyR_scale * ks * Y[21] * ( Y[1] - Y[2] ) * ( V_SR / Vc );

	// Background Sodium (I_bNa):
	// Ten Tusscher formulation
	double g_b_Na = GbNa_atrial * 0.00029 * 1.5 * x_scale_conductance[13];  // nS_per_pF (in i_b_Na)
	double i_b_Na = g_b_Na * ( Y[0]  - E_Na );

	// -------------------------------------------------------------------------------
	// Background Calcium (I_bCa):
	// Ten Tusscher formulation
	double g_b_Ca = GbCa_atrial * 0.000592 * 0.62 * x_scale_conductance[14];   // nS_per_pF (in i_b_Ca)
	double i_b_Ca = g_b_Ca * ( Y[0] - E_Ca );

	// -------------------------------------------------------------------------------
	// Calcium SL Pump (I_pCa):
	// Ten Tusscher formulation
	double g_PCa = GPCa_atrial * 0.025 * 10.5 * x_scale_conductance[15];  // pA_per_pF (in i_PCa)
	double KPCa = 0.0005;   // millimolar (in i_PCa)
	double i_PCa = g_PCa * Y[2] / ( Y[2] + KPCa );

	// -------------------------------------------------------------------------------
	// 2: CaSR (millimolar)
	// rapid equilibrium approximation equations -- not as formulated in ten Tusscher 2004 text
	double Buf_SR = Buffer_Scale * 10.0 * 1.2; // millimolar (in calcium_dynamics)
	double Kbuf_SR = 0.3; // millimolar (in calcium_dynamics)
	double Ca_SR_bufSR = ( 1.
	                       / ( 1.0
	                           + Buf_SR * Kbuf_SR
	                           / ( ( Y[1] + Kbuf_SR )
	                               * ( Y[1] + Kbuf_SR ) ) ) );

	dY[1] = Ca_SR_bufSR * Vc / V_SR * ( i_up - ( i_rel + i_leak ) );

	// -------------------------------------------------------------------------------
	// 3: Cai (millimolar)
	// rapid equilibrium approximation equations -- not as formulated in ten Tusscher 2004 text
	double Buf_C = 0.06; // millimolar (in calcium_dynamics)
	double Kbuf_C = 0.0006; // millimolar (in calcium_dynamics)
	double Cai_bufc = 1. / ( 1.0 + Buf_C * Kbuf_C / ( ( Y[2] + Kbuf_C ) * ( Y[2] + Kbuf_C ) ) );


	dY[5] = 0;// zoero here // https://github.com/ClancyLabUCD/IPSC-model/issues/1
	dY[2] = ( Cai_bufc ) * ( i_leak - i_up + i_rel - dY[5]
	                         - ( i_CaL_Ca + i_CaT + i_b_Ca + i_PCa - 2. * i_NaCa )
	                         * Cm / ( 2.0 * Vc * F ) );

	// -------------------------------------------------------------------------------
	// 4: Nai (millimolar) (in sodium_dynamics)
	dY[3] = -Cm * ( i_Na + i_b_Na + i_fNa + 3.0 * i_NaK + 3.0 * i_NaCa + i_CaL_Na ) / ( F * Vc );

	if ( stim_flag == 1 ) {
		dY[3] = 0;
	} // end

	//-------------------------------------------------------------------------------
	// 5: Ki (millimolar) (in potatssium_dynamics)
	dY[4] = -Cm * ( i_Kur + i_K1 + i_to + i_Kr + i_Ks + i_fK - 2. * i_NaK + i_CaL_K + i_SK + i_K2P + i_KAch) / ( F * Vc );
	// dY[4] = 0; // 12:34:04, Tue, 01-September-2020, By Haibo
	if ( stim_flag == 1 ) {
		dY[4] = 0;
	} // end
	dY[4] = 0;  // for GA, set to constant 
	// -------------------------------------------------------------------------------
	// 1: Vm (Membrane voltage)

	// I_stim:
	double time = t;
	double i_stim;

	if ( ( time >= i_stim_Start )
	        && ( time <= i_stim_End )
	        && ( fmod( time - i_stim_Start - 100, cyclelength ) < i_stim_PulseDuration ) ) {
		i_stim = stim_flag * i_stim_Amplitude;
	} else {
		i_stim = 0.0;
	} // end

	// Voltage Clamp:
	double v_clamp;
	if ( voltageclamp == 0 ) {
		v_clamp = Y[0] ; //set i_voltageclamp to 0
	} else if ( voltageclamp == 1 ) {  // train of square pulse:
		if ( fmod( time, cyclelength ) < cyclelength - steplength ) {
			v_clamp = v_clamp_rest;
		} else {
			v_clamp = v_clamp_step;
		} // end
	} //end

	double i_voltageclamp = ( v_clamp - Y[0]  ) / R_clamp;

	dY[0] = -( i_Kur + i_K1 + i_to + i_Kr + i_Ks + i_SK + i_K2P + i_KAch
	           + i_CaL + i_CaT + i_NaK + i_Na + i_NaCa + i_PCa
	           + i_f + i_b_Na + i_b_Ca - i_stim - i_voltageclamp );

	// currents = [i_K1, i_to, i_Kr, i_Ks, i_CaL, i_NaK, i_Na, i_NaCa, i_PCa, i_f, i_b_Na, i_b_Ca, i_rel, i_up, i_leak, i_stim, i_CaT];
	currents[0] = i_K1;
	currents[1] = i_to;
	currents[2] = i_Kr;
	currents[3] = i_Ks;
	currents[4] = i_CaL;
	currents[5] = i_NaK;
	currents[6] = i_Na;
	currents[7] = i_NaCa;
	currents[8] = i_PCa;
	currents[9] = i_f;
	currents[10] = i_b_Na;
	currents[11] = i_b_Ca;
	currents[12] = i_rel;
	currents[13] = i_up;
	currents[14] = i_leak;
	currents[15] = i_stim;  // 17
	currents[16] = i_CaT;  // 18
	currents[17] = i_Kur;
	currents[18] = i_K2P;  //20
	currents[19] = i_SK;
	currents[20] = i_KAch;  // 22

	// end
	//%===============================================================================
	//% End of file
	//%===============================================================================

}


