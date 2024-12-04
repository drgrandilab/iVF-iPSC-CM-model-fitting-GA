// Haibo Ni,
// new atrial ipsc cells



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




#include <fstream>
#include <iostream>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <cstdlib>
// #include <omp.h>
#include <string>

// #include "integrate_rk2.h"
// #include "ipsc_function.h"
#include "cvode_function.h"

// #define USE_V1

// #define OUT_CURRENTS

#define NEQ (23+3+2)
int main(int argc, char const *argv[])
{
    using namespace std;
    cout.precision(16);

    //% Main File to generate Baseline Model Figures 10-11
    //load ICs_baseline
    int sy0 = 23 + 5;
    double y0[sy0], y0n[sy0];

    char filename[] = "Input_data/ICs_baseline.txt";
    FILE *fp;
    fp = fopen( filename, "r" );
    for ( int idy = 0; idy < sy0; idy++ ) {
        fscanf( fp, "%lf", &y0[idy] );
        y0n[idy] = y0[idy];
        // cout << "Y_init[" << idy << "] = " << y0[idy] << endl;
    }
    fclose( fp );


    //load baseline_parameter_inputs
    int sp0 = 87;
    int sp1 = 87 + 18;
    double p0[sp1];
    fp = fopen( "Input_data/baseline_parameter_inputs.txt", "r" );
    for ( int idp = 0; idp < sp0; idp++ ) {
        fscanf( fp, "%lf", &p0[idp] );
        // cout << "P_init[" << idp << "] = " << p0[idp] << endl;
    }
    fclose( fp );

    for (int i = 0; i < 18; ++i)
    {
        p0[i + sp0] = 1; // (1:18);
        // std::cout <<  "P_init[" << i << "] = " << p0[i] << endl;
    }

    std::string thread_ID = "0";
    if (argc == 18) {
        // for (int i = 0; i < 18; ++i)
        // {
        //     p0[i+sp0] = atof(argv[i+1]); // (1:18);
        //     // std::cout <<  "P_init[" << i << "] = " << p0[i] << endl;
        // }
    } else     if (argc == 20) {
        for (int i = 0; i < 18; ++i)
        {
            p0[i + sp0] = p0[i] * atof(argv[i + 1]); // (1:18);
            // std::cout <<  "P_init[" << i << "] = " << p0[i] << endl;
        }
        thread_ID = argv[19];
    }




#ifdef DEBUG
    for ( int idp = 0; idp < sp1; idp++ ) {
        // fscanf( fp, "%lf", &p0[idp] );
        cout << "P_init[" << idp << "] = " << p0[idp] << endl;
    }
#endif

    double t, dt;
    dt = 1;
    t = 0;
    int sc = 30;
    double currents[sc] ;
    for ( int idc = 0; idc < sc; idc++ ) {
        currents[idc] = 0;
    }

    FILE *fp_currents, *fp_y;

    std::string filename_out = "ys.txt." + thread_ID;

    fp_y = fopen( filename_out.c_str(), "w" );
    fp_currents = fopen( "currents.txt", "w" );


    cell_parameters Cell;
    Cell.model_parameter_inputs = p0;
    Cell.currents = currents;
    Cell.ODE_NUM = NEQ;

    Cell.stim_file = new StimFromInputFile("stim.0.5Hz.iVF", false);

    cvode_solver cvode(NEQ, 0.2);
    // std::cout << "1";
    cvode.set_IC(y0);

    cvode.initialise_mem(f);

    cvode.set_user_data(&Cell);


    double total_time = 110e3;//90e3 +50e3;
    // std::cout << "1";


    int outfreq = 1;

    int counter = 0;

    for ( t = 0; t < total_time; t += dt ) {
        //%% Run iPSC_function
        //options = odeset('MaxStep',1,'InitialStep',2e-2);
        //run_time=3e3;
        //[Time, values] = ode15s(@ipsc_function,[0, run_time],Y_init, options, baseline_parameter_inputs);

        if ( counter % outfreq == 0  and (t >=98e3) ) {
            // fprintf( fp_y, "%f\t", t );
            printf( "%f\t", t );
            // for( int idy = 0; idy < sy0; idy++ ) {
            for ( int idy = 0; idy < 3; idy++ ) {
                // fprintf( fp_y, "%f\t", Ith(cvode.y, idy + 1) );
                printf( "%.8f\t", Ith(cvode.y, idy + 1) );
            }
            // fprintf( fp_y, "\n" );
            printf(  "\n" );


            #ifdef OUT_CURRENTS
            fprintf( fp_currents, "%16.14e\t", t );
            for ( int idc = 0; idc < sc; idc++ ) {
                fprintf( fp_currents, "%16.14e\t", Cell.currents[idc] );
            }
            fprintf( fp_currents, "\n" );
            #endif
        }

        double outputTime = t + dt;
        cvode.solve_single_step(outputTime); // 1 time step solution.
        // integrate_rk2( ipsc_function, &t, sy0, y0n, dt, p0 , currents );

        counter ++;
    }

    fclose( fp_y );
    fclose( fp_currents );

}
