// #include "Lsoda_function_wrap.h"
// #include <cvode/cvode.h>                          /* prototypes for CVODE fcts., consts.          */
#include "cvode_solver.hpp"

#include "atrial_ipsc_function.h"
#include "cell_parameters.h"


int f(realtype t, N_Vector y, N_Vector ydot, void *user_data) {

	cell_parameters *Data = (cell_parameters*) user_data;
	int NEQ = Data->ODE_NUM;

	// std::cout << NEQ<< std::endl;

	realtype Y[NEQ];
	realtype dY[NEQ];

	// reset to zoero
	for (int i = 0; i < NEQ; i++)
		dY[i] = 0;
	// double
	// #pragma novector
	for (int i = 0; i < NEQ; i++)
		Y[i] = Ith(y, i + 1);
	// Data->Master_ODE_update(t);

	atrial_ipsc_function(  t,  Y,  Data->model_parameter_inputs,  dY,  Data->currents );

	bool allow_stimulation = true;
	if (allow_stimulation) {
		if (Data->stim_file != NULL) {
			double i_stim = Data->stim_file-> ApplyStim(10.0, 3.0, 0, t); // duration 3.0ms; strength: 10.0
			dY[0] += i_stim; // add stimlus to dY0, which is Vm
			Data->currents[15] = i_stim;
		} else {
			std::cerr << " Data-> stim_file == NULL!!!!"<<std::endl;
		}
	}


	for (int i = 0; i < NEQ; i++)
		Ith(ydot, i + 1) = dY[i];


	// std::cout << Ith(ydot, 1)<< std::endl;

	// Y = Data->y;
	return 0;

}



