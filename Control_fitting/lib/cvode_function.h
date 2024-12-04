// #include "Lsoda_function_wrap.h"
// #include <cvode/cvode.h>                          /* prototypes for CVODE fcts., consts.          */




#ifndef CVODE_FUNCTION__H
#define CVODE_FUNCTION__H

#include "cvode_solver.hpp"
#include "cell_parameters.h"
int f(realtype t, N_Vector y, N_Vector ydot, void *user_data);


#endif