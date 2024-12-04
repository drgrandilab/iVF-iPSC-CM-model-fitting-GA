

#include "stimulus.h"




class cell_parameters
{
public:
    cell_parameters(){

       stim_file =  NULL;


    };
    ~cell_parameters(){ 

        if(stim_file != NULL)
            delete stim_file;

    };
    double *model_parameter_inputs;
    double *currents;
    int ODE_NUM;

    StimFromInputFile *stim_file;

};

