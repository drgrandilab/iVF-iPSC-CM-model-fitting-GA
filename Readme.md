


### Packages and Software

* Matlab R2022b with optimization toolbox
* Intel compiler for C++ (icpc version 2021.6.0 [gcc version 9.4.0 compatibility])

### Folders
* Control_fitting: source code for fitting to control traces
* iVF_fitting: source code for fitting to iVF traces


### Test runs
* Compile C++ code for iPSC-CMs
```sh

make clean;
make cvode_main_ipsc_baseline_05Hz
make cvode_main_ipsc_baseline_10Hz

```

* Matlab optimization using GA  
Run Matlab -> test_optimisation_ga.m

* Post analysis
Run Matlab -> matlab_population_analysis.m (Example output in res_lim_800)

