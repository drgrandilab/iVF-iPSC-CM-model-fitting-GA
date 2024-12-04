clear;
addpath('lib/')


% global data_t;

% parpool(4)


% load('C19911_T13_w3_RA1_filtered.mat')


% global Exp_Ca;

% Exp_Ca = align_RA_trace;


Exp_Ca_1Hz = load('exp_data/iVF_1Hz.dat');
Exp_Ca = load('exp_data/iVF_0.5Hz.dat');

% here, add Exp_Ca as a function parameter
fun = @(x)Cost_Function(x,Exp_Ca,Exp_Ca_1Hz,0);
% Cost_Function_whole_trace(zeros(18,1),1);

% pause
% break

for i = 1:50
run_ga(i,fun,Exp_Ca)
end

%% run_ga: function description
function [outputs] = run_ga(rng_number,fun,Exp_Ca)


    rng(rng_number, "twister");
    nval = 18;
    values = log(5);

    lb = ones(1,nval)*(-values);  % 0.4
    ub = ones(1,nval)*values;  % 2.5


    % initial_population = lhsdesign(100,19)*1.6 + 0.4; % within range of [0.4,2.0]
    initial_population = lhsdesign(300,nval)*2*values - values;  % from  exp^  -log(2.5) to log(2.5)
    % ub = [Inf; Inf; Inf; Inf; Inf; Inf; Inf; Inf; Inf; Inf; Inf; Inf; Inf; Inf; Inf; Inf];
    % initial_population(1,:)=0;
    
    A = [];
    b = [];
    Aeq = [];
    beq = [];
    nonlcon = [];
    % IntCon= 1:16;%[]
     
        options = optimoptions('ga','PlotFcn','gaplotbestf','UseParallel',true,'PopulationSize',300,'MaxGenerations',50, 'FitnessLimit',800, ...
            'InitialPopulationMatrix', initial_population);
        % Inputs_Final = fminsearchbnd_2(@ras, Inputs, lb, ub, options);
        [x,fval,exitflag,output, population, costs] = ga(fun,nval,A,b,Aeq,beq,lb,ub,nonlcon,options);

filename = sprintf('normal_population_seed_%i.mat', rng_number);

save (filename, 'population', 'costs', 'x');
% delete(gcp('nocreate'));
outputs = x;
end
