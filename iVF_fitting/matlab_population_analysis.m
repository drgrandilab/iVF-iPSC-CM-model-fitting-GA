clear; 
% close all;
lim = 1.0e3;
pop = [];
fitness=[]
Exp_Ca_1Hz = load('exp_data/iVF_1Hz.dat');
Exp_Ca = load('exp_data/iVF_0.5Hz.dat');

for i = [1:50]
    rng_number = i;

    filename = sprintf('res_lim_800/normal_population_seed_%i.mat', rng_number);

    d = load(filename);
    % x =d.x;
    [ft,index] = min(d.costs);
    tmp = d.population;

    solu_best = tmp(index,:);
    % d = readtable(filename);
    % x = d{:,2};
    % ft = mean(d{:,3});
    % pop = [pop; x'];
    


    pop = [pop; solu_best];
    fitness = [fitness;ft];

end

length(fitness)


X = ones(size(pop));

[XX, YY] = size(X);

for i = 1:XX
    X(i,:) = 1:YY;
end

tick_labels = {'Gto', 'GK1', 'GNaca', 'Gserca', 'GNaK', 'CaL_v_shift', 'G_CaL', 'GCat', 'Ca_Buffer', 'ec50SR', 'Kmf', 'Ina_shift', 'GbNa', 'GNa', 'Gf', 'GKr', 'GPCa', 'GbCa'}
index= 1:length(fitness)

figure(20), 
swarmchart(X,pop(index,:));

xticks(1:18)
xticklabels(tick_labels)
hold on
ylabel('log(Parameter)')



 
% load('C19911_T13_w3_RA1_filtered.mat');

% Exp_Ca = align_RA_trace;

data = {};



% figure; hold on;
% plot(dca_2,'Marker','o', 'color', [0.00,0.45,0.74], 'MarkerSize',5, 'MarkerFaceColor', [0.00,0.45,0.74]);





index = [];

sim_err =[];
saved = {};
for i = 1:1:XX

    % [err, output] = Cost_Function_whole_trace_debug(pop(i,:), Exp_Ca);
    [err, outputs] = Cost_Function(pop(i,:), Exp_Ca, Exp_Ca_1Hz, 0);
    sim_err(i) = err;

   
    if(err < lim)

            figure(1);
        output = outputs{1};
         sim_vm=output.AP;
        sim_ca=output.Ca;
        norm_sim_Ca=output.Ca_norm;
            subplot(3,1,1); plot(norm_sim_Ca,'color',[1,0.,0., 0.5], 'LineWidth',1.5); hold on
             ylabel('Norm. Ca') ; box off
            subplot(3,1,2);hold on;  plot(1e6*sim_ca,'color',[0.5,0.5,0.5, 0.5], 'LineWidth',1.5 ); ylabel('Ca (nM)');, box off

            subplot(3,1,3);hold on; plot(sim_vm,'color',[0.5,0.5,0.5, 0.5] , 'LineWidth',1.5); ylabel('V_m (mV'); xlabel('Time (ms)'), box off


            figure(2);
        output = outputs{2};
         sim_vm=output.AP;
        sim_ca=output.Ca;
        norm_sim_Ca=output.Ca_norm;
            subplot(3,1,1); plot(norm_sim_Ca,'color',[1,0.,0., 0.5], 'LineWidth',1.5); hold on
             ylabel('Norm. Ca') ; box off
            subplot(3,1,2);hold on;  plot(1e6*sim_ca,'color',[0.5,0.5,0.5, 0.5], 'LineWidth',1.5 ); ylabel('Ca (nM)');, box off

            subplot(3,1,3);hold on; plot(sim_vm,'color',[0.5,0.5,0.5, 0.5] , 'LineWidth',1.5); ylabel('V_m (mV'); xlabel('Time (ms)'), box off



        index(end+1) = i;

        output.err = err;
        saved{i} = outputs;

    end
end

figure(1);
subplot(3,1,1); hold on; plot( norm_percentile(Exp_Ca), 'color',[0,0,0, 1], 'LineWidth',2.5);

figure(2);
subplot(3,1,1); hold on; plot( norm_percentile(Exp_Ca_1Hz), 'color',[0,0,0, 1], 'LineWidth',2.5);





% save('population_analysis_data.mat');


% save('population_analysis_data.mat');



%% functionname: function description
function outputs = analyzed(data, lim)

% % return
% index = (data.costs<=lim);
% outputs = data.population(index,:);
%
% costs = data.costs(index);
%
% [costs, Ind] = sort(costs,'ascend');
% len = min(10, length(Ind));
% outputs = outputs(Ind(1:len),:);
%
% size(outputs);
if(min(data.costs) <=lim)
    outputs = data.x;
else
    outputs=[];
end
end











% evaluate your const function here.

% place targeting parameters, including APD, RMP, pleteau potential.






    function out = normalise(data)
    if (max(data) - min(data) ~= 0)
        out= (data - min(data)) ./ (max(data) - min(data));
    else
        data
    end
    end



    function out = normalise_exp_normal(data)
    if (max(data) - min(data) ~= 0)
        out=  (data - (-1.1)) / (3.85 - (-1.1));
    else
        data
    end
    end


    function out = normalise_exp_iVF(data)
    if (max(data) - min(data) ~= 0)
        out=  (data - (0.64)) / (4.25 - (0.64));
    else
        data
    end
    end
%% normalize to percentile data at [5, 95] instead of simple max/min
function [outputs] = norm_percentile(inputs)
out_min  = prctile(inputs,5);
out_max  = prctile(inputs,98);

outputs = (inputs - out_min) ./ (out_max - out_min);
end
