function [Total_Error, output] = Cost_Function(para,varargin)

%Parallel loop of protocols to simulate:
% SSA = steady state availability; ACT = steady state activation; RUDB =
% recovery from use-dependent block; TAU = tau to 50% decay; RFI = recovery
% from 1 pulse inactivation

debug = 0;
output={};


% output.AP=[];
% output.Ca=[];
% output.Ca_norm=[];


% global Exp_Ca


% load('C19911_T13_w7_RA1_filtered.mat')


% global Exp_Ca;


Exp_Ca = [];

if(length(varargin) == 1)
    % debug = varargin{2};
    Exp_Ca = varargin{1};
end



if(length(varargin) == 2)
    Exp_Ca = varargin{1};
    Exp_Ca_1Hz = varargin{2};
end

if(length(varargin) == 3)
    Exp_Ca = varargin{1};
    Exp_Ca_1Hz = varargin{2};
    debug = varargin{3};
end


% if(size(varargin) == 2)
%     % debug = varargin{2};
%     Exp_Ca = varargin{1};
%     debug = varargin{2};
% end


% para = ones(1,18)*0.1;
% para = para/10.0;

% np.sum(10.0 * np.abs(np.array(list(individual[:]))))
lasso_error = sum(10* abs(para));
para = exp(para);
spError = 0e4;
Total_Error=1e4;
Total_Error_05 = 2e4;
Total_Error_10 = 2e4;
Total_Error = 2e4;

try
    comand = ["./cvode_main_ipsc_baseline_05Hz", string(para), '1'];
    new_command = strjoin(comand, " ")
    [a,out] = system(new_command);
data = str2num(out);

[Total_Error_05, output_05] = get_error(data, Exp_Ca, debug);


comand = ["./cvode_main_ipsc_baseline_10Hz", string(para), '1'];
    new_command = strjoin(comand, " ")
    [a,out] = system(new_command);
data = str2num(out);

[Total_Error_10, output_10] = get_error(data, Exp_Ca_1Hz, debug);


output = {output_05, output_10};


Ca = output_05.Ca_norm;


dCa =  get_dVm(Ca);


[pos_pks_ca,pos_ca_ind] = findpeaks(dCa, 'MinPeakProminence',0.002,'MinPeakDistance',200, 'Annotate','extents'); % traces are normalized

if(length(pos_pks_ca)~=3) 
    Total_Error_05 = Total_Error_05 + 2000;
end

Total_Error = Total_Error_05 + Total_Error_10 + lasso_error;

catch
    disp('file name loaded incorrectely');
    Total_Error = 10e4;
end




% evaluate your const function here.

% place targeting parameters, including APD, RMP, pleteau potential.

fprintf ('Total_Error_05=%f; *******; Total_Error_10=%f; *******; Total_Error =%f\n', Total_Error_05,Total_Error_10,Total_Error);

end



function [Total_Error, output] = get_error(data, Exp_Ca, debug)


t_shift = 100e3;
    

    sim_t = data(:,1);
    sim_vm = data(:,2);
    sim_ca = data(:,4);

    APA = max(sim_vm) - min(sim_vm);
    AMP_Ca = max(sim_ca) - min(sim_ca);

    APA_Error_apr = 1000*calculate_error(APA, 125, 90, 105);

    CaAMP_Error_apr = 1000*calculate_error(AMP_Ca, 1e-3, 1e-4, 0.5e-3);



    sim_Ca_rise_index = get_median_up_time(normalise_percentile(sim_ca));
    sim_Ca_num = length(sim_Ca_rise_index);
    % sim_Ca_rise_index = sim_Ca_rise_index(1);

    Exp_Ca_rise_index = get_median_up_time(normalise_percentile(Exp_Ca));


    Exp_Ca_rise_index = Exp_Ca_rise_index(1);
    index_ca = Exp_Ca_rise_index;

    % dd_ca_2 = normalise(sim_ca(1:end-1)) - normalise(sim_ca(2:end));

    time_2_ca_vec = sim_t(sim_Ca_rise_index) ;
    [C,idx_ca] = min(abs((time_2_ca_vec-t_shift) - Exp_Ca_rise_index));
    [C,idx_ca_2] = min(abs(sim_t  - time_2_ca_vec(idx_ca)));
    in_ca_2 = idx_ca_2;
    length_ca = length(Exp_Ca);
    norm_Exp_Ca = normalise_percentile (Exp_Ca(index_ca-index_ca+1:index_ca+length_ca-index_ca-0));
    norm_sim_Ca = normalise(sim_ca(in_ca_2-index_ca+1:in_ca_2+length_ca-index_ca-0));

    Total_Error = sum(abs(norm_sim_Ca-norm_Exp_Ca)) + APA_Error_apr + CaAMP_Error_apr;

    output.AP = sim_vm(in_ca_2-index_ca+1:in_ca_2+length_ca-index_ca-0);
    output.Ca = sim_ca(in_ca_2-index_ca+1:in_ca_2+length_ca-index_ca-0);
    output.Ca_norm = norm_sim_Ca;
    output.Total_Error = Total_Error;

        % data = load('ys.txt.1');
    if(debug)
        figure; subplot(3,1,1); plot(norm_sim_Ca,'r','LineWidth',1.5); hold on
        plot(norm_Exp_Ca,'k','LineWidth',1.5); ylabel('Norm. Ca') ; 
        legend({'Sim.', 'Exp'}); box off
            % set(gca,"FontSize",20)
            % set(gca, 'LineWidth',1.5)
        subplot(3,1,2);  plot(1e6*sim_ca(in_ca_2-index_ca+1:in_ca_2+length_ca-index_ca-0),'r','LineWidth',1.5); ylabel('Ca (nM)');, box off
            % set(gca,"FontSize",20)
            % set(gca, 'LineWidth',1.5)
        subplot(3,1,3); plot(sim_vm(in_ca_2-index_ca+1:in_ca_2+length_ca-index_ca-0),'r','LineWidth',1.5); ylabel('V_m (mV'); xlabel('Time (ms)'), box off
            % set(gca,"FontSize",20)
            % set(gca, 'LineWidth',1.5)
        set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',15, 'LineWidth', 1.5); %'FontWeight','Bold',
    end
    % end

end





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
function [outputs] = normalise_percentile(inputs)
out_min  = prctile(inputs,5);
out_max  = prctile(inputs,98);

outputs = (inputs - out_min) ./ (out_max - out_min);
end


%% calculate_error: function description
function [outputs] = calculate_error(value, ub, lb, mean_v)
    if( value <= ub && value>=lb)
        outputs = 0;
    else 
        a = abs(value - ub)/abs(mean_v);
        b = abs(value - lb)/abs(mean_v);
        outputs = min(a,b);
    end
end




function [out] = get_median_up_time(Ca)

dCa =  get_dVm(Ca);


value = 0.5*(max(Ca) - min(Ca)) + min(Ca);



[pos_pks_ca,pos_ca_ind] = findpeaks(-abs(Ca-value), 'MinPeakProminence',0.2,'MinPeakDistance',300, 'Annotate','extents'); % traces are normalized

% out = pos_ca_ind;
out = [];
for i =1:length(pos_ca_ind)
    index=pos_ca_ind(i);
%     dCa(index)
    if(dCa(index) > 0) 
        out(end+1)=index;
    end
end

% out
% out


end


function dVm = get_dVm(Vm)


dVm =  Vm(2:end) - Vm(1:end-1) ;
dVm(end+1) = 0; % add 0 to keep length constant 

end