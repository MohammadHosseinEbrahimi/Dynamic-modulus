function [] = biomech_analysis_mohammad_dyn(~)
%%%%%%%%%%%%%%%%%%%%%%%
% This script is to the analysis Dynamic modulus and phase differnce of materials basede on 
% dynamic sinusoidal mechanical tests in indentation geometry.  
% 
%
% Script loads the dynamic testing data, fits a sinusoidal function to the
% data and further calculates the dynamic Youngs modulus and phase difference
% for used frequencies. It also modifies the values based on Hayes correction factor
%
% Parameters such as dynamic frequencies, indenter size, sample thickness, number of cycles, 
% units of load and displacement etc must be modified based on measurement protocol
% 
% Mohammadhossein Ebrahimi 14.2.2018
%
%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc

disp('Pick up the biomechanical testing data folder including both static and dynamic loading.')
folder = uigetdir(pwd, 'Pick up the biomechanical testing data folder including both static and dynamic loading.');

folder_stress_strain = strcat(folder, '\01_Stress_Strain');
folder_dynamic = strcat(folder, '\02_Dynamic');

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%   DYNAMIC MODULUS    %%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


cd(folder_dynamic)
d = dir(['*.dyn']);
files = {d(1:end,1).name}';
      
d = dir();
files = {d(1:end,1).name}';
files(ismember(files,{'.','..'})) = [];
filename1 = files(1);
filename = strrep(filename1, '_0.dyn', '');

r = (0.730e-3)/2;                           % indenter radius (m)
A = pi*(r^2);                               % indenter area (m^2)
freqs = [0.005 0.05 0.1 0.25 0.5 0.625 0.833 1 2];     % Indentation measurement frequencies [Hz]

cycs = 4;                                   % nro of cycles
steps = 4;

% thickness 
thickness = 2.2810; % (mm)
h = thickness*10^-3;

for kk = 1:size(files,1);

    
% Load data:
data2 = [];
tempFile = fopen(char(files(kk)),'r');
tempData = fread(tempFile, [1 inf]);
tempBeg = findstr('position, um', tempData);
fseek(tempFile, tempBeg,-1);
fgetl(tempFile);
datarow = sscanf(fgetl(tempFile),'%f %f %f')';
while isempty(datarow)==0
    data2 = [data2; datarow];
    datarow = sscanf(fgetl(tempFile),'%f %f %f')';
end

% Set new Offset:
displ_temp = (data2(:,1) - mean(data2(:,1)))*1e-6;
load_temp = data2(:,2) - mean(data2(:,2));

pts = length(displ_temp)/cycs;                             % measurement points in each cycle
    
% Define:
time_per_cycle = 1./freqs(kk);              % Time taken by one measurement (sine) cycle
tot_time = cycs*time_per_cycle;             % Total time taken by the measurement
tt = linspace(0,tot_time,cycs*pts)';        % Create time vector
Fs = cycs*pts/tot_time;                     % Sampling frequency



% ignore the first data points (usually prone to error)
displ = displ_temp(2:end);     % displacement 
load = load_temp(2:end);
t = tt(2:end);

% Convert load data to force, from gram to Newton
force = load*9.81*1e-3; % F = ma , where a = g = 9.81 (kgm/s^2), [load] = g, [F] = N

% Fit simple first order sinusoidal curve to the data:
[fitresult_disp, gof_disp] = createFit_displ(t, displ, kk);
[fitresult_force, gof_force] = createFit_force(t, force, kk);

fit_disp_R2(kk) = gof_disp.rsquare;
fit_force_R2(kk) = gof_force.rsquare;

disp_fit = fitresult_disp.a1*sin(fitresult_disp.b1*t+fitresult_disp.c1);
force_fit = fitresult_force.a1*sin(fitresult_force.b1*t+fitresult_force.c1);

fit_param_disp(kk,:) = [fitresult_disp.a1 fitresult_disp.b1 fitresult_disp.c1];
fit_param_force(kk,:) = [fitresult_force.a1 fitresult_force.b1 fitresult_force.c1];


% Calculate amplitudes and phase shifts in a subfunction "amplitude":
[delta_d(kk) , delta_F(kk), phase_lag(kk)] = amplitude(force,displ,t,freqs(kk) , files(kk) , pts);
[delta_d_fit(kk) , delta_F_fit(kk), phase_lag_fit(kk)] = amplitude(force_fit,disp_fit,t,freqs(kk) , files(kk) , pts);

phase_lag_deg(kk) = phase_lag(kk)*360/(2*pi);
phase_lag_fit_deg(kk) = phase_lag_fit(kk)*360/(2*pi);


Edyn(kk) = ((delta_F(kk)*1e-6)./A)./delta_d(kk);
Edyn_fit(kk) = ((delta_F_fit(kk)*1e-6)./A)./delta_d(kk);


% In comparison, calculate phase shift in the frequency domain in a subfunction:
phase_lag_freq(kk) = phase_shift(displ, force, Fs);         % rad
phase_lag_freq_deg(kk) = phase_lag_freq(kk)*360/(2*pi);     % degrees

% Phase shift OF THE SINE-FIT in the frequency domain
phase_lag_freq_fit(kk) = phase_shift(disp_fit, force_fit, Fs);  % rad
phase_lag_freq_fit_deg(kk) = phase_lag_freq_fit(kk)*360/(2*pi);     % degrees

% Phase shift between diplacement and measured load acquired from the sine-fit:
phase_disp(kk) = fitresult_disp.c1;                                % phase of displacement
phase_force(kk) = fitresult_force.c1;                              % phase of load (force)
delta_phase(kk) = abs((phase_disp(kk) - phase_force(kk)));         % phase difference between diplacement and load (force)
delta_phase_deg(kk) = delta_phase(kk)*360/(2*pi);


v = 0.5;        % Poisson's ratio, determined to be 0.5 in a dynamic loadingm due to non-compresible materials
kappa(kk) = kappa_value(r, h*0.95.^steps , v);
Edyn_cor(kk) = ((1-v.^2)*delta_F(kk)) ./ (2*kappa(kk)*r*delta_d(kk))/1e6;
Edyn_cor_fit(kk) = ((1-v.^2)*delta_F_fit(kk)) ./ (2*kappa(kk)*r*delta_d(kk))/1e6;

% Calculate storage and loss modulae to the dynamic loading:
Estorage(kk) = Edyn_cor(kk)*cos(phase_lag_freq(kk));
Eloss(kk) = Edyn_cor(kk)*sin(phase_lag_freq(kk));

% In comparison, calculate storage and loss modulae to the dynamic loading also to the sine-fitted data:
Estorage_fit(kk) = Edyn_cor(kk)*cos(phase_lag_fit(kk));
Eloss_fit(kk) = Edyn_cor(kk)*sin(phase_lag_fit(kk));

    
% Plot fits with datas:
h2 = figure('Name', ['Displacement frequency: ' , num2str(freqs(kk)), ' Hz'], ...
    'units','normalized','position',[0.5 0.046 0.5 0.88]);
% figure('Name', ['Displacement frequency: ' , num2str(freqs(kk)), ' Hz']); 
subplot 211, 
plot( fitresult_disp, tt, displ_temp );
title(['Fit (R^2): ', num2str(fit_disp_R2(kk))])
legend('Displ vs. t', 'Sine-fit', 'Location', 'Best' );
% Label axes
xlabel( 't [s]' );
ylabel( 'Displacement [?m]' );
subplot 212, 
plot( fitresult_force, t, force);
title(['Fit (R^2): ', num2str(fit_force_R2(kk))])
legend('Force vs. t', 'Sine-fit', 'Location', 'Best' );
% Label axes
xlabel( 't [s]' );
ylabel( 'Force [N]' );    
    
disp([' ']);
disp([' -------------------------------- ']);
disp(['File: ' , char(files(kk)) , ', Frequency: ' , num2str(freqs(kk))])
% disp(['maximum differences of measured values:']);
%     disp([tempNimi]);
% disp(['Displacement: ', num2str(delta_d(kk).*1e-3) ,' m'])
disp(['maximum difference of measured Force: ' , num2str(delta_F(kk)), ' N'  ]);
% disp(['Force: ', num2str(delta_F(kk)), ' N' ])
% disp(['Phase of displacement: ', num2str(phase_disp(kk)) , ' rad, phase of force: ' , num2str(phase_force(kk)) ,' rad'])
% disp(['Phase difference: ', num2str(delta_phase(kk)), ' rad'])

disp(['Phase lag (degrees): ' , num2str(phase_lag_deg(kk)) , ' ?'])
disp(['Phase lag of the sine-fit (degrees): ' , num2str(phase_lag_fit_deg(kk)) , ' ?'])
disp(['Phase lag calculated in freq. domain (degrees): ' , num2str(phase_lag_freq_deg(kk)), ' ?'])
disp(['Phase lag of SINE-FIT calculated in freq. domain (degrees): ' , num2str(phase_lag_freq_fit_deg(kk)), ' ?'])

% disp([num2str(delta_d(kk).*10^-3), ' m, ', num2str(delta_F(kk)), 'N']);
% disp(['Dynamic modulus without correction: ', num2str(Edyn(kk)), ' MPa']);
disp(['Dynamic modulus with correction: ', num2str(Edyn_cor(kk)), ' MPa']);

% disp(['Dynamic modulus from SINE-fit without correction: ', num2str(Edyn_fit(kk)), ' MPa']);
disp(['Dynamic modulus from SINE-fit with correction: ', num2str(Edyn_cor_fit(kk)), ' MPa']);



% close file:
fclose(tempFile);

end

h3 = figure('units','normalized','position',[0.5 0.046 0.5 0.88]);
subplot(4,2,[1 2])
plot(freqs,Edyn,'ro--','Linewidth',2.5), hold on, plot(freqs,Edyn_fit,'bo--','Linewidth',1.5), legend('Edyn - uncorrected','Edyn - sine fit - uncorrected','Location','Best')
xlabel('Frequency [Hz]'), ylabel('Dynamic modulus [MPa]')
subplot(4,2,[3 4]), 
plot(freqs, Edyn_cor,'ro--','Linewidth',2.5), hold on, plot(freqs, Edyn_cor_fit,'bo--','Linewidth',1.5), legend('Edyn - corrected','Edyn - sine fit - corrected','Location','Best')
xlabel('Frequency [Hz]'), ylabel('Dynamic modulus [MPa]')

subplot(4,2,5), 
plot(freqs, Estorage,'ro--','Linewidth',2.5), hold on, plot(freqs, Estorage_fit,'bo--','Linewidth',1.5), legend('E_s_t_o_r_a_g_e','E_s_t_o_r_a_g_e - sine fit','Location','Best')
xlabel('Frequency [Hz]'), ylabel('E_s_t_o_r_a_g_e [MPa]')

subplot(4,2,6), 
plot(freqs, Eloss,'ro--','Linewidth',2.5), hold on, plot(freqs, Eloss_fit,'bo--','Linewidth',1.5), legend('E_l_o_s_s','E_l_o_s_s - sine fit','Location','Best')
xlabel('Frequency [Hz]'), ylabel('E_l_o_s_s [MPa]')

% PHASE IN DEGREES: 
subplot(4,2,[7 8]), 


plot(freqs,phase_lag_freq_deg,'ro--','Linewidth',2), hold on, 
plot(freqs,phase_lag_freq_fit_deg,'bx--','Linewidth',2), 
legend('Freq. domain', 'Freq. domain FIT' ,'Location','Best')


xlabel('Frequency [Hz]'), ylabel('Phase shift \Delta\theta [degrees]')


%%%%% EXPORT DATA INTO .xlsx-file
caption2 = {};
caption2(1,1) = filename; 
caption2(2,1) = {'Frequency [Hz]'}; 
caption2(3,1) = {'Max difference of force '};
caption2(4,1) = {'Max difference of displacement '};

caption2(5,1) = {'Dynamic modulus of the sine-fit (MPa)'};
caption2(6,1) = {'Storage modulus (Es) - sine-fit (MPa)'};
caption2(7,1) = {'Loss modulus (El) - sine-fit (MPa)'};


caption2(8,1) = {'Dynamic modulus (MPa)'};
caption2(9,1) = {'Storage modulus (Es) (MPa)'};
caption2(10,1) = {'Loss modulus (El) (MPa)'};


caption2(11,1) = {'Phase lag calculated in freq. domain (?)'};
caption2(12,1) = {'Phase lag of SINE-FIT calculated in freq. domain (?)'};


caption2(13,1) = {'Phase lag '};
caption2(14,1) = {'Phase lag of the sine-fit '};

caption2(15,1) = {'Goodness of sine-fit, displacement R^2'};
caption2(16,1) = {'Goodness of sine-fit, force R^2'};

caption2(17,1) = {'Kappa-value'};

for mm = 1:length(freqs)
    caption2(2,mm+1) = {freqs(mm)};
    caption2(3,mm+1) = {delta_F(mm)};
    caption2(4,mm+1) = {delta_d(mm)};
    
    caption2(5,mm+1) = {Edyn_cor_fit(mm)};  
    caption2(6,mm+1) = {Estorage_fit(mm)};
    caption2(7,mm+1) = {Eloss_fit(mm)};
    
    caption2(8,mm+1) = {Edyn_cor(mm)};
    caption2(9,mm+1) = {Estorage(mm)};
    caption2(10,mm+1) = {Eloss(mm)};
    
    caption2(11,mm+1) = {phase_lag_freq_deg(mm)};
    caption2(12,mm+1) = {phase_lag_freq_fit_deg(mm)};
    
    
    caption2(13,mm+1) = {phase_lag(mm)};
    caption2(14,mm+1) = {phase_lag_fit(mm)};
    
    
    
    caption2(15,mm+1) = {fit_disp_R2(mm)};
    caption2(16,mm+1) = {fit_force_R2(mm)};
    
    caption2(17,mm+1) = {kappa(mm)};
end
    
caption2(18,1) = {['thickness (m): ']}; caption2(18,2) = {[num2str(h)]};

cd(folder)
xlswrite(strcat(char(filename), '.xlsx'), caption2 , 'Dyn. mod.') 

% Save last image:
saveas(h3 , char(filename) , 'fig')  % as a Matlab-file
saveas(h3 , char(filename) , 'png')  % as a regular image


% In addition, save data into struct-format for further analysis:
Biomech_data_dynamic.freqs = freqs;
Biomech_data_dynamic.delta_F = delta_F;
Biomech_data_dynamic.delta_d = delta_d;
Biomech_data_dynamic.phase_lag = phase_lag;
Biomech_data_dynamic.phase_lag_fit = phase_lag_fit;
Biomech_data_dynamic.phase_lag_freq_domain_rad = phase_lag_freq;
Biomech_data_dynamic.phase_lag_freq_domain_deg = phase_lag_freq_deg;

Biomech_data_dynamic.phase_lag_freq_domain_rad = phase_lag_freq_fit;
Biomech_data_dynamic.phase_lag_freq_domain_deg = phase_lag_freq_fit_deg;

Biomech_data_dynamic.fit_param_disp = fit_param_disp;
Biomech_data_dynamic.fit_param_force = fit_param_force;

Biomech_data_dynamic.fit_disp_R2 = fit_disp_R2;
Biomech_data_dynamic.fit_force_R2 = fit_force_R2;
Biomech_data_dynamic.Kappa_factor = kappa;
Biomech_data_dynamic.Young_cor = Edyn;
Biomech_data_dynamic.Young_cor_fit = Edyn_fit;
Biomech_data_dynamic.Estorage = Estorage;
Biomech_data_dynamic.Eloss = Eloss;
Biomech_data_dynamic.Estorage_fit = Estorage_fit;
Biomech_data_dynamic.Eloss_fit = Eloss_fit;


Biomech_data_dynamic.Thickness = h;
Biomech_data_dynamic.name = filename;

% save data and information to .mat-file
save(strcat(char(filename),'_dyn.mat'),'-struct','Biomech_data_dynamic');




end


%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % % % % % % % % % SEPARATE SUBFUNCTIONS % % % % % % % % % % % % %

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%  CALCULATE PHASE SHIFT IN FREQUENCY DOMAIN  %%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [ theta , amp_xx , amp_yy] = phase_shift( xx, yy , Fs)
%   phase_shift calculates the phase shift of two signals xx and yy.
% 
% 
%   Script calculates the phase shift between to measured signals: in this
%   case indenter displacement [m] and measured load [g].
%
%   Input arguments:
%   Signals xx (Displacement) and yy (Force), in addition sampling frequency Fs [Hz]
%
%   
% 
%   Written by Mohammadhossein Ebrahimi
%   
%%%%%%%%%

% take the FFT
x=fft(xx);
y=fft(yy);

npts = length(xx);


% Calculate the numberof unique points
NumUniquePts = ceil((npts+1)/2);

f = (0:NumUniquePts-1)*Fs/npts;

% Determine the max value and max point.
% This is where the sinusoidal
% is located. 
[mag_x idx_x] = max(abs(x));
[mag_y idx_y] = max(abs(y));

% determine the phase difference
% at the maximum point.
px = angle(x(idx_x));
py = angle(y(idx_y));
phase_lag = py - px;
% determine the amplitude scaling

theta = phase_lag;
amp_xx = mag_x;
amp_yy = mag_y;

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%  CALCULATE AMPLITUDES AND SHIFT   %%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [ delta_d, delta_F , theta_t] = amplitude( stress , strain , t , freqs , name, pts)
% 	Calculates the amplitudes of displacement and force applied to the
% 	a sample.
% 
% 
%   Input:
%   - applied force (stress)
%   - motor displacement (strain)
%   - time vector (t)
%   - sampling frequency (freq)
%   - name of the sample (name)
% 
% 
%   Output:
%   - displacement amplitude (delta_d)
%   - force amplitude (delta_F)
%   - phase shift between displacement and force (theta_t)
%

cycs = 4;       % cycles of sine wave used in the test
force = stress;
displ = strain;
kk = 1;

for ll = 1:cycs;
    
    if ll == 1;
    
        [max_d, ind_max_d] = max(displ(1+pts*(ll-1):pts*ll-1)); [min_d, ind_min_d] = min(displ(1+pts*(ll-1):pts*ll-1));
        [max_F, ind_max_F] = max(force(1+pts*(ll-1):pts*ll-1)); [min_F, ind_min_F] = min(force(1+pts*(ll-1):pts*ll-1));

        ma_d(ll) = max_d; ima_d(ll) = ind_max_d; mi_d(ll) = min_d; imi_d(ll) = ind_min_d;
        ma_F(ll) = max_F; ima_F(ll) = ind_max_F; mi_F(ll) = min_F; imi_F(ll) = ind_min_F;


    
    elseif ll <= 4;
    
        [max_d, ind_max_d] = max(displ(1+pts*(ll-1):pts*ll-1)); [min_d, ind_min_d] = min(displ(1+pts*(ll-1):pts*ll-1));
        [max_F, ind_max_F] = max(force(1+pts*(ll-1):pts*ll-1)); [min_F, ind_min_F] = min(force(1+pts*(ll-1):pts*ll-1));

        ma_d(ll) = max_d; ima_d(ll) = ind_max_d+pts*(ll-1); mi_d(ll) = min_d; imi_d(ll) = ind_min_d+pts*(ll-1);
        ma_F(ll) = max_F; ima_F(ll) = ind_max_F+pts*(ll-1); mi_F(ll) = min_F; imi_F(ll) = ind_min_F+pts*(ll-1);

    else
        disp('Something went wrong in the calculations of max/min')

    end
    
end

% Calculate amplitudes:
amp_d(kk,:) = abs(ma_d - mi_d); delta_d(kk) = mean(amp_d(kk,:));
amp_F(kk,:) = abs(ma_F - mi_F); delta_F(kk) = mean(amp_F(kk,:));

% Phase shift between load and displacement from the second cycle:
t_d = t(ima_d(2));
t_F = t(ima_F(2));
delta_t(kk) = abs(t_d - t_F);       % time shift (s)
theta_t(kk) = 2*pi*freqs(kk)*delta_t(kk);


% Plot raw data:
    h = figure('units','normalized','position',[0.005 0.04 0.492 0.88], 'name' , char(name));
    subplot(3,1,1);
    plot(t , displ); hold on, plot(t(ima_d), ma_d, 'ro'), plot(t(imi_d), mi_d, 'ro')
    ylabel('Displacement (m)'); xlabel('time (s)'); title(['Amplitude: ' , num2str(delta_d(kk)), ' m']);
    subplot(3,1,2);
    plot(t , force); hold on, plot(t(ima_F), ma_F, 'ro'), plot(t(imi_F), mi_F, 'ro')
    ylabel('Force (N)'); xlabel('time (s)'); title(['Amplitude: ' , num2str(delta_F(kk)) , ' N'])
    subplot(3,1,3);
    hold on
    [ax,p1,p2]=plotyy(t,displ,t,force);
    hold off
    ylabel(ax(1),'Displacement (?m)');ylabel(ax(2),'Force (N)'); xlabel('t (s)');
    title(['Phase shift: ' , num2str(theta_t) , ' rad'])
    set(p1,'Color','blue')
    set(p2,'Color','red')
    set(ax(1),'YColor','blue')
    set(ax(2),'YColor','red')
 
    
% Rearrange data such that, 
% First four rows, 1st column: maximum values, 2nd column: indexes of the data
% Second four rows, 1st column: minimum values, 2nd column: indexes of the data

Displacement = [ma_d' ima_d'; mi_d' imi_d'];        
Force = [ma_F' ima_F' ; mi_F' imi_F'];
results = [Displacement; Force];

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%  CALCULATE KAPPA VALUE %%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [ kappa ] = kappa_value( radius, thickness , v )
%   This function calculates tke KAPPA-value to the calculation of effective
%   Young's modulus with respect to the sample thickness [m].
%   Data to the calculations is taken from Zhang M. et al 1997. as 1st and 2nd
%   order polynomial fit is used to interpolate the needed KAPPA-value.
%
%   Effective Young's modulus:
%   Es = (1-v^2)*pi*a* E / 2*KAPPA*thickness
%   where, E = measured Young's modulus, a = radius of the indenter,
%   v = Poisson's ratio
%   thickness = sample thickness [m]
% 
%   Poisson's ratio is determined to be v = 0.1 in equilibrium modulus and 
%   v = 0.5 in dynamic modulus
%
%
%

a = (radius);                         % indenter radius [m]


if v == 0.3
    ah = [0 : 0.2: 2]';
    % Hayes et al, 1972
    poisson = [1000 1207 1472 1784 2124 2480 2845 3214 3586 3960 4336]'.*10^-3; 
   
elseif v == 0.5
    ah = [0 : 0.2: 2]';
    % Hayes et al, 1972
    poisson = [1000 1268 1645 2129 2704 3359 4085 4878 5737 6659 7644]'.*10^-3;
else
    disp('Poissons coefficient you gave, is not v = 0.1 or v = 0.5. Please, give one of those. ')
    
end

spline_ah = 0:0.001:2;
%interpolate cubic spline for kappa values presented in Hayes et al. 1976 (or so)
spline_poisson = spline(ah,poisson,spline_ah);
%if a/h value is not excatly in the cubic spline data points, use linear
%interpolation for the kappa value between cubic spline data points 
%(error here is negligible due to the dense sampling in spline)
interp_poisson = interp1(spline_ah, spline_poisson, a/thickness);

%cubic spline approximated kappa
kappa = interp_poisson;

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%  FIT SINE-WAVE TO THE DISPLACEMENT DATA %%% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fitresult, gof] = createFit_displ(t, displ , ind)
% Script is made to fit a sine function to the biomechanical dynamic
% loading data with loading frequenscies 
%
% [0.005, 0.05,0.1, 0.25, 0.5, 0.625, 0.833, 1, 2] [Hz]
%
% Separate fitting options has to be used in the different frequencies as
% well as in load and displacement datas. Script is based on Matlab's
% cftool-toolbox.
%
% Fitted sine-wave:
%   a1*sin(b1*x+c1)
% 
% 



if ind == 1;
    %% Fit: 'untitled fit 1'.
    
% tt = linspace(t(1), t(end), length(displ));
% tt = tt(:);

[xData, yData] = prepareCurveData( t, displ );


% Set up fittype and options.
ft = fittype( 'sin1' );
opts = fitoptions( ft );
opts.Display = 'Off';
opts.Lower = [-Inf 0 -Inf];
opts.StartPoint = [5.16243388021594e-05 0.0314946631938826 -0.042456161571354];
opts.Upper = [Inf Inf Inf];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );


elseif ind == 2;

    
% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( t, displ );

% Set up fittype and options.
ft = fittype( 'sin1' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [-Inf 0 -Inf];
opts.StartPoint = [5.21462854887703e-05 0.315147187577089 -0.0936428853429026];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );
   

elseif ind == 3;


[xData, yData] = prepareCurveData( t, displ );

% Set up fittype and options.
ft = fittype( 'sin1' );
opts = fitoptions( ft );
opts.Display = 'Off';
opts.Lower = [-Inf 0 -Inf];
opts.StartPoint = [6.30803378227838 0.636373896496394 -0.315063752780712];
opts.Upper = [Inf Inf Inf];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );


elseif ind == 4;


[xData, yData] = prepareCurveData( t, displ );

% Set up fittype and options.
ft = fittype( 'sin1' );
opts = fitoptions( ft );
opts.Display = 'Off';
opts.Lower = [-Inf 0 -Inf];
opts.StartPoint = [6.26972991827161 1.59093474124098 -0.40226863701901];
opts.Upper = [Inf Inf Inf];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

    
elseif ind == 5;

    
[xData, yData] = prepareCurveData( t, displ );

% Set up fittype and options.
ft = fittype( 'sin1' );
opts = fitoptions( ft );
opts.Display = 'Off';
opts.Lower = [-Inf 0 -Inf];
opts.StartPoint = [6.29972893415309 3.18186948248197 -0.482201885584438];
opts.Upper = [Inf Inf Inf];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

    
elseif ind == 6;


[xData, yData] = prepareCurveData( t, displ );

% Set up fittype and options.
ft = fittype( 'sin1' );
opts = fitoptions( ft );
opts.Display = 'Off';
opts.Lower = [-Inf 0 -Inf];
opts.StartPoint = [6.29908706745555 3.97733685310246 -0.533506899326938];
opts.Upper = [Inf Inf Inf];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

    
elseif ind == 7;    


[xData, yData] = prepareCurveData( t, displ );

% Set up fittype and options.
ft = fittype( 'sin1' );
opts = fitoptions( ft );
opts.Display = 'Off';
opts.Lower = [-Inf 0 -Inf];
opts.StartPoint = [6.29059276343008 5.30099455781496 -0.623661408410145];
opts.Upper = [Inf Inf Inf];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

    
elseif ind == 8;    

    
[xData, yData] = prepareCurveData( t, displ );

% Set up fittype and options.
ft = fittype( 'sin1' );
opts = fitoptions( ft );
opts.Display = 'Off';
opts.Lower = [-Inf 0 -Inf];
opts.StartPoint = [6.26524068681652 6.36373896496394 -0.669830014856764];
opts.Upper = [Inf Inf Inf];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

    
elseif ind == 9;    

    
[xData, yData] = prepareCurveData( t, displ );

% Set up fittype and options.
ft = fittype( 'sin1' );
opts = fitoptions( ft );
opts.Display = 'Off';
opts.Lower = [-Inf 0 -Inf];
opts.StartPoint = [3.20857855532133 12.7274779299279 -1.01774769335157];
opts.Upper = [Inf Inf Inf];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );
  
    
end


end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%  FIT SINE-WAVE TO THE FORCE DATA %%% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fitresult, gof] = createFit_force(t, force , ind)
% Script is made to fit a sine function to the biomechanical dynamic
% loading data with loading frequenscies 
%
% [0.005, 0.05,0.1, 0.25, 0.5, 0.625, 0.833, 1, 2] [Hz]
%
% Separate fitting options has to be used in the different frequencies as
% well as in load and displacement datas. Script is based on Matlab's
% cftool-toolbox.
%
% Fitted sine-wave:
%   a1*sin(b1*x+c1)
% 


if ind == 1;


[xData, yData] = prepareCurveData( t, force );

% Set up fittype and options.
ft = fittype( 'sin1' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [-Inf 0 -Inf];
opts.StartPoint = [0.249053424600289 0.0314946631938826 0.260611199559208];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );


elseif ind == 2;

% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( t, force );

% Set up fittype and options.
ft = fittype( 'sin1' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [-Inf 0 -Inf];
opts.StartPoint = [0.280552725144696 0.0790337774487998 0.0198032528209321];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

    
    
elseif ind == 3;

% Fit: ' Sine 1'.
[xData, yData] = prepareCurveData( t, force );

% Set up fittype and options.
ft = fittype( 'sin1' );
opts = fitoptions( ft );
opts.Display = 'Off';
opts.Lower = [-Inf 0 -Inf];
opts.StartPoint = [0.388042080225689 0.636373896496394 -0.0124630344066262];
opts.Upper = [Inf Inf Inf];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );



elseif ind == 4;
    
[xData, yData] = prepareCurveData( t, force );

% Set up fittype and options.
ft = fittype( 'sin1' );
opts = fitoptions( ft );
opts.Display = 'Off';
opts.Lower = [-Inf 0 -Inf];
opts.StartPoint = [0.399185087991594 1.59093474124098 -0.120392061543909];
opts.Upper = [Inf Inf Inf];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );



elseif ind == 5;
    
[xData, yData] = prepareCurveData( t, force );

% Set up fittype and options.
ft = fittype( 'sin1' );
opts = fitoptions( ft );
opts.Display = 'Off';
opts.Lower = [-Inf 0 -Inf];
opts.StartPoint = [0.417254646329386 3.18186948248197 -0.206292901031662];
opts.Upper = [Inf Inf Inf];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );


elseif ind == 6;
    
[xData, yData] = prepareCurveData( t, force );

% Set up fittype and options.
ft = fittype( 'sin1' );
opts = fitoptions( ft );
opts.Display = 'Off';
opts.Lower = [-Inf 0 -Inf];
opts.StartPoint = [0.422328699950688 3.97733685310246 -0.272940069015889];
opts.Upper = [Inf Inf Inf];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );



elseif ind == 7;   
    
[xData, yData] = prepareCurveData( t, force );

% Set up fittype and options.
ft = fittype( 'sin1' );
opts = fitoptions( ft );
opts.Display = 'Off';
opts.Lower = [-Inf 0 -Inf];
opts.StartPoint = [0.42666401652757 5.30099455781496 -0.371303320850403];
opts.Upper = [Inf Inf Inf];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

    
elseif ind == 8;    
    
[xData, yData] = prepareCurveData( t, force );

% Set up fittype and options.
ft = fittype( 'sin1' );
opts = fitoptions( ft );
opts.Display = 'Off';
opts.Lower = [-Inf 0 -Inf];
opts.StartPoint = [0.429049254624203 6.36373896496394 -0.415591730762127];
opts.Upper = [Inf Inf Inf];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

    
elseif ind == 9;    
    
[xData, yData] = prepareCurveData( t, force );

% Set up fittype and options.
ft = fittype( 'sin1' );
opts = fitoptions( ft );
opts.Display = 'Off';
opts.Lower = [-Inf 0 -Inf];
opts.StartPoint = [0.209646656634676 12.7274779299279 -0.586316231588214];
opts.Upper = [Inf Inf Inf];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );
   
end

end





