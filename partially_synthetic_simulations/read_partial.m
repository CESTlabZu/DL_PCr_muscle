%% Partially Synthetic Data Simulation
% 
% Author: Malvika Viswanathan
% Corresponding author: Dr. Zhongliang Zu
% Email: zhongliang.zu@vumc.org
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc;

%% ---------------------------------------------------------------------
% loading multiple pool lorentzian fit (mfit)  amine and MT curves here 
% for 0.5uT and 1uT
% load('...')
% load('...');

% loading sample mfit curves
load('sample_files_1p0uT.mat');

%% ---------------------------------------------------------------------
% Defining constants and ranges
% please make sure power of lorenztian fit files match the power used for
% simulations 
tt = 1 % B1 (uT)
tt_4p7T = [0.8, 0.9, 1, 1.1, 1.2]*0.5;   % B1 with shifts


offppm = 200;            % frequency of 4.7T (MHz)

% Frequency offsets of different pools
sep1 = 3.6*offppm; % APT
sep2 = 3.0*offppm; % amine
sep3 = 2.6*offppm; % PCr   
sep4 = 2*offppm; % Creatine
sep5 = -3.3*offppm; % NOE 


% relaxation constants required for simulation               
R1M   = 1/1.5;      % R1 MT (1/s)
R2S1 = 1/0.002;     % R2 APT (1/s)        
R2S5 = 1/0.0005;    % R2 NOE (1/s)

% constants required for simulation
fs5 = 0.015;        % fs NOE
ksw1 = 100;         % ksw APT
ksw5 = 10;          % ksw NOE

% Frequency offset 
offset = -1000:25:1000;         
k = [-2000,-1750,-1500,-1250,-1000:25:1000,1250,1500,1750,2000]; 
k_4p7T  = [k-100; k-50; k; k+50; k+100]; % B0 shifts

%% ---------------------------------------------------------------------
% measured components (MT and amine)
m_MT_n = sample_mt;

% extracting amine effect using polynomial fit
k_scaled   = [-2000,-1750,-1500,-1250,-1000:25:1000,1250,1500,1750,2000]/200;

mfit_arex  = sample_amine;

FitParam.PeakOffset = -3;                
arex_background = poly_fit(k_scaled, mfit_arex, FitParam);
amine_normal    = [mfit_arex(1:4)', arex_background, mfit_arex(38:end)'];

fm_n  = sample_fm;

%% ---------------------------------------------------------------------
% Define ranges for simulation

num_T1W  = 3;
num_T2W  = 3;
num_fs1  = 3;
num_fs2  = 3;
num_fs3  = 5;
num_fs4  = 3;
num_fm   = 3;
num_ksw3 = 6;
num_ksw4 = 3;
num_T2S3 = 3;
num_T2S4 = 3;
num_tt   = 5;
num_B0   = 5;

T1W_matrix   = 1.1:0.3:1.7;
T2W_matrix   = [20, 40, 60]*0.001;
T2S3_matrix  = [0.008, 0.010, 0.012];
T2S4_matrix  = [0.008, 0.010, 0.012];

fs1_matrix  = [0.0008, 0.001, 0.0012];
fs2_matrix  = [0.75, 1, 1.25];
fs3_matrix  = [0.0006, 0.0008, 0.001, 0.0012, 0.0014];
fs4_matrix  = [0.0016, 0.002, 0.0024];
fm_matrix = [0.5, 1, 1.5];

ksw3_matrix = [60, 80, 100, 120, 140, 160];
ksw4_matrix = [300, 500, 700];

%% ---------------------------------------------------------------------
% Main loop 

i = 1;
for ii_T1W = 1:num_T1W
    ii_T1W
    R1W_cal = 1./T1W_matrix(ii_T1W);
    for ii_T2W = 1:num_T2W
        ii_T2W
        R2W_cal = 1./T2W_matrix(ii_T2W);
        for ii_T2S3 = 1:num_T2S3
            R2S3_cal = 1./T2S3_matrix(ii_T2S3);
            for ii_T2S4 = 1:num_T2S4
                R2S4_cal = 1./T2S4_matrix(ii_T2S4);
                for ii_fs1 = 1:num_fs1
                    fs1 = fs1_matrix(ii_fs1);
                    for ii_fs2 = 1:num_fs2
                        fs2_cal = fs2_matrix(ii_fs2);
                        for ii_fs3 = 1:num_fs3
                            fs3_cal = fs3_matrix(ii_fs3);
                            for ii_fs4 = 1:num_fs4
                                fs4_cal = fs4_matrix(ii_fs4);
                                for ii_fm = 1:num_fm
                                    fm_cal = fm_matrix(ii_fm);
                                    for ii_ksw3 = 1:num_ksw3
                                        ksw3_cal = ksw3_matrix(ii_ksw3);
                                        for ii_ksw4 = 1:num_ksw4
                                            ksw4_cal = ksw4_matrix(ii_ksw4);
                                            for ii_tt = 1:num_tt
                                                tt_shift = tt_4p7T(ii_tt);
                                                for ii_kk = 1:num_B0
                                                    B0_shift = k_4p7T(ii_kk, :);

                                                    % Observed R1W 
                                                    R1W_cal_obs_n = (R1W_cal + fm_cal*fm_n*R1M) ./ (1 + fm_cal*fm_n);
        
                                                    % CEST and MT effects with B0 and B1 shifts
                                                    cal_Lorentzian1_cal  = (fs1.*ksw1.*(tt_shift.*42.6*2*pi).^2 ./ ((tt_shift.*42.6*2*pi).^2 + (R2S1+ksw1)*ksw1 + ksw1./(R2S1+ksw1).*((B0_shift+sep1)*2*pi).^2));
                                                    cal_Lorentzian2_n_cal = (interp1(k_4p7T(3,:), fs2_cal.*amine_normal, B0_shift) .* ((tt_shift)^2));
                                                    cal_Lorentzian3_n_cal = (fs3_cal.*ksw3_cal.*(tt_shift.*42.6*2*pi).^2 ./ ((tt_shift.*42.6*2*pi).^2 + (R2S3_cal+ksw3_cal)*ksw3_cal + ksw3_cal./(R2S3_cal+ksw3_cal).*((B0_shift+sep3)*2*pi).^2));
                                                    cal_Lorentzian4_n_cal = (fs4_cal.*ksw4_cal.*(tt_shift.*42.6*2*pi).^2 ./ ((tt_shift.*42.6*2*pi).^2 + (R2S4_cal+ksw4_cal)*ksw4_cal + ksw4_cal./(R2S4_cal+ksw4_cal).*((B0_shift+sep4)*2*pi).^2));
                                                    cal_Lorentzian5_n_cal = (fs5.*ksw5.*(tt_shift.*42.6*2*pi).^2 ./ ((tt_shift.*42.6*2*pi).^2 + (R2S5+ksw5)*ksw5 + ksw5./(R2S5+ksw5).*((B0_shift+sep5)*2*pi).^2));
                                                    cal_Lorentzian6_n_cal = (interp1(k_4p7T(3,:), fm_cal*m_MT_n, B0_shift)) .* ((tt_shift)^2);
                                                    
                                                    % Reff calculation
                                                    cal_eff_cal = R1W_cal_obs_n.*((B0_shift)*2*pi).^2./((tt_shift.*42.6*2*pi).^2 + ((B0_shift)*2*pi).^2) + R2W_cal.*(tt_shift.*42.6*2*pi).^2./((tt_shift.*42.6*2*pi).^2 + ((B0_shift)*2*pi).^2);
                                                    
                                                    % S/S_0 equation calculation = R1/(Reff + CEST + MT)
                                                    sscal_n = R1W_cal_obs_n ./ (cal_eff_cal + cal_Lorentzian1_cal./(1+fm_cal*fm_n) + cal_Lorentzian2_n_cal./(1+fm_cal*fm_n) + ...
                                                                               cal_Lorentzian3_n_cal./(1+fm_cal*fm_n) + cal_Lorentzian4_n_cal./(1+fm_cal*fm_n) + ...
                                                                               cal_Lorentzian5_n_cal./(1+fm_cal*fm_n) + cal_Lorentzian6_n_cal) .* (((B0_shift)*2*pi).^2 ./ ((tt_shift.*42.6*2*pi).^2 + ((B0_shift)*2*pi).^2));
        
                                                    SS_cal_n(:, ii_T1W, ii_T2W, ii_T2S3, ii_T2S4, ii_fs1, ii_fs2, ii_fs3, ii_fs4, ii_fm, ii_ksw3, ii_ksw4, ii_tt, ii_kk) = sscal_n;
        
                                                    SS_cal_value_ref_n = R1W_cal_obs_n ./ (cal_eff_cal + cal_Lorentzian1_cal./(1+fm_cal*fm_n) + cal_Lorentzian2_n_cal./(1+fm_cal*fm_n) + ...
                                                                                           0./(1+fm_cal*fm_n) + cal_Lorentzian4_n_cal./(1+fm_cal*fm_n) + ...
                                                                                           cal_Lorentzian5_n_cal./(1+fm_cal*fm_n) + cal_Lorentzian6_n_cal) .* ...
                                                                         (((B0_shift)*2*pi).^2 ./ ((tt_shift.*42.6*2*pi).^2 + ((B0_shift)*2*pi).^2));
        
                                                    SS_cal_value_ref_cr = R1W_cal_obs_n ./ (cal_eff_cal + cal_Lorentzian1_cal./(1+fm_cal*fm_n) + cal_Lorentzian2_n_cal./(1+fm_cal*fm_n) + ...
                                                                                            cal_Lorentzian3_n_cal./(1+fm_cal*fm_n) + 0./(1+fm_cal*fm_n) + ...
                                                                                            cal_Lorentzian5_n_cal./(1+fm_cal*fm_n) + cal_Lorentzian6_n_cal) .* ...
                                                                          (((B0_shift)*2*pi).^2 ./ ((tt_shift.*42.6*2*pi).^2 + ((B0_shift)*2*pi).^2));
        
                                                    cal_Lorentzian1_cal_ns  = (fs1.*ksw1.*(tt.*42.6*2*pi).^2 ./ ((tt.*42.6*2*pi).^2 + (R2S1+ksw1)*ksw1 + ksw1./(R2S1+ksw1).*((k+sep1)*2*pi).^2));
                                                    cal_Lorentzian2_n_cal_ns = amine_normal*fs2_cal;
                                                    cal_Lorentzian3_n_cal_ns = (fs3_cal.*ksw3_cal.*(tt.*42.6*2*pi).^2 ./ ((tt.*42.6*2*pi).^2 + (R2S3_cal+ksw3_cal)*ksw3_cal + ksw3_cal./(R2S3_cal+ksw3_cal).*((k+sep3)*2*pi).^2));
                                                    cal_Lorentzian4_n_cal_ns = (fs4_cal.*ksw4_cal.*(tt.*42.6*2*pi).^2 ./ ((tt.*42.6*2*pi).^2 + (R2S4_cal+ksw4_cal)*ksw4_cal + ksw4_cal./(R2S4_cal+ksw4_cal).*((k+sep4)*2*pi).^2));
                                                    cal_Lorentzian5_n_cal_ns = (fs5.*ksw5.*(tt.*42.6*2*pi).^2 ./ ((tt.*42.6*2*pi).^2 + (R2S5+ksw5)*ksw5 + ksw5./(R2S5+ksw5).*((k+sep5)*2*pi).^2));
                                                    cal_Lorentzian6_n_cal_ns = m_MT_n'*fm_cal;
        
                                                    cal_eff_cal_ns = R1W_cal_obs_n.*((k)*2*pi).^2./((tt.*42.6*2*pi).^2 + ((k)*2*pi).^2) + ...
                                                                     R2W_cal.*(tt.*42.6*2*pi).^2./((tt.*42.6*2*pi).^2 + ((k)*2*pi).^2);
        
                                                    sscal_n_ns = R1W_cal_obs_n ./ (cal_eff_cal_ns + cal_Lorentzian1_cal_ns./(1+fm_cal*fm_n) + cal_Lorentzian2_n_cal_ns./(1+fm_cal*fm_n) + ...
                                                                                   cal_Lorentzian3_n_cal_ns./(1+fm_cal*fm_n) + cal_Lorentzian4_n_cal_ns./(1+fm_cal*fm_n) + ...
                                                                                   cal_Lorentzian5_n_cal_ns./(1+fm_cal*fm_n) + cal_Lorentzian6_n_cal_ns) .* ...
                                                                  (((k)*2*pi).^2 ./ ((tt.*42.6*2*pi).^2 + ((k)*2*pi).^2));
        
                                                    SS_cal_n_ns(:,  ii_T1W, ii_T2W, ii_T2S3, ii_T2S4, ii_fs1, ii_fs2, ii_fs3, ii_fs4, ii_fm, ii_ksw3, ii_ksw4, ii_tt, ii_kk) = sscal_n_ns;
        
                                                    R1W_cal_matrix_n( ii_T1W, ii_T2W, ii_T2S3, ii_T2S4, ii_fs1, ii_fs2, ii_fs3, ii_fs4, ii_fm, ii_ksw3, ii_ksw4, ii_tt, ii_kk) = R1W_cal_obs_n;
                                                    fm_cal_matrix_n( ii_T1W, ii_T2W, ii_T2S3, ii_T2S4, ii_fs1, ii_fs2, ii_fs3, ii_fs4, ii_fm, ii_ksw3, ii_ksw4, ii_tt, ii_kk) = fm_cal*fm_n;
                                                    params_matrix(:,  ii_T1W, ii_T2W, ii_T2S3, ii_T2S4, ii_fs1, ii_fs2, ii_fs3, ii_fs4, ii_fm, ii_ksw3, ii_ksw4, ii_tt, ii_kk) = ...
                                                        [R1W_cal_obs_n; R2W_cal; R2S3_cal; R2S4_cal; fs1; fs2_cal; fs3_cal; fs4_cal; fm_cal; ksw3_cal; ksw4_cal; tt_shift; k(35)];
        
                                                    fs_pcr(1, ii_T1W, ii_T2W, ii_T2S3, ii_T2S4, ii_fs1, ii_fs2, ii_fs3, ii_fs4, ii_fm, ii_ksw3, ii_ksw4, ii_tt, ii_kk)  = fs3_cal;
                                                    ksw_pcr(1, ii_T1W, ii_T2W, ii_T2S3, ii_T2S4, ii_fs1, ii_fs2, ii_fs3, ii_fs4, ii_fm, ii_ksw3, ii_ksw4, ii_tt, ii_kk) = ksw3_cal;
        
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

%% ---------------------------------------------------------------------
% Post-processing (reshape to 2D matrices; identical to original intent)

% with shifts
matrix_input_all_n     = reshape(SS_cal_n,    [length(k')  num_T1W*num_T2W*num_T2S3*num_T2S4*num_fs1*num_fs2*num_fs3*num_fs4*num_fm*num_ksw3*num_ksw4*num_tt*num_B0]);

% without shifts
matrix_input_all_n_ns  = reshape(SS_cal_n_ns, [length(k')  num_T1W*num_T2W*num_T2S3*num_T2S4*num_fs1*num_fs2*num_fs3*num_fs4*num_fm*num_ksw3*num_ksw4*num_tt*num_B0]);

% requred parameters for fitting
R1W_cal_matrix_output_all_n = reshape(R1W_cal_matrix_n, [num_T1W*num_T2W*num_T2S3*num_T2S4*num_fs1*num_fs2*num_fs3*num_fs4*num_fm*num_ksw3*num_ksw4*num_tt*num_B0 1]);
fm_cal_matrix_output_all_n  = reshape(fm_cal_matrix_n,  [num_T1W*num_T2W*num_T2S3*num_T2S4*num_fs1*num_fs2*num_fs3*num_fs4*num_fm*num_ksw3*num_ksw4*num_tt*num_B0 1]);

% details of parameters
params_matrix = reshape(params_matrix, [13 num_T1W*num_T2W*num_T2S3*num_T2S4*num_fs1*num_fs2*num_fs3*num_fs4*num_fm*num_ksw3*num_ksw4*num_tt*num_B0]);

% target values
matrix_AREX_output_pcr(:,1) = reshape(fs_pcr,        [num_T1W*num_T2W*num_T2S3*num_T2S4*num_fs1*num_fs2*num_fs3*num_fs4*num_fm*num_ksw3*num_ksw4*num_tt*num_B0 1]);
matrix_AREX_output_pcr(:,2) = reshape(ksw_pcr,       [num_T1W*num_T2W*num_T2S3*num_T2S4*num_fs1*num_fs2*num_fs3*num_fs4*num_fm*num_ksw3*num_ksw4*num_tt*num_B0 1]);

% Plotting
figure(1);
plot(sample_amine);
hold on
plot(amine_normal);
ylim([0 0.5]);

figure(2);
plot(matrix_input_all_n_ns(:,1));