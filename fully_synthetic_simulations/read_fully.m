%% Partially Synthetic Data Simulation (Cleaned for GitHub/IEEE)
% Author: Malvika Viswanathan
% Corresponding author: Dr. Zhongliang Zu
% Email: zhongliang.zu@vumc.org
%
% NOTE: Logic preserved exactly; only formatting, comments, and missing
% semicolons were added for readability.

clear; clc;

%% ---------------------------------------------------------------------
% Defining constants and ranges
% Please ensure power of Lorentzian fit files match the power used for
% simulations.

tt = 0.5;                               % B1 (uT)
tt_4p7T = [0.8, 0.9, 1, 1.1, 1.2] * tt; % B1 with shifts

% Sequence parameters
pulseduration = 5;
gauss         = 100;


% Flip angle of saturation pulse
satangle = tt_4p7T * 42.6 * 360 * pulseduration;

offppm = 200; % frequency of 4.7T (MHz)

% Frequency offsets of different pools
sep1_4p7T = 3.6 * offppm;
sep2_4p7T = 3.0 * offppm;
sep3_4p7T = 2.6 * offppm;  % (added missing semicolon)
sep4_4p7T = 2.0 * offppm;
sep5_4p7T = -3.3 * offppm;

% Required relaxations
R1S  = 1/1.5;      % R1 of all pools (1/s)
R2S1 = 1/0.002;    % R2 APT (1/s)
R2S2 = 1/0.015;    % R2 amine CEST (1/s)
R2S5 = 1/0.0005;   % R2 NOE (1/s)
R1M  = 1/1.5;      % R1 MT (1/s)
R2M  = 1/0.00005;  % R2 MT (1/s)

% Frequency offset grid
maxv  = 1000;
step  = 25;
offset = -maxv:step:maxv;
k      = [-2000, -1750, -1500, -1250, offset, 1250, 1500, 1750, 2000];
k_4p7T = [k-100; k-50; k; k+50; k+100]'; % B0 shifts

% Constants required for simulation
fs5 = 0.015;
ksw1 = 100;
ksw5 = 10;
kmw  = 25;

%% ---------------------------------------------------------------------
% Define ranges for simulation

num_T1W  = 3;
num_T2W  = 3;
num_fs1  = 3;
num_fs2  = 3;
num_fs3  = 4;
num_fs4  = 3;
num_fm   = 3;
num_ksw2 = 3;
num_ksw3 = 6;
num_ksw4 = 3;
num_T2S3 = 3;
num_T2S4 = 3;
num_tt   = 3;
num_B0   = 3;

T1W_matrix  = 1.1:0.3:1.7;
T2W_matrix  = [20, 40, 60] * 0.001;
T2S3_matrix = [0.008, 0.010, 0.012];
T2S4_matrix = [0.008, 0.010, 0.012];

fs1_matrix = [0.0008, 0.0010, 0.0012];
fs2_matrix = [0.0018, 0.0030, 0.0042];
fs3_matrix = [0.0006, 0.0008, 0.0010, 0.0012, 0.0014];
fs4_matrix = [0.0016, 0.0020, 0.0024];
fm_matrix  = [0.075, 0.150, 0.225];

ksw2_matrix = [3000, 5000, 7000];
ksw3_matrix = [60, 80, 100, 120, 140, 160];
ksw4_matrix = [300, 500, 700];

i = 1;  % loop counter

%% ---------------------------------------------------------------------
% Main nested loops (PRESERVED exactly)
for ii_T1W = 1:num_T1W
    ii_T1W
    R1W = 1 ./ T1W_matrix(ii_T1W);
    for ii_T2W = 1:num_T2W
        ii_T2W
        R2W = 1 ./ T2W_matrix(ii_T2W);
        for ii_T2S3 = 1:num_T2S3
            R2S3 = 1 ./ T2S3_matrix(ii_T2S3);
            for ii_T2S4 = 1:num_T2S4
                R2S4 = 1 ./ T2S4_matrix(ii_T2S4);
                for ii_fs1 = 1:num_fs1
                    fs1 = fs1_matrix(ii_fs1);
                    for ii_fs2 = 1:num_fs2
                        fs2 = fs2_matrix(ii_fs2);
                        for ii_fs3 = 1:num_fs3
                            fs3 = fs3_matrix(ii_fs3);
                            for ii_fs4 = 1:num_fs4
                                fs4 = fs4_matrix(ii_fs4);
                                for ii_fm = 1:num_fm
                                    fm = fm_matrix(ii_fm);
                                    for ii_ksw2 = 1:num_ksw2
                                        ksw2 = ksw2_matrix(ii_ksw2);
                                        for ii_ksw3 = 1:num_ksw3
                                            ksw3 = ksw3_matrix(ii_ksw3);
                                            for ii_ksw4 = 1:num_ksw4
                                                ksw4 = ksw4_matrix(ii_ksw4);
                                                for ii_tt = 1:num_tt
                                                    tt_shift = satangle(ii_tt);
                                                    for ii_kk = 1:num_B0
                                                        B0_shift = k_4p7T(:, ii_kk);

                                                        R1W_cal_obs = (R1W + (fm * R1M)) ./ (1 + fm);
                                                        R1W_cal_matrix(ii_T1W, ii_T2W, ii_T2S3, ii_T2S4, ii_fs1, ii_fs2, ii_fs3, ii_fs4, ii_fm, ii_ksw2, ii_ksw3, ii_ksw4, ii_tt, ii_kk) = R1W_cal_obs;

                                                        a25mspulse = runsteadysimgauss(ksw1, ksw2, ksw3, ksw4, ksw5, kmw, fs1, fs2, fs3, fs4, fs5, 1, fm, R1S, R2S1, R2S2, R2S3, R2S4, R2S5, R1W_cal_obs, R2W, R1M, R2M, sep1_4p7T*2*pi, sep2_4p7T*2*pi, sep3_4p7T*2*pi, sep4_4p7T*2*pi, sep5_4p7T*2*pi, pulseduration, gauss, tt_shift, 1, 2, 1, .00, 1, 1, B0_shift*2*pi, 1);
                                                        a25mspulse_ref_pcr = runsteadysimgauss(ksw1, ksw2, ksw3, ksw4, ksw5, kmw, fs1, fs2, 0, fs4, fs5, 1, fm, R1S, R2S1, R2S2, R2S3, R2S4, R2S5, R1W_cal_obs, R2W, R1M, R2M, sep1_4p7T*2*pi, sep2_4p7T*2*pi, sep3_4p7T*2*pi, sep4_4p7T*2*pi, sep5_4p7T*2*pi, pulseduration, gauss, tt_shift, 1, 2, 1, .00, 1, 1, B0_shift*2*pi, 1);
                                                        a25mspulse_ref_cr = runsteadysimgauss(ksw1, ksw2, ksw3, ksw4, ksw5, kmw, fs1, fs2, fs3, 0, fs5, 1, fm, R1S, R2S1, R2S2, R2S3, R2S4, R2S5, R1W_cal_obs, R2W, R1M, R2M, sep1_4p7T*2*pi, sep2_4p7T*2*pi, sep3_4p7T*2*pi, sep4_4p7T*2*pi, sep5_4p7T*2*pi, pulseduration, gauss, tt_shift, 1, 2, 1, .00, 1, 1, B0_shift*2*pi, 1);

                                                        a25mspulse_ns = runsteadysimgauss(ksw1, ksw2, ksw3, ksw4, ksw5, kmw, fs1, fs2, fs3, fs4, fs5, 1, fm, R1S, R2S1, R2S2, R2S3, R2S4, R2S5, R1W_cal_obs, R2W, R1M, R2M, sep1_4p7T*2*pi, sep2_4p7T*2*pi, sep3_4p7T*2*pi, sep4_4p7T*2*pi, sep5_4p7T*2*pi, pulseduration, gauss, satangle(2), 1, 2, 1, .00, 1, 1, k_4p7T(:, 2)*2*pi, 1);
                                                        a25mspulse_ref_pcr_ns = runsteadysimgauss(ksw1, ksw2, ksw3, ksw4, ksw5, kmw, fs1, fs2, 0, fs4, fs5, 1, fm, R1S, R2S1, R2S2, R2S3, R2S4, R2S5, R1W_cal_obs, R2W, R1M, R2M, sep1_4p7T*2*pi, sep2_4p7T*2*pi, sep3_4p7T*2*pi, sep4_4p7T*2*pi, sep5_4p7T*2*pi, pulseduration, gauss, satangle(2), 1, 2, 1, .00, 1, 1, k_4p7T(:, 2)*2*pi, 1);
                                                        a25mspulse_ref_cr_ns = runsteadysimgauss(ksw1, ksw2, ksw3, ksw4, ksw5, kmw, fs1, fs2, fs3, 0, fs5, 1, fm, R1S, R2S1, R2S2, R2S3, R2S4, R2S5, R1W_cal_obs, R2W, R1M, R2M, sep1_4p7T*2*pi, sep2_4p7T*2*pi, sep3_4p7T*2*pi, sep4_4p7T*2*pi, sep5_4p7T*2*pi, pulseduration, gauss, satangle(2), 1, 2, 1, .00, 1, 1, k_4p7T(:, 2)*2*pi, 1);
                                                        a25mspulse_ns_amine = runsteadysimgauss(ksw1, ksw2, ksw3, ksw4, ksw5, kmw, 0, fs2, 0, 0, 0, 1, fm, R1S, R2S1, R2S2, R2S3, R2S4, R2S5, R1W_cal_obs, R2W, R1M, R2M, sep1_4p7T*2*pi, sep2_4p7T*2*pi, sep3_4p7T*2*pi, sep4_4p7T*2*pi, sep5_4p7T*2*pi, pulseduration, gauss, satangle(2), 1, 2, 1, .00, 1, 1, k_4p7T(:, 2)*2*pi, 1);

                                                        disp("step1")
                                                        aa_4p7T(:, ii_T1W, ii_T2W, ii_T2S3, ii_T2S4, ii_fs1, ii_fs2, ii_fs3, ii_fs4, ii_fm, ii_ksw2, ii_ksw3, ii_ksw4, ii_tt, ii_kk) = a25mspulse(:, 6);
                                                        aa_4p7T_ns(:, ii_T1W, ii_T2W, ii_T2S3, ii_T2S4, ii_fs1, ii_fs2, ii_fs3, ii_fs4, ii_fm, ii_ksw2, ii_ksw3, ii_ksw4, ii_tt, ii_kk) = a25mspulse_ns(:, 6);
                                                        aa_4p7T_ref(:, ii_T1W, ii_T2W, ii_T2S3, ii_T2S4, ii_fs1, ii_fs2, ii_fs3, ii_fs4, ii_fm, ii_ksw2, ii_ksw3, ii_ksw4, ii_tt, ii_kk) = a25mspulse_ref_pcr(:, 6);

                                                        disp("Zspec")

                                                        Slab     = a25mspulse(:, 6);
                                                        Sref_cr  = a25mspulse_ref_cr(:, 6);
                                                        Sref_Pcr = a25mspulse_ref_pcr(:, 6);

                                                        Slab_ns    = a25mspulse_ns(:, 6);
                                                        Sref_ns_cr = a25mspulse_ref_cr_ns(:, 6);
                                                        Sref_ns_pcr= a25mspulse_ref_pcr_ns(:, 6);

                                                        S0 = 1;

                                                        fm_cal_matrix(ii_T1W, ii_T2W, ii_T2S3, ii_T2S4, ii_fs1, ii_fs2, ii_fs3, ii_fs4, ii_fm, ii_ksw2, ii_ksw3, ii_ksw4, ii_tt, ii_kk) = fm;
                                                        params_matrix(:, ii_T1W, ii_T2W, ii_T2S3, ii_T2S4, ii_fs1, ii_fs2, ii_fs3, ii_fs4, ii_fm, ii_ksw2, ii_ksw3, ii_ksw4, ii_tt, ii_kk) = ...
                                                            [R1W_cal_obs; R2W; R2S3; R2S4; fs1; fs2; fs3; fs4; fs5; fm; ksw3; ksw4; tt_4p7T(ii_tt); k_4p7T(35, ii_kk)];

                                                        % PCr
                                                        fs_pcr(ii_T1W, ii_T2W, ii_T2S3, ii_T2S4, ii_fs1, ii_fs2, ii_fs3, ii_fs4, ii_fm, ii_ksw2, ii_ksw3, ii_ksw4, ii_tt, ii_kk) = fs3;
                                                        ksw_pcr(ii_T1W, ii_T2W, ii_T2S3, ii_T2S4, ii_fs1, ii_fs2, ii_fs3, ii_fs4, ii_fm, ii_ksw2, ii_ksw3, ii_ksw4, ii_tt, ii_kk) = ksw3;

                                                        X = sprintf('-------------------------------------%d', i);
                                                        disp(X)
                                                        i = i + 1;
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
end

%% ---------------------------------------------------------------------
% Reshape outputs

matrix_output1(:, 1) = reshape(fs_pcr,  [num_T1W*num_T2W*num_T2S3*num_T2S4*num_fs1*num_fs2*num_fs3*num_fs4*num_fm*num_ksw2*num_ksw3*num_ksw4*num_tt*num_B0, 1]);
matrix_output1(:, 2) = reshape(ksw_pcr, [num_T1W*num_T2W*num_T2S3*num_T2S4*num_fs1*num_fs2*num_fs3*num_fs4*num_fm*num_ksw2*num_ksw3*num_ksw4*num_tt*num_B0, 1]);

matrix_input_all = reshape(aa_4p7T(:,:), [89, num_T1W*num_T2W*num_T2S3*num_T2S4*num_fs1*num_fs2*num_fs3*num_fs4*num_fm*num_ksw2*num_ksw3*num_ksw4*num_tt*num_B0]);

R1W_cal_matrix_output_all = reshape(R1W_cal_matrix, [num_T1W*num_T2W*num_T2S3*num_T2S4*num_fs1*num_fs2*num_fs3*num_fs4*num_fm*num_ksw2*num_ksw3*num_ksw4*num_tt*num_B0, 1]);
fm_cal_matrix_output_all  = reshape(fm_cal_matrix,  [num_T1W*num_T2W*num_T2S3*num_T2S4*num_fs1*num_fs2*num_fs3*num_fs4*num_fm*num_ksw2*num_ksw3*num_ksw4*num_tt*num_B0, 1]);
params_matrix = reshape(params_matrix, [14, num_T1W*num_T2W*num_T2S3*num_T2S4*num_fs1*num_fs2*num_fs3*num_fs4*num_fm*num_ksw2*num_ksw3*num_ksw4*num_tt*num_B0]);
