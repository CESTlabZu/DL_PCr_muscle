%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Multiple Pool Lorenztian Fit (mfit)
% Please uncomment required lines of text
%
% Authors: Malvika Viswanathan, Zhongliang Zu
%
% Correspondance: zhongliang.zu@vumc.org 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
% load CEST Z spectrum, R1W and fm values from simulations or measured
% data.
% load('...')

% required initial parameters
maxv  = 1000;
step  = 25;
offset = -maxv:step:maxv;
k = [-2000, -1750, -1500, -1250, offset, 1250, 1500, 1750, 2000];
k_4p7T=k';

 for i=1:length(matrix_input_all)

    sig=(1-matrix_input_all(:,i)); 
    R1W_AREX=R1W_cal_matrix_output_all(i); % loaded R1W values
    fm_AREX=fm_cal_matrix_output_all(i); % loaded PSR values
   
    x =k_4p7T(index);
    % Pools: Water, APT, PCr, amines/guanidine, NOE, MT
    beta0= [0.9,   0,    280,   0.025, -700, 100,      0.01, -520,  200,    0.01,-400, 200,  0.02, 700,  600,       0.1,    0, 5000]; % initial test
    lb=[0.02,-200,     20,       0, -800,  80,     0,   -600,  100,        0, -500,  100,   0, 500,  200,         0, -800, 2000]; % lower bound
    ub=[ 1,  200,   2000,     0.2, -600, 600,     0.2, -500, 500,      0.2, -300, 400,      1, 900, 1000,         1,  800, 20000]; % upper bound
                    
    Delta=[1]; 
    options=optimset('lsqcurvefit') ; 
    options = optimset('Display','off','TolFun',1e-8,'TolX',1e-8,'MaxFunEvals',5e4*length(x),'MaxIter',2e5) ;
    
    [beta,resnorm,residual,exitflag,output,lambda,jacobian] = ...
        lsqcurvefit(@matsolv, beta0, x, sig, lb, ub, options, Delta) ;
    
    
    
    % amide
    beta_amide=beta;
    sig_simur_amide=matsolv(beta_amide,x,Delta);
    mor_AREX_amide2(:,i) = spectrum(beta_amide(4), beta_amide(6), -abs(beta_amide(5)));

    beta_amide(4)=0;
    sig_simur_ref_amide=matsolv(beta_amide,x,Delta);
    % save required values!
    mor_MTR_amide(:,i)=(sig_simur_amide-sig_simur_ref_amide);
    mor_AREX_amide(:,i)=(1./(1-sig_simur_amide)-1./(1-sig_simur_ref_amide))*R1W_AREX*(1+fm_AREX);
    
    % PCr
    beta_PCr=beta;
    sig_simur_PCr=matsolv(beta_PCr,x,Delta);
    beta_PCr(7)=0;
    sig_simur_ref_PCr=matsolv(beta_PCr,x,Delta);
    
    mor_PCr=(1./(1-sig_simur_PCr)-1./(1-sig_simur_ref_PCr))*R1W_AREX*(1+fm_AREX);
    mor_dir_PCr=(sig_simur_PCr-sig_simur_ref_PCr);

    % amines/ guanidine
    beta_amine=beta;
    sig_simur_amine=matsolv(beta_amine,x,Delta);
    beta_amine(10)=0;
    sig_simur_ref_amine=matsolv(beta_amine,x,Delta);

    mor_AREX_amine(:,i)=(1./(1-sig_simur_amine)-1./(1-sig_simur_ref_amine))*R1W_AREX*(1+fm_AREX);
    mor_MTR_amine(:,i)=(sig_simur_amine-sig_simur_ref_amine);

    % NOE
    beta_NOE=beta;
    sig_simur_NOE=matsolv(beta_NOE,x,Delta);
    beta_NOE(13)=0;
    sig_simur_ref_NOE=matsolv(beta_NOE,x,Delta);
    
    mor_NOE(:,i)=(1./(1-sig_simur_NOE)-1./(1-sig_simur_ref_NOE))*R1W_AREX*(1+fm_AREX);
    mor_dir_NOE(:,i)=(sig_simur_NOE-sig_simur_ref_NOE);

    % MT
    beta_MT=beta;
    sig_simur_MT=matsolv(beta_MT,x,Delta);
    beta_MT(16)=0;
    sig_simur_ref_MT=matsolv(beta_MT,x,Delta);
    
    mor_MT=(sig_simur_MT-sig_simur_ref_MT);
    mor_AREX_MT(:,i) = (mor_MT./(1-mor_MT))*R1W_AREX;
    mor_dir_MT(:,i)=(sig_simur_MT-sig_simur_ref_MT);

    
   
    sprintf("----------------------- %d",i)
        
 end

 