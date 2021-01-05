%% Working on some simulation for the QuiQI project 
% Main project lead by A. Lutti
%__________________________________________________________________________
% Copyright (C) 2020 Cyclotron Research Centre

% Written by C. Phillips, 2020.
% GIGA Institute, University of Liege, Belgium

% Some parameters
% ===============
Nsub = 20; % Number of subjects
sig_intens = [100 120]; % signal intensity for 1 or 2 groups
Qind_range = [2 10]; % quality index range
Qind_weight = 1; % weight of the quality index 
noise_std = [5 10]; % background noise std per group
age_range = [20 70];
age_weight = .2;
gr_ratio = [1 3]/4; % #sub group ratio for 2s t-test

reml_sc = 0; % positivity constraints on covariance parameters

% random ages and Qindex in interval
age_sub = age_range(1) + rand(Nsub,1)*diff(age_range);
Qind_sub = Qind_range(1) + rand(Nsub,1)*diff(Qind_range);

% 1-sample t-test model
% =====================
% design matrix
X = [ones(Nsub,1) age_sub-mean(age_sub)] ; 

% create noisy signal
Y = sig_intens(1) ... % baseline noise free signal
    + age_weight * (age_sub-min(age_sub)) ... % age effect
    + randn(Nsub,1)*noise_std(1) ... % adding background noise
    + Qind_weight*(Qind_sub.*randn(Nsub,1)); % adding QInd noise

% Simple OLS solution
% -------------------
beta = X\Y;
res = Y - X*beta;
fprintf('\nSimple OLS solution:')
fprintf('\n\tBeta''s (mean + age) : %.3f + %.3f x age', beta)
fprintf('\n\tMean and std of res  : %.3f +/- %.3f', mean(res), std(res))
fprintf('\n')

% Trying out ReML estimate with different Q's
% -------------------------------------------

% Create data cov matrix + adjusted Qind values
YY = Y*Y'; % assuming a single voxel
Qind_sub_m = Qind_sub-min(Qind_sub); % Qind min value to 0
Qind_sub_mm = (Qind_sub-min(Qind_sub))/max(Qind_sub); % Q ind in [0 1]

% Q1, Using just identity matrix, i.e. OLS
Q_lab{1} = 'Using just identity matrix, i.e. OLS';
Q1 = {eye(Nsub)};

% Q2, Using identity matrix + Qind_sub
Q_lab{2} = 'Using identity matrix + Qind_sub';
Q2 = {eye(Nsub), diag(Qind_sub)}; 

% Q3, Using identity matrix + Qind_sub with min=0
Q_lab{3} = 'Using identity matrix + Qind_sub with min=0';
Q3 = {eye(Nsub), diag(Qind_sub_m)}; 

% Q4, Using identity matrix + Qind_sub in [0 1]
Q_lab{4} = 'Using identity matrix + Qind_sub in [0 1]';
Q4 = {eye(Nsub) ,diag(Qind_sub_mm)}; 

% Q5, Using id matrix + 1st, 2nd and 3rd order for Qind in [0 1]
Q_lab{5} = 'Using id matrix + 1st, 2nd and 3rd order for Qind in [0 1]';
Q5 = {eye(Nsub), diag(Qind_sub_mm), diag(Qind_sub_mm.^2), diag(Qind_sub_mm.^3)}; 

% Q6, Using id matrix + 1st, 2nd and 3rd order for original Qind
Q_lab{6} = 'Using id matrix + 1st, 2nd and 3rd order for original Qind';
Q6 = {eye(Nsub), diag(Qind_sub), diag(Qind_sub.^2), diag(Qind_sub.^3)}; % Using Qind_sub

% Collecting all Q's
Qall = {Q1, Q2, Q3, Q4, Q5, Q6};
num_Q = numel(Qall);
beta_s = cell(1,num_Q);
res_s = cell(1,num_Q);
h = cell(1,num_Q);

% Try them all & collect results
for iQ = 1:num_Q
    fprintf('\nQ model #%d : %s\n', iQ, Q_lab{iQ});
    switch reml_sc
        case 0, % Classic ReML
            [V,h{iQ,1},Ph,F,Fa,Fc] = spm_reml(YY,X,Qall{iQ});
        case 1, % ReML with positivity constraint
            [V,h{iQ,1},Ph,F,Fa,Fc] = spm_reml_sc(YY,X,Qall{iQ});
    end
%     beta_ss = inv(X'*inv(V)*X)*X'*inv(V)*Y %#ok<*MINV>
    W = sqrt(inv(V));
    Ys = W*Y;
    Xs = W*X;
    beta_s{iQ} = Xs\Ys;
    res_s{iQ} = Ys - Xs*beta_s{iQ};
end

% Display results
fprintf('\nSimple OLS solution:')
fprintf('\n\tBeta''s (mean + age) : %.3f + %.3f x age', beta)
fprintf('\n\tMean and std of res  : %.3f +/- %.3f', mean(res), std(res))
fprintf('\n')
fprintf('\nWLS with ReML solution:')
for iQ = 1:num_Q
    fprintf('\nQ model #%d : %s', iQ, Q_lab{iQ});
    fprintf('\n\tHyper-parameter(s) : ')
    for iQ_h = 1:numel(Qall{iQ})
        fprintf('%.3f ',full(h{iQ}(iQ_h)))
    end
    fprintf('\n\tBeta''s (mean + age) : %.3f + %.3f x age', beta_s{iQ})
    fprintf('\n\tMean and std of W-residual : %.3f +/- %.3f', ...
        mean(res_s{iQ}), std(res_s{iQ}))
    fprintf('\n')
end
fprintf('\n')

% 2-sample t-test model
% =====================
% design matrix
Nsub_gr12(1) = round(Nsub*gr_ratio(1));
Nsub_gr12(2) = Nsub-Nsub_gr12(1);

X = [ones(1,Nsub_gr12(1)) zeros(1,Nsub_gr12(2)) ;  ...
     zeros(1,Nsub_gr12(1)) ones(1,Nsub_gr12(2))]' ; 

% create noisy signal
Y = X*sig_intens' ... % baseline noise free signal
    + [randn(Nsub_gr12(1),1)*noise_std(1) ; randn(Nsub_gr12(2),1)*noise_std(2)]... % adding background noise
    + Qind_weight*(Qind_sub.*randn(Nsub,1)); % adding QInd noise

% Simple OLS solution
% -------------------
beta = X\Y;
res = Y - X*beta;
fprintf('\nSimple OLS solution:')
fprintf('\n\tBeta''s (gr1 + gr2) : %.3f , %.3f', beta)
fprintf('\n\tMean and std of res  : %.3f +/- %.3f', mean(res), std(res))
fprintf('\n')

% Trying out ReML estimate with different Q's
% -------------------------------------------

% Create data cov matrix + adjusted Qind values
YY = Y*Y'; % assuming a single voxel
Qind_sub_m = Qind_sub-min(Qind_sub); % Qind min value to 0
Qind_sub_mm = (Qind_sub-min(Qind_sub))/max(Qind_sub); % Q ind in [0 1]
Qgroup = {diag(X(:,1)), diag(X(:,2))};

% Q1, Using just non-indentical variance
Q_lab{1} = 'Using just non-indentical variance';
Q1 = Qgroup;

% Q2, Using non-ident var + Qind_sub
Q_lab{2} = 'Using non-ident var + Qind_sub';
Q2 = {Qgroup{:}, diag(Qind_sub)}; 

% Q3, Using non-ident var + Qind_sub with min=0
Q_lab{3} = 'Using non-ident var + Qind_sub with min=0';
Q3 = {Qgroup{:}, diag(Qind_sub_m)}; 

% Q4, Using non-ident var + Qind_sub in [0 1]
Q_lab{4} = 'Using non-ident var + Qind_sub in [0 1]';
Q4 = {Qgroup{:} ,diag(Qind_sub_mm)}; 

% Q5, Using non-ident var + 1st, 2nd and 3rd order for Qind in [0 1]
Q_lab{5} = 'Using non-ident var + 1st, 2nd and 3rd order for Qind in [0 1]';
Q5 = {Qgroup{:}, diag(Qind_sub_mm), diag(Qind_sub_mm.^2), diag(Qind_sub_mm.^3)}; 

% Q6, Using non-ident var + 1st, 2nd and 3rd order for original Qind
Q_lab{6} = 'Using non-ident var + 1st, 2nd and 3rd order for original Qind';
Q6 = {Qgroup{:}, diag(Qind_sub), diag(Qind_sub.^2), diag(Qind_sub.^3)}; % Using Qind_sub

% Collecting all Q's
Qall = {Q1, Q2, Q3, Q4, Q5, Q6};
num_Q = numel(Qall);
beta_s = cell(1,num_Q);
res_s = cell(1,num_Q);
h = cell(1,num_Q);

% Try them all & collect results
for iQ = 1:num_Q
    fprintf('\nQ model #%d : %s\n', iQ, Q_lab{iQ});
    switch reml_sc
        case 0, % Classic ReML
            [V,h{iQ,1},Ph,F,Fa,Fc] = spm_reml(YY,X,Qall{iQ});
        case 1, % ReML with positivity constraint
            [V,h{iQ,1},Ph,F,Fa,Fc] = spm_reml_sc(YY,X,Qall{iQ});
    end
%     beta_ss = inv(X'*inv(V)*X)*X'*inv(V)*Y %#ok<*MINV>
    W = sqrt(inv(V));
    Ys = W*Y;
    Xs = W*X;
    beta_s{iQ} = Xs\Ys;
    res_s{iQ} = Ys - Xs*beta_s{iQ};
end

% Display results
fprintf('\nSimple OLS solution:')
fprintf('\n\tBeta''s (gr1 + gr2) : %.3f , %.3f', beta)
fprintf('\n\tMean and std of res  : %.3f +/- %.3f', mean(res), std(res))
fprintf('\n')
fprintf('\nWLS with ReML solution:')
for iQ = 1:num_Q
    fprintf('\nQ model #%d : %s', iQ, Q_lab{iQ});
    fprintf('\n\tHyper-parameter(s) : ')
    for iQ_h = 1:numel(Qall{iQ})
        fprintf('%.3f ',full(h{iQ}(iQ_h)))
    end
    fprintf('\n\tBeta''s (gr1 + gr2) : %.3f , %.3f', beta_s{iQ})
    fprintf('\n\tMean and std of W-residual : %.3f +/- %.3f', ...
        mean(res_s{iQ}), std(res_s{iQ}))
    fprintf('\n')
end
fprintf('\n')
