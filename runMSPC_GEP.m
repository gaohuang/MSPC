clc; clear; close all;
format compact;

addpath('data')
addpath('functions')
% mex functions/find_opt_label.c

% load ionosphere
load breast; 
% load australian
% load diabetes
% load letter;  
% load sate_scale
% load spam_scale

% train MSPC-GEP
rho_set=[10.^(4:-1:-4),1e-10];
for i=1:length(rho_set)
    paras.verbose=1;
    paras.rho=rho_set(i)*var(X);
    paras.y=y;    
    model=mspc_gep(X,paras);    
    Acc(i)=accuracy(model.pred,y);
end

disp(['Best error rate of MSPC-GEP: ',num2str(100-max(Acc)),'%'])


