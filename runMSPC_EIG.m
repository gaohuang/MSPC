clc; clear; format compact;

addpath('data')
addpath('functions')

% load ionosphere
load breast; 
% load australian
% load diabetes
% load letter; 
% load sate_scale
% load spam_scale


% train MSPC-EIG
rho_set=[10.^(4:-1:-4),1e-10];
for i=1:length(rho_set)
    paras.rho=rho_set(i);
    paras.y=y;
    model=mspc_eig(X,paras);
    Acc(i)=model.acc_w
end

disp(['Best error rate of MSPC-EIG: ',num2str(100-max(Acc)),'%'])

