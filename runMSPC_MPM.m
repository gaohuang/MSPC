clc;clear; format compact;

addpath('data')
addpath('functions')

% load ionosphere
load breast; 
% load australian
% load diabetes
% load letter;
% load sate_scale
% load spam_scale

% train MSPC-MPM
rho_set=[10.^(4:-1:-4),1e-10];
for i=1:length(rho_set)
    paras.verbose=1;
    paras.rho=rho_set(i);
    paras.y=y;
    model=mspc_mpm(X,paras);
    Acc(i)=accuracy(model.pred,y);
    Err_MPC(i)=100-Acc(i);
    Err_KM(i)=100-model.acc_km;
end

disp(['Best error rate of k-means: ',num2str(min(Err_KM)),'%'])
disp(['Best error rate of MSPC-MPM: ',num2str(min(Err_MPC)),'%'])
