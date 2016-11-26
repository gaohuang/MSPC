function model=mspc_mpm(X,paras)

if ~isfield(paras,'verbose')
    paras.verbose=1;
end
[N,d]=size(X);

%% initiate with kmeans
acc_km=0;
for i=1:1
    [label0, center] = litekmeans(X, 2, 'MaxIter',200);
    acc_tmp(i)=accuracy(paras.y,label0);
    if acc_tmp(i)>acc_km
        acc_km=acc_tmp(i);
        label=label0;
    end
end

if paras.verbose==1
    disp(['Initial clustering error using k-means: ',num2str(100-max(acc_km)),'%']);
end

%%
if ~isfield(paras,'MaxIter')
    paras.MaxIter=50;
end

Miu1=mean(X(label==1,:),1);
Miu2=mean(X(label==2,:),1);
Cov1=cov(X(label==1,:),1);
Cov2=cov(X(label==2,:),1);
label0=zeros(size(label));

acc=acc_km;
rho=paras.rho*var(X,1);
for iter=1:paras.MaxIter
    
    if iter>1
        [Miu1,Miu2,Cov1,Cov2]=update_mean_cov(Miu1,Miu2,Cov1,Cov2,X,label0,label);
    end
    
    label0=label;
    [alfa,a,b]=mpm_linear(Miu1',Miu2',Cov1,Cov2,0,rho,rho,0,1e-6,1e-6,50);
    out=X*a-b;
    label=-(sign(out)-3)/2;
    
    if norm(label-label0)<1e-6
        if paras.verbose==1
            disp(['Algorithm converges at iteration ', num2str(iter-1),'!'])
        end
        break;
    end
    
    acc(iter)=accuracy(paras.y,label);
    if paras.verbose==1
        disp(['Accuracy at iteration ',num2str(iter),': ',num2str(acc(iter))])
    end
end

% record p and accuracy
Acc=acc;

kappa=abs(Miu1*a-Miu2*a)/(sqrt(a'*Cov1*a)+sqrt(a'*Cov2*a));

p=kappa^2/(1+kappa^2);
kappa0=abs(Miu1*a-Miu2*a)/(sqrt(a'*(Cov1+diag(rho))*a)+sqrt(a'*(Cov2+diag(rho))*a));
p0=kappa0^2/(1+kappa0^2);


model.a=a;
model.b=b;
model.alfa=alfa;
model.pred=label;
model.iteration=iter;
model.p=100*p;
model.p0=100*p0;
model.kappa=kappa;
model.kappa0=kappa0;
model.acc=Acc;
model.acc_km=acc_km;
model.center=[Miu1;Miu2];


function [Miu1,Miu2,Cov1,Cov2]=update_mean_cov(Miu1,Miu2,Cov1,Cov2,X,label0,label)

[N,d]=size(X);

idx1to2=(label0==1)&(label==2);
idx2to1=(label0==2)&(label==1);

N1to2=sum(idx1to2);
N2to1=sum(idx2to1);
N1=sum(label0==1);
N2=sum(label0==2);
N1r=N1-N1to2;
N2r=N2-N2to1;

if N1r==0|N2r==0
    
else
    if N1to2==0
        Miu1r=Miu1;
        Cov1r=Cov1;
        Cov1to2=zeros(d);
        Miu1to2=zeros(1,d);
    else
        Miu1to2=mean(X(idx1to2,:),1);
        Miu1r=(N1/N1r)*Miu1-(N1to2/N1r)*Miu1to2;
        if N1to2==1
            Cov1to2=zeros(d);
        else
            Cov1to2=cov(X(idx1to2,:),1);
        end
        Cov1r=(N1/N1r)*(Cov1-(N1to2/N1)*Cov1to2-(N1to2*N1r/N1^2)*(Miu1r-Miu1to2)'*(Miu1r-Miu1to2));
    end
    
    if N2to1==0
        Miu2r=Miu2;
        Cov2r=Cov2;
        Cov2to1=zeros(d);
        Miu2to1=zeros(1,d);
    else
        Miu2to1=mean(X(idx2to1,:),1);
        Miu2r=(N2/N2r)*Miu2-(N2to1/N2r)*Miu2to1;
        if N2to1==1
            Cov2to1=zeros(d);
        else
            Cov2to1=cov(X(idx2to1,:),1);
        end
        Cov2r=(N2/N2r)*(Cov2-(N2to1/N2)*Cov2to1-(N2to1*N2r/N2^2)*(Miu2r-Miu2to1)'*(Miu2r-Miu2to1));
    end
    
    Miu1=(N1r/(N1r+N2to1))*Miu1r+(N2to1/(N1r+N2to1))*Miu2to1;
    Miu2=(N2r/(N2r+N1to2))*Miu2r+(N1to2/(N2r+N1to2))*Miu1to2;
    Cov1=(N1r/(N1r+N2to1))*Cov1r+(N2to1/(N1r+N2to1))*Cov2to1+(N1r*N2to1/(N1r+N2to1)^2)*(Miu2to1-Miu1r)'*(Miu2to1-Miu1r);
    Cov2=(N2r/(N2r+N1to2))*Cov2r+(N1to2/(N2r+N1to2))*Cov1to2+(N2r*N1to2/(N2r+N1to2)^2)*(Miu1to2-Miu2r)'*(Miu1to2-Miu2r);   
end
if isnan(Miu1)|isnan(Miu2)
    keyboard
end








