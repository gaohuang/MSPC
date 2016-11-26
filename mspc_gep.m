function model=mspc_gep(X,paras)

if ~isfield(paras,'verbose')
    paras.verbose=1;
end
[N,d]=size(X);

%
%% initiate with kmeans
[label, center] = litekmeans(X, 2, 'MaxIter',200);
acc=accuracy(paras.y,label);
if paras.verbose==1
    disp(['Initial accuracy using kmeans: ',num2str(acc),'%'])
end

covZ=(N-1)/N*cov(X)+diag(paras.rho)+1e-6*eye(d);
%%
for iter=1:50
    Nx=sum(label==1);
    Ny=N-Nx;
    if iter==1
        mX=mean(X(label==1,:));
        mY=mean(X(label==2,:));
    else
        [mX,mY]=update_mean(mX,mY,X,label0,label);
    end
    
    
    A=(mX-mY)'*(mX-mY);
    opts.tol = 1e-9;
    opts.issym=1;
    opts.disp = 0;
    [w,v] = eigs(A,covZ,1,'lm',opts);
    
    out=X*w;
    flag=mean(sign(label-3/2)==sign(out))>0.5;
    out=out*(flag*2-1); % Avoid switching majority samples
    
    [val,idx]=sort(out,'descend');
    Num=w'*covZ*w;
    Sum=sum(val);
    obj=find_opt_label(val,Sum,Num);
    [Cost(iter),split]=min(obj);
        
    label0=label;
    label(idx(1:split))=2;
    label(idx(split+1:end))=1;

    
    if label==label0
        if paras.verbose==1
            disp(['Algorithm converges at iteration ', num2str(iter),'!'])
        end
        break;
    end
    acc(iter)=accuracy(paras.y,label);
end

% record p and accuracy
Acc=acc;

covX=(Nx-1)/Nx*cov(X(label==1,:));
covY=(Ny-1)/Ny*cov(X(label==2,:));
kappa=abs(mX*w-mY*w)/(sqrt(w'*covX*w)+(w'*covY*w));
p=kappa^2/(1+kappa^2);



model.w=w;
model.pred=label;
model.iteration=iter;
model.p=100*p;
model.acc=Acc;
model.Cost=Cost;


function [Miu1,Miu2]=update_mean(Miu1,Miu2,X,label0,label)

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
    Miu1=mean(X(label==1,:));
    Miu2=mean(X(label==2,:));
else
    if N1to2==0
        Miu1r=Miu1;
        Miu1to2=zeros(1,d);
    else
        Miu1to2=mean(X(idx1to2,:),1);
        Miu1r=(N1/N1r)*Miu1-(N1to2/N1r)*Miu1to2;
    end

    if N2to1==0
        Miu2r=Miu2;
        Miu2to1=zeros(1,d);
    else
        Miu2to1=mean(X(idx2to1,:),1);
        Miu2r=(N2/N2r)*Miu2-(N2to1/N2r)*Miu2to1;
    end
    
    
    Miu1=(N1r/(N1r+N2to1))*Miu1r+(N2to1/(N1r+N2to1))*Miu2to1;
    Miu2=(N2r/(N2r+N1to2))*Miu2r+(N1to2/(N2r+N1to2))*Miu1to2;   
end
if isnan(Miu1)|isnan(Miu2)
    keyboard
end




