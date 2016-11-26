function model=mspc_eig(X,paras)

[N,d]=size(X);

Miu=sum(X,1)/N;
X=bsxfun(@minus,X,Miu);

reg=paras.rho*var(X)+1e-10;

t_start=tic;
if N<d
    Xbar=bsxfun(@times,X,1./sqrt(reg));
    G=Xbar*Xbar';
    A=eye(N)-inv(eye(N)+G);
    opts.tol = 1e-9;
    opts.issym=1;
    opts.disp = 0;
    [q,eval]=eigs(A,1,'lm',opts);
else
    C=X'*X;
    A=X*((C+diag(reg))\X');
    A=(A'+A)/2;
    opts.tol = 1e-9;
    opts.issym=1;
    opts.disp = 0;
    [q,eval]=eigs(A,1,'lm',opts);
end
model.T=toc(t_start);

model.q=(sign(q)+3)/2;
acc_sign=accuracy(model.q,paras.y);
w=X'*q;
pred=(sign(X*w)+3)/2;
acc_w=accuracy(pred,paras.y);

model.w=w;
model.acc_sign=acc_sign;
model.acc_w=acc_w;
model.pred=pred;
