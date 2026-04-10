rng(0);
n=20000;
l=50;
K=2000;
U=rand(n,l);
V=rand(l,l);
S=zeros(l,l);
for i=1:l
    S(i,i)=1/exp(i);
end
U=orth(U);
V=orth(V);
A=U*S*V';
[n,l]=size(A);
if n<l
    A=A';
    [n,l]=size(A);
end
s=min(l,8);
P=zeros(n,1);
sig=zeros(n,1);
for i=1:n
    j=1;
    while j<=s
        P(i,j)=ceil(rand*K);
        tp=1;
        for o=1:j-1
            if P(i,j)==P(i,o)
                tp=0;
                break;
            end
        end
        if tp==1
            j=j+1;
        end
    end
end
for i=1:n
    for j=1:s
        sig(i,j)=ceil(rand*2);
        if(sig(i,j)>1)
            sig(i,j)=-1;
        end
    end
end

[Q,R]=QR(A);
err_QR=norm(A-Q*R);
cond_QR=cond(Q);

[Q,R]=rQR_t(A,P,sig,K);
err_rQRt=norm(A-Q*R);
cond_rQRt=cond(Q);

[Q,R]=shiftedCholeskyQR3(A);
err_shiftedCholeskyQR3=norm(A-Q*R);
cond_shiftedCholeskyQR3=cond(Q);
