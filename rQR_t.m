function [Q,R]=rQR_t(A,P,sig,K)
[n,l]=size(A);
k=K;
Q=A;
s=1;
R=zeros(l,l);
S=zeros(k,l);
QS=zeros(k,l);
RS=zeros(l,l);
omega=zeros(1,l);
c=zeros(l,l);
for i=1:n
    for j=1:s
        S(P(i,j),1)=S(P(i,j),1)+Q(i,1).*sig(i,j)./sqrt(s);
    end
end
R(1,1)=norm(Q(:,1));
Q(:,1)=Q(:,1)/R(1,1);
S(:,1)=S(:,1)/R(1,1);
RS(1,1)=norm(S(:,1));
QS(:,1)=S(:,1)/RS(1,1);
c(1,1)=S(:,1)'*S(:,1);
vis=zeros(1,l);
for o=2:l
    for i=1:n
        for j=1:s
            S(P(i,j),o)=S(P(i,j),o)+Q(i,o).*sig(i,j)./sqrt(s);
        end
    end
    R(1:o-1,o)=RS(1:o-1,1:o-1)\(QS(:,1:o-1)'*S(:,o));
    omega(1:o-1)=S(:,o)'*S(:,1:o-1);
    for i=1:o-1
        vis(i)=0;
    end
    %t=cos(80/360*2*pi);
    t=0.4;
    getpos=[];
    for i=1:o-1
        abs_omega=abs(omega(1:o-1));
        abs_omega(vis(1:o-1)==1)=-1;
        [~,mxpos]=max(abs_omega);
        omega(1:o-1)=omega(1:o-1)-R(mxpos,o)*c(mxpos,1:o-1);
        S(:,o)=S(:,o)-R(mxpos,o)*S(:,mxpos);
        Q(:,o)=Q(:,o)-R(mxpos,o)*Q(:,mxpos);
        getpos=[getpos,mxpos];
        vis(mxpos)=1;
        mx=norm(omega,'inf');
        if mx/norm(Q(:,o))<=t
            break;
        end
    end
    for i=1:o-1
        if vis(i)==0
            R(i,o)=0;
            continue;
        end
    end
    %Q(:,o)=Q(:,o)-Q(:,getpos)*R(getpos,o);
    R(o,o)=norm(Q(:,o));
    S(:,o)=S(:,o)/R(o,o);
    Q(:,o)=Q(:,o)/R(o,o);
    QS(:,o)=S(:,o);
    for i=1:o-1
        RS(i,o)=(QS(:,o)')*QS(:,i);
        QS(:,o)=QS(:,o)-RS(i,o)*QS(:,i);
    end
    RS(o,o)=norm(QS(:,o));
    QS(:,o)=QS(:,o)/RS(o,o);
    for i=1:o
        c(i,o)=S(:,i)'*S(:,o);
        c(o,i)=c(i,o);
    end
end
