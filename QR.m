function [Q,R]=QR(A)
rng(0);
l=size(A,2);
R=zeros(l,l);
Q=A;
R(1,1)=norm(Q(:,1));
Q(:,1)=Q(:,1)/R(1,1);
for o=2:l
    for j=1:o-1
        R(j,o)=(Q(:,o)')*Q(:,j);
        Q(:,o)=Q(:,o)-R(j,o)*Q(:,j);
    end
    R(o,o)=norm(Q(:,o));
    Q(:,o)=Q(:,o)/R(o,o);
end
