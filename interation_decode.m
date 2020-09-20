N=10;
B=round(rand(N,4));
G=[1 0 0 0 1 0 1;0 1 0 0 1 1 1;0 0 1 0 1 1 0;0 0 0 1 0 1 1];
C=mod(B*G,2);
C=[C C(:,7)];
%BSC system
p=0.2;
R=bsc(C,p);
pri=1;
D=zeros(N,4);
for h=1:N
    [m51,m52]=state(p,R(h,5),R(h,3));
    [m61,m62]=state(p,R(h,6),R(h,4));
    [m71,m62]=state(p,R(h,7),R(h,1));
    D(h,3)=R(h,3);
    for i=1:10
        if D(h,3)==R(h,2)
            u1=m61;
        else
            u1=m62;
        end
        
        