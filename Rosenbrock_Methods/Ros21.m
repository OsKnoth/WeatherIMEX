function Ros21
global A Ahat b bhat c chat r
ns=3;
r=3;
g11=1+1/sqrt(2);

a31=1/2;
a32=3/2;
a21=1/(2*a32);
ARos=[0 0 0 
   a21 0 0 
   a31 a32 0];
g21=-g11/a32;
g22=g11;
Gamma=[g11 0
       g21 g22];

AL=ARos(2:end,1:end-1);
A=ARos;
AIU=AL*Gamma*inv(AL);
AI=[zeros(1,3)
    zeros(2,1) AIU];
%Ahat=AI+A;
Ahat=AI;
b=A(r,:)';
c=A*ones(r,1);
bhat=Ahat(r,:)';
chat=Ahat*ones(r,1);

end