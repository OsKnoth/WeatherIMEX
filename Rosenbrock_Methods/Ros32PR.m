function Ros32PR
global A Ahat b bhat c chat r
r=4;
gg = 7.8867513459481287e-01;
ARos=[0 0 0 0
   2.3660254037844388e+00 0 0 0
   0                      1 0 0
   2.9266384402395124e-01 -8.1338978618764143e-02 7.8867513459481287e-01 0];
Gamma=[gg 0 0 
       -2.3660254037844388e+00 gg 0
       -2.8468642516567449e-01 -1.0813389786187642e+00 gg];


AL=ARos(2:end,1:end-1);
A=ARos;
AIU=AL*Gamma*inv(AL);
AI=[zeros(1,r)
    zeros(r-1,1) AIU];
Ahat=AI;
b=A(r,:)';
c=A*ones(r,1);
bhat=Ahat(r,:)';
chat=Ahat*ones(r,1);

end