function Ros2
global A Ahat b bhat c chat r
ns=3;
r=3;
g11=1/2+1/6*sqrt(3);

a21=2/3;
a31=1/4
a32=3/4
ARos=[0 0 0 
   a21 0 0 
   a31 a32 0];
g21=-g11/a32
g22=g11
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


%     r = 6; coefa1 = 1/5; coefa2 = 1/5; coefa3 = 1/3; coefa4 = 1/2; coefa5 = 1;
%     A = [0 0 0 0 0 0 ; ...
%          coefa1 0 0 0 0 0 ; ...
%          0 coefa2 0 0 0 0 ; ...
%          0 0 coefa3 0 0 0 ; ...
%          0 0 0 coefa4 0 0 ; ...
%          0 0 0 0 coefa5 0];
%     c = A*ones(r,1); b = A(r,:)';
%     
%     ahat1 = 0; 
%     dhat1 = coefa1;
%     
%     ahat2 = 0; 
%     dhat2 = coefa2;
%     
%     ahat3 = 0; 
%     dhat3 = coefa3;
%     
%     ahat4 = 0; 
%     dhat4 = coefa4;
% 
%     eps = 1/18;
%     ahat51 = 1/2-4*eps; 
%     ahat52 = 5*eps;
%     dhat5 = 1/2 - eps;
% 
%     Ahat = [ 0 0 0 0 0 0 ; ahat1 dhat1 0 0 0 0 ; ahat2 0 dhat2 0 0 0 ; ...
%              ahat3 0 0 dhat3 0 0 ; ahat4 0 0 0 dhat4 0 ; ...
%              ahat51 ahat52 0 0 0 dhat5];
%     chat = Ahat*ones(r,1); bhat = Ahat(r,:)';
end