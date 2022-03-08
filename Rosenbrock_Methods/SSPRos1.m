function SSPRos1
global A Ahat b bhat c chat r
r=5;
ARos = [0 0 0 0 0
        1./2. 0 0 0 0
        1./2. 1./2. 0 0 0
        1./6. 1./6. 1./6. 0 0
        1./6. 1./6. 1./6. 1./2. 0];

Gamma= [1./2. 0 0 0
        0.0 3./4. 0 0
        -2./3. -23./9. 2./9. 0
        1./18. 65./108. -2./27 0];


AL=ARos(2:end,1:end-1);
A=ARos;
AIU=AL*Gamma*inv(AL);
AI=[zeros(1,r)
    zeros(r-1,1) AIU];
%Ahat=AI+A;
Ahat=AI;
b=A(r,:)';
c=A*ones(r,1);
bhat=Ahat(r,:)';
chat=Ahat*ones(r,1);
end