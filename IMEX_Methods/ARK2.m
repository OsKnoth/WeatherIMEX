function ARK2
global A Ahat b bhat c chat r
%An efficient IMEX-DG solver for the compressible
%Navier-Stokes equations with a general equation
%of state
%Giuseppe Orlando(1)
%Paolo Francesco Barbante(1), Luca Bonaventura(1)
%November 29, 2021

r = 3;
gamma=2-sqrt(2);
alpha=1/2;
    
A=zeros(r,r);
b=zeros(r,1);
A(2,1)=gamma;
A(3,1)=1-alpha;
A(3,2)=alpha;
b(1)=1/2-gamma/4;
b(2)=1/2-gamma/4;
b(3)=gamma/2;

Ahat=zeros(r,r);
bhat=zeros(r,1);
Ahat(2,1)=gamma/2;
Ahat(2,2)=gamma/2;
Ahat(3,1)=1/(2*sqrt(2));
Ahat(3,2)=1/(2*sqrt(2));
Ahat(3,3)=1-1/sqrt(2);
bhat(1)=1/2-gamma/4;
bhat(2)=1/2-gamma/4;
bhat(3)=gamma/2;

c=A*ones(r,1); 
    
chat=Ahat*ones(r,1);
end