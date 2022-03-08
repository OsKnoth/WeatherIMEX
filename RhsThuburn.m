function [f]=RhsThuburn(y,Param)
% rho u v w Th
n=5*Param.nz-1;
Rho=y(1:Param.nz);
u=y(Param.nz+1:2*Param.nz);
v=y(2*Param.nz+1:3*Param.nz);
w=y(3*Param.nz+1:4*Param.nz);
Th=y(4*Param.nz:end);
% Rho
fRho=zeros(Param.nz,1);
% Advection
for iz=1:Param.nz
  fRho(iz)=-sqrt(-1)*Param.rho_ref(iz)*(Param.k_x*u(iz)...
    +Param.k_y*v(iz))...
    -sqrt(-1)*(Param.U*Param.k_x+Param.V*Param.k_y)*Rho(iz);
end
for iz=1:Param.nz-1
  Flux=0.5*(Param.rho_ref(iz)+Param.rho_ref(iz+1))*w(iz)/Param.dz;
  fRho(iz)=fRho(iz)-Flux;
  fRho(iz+1)=fRho(iz+1)+Flux;
end
% u
fu=zeros(Param.nz,1);
for iz=1:Param.nz
  pC=1/(1-Param.kappa)*Param.p_ref(iz)...
    *(Th(iz)/Param.Th_ref(iz)+Rho(iz)/Param.rho_ref(iz));
  fu(iz)=-sqrt(-1)*Param.k_x*pC/Param.rho_ref(iz)...
    -sqrt(-1)*(Param.U*Param.k_x+Param.V*Param.k_y)*u(iz)...
    -Param.D_mom*(Param.k_x^2+Param.k_y^2)^2*u(iz);
end
% v
fv=zeros(Param.nz,1);
for iz=1:Param.nz
  pC=1/(1-Param.kappa)*Param.p_ref(iz)...
    *(Th(iz)/Param.Th_ref(iz)+Rho(iz)/Param.rho_ref(iz));
  fv(iz)=-sqrt(-1)*Param.k_y*pC/Param.rho_ref(iz)...
    -sqrt(-1)*(Param.U*Param.k_x+Param.V*Param.k_y)*v(iz)...
    -Param.D_mom*(Param.k_x^2+Param.k_y^2)^2*v(iz);
end
% w
fw=zeros(Param.nz-1,1);
for iz=1:Param.nz-1
  pM=1/(1-Param.kappa)*Param.p_ref(iz)...
    *(Th(iz)/Param.Th_ref(iz)+Rho(iz)/Param.rho_ref(iz));
  pP=1/(1-Param.kappa)*Param.p_ref(iz+1)...
    *(Th(iz+1)/Param.Th_ref(iz+1)+Rho(iz+1)/Param.rho_ref(iz+1));
  RhoInvF=0.5*(1/Param.rho_ref(iz)+1/Param.rho_ref(iz+1));
  RhoF=0.5*(Rho(iz)+Rho(iz+1));
  fw(iz)=-RhoInvF*((pP-pM)/Param.dz+Param.Grav1*RhoF)...
    -sqrt(-1)*(Param.U*Param.k_x+Param.V*Param.k_y)*w(iz);
end
% Th
%dTh/dt = -w*dTh_ref/dz=-d(w*Th_ref)/dz+Th_ref*dw/dz
fTh=zeros(Param.nz,1);
for iz=1:Param.nz-1
  Flux=0.5*(Param.Th_ref(iz)+Param.Th_ref(iz+1))*w(iz)/Param.dz;
  fTh(iz)=fTh(iz)-Flux;
  fTh(iz+1)=fTh(iz+1)+Flux;
  fTh(iz)=fTh(iz)+w(iz)*Param.Th_ref(iz)/Param.dz;
  fTh(iz+1)=fTh(iz+1)-w(iz)*Param.Th_ref(iz+1)/Param.dz;   
end
for iz=1:Param.nz
  fTh(iz)=fTh(iz)-sqrt(-1)*(Param.U*Param.k_x+Param.V*Param.k_y)*Th(iz)...
    -Param.D_Th*(Param.k_x^2+Param.k_y^2)^2*Th(iz);
end
f=zeros(n,1);
f(1:Param.nz,1)=fRho;
f(Param.nz+1:2*Param.nz,1)=fu;
f(2*Param.nz+1:3*Param.nz,1)=fv;
f(3*Param.nz+1:4*Param.nz-1,1)=fw;
f(4*Param.nz:end,1)=fTh;
end





