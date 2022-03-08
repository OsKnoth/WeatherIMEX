function LinearThuburn()
global A Ahat b bhat c chat r
close all
% du/dt = -1/rho_ref*dp/dx + f*v - Dd4u/dx4 
% dv/dt = -1/rho_ref*dp/dy - f*u - Dd4v/dx4 
% dw/dt = -1/rho_ref*(dp/dz + rho*g) - Dd4w/dx4
% drho/dt = -rho_ref*(du/dx+dv/dy+dw/dz)-w*drho_ref/dz
% dTh/dt = -w*dTh_ref/dz - Dd4Th/dx4
% Linearized equation of state
%
% p/p_ref=T/T_ref + rho/rho_ref
%
% Th/Th_ref = T/T_ref-kappa*p/p_ref
%
%Th/Th_ref = (1-kappa)*p/p_ref - rho/rho_ref 
%
%p = 1/(1-kappa)*p_ref*(Th/Th_ref + rho/rho_ref)
%
% Fourier Ansatz \psi= \hat{psi}*exp(i(k_x*x+k_y*y))
%
% drho/dt = -i*rho_ref*(k_x*u+k_y*v)-rho_ref*dw/dz-w*drho_ref/dz
% du/dt = -ik_x/rho_ref/(1-kappa)*p_ref*(Th/Th_ref+rho/rho_ref)
% dv/dt = -ik_y/rho_ref/(1-kappa)*p_ref*(Th/Th_ref+rho/rho_ref)
% dw/dt = -1/rho_ref/(1-kappa)*p_ref/Th_ref*dTh/dz...
%         -1/rho_ref/(1-kappa)*p_ref/rho_ref*drho/dz- 1/rho_ref*rho*g
% dTh/dt = -w*dTh_ref/dz
% Height in vertical direction D = 10^5 m
% Meridional wavenumber k_y = 0
% Temperature of reference state T0 = 250 K
% Gravitational acceleration g = 9.80616 m s^?2
% Thermodynamic constants R = 287.05 J kg^?1 K^?1 
% and cp = 1005.0 J kg^?1 K^?1

% Vertical discretization
T_x=220*10^3;
T_y=220*10^3;
Param.Profile='Isentropic';
Param.R=287.05;
Param.Grav=9.80616;
Param.Grav1=9.80616;
Param.T0=250;
Param.Cp=1005.0;
Param.kappa=287.05/Param.Cp;
Param.p_0=1.e5;
Param.k_x=2*pi/T_x;
Param.k_y=2*pi/T_y;
Param.U=30; 
Param.V=30;
Param.D_Th=2.e12;
Param.D_mom=2.e12;
nz=10;
H=30*10^3; %30*10^3;
dz=H/nz;
Param.numVar=5*nz-1;
zP=zeros(nz+1,1);
zM=zeros(nz,1);
for iz=1:nz
  zP(iz+1)=zP(iz)+dz;
  zM(iz)=0.5*(zP(iz)+zP(iz+1));
end
Param.nz=nz;
Param.H=H;
Param.dz=dz;
Param.zM=zM;
Param.zP=zP;
u=zeros(nz,1);
rho=zeros(nz,1);
Th=zeros(nz,1);
w=zeros(nz-1,1);
Param.rho_ref=zeros(nz,1);
Param.p_ref=zeros(nz,1);
Param.Th_ref=zeros(nz,1);
for iz=1:nz
  Param.rho_ref(iz)=RhoRef(zM(iz),Param);
  Param.p_ref(iz)=pRef(zM(iz),Param);
  Param.Th_ref(iz)=ThRef(zM(iz),Param);
end


%chooseIMEX(10030)
%chooseIMEX(102321);
%chooseIMEX(10101)
chooseIMEX('Rokhzadi') % Rokhzadi
%SSPRos1()
%ARKODE_ARK437L2SA_ERK_7_3_4()
%ARK2()
locdt=1;
maxeig=zeros(300,2);
dt=zeros(300,1);
[L,N]=JacThuburn(Param);
for i=1:300
  imexA=IMEXRKstabmat(N,L,Param.numVar,locdt,A,Ahat,b,bhat,r);
  kk=eig(imexA);
  maxeig(i,1)=max(abs(kk));
  dt(i)=locdt;
  locdt=locdt+1.0;
end
figure(1)
plot(dt,maxeig(:,1)-1)

T_x=110*10^3;
Param.k_x=2*pi/T_x;
locdt=1;
[L,N]=JacThuburn(Param);
for i=1:300
  imexA=IMEXRKstabmat(N,L,Param.numVar,locdt,A,Ahat,b,bhat,r);
  kk=eig(imexA);
  maxeig(i,2)=max(abs(kk));
  dt(i)=locdt;
  locdt=locdt+1.0;
end
figure(1)
plot(dt,maxeig(:,1)-1)
figure(2)
plot(dt,maxeig(:,2)-1)

end

