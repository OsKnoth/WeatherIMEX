function [L,N]=JacThuburn(Param)
% rho u v w Th
n=5*Param.nz-1;
RowIndL=zeros(10*2*Param.nz,1);
ColIndL=zeros(10*2*Param.nz,1);
ValL=zeros(10*2*Param.nz,1);
RowIndN=zeros(10*2*Param.nz,1);
ColIndN=zeros(10*2*Param.nz,1);
ValN=zeros(10*2*Param.nz,1);
iSN=1;
iSL=1;

%rho rho
irho=1;
for iz=1:Param.nz
  RowIndN(iSN)=irho;
  ColIndN(iSN)=irho;
  ValN(iSN)=-sqrt(-1)*(Param.U*Param.k_x+Param.V*Param.k_y);
  iSN=iSN+1;
  irho=irho+1;
end
%rho u
irho=1;
iu=Param.nz+1;
for iz=1:Param.nz
  RowIndN(iSN)=irho;
  ColIndN(iSN)=iu;
  ValN(iSN)=-sqrt(-1)*Param.rho_ref(iz)*Param.k_x;
  iSN=iSN+1;
  irho=irho+1;
  iu=iu+1;
end

%rho v
irho=1;
iv=2*Param.nz+1;
for iz=1:Param.nz
  RowIndN(iSN)=irho;
  ColIndN(iSN)=iv;
  ValN(iSN)=-sqrt(-1)*Param.rho_ref(iz)*Param.k_y;
  iSN=iSN+1;
  irho=irho+1;
  iv=iv+1;
end

%rho w
irho=1;
iw=3*Param.nz+1;
for iz=1:Param.nz-1
  rho_ref_F=0.5*(Param.rho_ref(iz)+Param.rho_ref(iz+1));
  RowIndL(iSL)=irho;
  ColIndL(iSL)=iw;
  ValL(iSL)=-rho_ref_F/Param.dz;
  iSL=iSL+1;
  RowIndL(iSL)=irho+1;
  ColIndL(iSL)=iw;
  ValL(iSL)=+rho_ref_F/Param.dz;
  iSL=iSL+1;
  irho=irho+1;
  iw=iw+1;
end

% u rho
irho=1;
iu=Param.nz+1;
for iz=1:Param.nz
  RowIndN(iSN)=iu;
  ColIndN(iSN)=irho;
  ValN(iSN)=-sqrt(-1)*Param.k_x/Param.rho_ref(iz)/(1-Param.kappa)...
    *Param.p_ref(iz)/Param.rho_ref(iz);
  iSN=iSN+1;
  irho=irho+1;
  iu=iu+1;
end

% u u
iu=Param.nz+1;
for iz=1:Param.nz
  RowIndN(iSN)=iu;
  ColIndN(iSN)=iu;
  ValN(iSN)=-Param.D_mom*(Param.k_x^2+Param.k_y^2)^2....
    -sqrt(-1)*(Param.U*Param.k_x+Param.V*Param.k_y);
  iSN=iSN+1;
  iu=iu+1;
end

% u Th
iTh=4*Param.nz;
iu=Param.nz+1;
for iz=1:Param.nz
  RowIndN(iSN)=iu;
  ColIndN(iSN)=iTh;
  ValN(iSN)=-sqrt(-1)*Param.k_x/Param.rho_ref(iz)/(1-Param.kappa)...
    *Param.p_ref(iz)/Param.Th_ref(iz);
  iSN=iSN+1;
  iTh=iTh+1;
  iu=iu+1;
end

% v rho
irho=1;
iv=2*Param.nz+1;
for iz=1:Param.nz
  RowIndN(iSN)=iv;
  ColIndN(iSN)=irho;
  ValN(iSN)=-sqrt(-1)*Param.k_y/Param.rho_ref(iz)/(1-Param.kappa)...
    *Param.p_ref(iz)/Param.rho_ref(iz);
  iSN=iSN+1;
  irho=irho+1;
  iv=iv+1;
end

% v v
iv=2*Param.nz+1;
for iz=1:Param.nz
  RowIndN(iSN)=iv;
  ColIndN(iSN)=iv;
  ValN(iSN)=-Param.D_mom*(Param.k_x^2+Param.k_y^2)^2....
    -sqrt(-1)*(Param.U*Param.k_x+Param.V*Param.k_y);
  iSN=iSN+1;
  iv=iv+1;
end



% v Th
iTh=4*Param.nz;
iv=2*Param.nz+1;
for iz=1:Param.nz
  RowIndN(iSN)=iv;
  ColIndN(iSN)=iTh;
  ValN(iSN)=-sqrt(-1)*Param.k_y/Param.rho_ref(iz)/(1-Param.kappa)...
    *Param.p_ref(iz)/Param.Th_ref(iz);
  iSN=iSN+1;
  iTh=iTh+1;
  iv=iv+1;
end
%w rho
irho=1;
iw=3*Param.nz+1;
for iz=1:Param.nz-1
  FakRho=0.5*(1/Param.rho_ref(iz)+1/Param.rho_ref(iz+1)); 
  RowIndL(iSL)=iw;
  ColIndL(iSL)=irho;
  ValL(iSL)=FakRho/(1-Param.kappa)*Param.p_ref(iz)...
          /Param.rho_ref(iz)/Param.dz;
  iSL=iSL+1;
  RowIndL(iSL)=iw;
  ColIndL(iSL)=irho+1;
  ValL(iSL)=-FakRho/(1-Param.kappa)*Param.p_ref(iz+1)...
          /Param.rho_ref(iz+1)/Param.dz;
  iSL=iSL+1;
  RowIndL(iSL)=iw;
  ColIndL(iSL)=irho;
  ValL(iSL)=-0.5*FakRho*Param.Grav1;
  iSL=iSL+1;
  RowIndL(iSL)=iw;
  ColIndL(iSL)=irho+1;
  ValL(iSL)=-0.5*FakRho*Param.Grav1;
  iSL=iSL+1;
  iw=iw+1;
  irho=irho+1;
end

% w w
iw=3*Param.nz+1;
for iz=1:Param.nz-1
  RowIndN(iSN)=iw;
  ColIndN(iSN)=iw;
  ValN(iSN)=-sqrt(-1)*(Param.U*Param.k_x+Param.V*Param.k_y);
  iSN=iSN+1;
  iw=iw+1;
end

%w Th
iTh=4*Param.nz;
iw=3*Param.nz+1;
for iz=1:Param.nz-1
%   pM=1/(1-Param.kappa)*Param.p_ref(iz)...
%     *(Th(iz)/Param.Th_ref(iz)+Rho(iz)/Param.rho_ref(iz));
%   pP=1/(1-Param.kappa)*Param.p_ref(iz+1)...
%     *(Th(iz+1)/Param.Th_ref(iz+1)+Rho(iz+1)/Param.rho_ref(iz+1));
%   RhoInvF=0.5*(1/Param.rho_ref(iz)+1/Param.rho_ref(iz+1));
%   RhoF=0.5*(Rho(iz)+Rho(iz+1));
%   fw(iz)=-RhoInvF*((pP-pM)/Param.dz+Param.Grav1*RhoF);
  FakRho=0.5*(1/Param.rho_ref(iz)+1/Param.rho_ref(iz+1)); 
  RowIndL(iSL)=iw;
  ColIndL(iSL)=iTh;
  ValL(iSL)=FakRho/(1-Param.kappa)*Param.p_ref(iz)...
          /Param.Th_ref(iz)/Param.dz;
  iSL=iSL+1;
  RowIndL(iSL)=iw;
  ColIndL(iSL)=iTh+1;
  ValL(iSL)=-FakRho/(1-Param.kappa)*Param.p_ref(iz+1)...
          /Param.Th_ref(iz+1)/Param.dz;
  iSL=iSL+1;
  iw=iw+1;
  iTh=iTh+1;
end
%Th w
iw=3*Param.nz+1;
iTh=4*Param.nz;
for iz=1:Param.nz-1
%   Flux=0.5*(Param.Th_ref(iz)+Param.Th_ref(iz+1))*w(iz)/Param.dz;
%   fTh(iz)=fTh(iz)-Flux;
%   fTh(iz+1)=fTh(iz+1)+Flux;
%   fTh(iz)=fTh(iz)+w(iz)*Param.Th_ref(iz)/Param.dz;
%   fTh(iz+1)=fTh(iz+1)-w(iz)*Param.Th_ref(iz+1)/Param.dz;
  Fak=0.5*(Param.Th_ref(iz+1)...
          +Param.Th_ref(iz));
  RowIndL(iSL)=iTh;
  ColIndL(iSL)=iw;
  ValL(iSL)=-Fak/Param.dz;
  iSL=iSL+1;
  RowIndL(iSL)=iTh+1;
  ColIndL(iSL)=iw;
  ValL(iSL)=Fak/Param.dz;
  iSL=iSL+1;
%   fTh(iz)=fTh(iz)+w(iz)*Param.Th_ref(iz)/Param.dz;
%   fTh(iz+1)=fTh(iz+1)-w(iz)*Param.Th_ref(iz+1)/Param.dz;
  RowIndL(iSL)=iTh;
  ColIndL(iSL)=iw;
  ValL(iSL)=Param.Th_ref(iz)/Param.dz;
  iSL=iSL+1;
  RowIndL(iSL)=iTh+1;
  ColIndL(iSL)=iw;
  ValL(iSL)=-Param.Th_ref(iz+1)/Param.dz;
  iSL=iSL+1;
  iTh=iTh+1;
  iw=iw+1;
end

% Th Th
iTh=4*Param.nz;
for iz=1:Param.nz
  RowIndN(iSN)=iTh;
  ColIndN(iSN)=iTh;
  ValN(iSN)=-Param.D_Th*(Param.k_x^2+Param.k_y^2)^2....
    -sqrt(-1)*(Param.U*Param.k_x+Param.V*Param.k_y);
  iSN=iSN+1;
  iTh=iTh+1;
end

L=sparse(RowIndL(1:iSL-1),ColIndL(1:iSL-1),ValL(1:iSL-1),n,n);
N=sparse(RowIndN(1:iSN-1),ColIndN(1:iSN-1),ValN(1:iSN-1),n,n);
end





