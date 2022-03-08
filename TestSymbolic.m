% wP
% rhoP ThP
% w
% rhoM ThM
% wM
%p/pRef=T/TRef + Rho/RhoRef;
%Th/ThRef=T/TRef-kappa*p/pRef
%T/TRef=Th/ThRef+kappa*p/pRef
%p/Pref=Th/ThRef+kappa*p/pRef+Rho/RhoRef;
%(1-kappa)*p/pRef=Th/ThRef+Rho/RhoRef;
%p=1/(1-kappa)*pRef*(Th/ThRef+Rho/RhoRef);
%dp/dt=1/(1-kappa)*pRef*(dTh/dt/ThRef+dRho/dt/RhoRef);
%     =1/(1-kappa)*pRef*(-w*dThRef/dz/Thref+d(w*Rho)/dz/RhoRef)

Param.R=287.05;
Param.T0=250;
Param.Cp=1005.0;
Param.Cv=Param.Cp-Param.R;
Param.Gamma=Param.Cp/Param.Cv;
kappa=Param.R/Param.Cp;
Param.Profile='Isentropic';
Param.T0=250;
Param.Grav=9.80616;
Param.p_0=1.e5;
Param.kappa=Param.R/Param.Cp;

zP=15;
zM=5;
pRefP=pRef(zP,Param);
pRefM=pRef(zM,Param);
RhoRefP=RhoRef(zP,Param);
RhoRefM=RhoRef(zM,Param);
ThRefP=ThRef(zP,Param);
ThRefM=ThRef(zM,Param);

pRef1P=(Param.R*RhoRefP*ThRefP/Param.p_0^Param.kappa)^(1/(1-Param.kappa));
dpdRho=1/(1-Param.kappa)...
  *(Param.R*RhoRefP*ThRefP/Param.p_0^Param.kappa)^(1/(1-Param.kappa))...
  *(Param.R*ThRefP/Param.p_0^Param.kappa)...
  /(Param.R*RhoRefP*ThRefP/Param.p_0^Param.kappa);
%pM=1/(1-kappa)*pRefM*(ThM/ThRefM+RhoM/RhoRefM);
%pP=1/(1-kappa)*pRefP*(ThP/ThRefP+RhoP/RhoRefP);

RhoRefPM=0.5*(RhoRefP+RhoRefM);
ThRefPM=0.5*(ThRefP+ThRefM);
%RhoP=-(wP*RhoRefPP - w*RhoRefPM);
%RhoM=-( w*RhoRefPM -wM*RhoRefMM);
%w = -1/RhoRefPM*(pP-pM);
%Th=-w*dThRef/dz = -d(w*ThRef)/dz+ThRef*dw/dz
%ThP=-(wP*ThRefPP-w*ThRefPM)+(wP-w)*ThRefP, w*(ThRefPM-ThRefP)
%ThM=-(w*ThRefPM-wM*ThRefMM)+(w-wM)*ThRefM, w*(ThRefM-ThRefPM)
A=[         0 0 -RhoRefPM 0 0
            0 0  RhoRefPM 0 0
            0 0     0     0 0
            0 0  (ThRefPM-ThRefP)     0 0
            0 0  (ThRefM-ThRefPM)     0 0];
A(3,1)=1/RhoRefPM/(1-kappa)*pRefM/RhoRefM;
A(3,2)=-1/RhoRefPM/(1-kappa)*pRefP/RhoRefP;
A(3,4)=1/RhoRefPM/(1-kappa)*pRefM/ThRefM;
A(3,5)=-1/RhoRefPM/(1-kappa)*pRefP/ThRefP;
cS=sqrt(abs(A(3,:)*A(:,3)));
%c=sqrt(Gamma*R*T);
cSE=sqrt(Param.Gamma*Param.R*Param.T0);
cSEN=sqrt(Param.Gamma*Param.R*293.15);
aa=3;

