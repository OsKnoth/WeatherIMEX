function p=pRef(z,Param)
switch Param.Profile
  case 'Isentropic'
    % T(z)=T_0
    % dp/dz = -rho*g = -p/(R*T)*g = -g/(R*T_0)*p
    % ln(p) = -g/(R*T_0)*z
    % p = p_0*exp(-g/(R*T_0)*z)
    p=Param.p_0*exp(-Param.Grav/(Param.R*Param.T0)*z);
end
end