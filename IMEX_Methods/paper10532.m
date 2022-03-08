function Rokhzadi
global A Ahat b bhat c chat r
% An Optimally Stable and Accurate Second-Order SSPRunge-Kutta IMEX Scheme 
%for Atmospheric Applications
% Arman Rokhzadi, Abdolmajid Mohammadian, and Martin Charron
r = 3;

A = [0 0 0;...
  0.711664700366941 0 0;...
  0.077338168947683 0.917273367886007 0];
b = [0.398930808264688 0.345755244189623 0.255313947545689]';

c = A*ones(r,1);

Ahat = [0 0 0 ;...
    0.353842865099275 0.353842865099275 0;...
    0.398930808264689 0.345755244189622 0.255313947545689];

bhat = Ahat(r,:)';
chat = Ahat*ones(r,1);
    
end