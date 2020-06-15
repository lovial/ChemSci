function [Residuals] = Function(G0,Constants,Hi)
%% Calculate the concentration of unbound species as a function of association constant values.

options = optimset('Display','off');

x0 = [0,0,0,0,0];
x = fsolve(@(Res) Solver(Res, Constants, G0,Hi),x0,options);

Residuals = real(x);
end
