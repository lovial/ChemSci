function [Result] = Avancement(G0,K,Hi,Dmax)
%% For a set of parameters, calculate the expected 1H NMR signal displacement for each species.


Nb_Pts = size(G0,1);
Residuals = zeros(Nb_Pts,5); 
Result = zeros(Nb_Pts,5);
Constants = K(1:5); % Association constants
DHG = K(6:11); %Ratio of induced 1H NMR signal displacement upon second binding event by that of first binding event.

for i = 1:Nb_Pts
    X = Function (G0(i),Constants,Hi);
    Residuals(i,1:5)= X(1:5); % Concentration of unbound species.
%%  Calculated 1H NMR signal displacement of (pR)3/(pS)3
    Result(i,1) = Dmax(1)*(Hi(2)-Residuals(i,2))/Hi(2);
%%  Calculated displacement of the four 1H NMR Signals of (pR)3pS / (pS)3pR
    Result(i,2) = Dmax(2)*(Constants(3)*Constants(2)*Residuals(i,3)*Residuals(i,5)^2)/Hi(3)+(DHG(1)*Constants(2)*Residuals(i,3)*Residuals(i,5))/Hi(3);
    Result(i,3) = Dmax(3)*(Constants(3)*Constants(2)*Residuals(i,3)*Residuals(i,5)^2)/Hi(3)+(DHG(2)*Constants(2)*Residuals(i,3)*Residuals(i,5))/Hi(3);
    Result(i,4) = Dmax(4)*(Constants(3)*Constants(2)*Residuals(i,3)*Residuals(i,5)^2)/Hi(3)+(DHG(3)*Constants(2)*Residuals(i,3)*Residuals(i,5))/Hi(3);
    Result(i,5) = Dmax(5)*(Constants(3)*Constants(2)*Residuals(i,3)*Residuals(i,5)^2)/Hi(3)+(DHG(4)*Constants(2)*Residuals(i,3)*Residuals(i,5))/Hi(3);
%%  Calculated 1H NMR signal displacement of(pR)2(pS)2 
    Result(i,6) = Dmax(6)*(Constants(4)*Constants(5)*Residuals(i,4)*Residuals(i,5)^2)/Hi(4)+(DHG(5)*Constants(4)*Residuals(i,4)*Residuals(i,5))/Hi(4);
    Result(i,7) = Dmax(7)*(Constants(4)*Constants(5)*Residuals(i,4)*Residuals(i,5)^2)/Hi(4)+(DHG(6)*Constants(4)*Residuals(i,4)*Residuals(i,5))/Hi(4);
end


end

