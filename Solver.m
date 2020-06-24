function [F] = Solver(GH, Constants, G0,Hi)
%% Mass Balance Equations. 

KA41= 2.75*10^7; %Association constant of (pR)4/(pS)4 with cadaverine measured by ITC.
KA31= Constants(1); %Association constant of (pR)3/(pS)3 with cadaverine.
KA4C11= Constants(2); %First Association constant of (pR)3pS / (pS)3pR with cadaverine.
KA4C12= Constants(3); %Second Association constant of (pR)3pS / (pS)3pR with cadaverine.
KA4H1= Constants(4); %First Association constant of (pR)2(pS)2 with cadaverine.
KA4H2 = Constants(5); %Second Association constant of (pR)2(pS)2 with cadaverine.

%% Mass Balances on Hosts
F(1) = Hi(1) - GH(1) - KA41*GH(1)*GH(5);
F(2) = Hi(2) - GH(2) - KA31*GH(2)*GH(5);
F(3) = Hi(3) - GH(3) - KA4C11*GH(3)*GH(5) - KA4C11*KA4C12*GH(3)*GH(5)^2;
F(4) = Hi(4) - GH(4) - KA4H1*GH(4)*GH(5)- KA4H1*KA4H2*GH(4)*GH(5)^2;
%% Mass Balance on Guest
F(5) = G0 - GH(5) -  KA41*GH(1)*GH(5) - KA31*GH(2)*GH(5) - KA4C11*GH(3)*GH(5) - 2*KA4C11*KA4C12*GH(3)*GH(5)^2 - KA4H1*GH(4)*GH(5) - 2*KA4H1*KA4H2*GH(4)*GH(5)^2;

end

