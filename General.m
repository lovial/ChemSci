%% General script for the concomitant fitting of association constants, as 
%%in "The dark side of disulfide-based dynamic combinatorial chemistry" 
%%from L. Vial et al.
%%Example is given for 1:1 associations stoechiometries for (pR)4/(pS)4 and
%%(pR)3/(pS)3 pairing with cadaverine and 1:2 pairing of (pR)3pS / (pS)3pR and (pR)2(pS)2 with cadaverine
%%


%% Let's get some data
G0 = DATA(:,1); % matrix ov Guest concentration (M)
HData = DATA(:,2:15); %matrix of 1H NMR signal displacement (ppm)
H0 = 9.13 * 10^-3 ; %M , Initial Host concentration
Repartition= [30;40;20;10]; % Molar fraction of (pR)4/(pS)4, (pR)3/(pS)3, (pR)3pS / (pS)3pR, (pR)2(pS)2 in Host pool respectively

Hi = Repartition (:,1)*H0/100 ;% Initial concentration in each host species

for i = 1:size(HData,2)
    Dmax(i) = max(HData(:,i)); %Maximal 1H NMR signal displacement
end
HDataNorm = HData./Dmax;
Kres(1:5)
%% Parameter Fitting

fun = @(K,G0) Avancement(G0,K,Hi);
Constants0 = [ 10^2; 10^2;10^2; 10^2;10^2;1;1;1;1;1;1;1;1;1;1;1;1;1;1]; %Initiation of parameters. Here, the first five are association costants, the subsequents are the ratio of induced 1H NMR signal displacement upon second binding event by that of first binding event.
LB = [0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0]; %Lower boundaries of parameters
[Kres, resnorm] = lsqcurvefit(fun,Constants0,G0,HDataNorm, LB); % Kres contains the fitting parameters after fitting
HResult = fun(Kres,G0).*Dmax;
%% Calcul of the coefficient of determination
Coeff_Deter = zeros(size(HData,2),1);
for j =1:size(HData,2)
    Somme = 0;
    e = 0;
    Mean = mean(HData(:,j));
    for i = 1:size(HData,1)
        e = e + (HResult(i,j)-HData(i,j))^2;
        Somme = Somme + (HData(i,j)-MoyenneA3)^2;
    end
    Coeff_Deter(j)=1-e/Somme;
   
end
Coeff_Deter


%% Plotting

figure
for i = 1:size(HData,2)
    subplot(4,4,i)
    plot(G0,HData(:,i),'+')
    hold on
    plot(G0,HResult(:,i),'-')
end
figure

plot(G0,HData,'+')
hold on
set(gca,'ColorOrderIndex',1)
plot(G0,HResult,'-')


