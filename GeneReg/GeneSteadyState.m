function deltamRNA=GeneSteadyState(t,x,ISGF3,params)

%params(1)= kact, maximal expression level of the promoter
%params(2)= Ka, activation coeff.(concentration at half-maximal expression)
%params(3)= kdeg, mRNA degradation rate
%params(4)= n, Hill coeff. (governs steepness of input function)

kact=params(1); %units: min^-1
Ka=params(2); %units: nM
kdeg=params(3); %units: min^-1
n=params(4); %units:
%ISGF3=;


productionHillFunc=(kact*(ISGF3.^n))./((Ka.^n)+(ISGF3.^n));
degrRate=kdeg*x;

deltamRNA=productionHillFunc-degrRate;

end