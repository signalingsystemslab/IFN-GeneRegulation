function deltamRNA = ISGF3GeneReg(t,x,ISGF3makimaFit,params)

%Parameter definition
%params(1)= kact, maximal expression level of the promoter
%params(2)= Ka, activation coeff.(concentration at half-maximal expression)
%params(3)= kdeg, mRNA degradation rate
%params(4)= n, Hill coeff. (governs steepness of input function)

kact=params(1); %units: min^-1
Ka=params(2); %units: nM
kdeg=params(3); %units: min^-1
n=params(4); %units:
ISGF3=ISGF3makimaFit;

%Since ISGF3 is a function of time, need to interpolate the curve to find
%the ISGF3 amount at each time point
if t<=1
    ISGF3_t=ISGF3(1);
elseif t>=length(ISGF3)
    ISGF3_t=ISGF3(end);
else
    ISGF3_tPrev=floor(t); %rounds the toward negative infinity to determine the previous time
    ISGF3_tNext=floor(t+1); %rounds the next time toward negative infinity
    ISGF3_rise=ISGF3(ISGF3_tNext)-ISGF3(ISGF3_tPrev);
    t_run=ISGF3_tNext-ISGF3_tPrev;
    slope=ISGF3_rise./t_run;
    tdiff=t-floor(t);
    ISGF3_t=ISGF3(ISGF3_tPrev)+(slope*tdiff); %using the point-slope equation to calculate ISGF3 at t
end


productionHillFunc=(kact*(ISGF3_t.^n))./((Ka.^n)+(ISGF3_t.^n));
degrRate=kdeg*x;

deltamRNA=productionHillFunc-degrRate;

end