function [newParams,minRMSD]=RNALocalOptim(params,BetamakimaFit,LambdamakimaFit,GeneTime,BetaRNA,LambdaRNA,avgBasalISGF3,avgBasalRNA)


%Only fitting kdeg, Ka, and n parameters since the kact does not appear
%sensitive to achieving the IFN specific responses
kact=params(1);
params=params(2:4);

%options for fminsearch
options=optimset('PlotFcns',@optimplotfval,'MaxIter',300*length(params));


%use when fitting all parameters except kact with constraints
LB = [0.1,0.0006931,0]; %define bounds for lower constraints
UB = [10,0.6931,10]; %define bounds for upper constraints

%run parameter fitting
[optim_param,fval,exitflag]=fminsearchbnd(@(params)...
                                mRNACostFunction(BetamakimaFit,...
                                LambdamakimaFit,10.^params,10.^kact, ...
                                GeneTime,BetaRNA,LambdaRNA, ...
                                avgBasalISGF3,avgBasalRNA), ...
                                log10(params),log10(LB),log10(UB),options);

newParams=10.^optim_param;
minRMSD=fval;

%Determine why optimization was stopped
    if exitflag==1
        stop="Solution Found";
    
        elseif exitflag==0
            stop="Exceeded Iterations Allowed";
            
    end

%Redefining params and new params to include kact
params=[kact,params];
newParams=[kact,newParams];

fprintf('OldParams values:')
disp(params)

fprintf('NewParams values:')
disp(newParams)

fprintf('Minimum RMSD:')
disp(minRMSD)

fprintf('Stop Value: ')
disp(stop)


end
