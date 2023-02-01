function RMSD=mRNACostFunction(BetamakimaFit,LambdamakimaFit,params,kact,GeneTime,BetaRNA,LambdaRNA,avgBasalISGF3,avgBasalRNA)

x=avgBasalRNA;
time=[0:10:84000];
ISGF3=avgBasalISGF3;
kact=log10(kact);
params=[kact,params];

%Run model to steady state for initial values
[~,y_ss_opt]=ode15s(@(t,x) GeneSteadyState(t,x,ISGF3,params),time,x);

SSInitial=y_ss_opt(end);
t=[1:1:2500]; %time

%Using the Modified Akima cubic Hermite interpolation as the input function
[~,y_beta]=ode15s(@(t,x) ISGF3GeneReg(t,x,BetamakimaFit,params), ...
                        t,SSInitial);

[~,y_lambda]=ode15s(@(t,x) ISGF3GeneReg(t,x,LambdamakimaFit,params), ...
                        t,SSInitial);

%normalizing all simulated data to be scaled together
allGenes_sim=[y_beta,y_lambda];
minGene_sim=min(allGenes_sim,[],'all');
allGenes_sim=allGenes_sim-minGene_sim;
maxGene_sim=max(allGenes_sim,[],'all');

%Only take time points for RNA at experimental RNAseq time points
SimBeta=allGenes_sim(1+GeneTime,1);
SimNormBeta=SimBeta./maxGene_sim;

SimLambda=allGenes_sim(1+GeneTime,2);
SimNormLambda=SimLambda./maxGene_sim;
    
%Calculate RMSD between RNAseq IFN-beta induced expression data and model
error_beta=(BetaRNA-SimNormBeta).^2;
SSE_beta=sum(error_beta);
n=size(GeneTime,1);
RMSD_beta=sqrt((1/n)*SSE_beta);

%Calculate RMSD between RNAseq IFN-lambda induced expression data and model
error_lambda=(LambdaRNA-SimNormLambda).^2;
SSE_lambda=sum(error_lambda);
n=size(GeneTime,1);
RMSD_lambda=sqrt((1/n)*SSE_lambda);
 

RMSD=plus(RMSD_beta,RMSD_lambda);

fprintf('RMSD_beta:%1$.4f\n RMSD_lambda:%2$.4f\n RMSD_total:%3$.4f\n', ...
            RMSD_beta,RMSD_lambda,RMSD);
