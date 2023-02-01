function RMSD=mRNACostFunction_ParamSensitivity(allGenes_sim,maxGene_sim,GeneTime,BetaRNA,LambdaRNA)


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

fprintf('RMSD_beta:%1$.4f\n RMSD_lambda:%2$.4f\n RMSD_total:%3$.4f\n',...
             RMSD_beta,RMSD_lambda,RMSD);
