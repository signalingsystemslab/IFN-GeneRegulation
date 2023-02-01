%This script uses an ODE model to simulate changes in mRNA abundance as a
% function of mRNA synthesis activated by the ISGF3 transcription factor
% with an activation constant (KA), Hill coefficient (n), and an mRNA
% degradation rate constant (kdeg). The parameter values are optimized
% using a local optimizer with multiple starting seeds and
% parameter sensitivity is assessed.


% add path for data
dir_data = './Data';
addpath(dir_data)

%Load experimental ISGF3 activity and mRNA data
load([dir_data,'/BetaISGF3.mat']);
load([dir_data,'/LambdaISGF3.mat']);
load([dir_data,'/BetaTime.mat']);
load([dir_data,'/BetaCommonGenes.mat']);
load([dir_data,'/LambdaCommonGenes.mat']);
load([dir_data,'/GeneTime.mat']);
load([dir_data,'/BetaBetaGenes.mat']);
load([dir_data,'/LambdaBetaGenes.mat']);

%Use these two line commands if running script for IFN-beta specific genes
%BetaCommonGenes=BetaBetaGenes;
%LambdaCommonGenes=LambdaBetaGenes;

%Defined parameters
%params(1)= kact, maximal expression level of the promoter
%params(2)= Ka, activation coeff.(concentration at half-maximal expression)
%params(3)= kdeg, mRNA degradation rate
%params(4)= n, Hill coeff. (governs steepness of input function)

%Use this set of parameters for the IFN common genes
%the best fit parameters from 5 sets of 50 constrained optimization for the
% common gene cluster
params=[0.0015,6.94,0.0061,1.04]; 

%Use this set of parameters for the IFN-beta specific genes
%the best fit parameters from 5 sets of 50 constrained optimization for the
% IFNbeta gene cluster
%params=[0.0015,6.43,0.0094,2.16]; 

%% Normalize ISGF3 and RNA experimental data

totalISGF3=1; %max nuclear ISGF3
maxPercentage=1;
minPercentage=0.0025;

%scale EMSA proportional to amount of ISGF3 in nucleus (0.25% basal, %max)
allEMSA=[BetaISGF3,LambdaISGF3];
minEMSA=min(allEMSA,[],'all'); 
allEMSA=allEMSA-minEMSA;
maxEMSA=max(allEMSA,[],'all');
EMSAScaled=(allEMSA./maxEMSA)*(maxPercentage*totalISGF3);
EMSAScaled=EMSAScaled+(minPercentage*totalISGF3);

EMSABetaScaled=EMSAScaled(:,1);
EMSALambdaScaled=EMSAScaled(:,2);

%Normalize RNA    
allCommonGenes_data=[BetaCommonGenes,LambdaCommonGenes];
minGene_data=min(allCommonGenes_data,[],'all');
allCommonGenes_data=allCommonGenes_data-minGene_data;
maxGene_data=max(allCommonGenes_data,[],'all'); 

BetaRNA=allCommonGenes_data(:,1)./maxGene_data;
BetaRNA(BetaRNA<0)=0;

LambdaRNA=allCommonGenes_data(:,2)./maxGene_data;
LambdaRNA(LambdaRNA<0)=0;

%Calculating basal RNA concentration by taking average of Beta and IFN
% Lambda basal conditions
avgBasalRNA=mean([BetaRNA(1),LambdaRNA(1)]);

%Calculating basal ISGF3 concentration by taking average of IFN Beta and
% Lambda basal conditions
avgBasalISGF3=mean([EMSABetaScaled(1),EMSALambdaScaled(1)]);

%% Interpolate and plot ISGF3 data

%Interpolation of IFNBeta-induced ISGF3 data
BetamakimaFit=interp1(BetaTime,EMSABetaScaled,[0:800],'makima');

%Interpolation of IFNLambda-induced ISGF3 data
LambdamakimaFit=interp1(BetaTime,EMSALambdaScaled,[0:800],'makima');

%Plot figures of experimental IFNbeta- and IFNlambda-induced ISGF3 data 
% with Modified Akima Cubic fit
figure

plot([0:800],BetamakimaFit,'b-',BetaTime,EMSABetaScaled,'bo',[0:800],...
        LambdamakimaFit,'r:',BetaTime,EMSALambdaScaled,'rx','LineWidth',5)

legend({'Beta-ModifiedAkimaCubicFit','BetaData',...
        'Lambda-ModifiedAkimaCubicFit','LambdaData'},...
        'FontSize',18,'FontWeight','bold')
 
xlabel ('Time(minutes)','FontSize',18,'FontWeight','bold')
ylabel('Scaled','FontSize',18,'FontWeight','bold')
title('Active ISGF3 Input Curves','FontSize',18,'FontWeight','bold');
ax=gca;
ax.YLim=[0,1.2];


%% Run model to steady state to get initial values

% add path for scripts
dir_GeneReg = './GeneReg';
addpath(dir_GeneReg)

x=avgBasalRNA;
time=[0:10:84000];
ISGF3=avgBasalISGF3;

[t_ss,y_ss]=ode15s(@(t,x) GeneSteadyState(t,x,ISGF3,params),time,x);


 %Create figure with plots for the steady state model.
    figure
    plot(t_ss,y_ss,'k','LineWidth',1);
    legend({'RNA SS model'},'Fontsize',18,'FontWeight','bold',...
        'Location','southeast')
    xlabel ('Time(minutes)','FontSize',18,'FontWeight','bold')
    ylabel('RNA','FontSize',18,'FontWeight','bold')
    ax=gca;
    autoY=get(gca,'YLim');
    ax.YLim=[0,autoY(2)];
    hold on

%% Run model with initial values from steady state model and plot

SSInitial=y_ss(end);
t=[1:1:2500]; %time


[t_b,y_b]=ode15s(@(t,x) ISGF3GeneReg(t,x,BetamakimaFit,params),...
            t,SSInitial);

[t_l,y_l]=ode15s(@(t,x) ISGF3GeneReg(t,x,LambdamakimaFit,params),...
            t,SSInitial);

%scale simulated RNA data in both conditions
allRNA_sim=[y_b,y_l];
minRNA=min(allRNA_sim,[],'all');
allRNA_sim=allRNA_sim-minRNA;
maxRNA_sim=max(allRNA_sim,[],'all');

%Calculating percent max of model RNA output
normmRNA_b=allRNA_sim(:,1)./maxRNA_sim;

normmRNA_l=allRNA_sim(:,2)./maxRNA_sim;

%Calculate RMSD to evaluate curve fit
RMSD0=mRNACostFunction_ParamSensitivity(allRNA_sim,maxRNA_sim,...
        GeneTime,BetaRNA,LambdaRNA);

  figure
    subplot(2,1,1)
        plot(GeneTime,BetaRNA,'bx',t_b,normmRNA_b,'-b',...
            'LineWidth',2,'MarkerSize',12)
        sgtitle('RNA Plots with Initial Parameters',...
            'FontSize',18,'FontWeight','bold')
        xlabel('Time(minutes)','FontSize',18,'FontWeight','bold')
        ylabel('Percent Max','FontSize',18,'FontWeight','bold')
        legend({'BetamRNA-data','BetamRNA-model'},'FontSize',18,...
            'FontWeight','bold')
        ax=gca;
        autoY=get(gca,'YLim');
        ax.YLim=[0,1];     
        ax.XLim=[0,1560];
        ax.XTick=[0:120:1560];
        hold on
        
    subplot(2,1,2)
       plot(GeneTime,LambdaRNA,'ro',t_l, normmRNA_l,'-r',...
           'LineWidth',2,'MarkerSize',12)
       xlabel('Time(minutes)','FontSize',18,'FontWeight','bold')
       ylabel('Percent Max','FontSize',18,'FontWeight','bold')
       legend({'LambdamRNA-data','LambdamRNA-model'},...
           'Location','southeast','FontSize',18,'FontWeight','bold')
       ax=gca;
       autoY=get(gca,'YLim');
       ax.YLim=[0,1];    
       ax.XLim=[0,1560];
       ax.XTick=[0:120:1560];
       hold on
       
%% Parameter sensititivity scan

%Run the following function for a parameter sensitivity scan
evaluateParamSensitivity(params,GeneTime,BetaRNA,LambdaRNA,t,...
                            BetamakimaFit,LambdamakimaFit,...
                            avgBasalRNA,avgBasalISGF3)

%Calculate the RMSD to evaluate curve fit
RMSD=mRNACostFunction_ParamSensitivity(allRNA_sim,maxRNA_sim,...
                                        GeneTime,BetaRNA,LambdaRNA);
  
%% Optimization with multiple random seeds for initial parameter values

fprintf('Starting optimization\n')

% add path for scripts
dir_optimization = './Optimization';
addpath(dir_optimization)

%Running a optimization using randomized starting parameter 
% values at the beginning of n number of loops
iterNum=50;
savedParams=zeros(iterNum,length(params));
savedMinRMSD=zeros(iterNum,1);
kact=params(1);

%Optimization will loop with randomized initial parameter values except
%kact will remain constant since it's not sensitive
for loop=1:iterNum
    
    % constrain the search range for Ka between 0.1 and 10
    randParams1 = (rand(1)*(10-0.1)) + 0.1; 

    % constrain the search range for kdeg between 0.0006931 and 0.6931
    randParams2 = (rand(1)*(0.6931-0.0006931)) + 0.0006931;
    
    % constrain the search range for n between 0 and 10
    randParams3 = (rand(1)*(10-0)) + 0; 

    paramsLoop=[kact,randParams1,randParams2,randParams3];

    [newParams,minRMSD]=RNALocalOptim(paramsLoop,BetamakimaFit,...
                            LambdamakimaFit,GeneTime,BetaRNA,LambdaRNA,...
                            avgBasalISGF3,avgBasalRNA);
  
    fprintf('Iteration: ')
    disp(loop)

    savedParams(loop,:)=newParams;
    savedMinRMSD(loop,:)=minRMSD;
    
end
%% Run model with best fit parameter set from optimization and plot

save('newParams.mat','savedParams') 
save('RMSD.mat','savedMinRMSD') 


%Find the minimum RMSD value from the RMSD loop matrix
[optRMSD,Index]=min(savedMinRMSD);
fprintf('RMSD_Optimal:%1$.4f\n',optRMSD);

%Find the parameter set that produce the minimum RMSD
optParams=savedParams(Index,:);

%Plot histogram of all of the RMSD values
figure
histogram(savedMinRMSD(:,1),[-2 -1:0.1:1 2])
xlabel ('RMSD Value','FontSize',18,'FontWeight','bold')
ylabel('Number of Optimizations','FontSize',18,'FontWeight','bold')

x=avgBasalRNA;
time=[0:10:84000];
ISGF3=avgBasalISGF3;

%Run steady state model with optimized parameter set
[t_ss_new,y_ss_new]=ode15s(@(t,x) GeneSteadyState(t,x,ISGF3,optParams), ...
                                time,x);
                                                    

 %Create figure with plots for the steady state model.
    figure
    plot(t_ss_new,y_ss_new,'k','LineWidth',1);
    legend({'RNA SS'},'Fontsize',18,'FontWeight','bold', ...
                'Location','southeast')
    xlabel ('Time(minutes)','FontSize',18,'FontWeight','bold')
    ylabel('RNA','FontSize',18,'FontWeight','bold')
    ax=gca;
    autoY=get(gca,'YLim');
    ax.YLim=[0,autoY(2)];
    hold on

SSInitial_new=y_ss_new(end);   

%Run ODE with best fit parameter set
[t_beta_new,y_beta_new]=ode15s(@(t,x) ISGF3GeneReg(t,x, ...
                                    BetamakimaFit,optParams),t, ...
                                    SSInitial_new);

[t_lambda_new,y_lambda_new]=ode15s(@(t,x) ISGF3GeneReg(t,x, ...
                                        LambdamakimaFit,optParams),t, ...
                                        SSInitial_new);

%Scale model proportional to amount of RNA
allRNA_sim=[y_beta_new,y_lambda_new];
minRNA=min(allRNA_sim,[],'all');
allRNA_sim=allRNA_sim-minRNA;
maxRNA_sim=max(allRNA_sim,[],'all');

%Calculating percent max of model RNA output
normmRNA_beta_new=allRNA_sim(:,1)./maxRNA_sim;

normmRNA_lambda_new=allRNA_sim(:,2)./maxRNA_sim;

  figure
    subplot(2,2,1)
        plot(GeneTime,BetaRNA,'bx',t_beta_new,normmRNA_beta_new, ...
                ':b','LineWidth',2.5,'MarkerSize',12)
        title('Plots with Optimized Parameters','FontSize',18, ...
                'FontWeight','bold')
        xlabel('Time(minutes)','FontSize',18,'FontWeight','bold')
        ylabel('Percent Max','FontSize',18,'FontWeight','bold')
        legend({'BetamRNA-data','BetamRNA-model'},'FontSize',18, ...
                'FontWeight','bold')
        ax=gca;
        %autoY=get(gca,'YLim');
        ax.YLim=[0,1.2]; 
        hold on

    subplot(2,2,2)
       plot(GeneTime,LambdaRNA,'ro',t_lambda_new, normmRNA_lambda_new, ...
               ':r','LineWidth',2.5,'MarkerSize',12)
       xlabel('Time(minutes)','FontSize',18,'FontWeight','bold')
       ylabel('Percent Max','FontSize',18,'FontWeight','bold')
       legend({'LambdamRNA-data','LambdamRNA-model'}, ...
               'Location','southeast','FontSize',18,'FontWeight','bold')
       ax=gca;
       %autoY=get(gca,'YLim');
       ax.YLim=[0,1.2]; 
       hold on

     subplot(2,2,3)
        plot(t_beta_new,y_beta_new,':b','LineWidth',2.5)
        title('IFNBeta Model-not normalized','FontSize',12, ...
                'FontWeight','bold')
        xlabel('Time(minutes)','FontSize',18,'FontWeight','bold')
        ylabel('RNA (A.U.)','FontSize',18,'FontWeight','bold')
        ax=gca;
        autoY=get(gca,'YLim');
        ax.YLim=autoY;
        hold on

    subplot(2,2,4)
        plot(t_lambda_new,y_lambda_new,':r','LineWidth',2.5)    
        title('IFNLambda Model-not normalized','FontSize',12, ...
                'FontWeight','bold')
        xlabel('Time(minutes)','FontSize',18,'FontWeight','bold')
        ylabel('RNA (A.U.)','FontSize',18,'FontWeight','bold')
        ax=gca;
        ax.YLim=autoY;
        hold on
