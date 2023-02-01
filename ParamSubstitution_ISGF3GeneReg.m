% This script runs the mRNA expression model and substitutes the parameters
% to determine which substitution improves the model performance 

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

%Defined parameters
%params(1)= kact, maximal expression level of the promoter
%params(2)= Ka, activation coeff.(concentration at half-maximal expression)
%params(3)= kdeg, mRNA degradation rate
%params(4)= n, Hill coeff. (governs steepness of input function)


%the best fit parameters from 5 sets of 50 constrained optimization for the
% IFNbeta gene cluster
params=[0.0015,6.43,0.0094,2.16]; 

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

%% Interpolate ISGF3 data

%Interpolation of IFNBeta-induced ISGF3 data
BetamakimaFit=interp1(BetaTime,EMSABetaScaled,[0:800],'makima');

%Interpolation of IFNLambda-induced ISGF3 data
LambdamakimaFit=interp1(BetaTime,EMSALambdaScaled,[0:800],'makima');

%% Run model to steady state to get initial values

% add path for scripts
dir_GeneReg = './GeneReg';
addpath(dir_GeneReg)

x=avgBasalRNA;
time=[0:10:84000];
ISGF3=avgBasalISGF3;

[t_ss,y_ss]=ode15s(@(t,x) GeneSteadyState(t,x,ISGF3,params),time,x);


%% Run model with initial parameter set from steady state model and plot

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
        sgtitle('RNA Plots',...
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

%% 
%Run the model with best fit parameters from the Beta gene cluster but
%replace parameters with the Common paramaters one by one to
%determine how sensitive they are

%the best fit parameters from 5 sets of 50 constrained optimization for
% the common gene cluster
paramsC=[0.0015,6.94,0.0061,1.04]; 
%the best fit parameters from 5 sets of 50 constrained optimization for
% the IFNbeta gene cluster
paramsB=[0.0015,6.43,0.0094,2.16]; 

RMSDstored=[];

% add path for colormap files
dir_cbrewer = './cbrewer'
addpath(dir_cbrewer)

c=cbrewer('seq','Greys',5);
c=flip(c); %flip color map so the darker shades will be selected

x=avgBasalRNA;
time=[0:10:84000];
ISGF3=avgBasalISGF3;


RMSDstored2=[];

for  a=1:4

    if a==1
        
        params2=[paramsB(1),paramsC(2),paramsB([3,4])]; % replace Ka
        name='Cluster 2 Ka';
        
    elseif a==2
        
        params2=[paramsB([1,2]),paramsC(3),paramsB(4)]; % replace kdeg
        name='Cluster 2 kdeg';
        
    elseif a==3
        
        params2=[paramsB([1:3]),paramsC(4)]; % replace n
        name='Cluster 2 n';
        
    else
        params2=[paramsB([1,2]),paramsC([3,4])]; % replace Ka and kdeg
        name='Cluster 2 n and kdeg';
                
    end

    [t_ss,y_ss]=ode15s(@(t,x) GeneSteadyState(t,x,ISGF3,params2),time,x);

    %Run model with input from steady state model
    SSInitial=y_ss(end);
    t=[1:1:2500]; %time


    [t_b,y_b]=ode15s(@(t,x) ISGF3GeneReg(t,x,BetamakimaFit,params2),...
                t,SSInitial);

    [t_l,y_l]=ode15s(@(t,x) ISGF3GeneReg(t,x,LambdamakimaFit,params2),...
                t,SSInitial);

    %scale model proportional to amount of RNA
    allRNA_sim=[y_b,y_l];
    minRNA=min(allRNA_sim,[],'all');
    allRNA_sim=allRNA_sim-minRNA;
    maxRNA_sim=max(allRNA_sim,[],'all');

    %Calculating percent max of model RNA output
    normmRNA_b=allRNA_sim(:,1)./maxRNA_sim;

    normmRNA_l=allRNA_sim(:,2)./maxRNA_sim;

    %Adding line graphs to previous plots
    subplot(2,1,1)
    plot(normmRNA_b,'--','Color',c(a,:),'LineWidth',2,...
            'DisplayName',name)
    hold on

    subplot(2,1,2)
    plot(t_l, normmRNA_l,'--','Color',c(a,:),'LineWidth',2,...
            'DisplayName',name)
    hold on


    %Calculate the RMSD to evaluate curve fit
    RMSD=mRNACostFunction_ParamSensitivity(allRNA_sim,maxRNA_sim,...
            GeneTime,BetaRNA,LambdaRNA);
    RMSDstored2=[RMSDstored2,RMSD];
end