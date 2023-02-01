% This script model plots ISGF3 dose-response curves for genes
% in cluster 1 (beta-specific) and 2 (common) using the Hill function


%params(1)= kact, maximal expression level of the promoter
%params(2)= Ka, activation coeff.(concentration at half-maximal expression)
%params(3)= kdeg, mRNA degradation rate
%params(4)= n, Hill coeff. (governs steepness of input function)

%make empty vector for storing data
storeB=[];
storeC=[];

%define the dose range for ISGF3
ISGF3Dose=[0:0.1:2.5];
ISGF3Dose=transpose(ISGF3Dose);

% add path for script
dir_doseresponse = './Dose Response';
addpath(dir_doseresponse)

% add path for colormap files
dir_cbrewer = './cbrewer';
addpath(dir_cbrewer)

%defines list of colors from a color map
c=cbrewer('seq','Reds',15);

%% Run model simulations at steady state at each dose

 for ISGF3=ISGF3Dose
     
        %the best fit parameters from 5 sets of 50 constrained 
        % optimization for the IFNbeta gene cluster
        paramsB=[0.0015,6.43,0.0094,2.16];

        %the best fit parameters from 5 sets of 50 constrained
        % optimization for the common gene cluster
        paramsC=[0.0015,6.94,0.0061,1.04]; 
                
        %Calculate RNA at steady state    
        mRNA_B=ISGF3GeneReg_SteadyState(ISGF3,paramsB);
        mRNA_C=ISGF3GeneReg_SteadyState(ISGF3,paramsC);
        
        %Store RNA value
        storeB=[storeB;mRNA_B];
        storeC=[storeC;mRNA_C];
 end
 
 normB=storeB./storeB(end);
 normC=storeC./storeC(end);
 
%% Make dose curve plot 

figure
 plot(ISGF3Dose,normB,'-k',ISGF3Dose,normC,':k','LineWidth',3.5)
        title('RNA Model dose response curve','FontSize',12, ...
                'FontWeight','bold')
        xlabel('ISGF3 Dose','FontSize',18,'FontWeight','bold')
        ylabel('mRNA Abundance(Percent Max of Each)','FontSize',18, ...
                'FontWeight','bold')
        legend({'Cluster1','Cluster2'},'Location','northeastoutside', ...
                'FontSize',18,'FontWeight','bold')
        ax=gca;
        ax.XTick = [0:0.5:1.0,1.4,1.5:0.5:max(ISGF3Dose)];
        hold on
 
             