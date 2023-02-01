function evaluateParamSensitivity(params,GeneTime,BetaRNA,LambdaRNA,t,BetamakimaFit,LambdamakimaFit,avgBasalRNA,avgBasalISGF3)

%Make an empty table for storing data
RMSDtable=[]; 

for a=2:length(params)
    
    %define an array of multipliers that I want to scan some parameter over
    multipliers=[0.25,0.35,0.5,0.70,1,1.41,2,2.82,4];

    %defines list of colors from a color map
    c=colormap(copper(length(multipliers)));
    
    %Making array of parameter names to use for figure titles
    paramName={'kact','Ka','kdeg','n'};
    name=strcat('Median mRNA Response of Common IFN Induced Genes - ', ...
                    paramName(a),' parameter scan');
        
    RMSDstored=[];
    
        for b=1:length(multipliers)
            %make a clean copy of updateParams unaffected by previous loops
            updateParams=params;
            
            updateParams(a)=params(a)*multipliers(b);

            %Run steady state model for new parameter set
            x=avgBasalRNA;
            time=[0:10:84000];
            ISGF3=avgBasalISGF3;
            [~,y_ss]=ode15s(@(t,x) GeneSteadyState(t,x, ...
                                ISGF3,updateParams),time,x);
            
            SSInitial=y_ss(end);
            
            %Using the Modified Akima cubic Hermite interpolation 
            % as the input function
            [t_beta,y_beta]=ode15s(@(t,x) ISGF3GeneReg(t,x,BetamakimaFit, ...
                                        updateParams),t,SSInitial);

            [t_lambda,y_lambda]=ode15s(@(t,x) ISGF3GeneReg(t,x, ...
                                            LambdamakimaFit,updateParams), ...
                                            t,SSInitial);
            
            %scale model proportional to amount of RNA
            allRNA_sim=[y_beta,y_lambda];
            minRNA=min(allRNA_sim,[],'all');
            allRNA_sim=allRNA_sim-minRNA;
            maxRNA_sim=max(allRNA_sim,[],'all');
            
            %Calculate the RMSD to evaluate curve fit
            RMSD=mRNACostFunction_ParamSensitivity(allRNA_sim, ...
                                maxRNA_sim,GeneTime,BetaRNA,LambdaRNA);
            
            RMSDstored=[RMSDstored,RMSD];         
                           
        end
              
        %Make a table for the RMSD of all the parameter scans for each
        %parameter in order to make a heatmap
        RMSDtable=[RMSDtable;RMSDstored];        
                
end


%Since reporting the Ka, which is the inverse of Kd, need to flip the table
%so that the multipliers match the correct term
 RMSDtable(1,:)=flip(RMSDtable(1,:),2);

figure
Mult = {'0.25x','0.35x','0.50x','0.70x','1x','1.41x','2x','2.82x','4x'};
h=heatmap(Mult,paramName(2:end),RMSDtable);

% add path for colormap files
dir_cbrewer = './cbrewer';
addpath(dir_cbrewer)

c=abs(cbrewer('seq','YlOrRd',50));
h.Colormap = c;

h.Title = 'RMSD - Parameter Sensitivity Analysis';
h.XLabel = 'Multiplier';
h.YLabel = 'Parameter';
