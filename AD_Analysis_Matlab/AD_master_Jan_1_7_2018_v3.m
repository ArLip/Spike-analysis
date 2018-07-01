% Main Analysis Script of neuronal Firing Rates by genotype,neuron type and treatment
%with Levetiracetam Arto Lipponen Jan Klee 15.11.17

clear all

%% INPUT PARAMETERS
%DataFolder
DataFolder='E:\analysis_Jan\';

%Sessions to include: Add newly clustered sessions here in no particular
%order
sessions=  {'32397_Lev_lM2-1000_2017-10-04_12-33-44',...
            '32397_saline_lM2-1000_2017-10-05_11-21-04',...            
            '32398_saline_M2_1000down_2017-10-03_15-37-39',...
            '32398_saline_M2_1000down_2017-10-03_15-37-39',...
            '32412_saline_lPFC-1000_2017-10-05_14-32-14',...
            '32412_saline_lPFC2000down+_2017-10-03_16-37-59',...
            '32412_saline_M2_1300down_2017-10-03_16-54-10',...
            '32423_Lev_lM2-1000_2017-10-10_12-18-55',...
            '32423_Lev_lPFC-2000_2017-10-09_12-44-36',...
            '32423_Saline_lM2-1000_2017-10-08_11-01-42',...
            '32424_Lev_lM2-1000_2017-10-09_16-01-30',...
            '32424_Saline_lPFC-2000_2017-10-08_14-18-26',...
            '32424_Lev_lM2-1000_2017-10-07_15-53-18',...
            '32424_Lev_rM2-1000_2017-10-10_15-53-53',...
            '32425_Lev_lM2-1000_2017-10-07_15-07-58',...
            '32425_Lev_lM2-1000_2017-10-09_15-28-05',...
            '32425_Lev_lM2-1000_2017-10-10_14-41-51',...
            '32425_Saline_lPFC-1000_2017-10-08_13-39-37',...
            '38606_Lev_lM-1000_2017-10-09_13-32-48',...
            '38606_Saline_lM2-1000_2017-10-08_12-04-28'            
            }
            
 %'32425_Lev_lM2-1000_2017-10-10_14-41-51',...
 %'32424_Saline_lPFC-2000_2017-10-08_14-18-26',...
 %'32424_Lev_lM2-1000_2017-10-07_15-53-18',...
 %'32424_Lev_rM2-1000_2017-10-10_15-53-53',...
%list of transgenicAnimals
transgenicAnimals={'32397';'32412';'32423';'38606'};
batch1=    {'32397';'32412';'32398';'38624'};   

%FiringRates 
bins=60*30000;  %bins in seconds for which firing rates are calculated, default is 1m
binsInc=1:5;  %which bins to inlcude in analysis, must be smaller than the length of the shortes session
spikes=2000; %number of spikes to inlcude to form waveform avarages for waveform analysis

%% main data collection loop
FRall5m=[];
FRmeans=[];
FRstd=[];
wfs=[];
includeWF=[];
SpikeV2P=[];
spikeWidth=[];
Depth=[];

cellcounter=1;
for s=1:size(sessions,2)
    s
    cd([DataFolder,sessions{s}])
    % extract Firing Rates with my FiringRates function
    [FRperiods,FRmeansS,FRstdS]=AD_FiringRates(bins);
    FRall5m=cat(1,FRall5m,FRperiods(:,binsInc));
    FRmeans=cat(2,FRmeans,FRmeansS);
    FRstd=cat(2,FRstd,FRstdS);
      
    % extract Waveforms and Waveform params with my waveFormAnlysis function
    [wfsS,includeWFS,SpikeV2PS,spikeWidthS,DepthS]=AD_waveFormAnalysis(spikes,0); %AD_master_Jan_1_7_2018_v2 [wfsS,includeWFS,SpikeV2PS,spikeWidthS,DepthS,wfsRawS]=AD_waveFormAnalysis(spikes,0); 
    wfs=cat(1,wfs,wfsS);
    includeWF=cat(2,includeWF,includeWFS);
    SpikeV2P=cat(1,SpikeV2P,SpikeV2PS);
    spikeWidth=cat(2,spikeWidth,spikeWidthS);
    Depth=cat(2,Depth,DepthS)  %AD_master_Jan_1_7_2018_v2  wfsRaw=cat(1,wfsRaw,wfsRawS);

    % checks and saves if session was drug (1) or not (0)
    drug(cellcounter:cellcounter+length(FRmeansS)-1)= contains(sessions{s},'Lev');
    % checks and saves if session was from a transgenic (1) or not (0)
    genotype(cellcounter:cellcounter+length(FRmeansS)-1)=contains(sessions{s},transgenicAnimals);
    % keeps count of batch
    batch(cellcounter:cellcounter+length(FRmeansS)-1)=contains(sessions{s},batch1);
    % keeps count of animal ID
    animals(cellcounter:cellcounter+length(FRmeansS)-1)=str2num(sessions{s}(1:5));
    
    %area optioal
    %area=(cellcounter:cellcounter+length(FRmeansS)-1)= contains(sessions{s},'M2');
    
    % change cell counter every iteration
    cellcounter=cellcounter+length(FRmeansS);
    
end

%% Waveform Analysis on all extracted waveforms 
% determine NeuronType by cluster analysis of all cells
good=find(includeWF==1);
% X=horzcat(spikeWidth(good)',SpikeV2P(good));
% [idxInt,C] = kmeans(X,2);
% neuronType(1:length(spikeWidth))=0;
% neuronType(good(idxInt==1))=1;
% neuronType(good(idxInt==2))=2;


neuronType(1:length(spikeWidth))=0;
neuronType(intersect(good,find(SpikeV2P<19)))=1;
neuronType(intersect(good,find(SpikeV2P>19)))=2;
neuronType(intersect(good,find(spikeWidth>25)))=3;


figure()
subplot(1,2,1)
for i=1:length(neuronType)
    if neuronType(i)==1
    plot(wfs(i,:),'r')
    hold on
    elseif neuronType(i)==2
    plot(wfs(i,:),'b')
    hold on
    else
        continue
    end
end

subplot(1,2,2)
scatter(spikeWidth(neuronType==1),SpikeV2P(neuronType==1),'r')
hold on
scatter(spikeWidth(neuronType==2),SpikeV2P(neuronType==2),'b')


% subplot(1,2,2)
% scatter(spikeWidth(good(idxInt==1)),SpikeV2P(good(idxInt==1)),'r')
% hold on
% scatter(spikeWidth(good(idxInt==2)),SpikeV2P(good(idxInt==2)),'b')
%set manual threshold
% manualTresh=find(SpikeV2P(good)>19);
% idxInt(manualTresh)=1;

%% Save
cd('E:\analysis_Jan\')
save AD_FiringRates_010718_v1     

%%% analysis scatch
% stats 

clear all
cd('Y:\data\2016-0098-003\AD_Lev Experiment_final\analysis_Jan')
load AD_FiringRates_002

% I think 3-way anova (genotype,drug,neuronType) is the way to go here
% just ttest would not be enough because no correction for multiple
% comparisons
% if significant you can go on with individual ttests
[p,tbl,stats] = anovan(FRmeans,{genotype,drug,neuronType},'model','interaction');

%%ttests
%genotype
[hGT,pGT]=ttest2(FRmeans(genotype==1),FRmeans(genotype==0))
%Drug
[hDrug,pDrug]=ttest2(FRmeans(drug==1),FRmeans(drug==0));
%Neurontype
[hneurontype,pneurontype]=ttest2(FRmeans(neuronType==1),FRmeans(neuronType==2));

%WT drug vs no drug
[hWTd,pWTd]=ttest2(FRmeans(genotype==0&drug==1),FRmeans(genotype==0&drug==0))
%TG drug vs no drug
[hTGd,pTGd]=ttest2(FRmeans(genotype==1&drug==1),FRmeans(genotype==1&drug==0))

%WT Int vs Pyr
[hWTip,pWTip]=ttest2(FRmeans(genotype==0&neuronType==1),FRmeans(genotype==0&neuronType==2));
%TG Int vs Pyr
[hTGip,pTGip]=ttest2(FRmeans(genotype==1&neuronType==1),FRmeans(genotype==1&neuronType==2));

%Drug Int vs Pyr
[hDip,pDip]=ttest2(FRmeans(drug==1&neuronType==1),FRmeans(drug==1&neuronType==2));

%Saline Int vs Pyr
[hSip,pSip]=ttest2(FRmeans(drug==0&neuronType==1),FRmeans(drug==0&neuronType==2));

%Drug Int vs Sal Int in wt
[hDip,pDip]=ttest2(FRmeans(drug==1&genotype==0&neuronType==2),FRmeans(drug==0&genotype==0&neuronType==2));

%Wt + leve + pyr/int 
[hDip,pDip]=ttest2(FRmeans(drug==1&genotype==0&neuronType==1),FRmeans(drug==1&genotype==0&neuronType==2));

%Wt + saline + pyr/int 
[hDip,pDip]=ttest2(FRmeans(drug==0&genotype==0&neuronType==1),FRmeans(drug==0&genotype==0&neuronType==2));

%TG + leve + pyr/int 
[hDip,pDip]=ttest2(FRmeans(drug==1&genotype==1&neuronType==1),FRmeans(drug==1&genotype==1&neuronType==2));


%TG + int + leve/sal 
[hDip,pDip]=ttest2(FRmeans(drug==1&genotype==1&neuronType==1),FRmeans(drug==0&genotype==1&neuronType==1));

%TG + pyr + leve/sal 
[hDip,pDip]=ttest2(FRmeans(drug==1&genotype==1&neuronType==2),FRmeans(drug==0&genotype==1&neuronType==2));

%wt + int + leve/sal 
[hDip,pDip]=ttest2(FRmeans(drug==1&genotype==0&neuronType==1),FRmeans(drug==0&genotype==0&neuronType==1));

%wt + pyr + leve/sal 
[hDip,pDip]=ttest2(FRmeans(drug==1&genotype==0&neuronType==2),FRmeans(drug==0&genotype==0&neuronType==2));

%TG + saline + pyr/int 
[hDip,pDip]=ttest2(FRmeans(drug==0&genotype==1&neuronType==1),FRmeans(drug==0&genotype==1&neuronType==2));


%Pyr + tg/wt
[hDip,pDip]=ttest2(FRmeans(genotype==1&neuronType==2),FRmeans(genotype==0&neuronType==2));

%Int + tg/wt
[hDip,pDip]=ttest2(FRmeans(genotype==1&neuronType==1),FRmeans(genotype==0&neuronType==1));

%Genotype saline
[hDip,pDip]=ttest2(FRmeans(drug==0&genotype==1&neuronType==2),FRmeans(drug==0&genotype==0&neuronType==2));

%% rank sum tests

%genotype
[hGT,pGT]=ranksum(FRmeans(genotype==1),FRmeans(genotype==0));
%Drug
[hDrug,pDrug]=ranksum(FRmeans(drug==1),FRmeans(drug==0));
%Neurontype
[hneurontype,pneurontype]=ranksum(FRmeans(neuronType==1),FRmeans(neuronType==2));

%WT drug vs no drug
[hWTd,pWTd]=ranksum(FRmeans(genotype==0&drug==1),FRmeans(genotype==0&drug==0));
%TG drug vs no drug
[hTGd,pTGd]=ranksum(FRmeans(genotype==1&drug==1),FRmeans(genotype==1&drug==0));

%WT Int vs Pyr
[hWTip,pWTip]=ranksum(FRmeans(genotype==0&neuronType==1),FRmeans(genotype==0&neuronType==2));
%TG Int vs Pyr
[hTGip,pTGip]=ranksum(FRmeans(genotype==1&neuronType==1),FRmeans(genotype==1&neuronType==2));

%Drug Int vs Pyr
[hDip,pDip]=ranksum(FRmeans(drug==1&neuronType==1),FRmeans(drug==1&neuronType==2));

%Saline Int vs Pyr
[hSip,pSip]=ranksum(FRmeans(drug==0&neuronType==1),FRmeans(drug==0&neuronType==2));

%Drug Int vs Sal Int in wt
[hDip,pDip]=ranksum(FRmeans(drug==1&genotype==0&neuronType==2),FRmeans(drug==0&genotype==0&neuronType==2));

%Wt + leve + pyr/int 
[hDip,pDip]=ranksum(FRmeans(drug==1&genotype==0&neuronType==1),FRmeans(drug==1&genotype==0&neuronType==2));

%Wt + saline + pyr/int 
[hDip,pDip]=ranksum(FRmeans(drug==0&genotype==0&neuronType==1),FRmeans(drug==0&genotype==0&neuronType==2));

%TG + leve + pyr/int 
[hDip,pDip]=ranksum(FRmeans(drug==1&genotype==1&neuronType==1),FRmeans(drug==1&genotype==1&neuronType==2));


%TG + int + leve/sal 
[hDip,pDip]=ranksum(FRmeans(drug==1&genotype==1&neuronType==1),FRmeans(drug==0&genotype==1&neuronType==1));

%TG + pyr + leve/sal 
[hDip,pDip]=ranksum(FRmeans(drug==1&genotype==1&neuronType==2),FRmeans(drug==0&genotype==1&neuronType==2));

%wt + int + leve/sal 
[hDip,pDip]=ranksum(FRmeans(drug==1&genotype==0&neuronType==1),FRmeans(drug==0&genotype==0&neuronType==1));

%wt + pyr + leve/sal 
[hDip,pDip]=ranksum(FRmeans(drug==1&genotype==0&neuronType==2),FRmeans(drug==0&genotype==0&neuronType==2));

%TG + saline + pyr/int 
[hDip,pDip]=ranksum(FRmeans(drug==0&genotype==1&neuronType==1),FRmeans(drug==0&genotype==1&neuronType==2));


%Pyr + tg/wt
[hDip,pDip]=ranksum(FRmeans(genotype==1&neuronType==2),FRmeans(genotype==0&neuronType==2));

%Int + tg/wt
[hDip,pDip]=ranksum(FRmeans(genotype==1&neuronType==1),FRmeans(genotype==0&neuronType==1));

%Genotype saline
[hDip,pDip]=ranksum(FRmeans(drug==0&genotype==1&neuronType==2),FRmeans(drug==0&genotype==0&neuronType==2));




%% plots

%tg vs wt

names={'tg','wt'};
ylab={'Firing Rates Hz'};
ttl={'Firing Rates for tg / wt'};

data{1}=FRmeans(genotype==1);
data{2}=FRmeans(genotype==0);

AD_BeehivePlot(data,names,ylab,ttl)

%leve vs sal

names={'leve','sal'};
ylab={'Firing Rates Hz'};
ttl={'Firing Rates for leve vs sal'};

data{1}=FRmeans(genotype==1);
data{2}=FRmeans(genotype==0);

AD_BeehivePlot(data,names,ylab,ttl)


%pyr vs int

names={'pyr','int'};
ylab={'Firing Rates Hz'};
ttl={'Firing Rates for pyr - int'};

data{1}=FRmeans(neuronType==2);
data{2}=FRmeans(genotype==1);

AD_BeehivePlot(data,names,ylab,ttl)

% pyr - tg vs wt
names={'pyr + tg','pyr + wt'};
ylab={'Firing Rates Hz'};
ttl={'Firing Rates for pyr - tg vs wt'};

data{1}=FRmeans(neuronType==2&genotype==1);
data{2}=FRmeans(neuronType==2&genotype==0);

AD_BeehivePlot(data,names,ylab,ttl)

% int - tg vs wt
names={'int + tg','int + wt'};
ylab={'Firing Rates Hz'};
ttl={'Firing Rates for int - tg vs wt'};

data{1}=FRmeans(neuronType==1&genotype==1);
data{2}=FRmeans(neuronType==1&genotype==0);

AD_BeehivePlot(data,names,ylab,ttl)

% tg - sal vs leve
names={'tg + saline ','tg + leve'};
ylab={'Firing Rates Hz'};
ttl={'Firing Rates for tg - sal vs leve'};

data{1}=FRmeans(genotype==1&drug==0);
data{2}=FRmeans(genotype==1&drug==1);

AD_BeehivePlot(data,names,ylab,ttl)

% wt - sal vs leve
names={'wt + sal ','wt + leve'};
ylab={'Firing Rates Hz'};
ttl={'Firing Rates for wt - sal vs leve'};

data{1}=FRmeans(genotype==0&drug==0);
data{2}=FRmeans(genotype==0&drug==1);

AD_BeehivePlot(data,names,ylab,ttl)

% wt - int vs pyr
names={'wt + int ','wt + pyr'};
ylab={'Firing Rates Hz'};
ttl={'Firing Rates for wt - int vs pyr'};

data{1}=FRmeans(genotype==0&neuronType==1);
data{2}=FRmeans(genotype==0&neuronType==2);

AD_BeehivePlot(data,names,ylab,ttl)

% tg - int/pyr
names={'tg + int ','tg + pyr'};
ylab={'Firing Rates Hz'};
ttl={'Firing Rates for tg - int/pyr'};

data{1}=FRmeans(genotype==1&neuronType==1);
data{2}=FRmeans(genotype==1&neuronType==2);

AD_BeehivePlot(data,names,ylab,ttl)

% wt - leve - int/pyr

names={'wt - leve - int','wt - leve - pyr'};
ylab={'Firing Rates Hz'};
ttl={'Firing Rates for wt - leve - int/pyr'};

data{1}=FRmeans(genotype==0&neuronType==1&drug==1);
data{2}=FRmeans(genotype==0&neuronType==2&drug==1);

AD_BeehivePlot(data,names,ylab,ttl)

% tg - leve - int/pyr

names={'tg - leve - int','tg - leve - pyr'};
ylab={'Firing Rates Hz'};
ttl={'Firing Rates for tg - leve - int/pyr'};

data{1}=FRmeans(genotype==1&neuronType==1&drug==1);
data{2}=FRmeans(genotype==1&neuronType==2&drug==1);

AD_BeehivePlot(data,names,ylab,ttl)

% wt - sal - int/pyr

names={'wt - sal - int','wt - sal - pyr'};
ylab={'Firing Rates Hz'};
ttl={'Firing Rates for wt - sal - int/pyr'};

data{1}=FRmeans(genotype==0&neuronType==1&drug==0);
data{2}=FRmeans(genotype==0&neuronType==2&drug==0);

AD_BeehivePlot(data,names,ylab,ttl)

% tg - sal - int/pyr

names={'tg - sal - int','tg - sal - pyr'};
ylab={'Firing Rates Hz'};
ttl={'Firing Rates for tg - sal - int/pyr'};

data{1}=FRmeans(genotype==1&neuronType==1&drug==0);
data{2}=FRmeans(genotype==1&neuronType==2&drug==0);

AD_BeehivePlot(data,names,ylab,ttl)

% tg - int - leve/sal

names={'tg - int - leve','tg - int - sal'};
ylab={'Firing Rates Hz'};
ttl={'Firing Rates for tg - int - leve/sal'};

data{1}=FRmeans(genotype==1&neuronType==1&drug==1);
data{2}=FRmeans(genotype==1&neuronType==1&drug==0);

AD_BeehivePlot(data,names,ylab,ttl)


% tg - pyr - leve/sal

names={'tg - pyr - leve','tg - pyr - sal'};
ylab={'Firing Rates Hz'};
ttl={'Firing Rates for tg - pyr - leve/sal'};

data{1}=FRmeans(genotype==1&neuronType==2&drug==1);
data{2}=FRmeans(genotype==1&neuronType==2&drug==0);

AD_BeehivePlot(data,names,ylab,ttl)

% wt - int - leve/sal

names={'wt - int - leve','wt - int - sal'};
ylab={'Firing Rates Hz'};
ttl={'Firing Rates for wt - int - leve/sal'};

data{1}=FRmeans(genotype==0&neuronType==1&drug==1);
data{2}=FRmeans(genotype==0&neuronType==1&drug==0);

AD_BeehivePlot(data,names,ylab,ttl)


% wt - pyr - leve/sal

names={'wt - pyr - leve','wt - pyr - sal'};
ylab={'Firing Rates Hz'};
ttl={'Firing Rates for wt - pyr - leve/sal'};

data{1}=FRmeans(genotype==0&neuronType==2&drug==1);
data{2}=FRmeans(genotype==0&neuronType==2&drug==0);

AD_BeehivePlot(data,names,ylab,ttl)


%y=find (idxInt==2)

