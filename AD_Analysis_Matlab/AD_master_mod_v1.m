% Main Analysis Script of neuronal Firing Rates by genotype,neuron type and treatment
%with Levetiracetam Arto Lipponen Jan Klee 15.11.17

clear all

%% INPUT PARAMETERS
%DataFolder
DataFolder='T:\arto\data\2016-0098-003\AD_Lev Experiment_final\analysis_Jan\';

%Sessions to include: Add newly clustered sessions here in no particular
%order
sessions=  {'32397_Lev_lM2-1000_2017-10-04_12-33-44',...
            '32397_saline_lM2-1000_2017-10-05_11-21-04',...
            '32412_saline_lPFC2000down+_2017-10-03_16-37-59',...
            '32412_saline_lPFC-1000_2017-10-05_14-32-14',...
            '32412_saline_M2_1300down_2017-10-03_16-54-10'
            }

%list of transgenicAnimals
transgenicAnimals={'32397';'32412';'32423';'38606'};
batch1=    {'32397';'32412';'32398';'38624'};   

%FiringRates 
bins=60*30000;  %bins in seconds for which firing rates are calculated, default is 1m
binsInc=1:5;  %which bins to inlcude in analysis, must be smaller than the length of the shortes session
spikes=100; %number of spikes to inlcude to form waveform avarages for waveform analysis

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
    
    cd([DataFolder,sessions{s}])

    % extract Firing Rates with my FiringRates function
    [FRperiods,FRmeansS,FRstdS]=AD_FiringRates(bins);
    FRall5m=cat(1,FRall5m,FRperiods(:,binsInc));
    FRmeans=cat(2,FRmeans,FRmeansS);
    FRstd=cat(2,FRstd,FRstdS);
      
    % extract Waveforms and Waveform params with my waveFormAnlysis function
    [wfsS,includeWFS,SpikeV2PS,spikeWidthS,DepthS]=AD_waveFormAnalysis(spikes,0);
    wfs=cat(1,wfs,wfsS);
    includeWF=cat(2,includeWF,includeWFS);
    SpikeV2P=cat(1,SpikeV2P,SpikeV2PS);
    spikeWidth=cat(2,spikeWidth,spikeWidthS);
    Depth=cat(2,Depth,DepthS)

    % checks and saves if session was drug (1) or not (0)
    drug(cellcounter:cellcounter+length(FRmeansS)-1)= contains(sessions{s},'Lev');
    % checks and saves if session was from a transgenic (1) or not (0)
    genotype(cellcounter:cellcounter+length(FRmeansS)-1)=contains(sessions{1},transgenicAnimals);
    % keeps count of batch
    batch(cellcounter:cellcounter+length(FRmeansS)-1)=contains(sessions{1},batch1);
    % keeps count of animal ID
    animals(cellcounter:cellcounter+length(FRmeansS)-1)=str2num(sessions{1}(1:5));
    
    %area optioal
    %area=(cellcounter:cellcounter+length(FRmeansS)-1)= contains(sessions{s},'M2');
    
    % change cell counter every iteration
    cellcounter=cellcounter+length(FRmeansS);
    
end

%% Waveform Analysis on all extracted waveforms 
% determine NeuronType by cluster analysis of all cells
good=find(includeWF==1);
X=horzcat(spikeWidth(good)',SpikeV2P(good));
[idxInt,C] = kmeans(X,2);
neuronType(1:length(spikeWidth))=0;
neuronType(good(idxInt==1))=1;
neuronType(good(idxInt==2))=2;

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
scatter(spikeWidth(good(idxInt==1)),SpikeV2P(good(idxInt==1)),'r')
hold on
scatter(spikeWidth(good(idxInt==2)),SpikeV2P(good(idxInt==2)),'b')

%% Save
cd('T:\arto\data\2016-0098-003\AD_Lev Experiment_final\analysis_Jan')
save AD_FiringRates     

%% analysis scatch
% stats 

clear all
cd('T:\arto\data\2016-0098-003\AD_Lev Experiment_final\analysis_Jan')
load AD_FiringRates

% I think 3-way anova (genotype,drug,neuronType) is the way to go here
% just ttest would not be enough because no correction for multiple
% comparisons
% if significant you can go on with individual ttests
[p,tbl,stats] = anovan(FRmeans,{genotype,drug,neuronType});

%ttests
%genotype
[hGT,pGT]=ttest2(FRmeans(genotype==1),FRmeans(genotype==0));
%Drug
[hDrug,pDrug]=ttest2(FRmeans(drug==1),FRmeans(drug==0));
%Neurontype
[hneurontype,pneurontype]=ttest2(FRmeans(neuronType==1),FRmeans(neuronType==02));

%WT drug vs no drug
[hWTd,pWTd]=ttest2(FRmeans(genotype==0&&drug==1),FRmeans(genotype==0&&drug==0));
%TG drug vs no drug
[hTGd,pTGd]=ttest2(FRmeans(genotype==1&&drug==1),FRmeans(genotype==1&&drug==0));

%WT Int vs Pyr
[hWTip,pWTip]=ttest2(FRmeans(genotype==0&&neuronType==1),FRmeans(genotype==0&&dneuronType==2));
%TG Int vs Pyr
[hTGip,pTGip]=ttest2(FRmeans(genotype==1&&neuronType==1),FRmeans(genotype==1&&dneuronType==2));

%Drug Int vs Pyr
[hDip,pDip]=ttest2(FRmeans(drug==1&&neuronType==1),FRmeans(drug==1&&neuronType==2));

%Saline Int vs Pyr
[hSip,pSip]=ttest2(FRmeans(drug==0&&neuronType==1),FRmeans(drug==0&&neuronType==2));

%% plots

%names={'Interneurons','Pyramidal Cells'};
names={'levetiracetam','vehicle'
ylab={'Firing Rates Hz'};

ttl={'Firing Rates for Interneurons in control and drug groups in Frontal Cortex of mice'};
%ttl={'Firing Rates for Different Neuron Types in Frontal Cortex of mice'};

data{1}=FRmeans(neuronType==1&drug==1);
data{2}=FRmeans(neuronType==1&drug==0);

AD_BeehivePlot(data,names,ylab,ttl)







