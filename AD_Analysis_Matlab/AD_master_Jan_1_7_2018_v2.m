%Main Analysis Script of neuronal Firing Rates by genotype,neuron type and treatment
%with Levetiracetam Arto Lipponen Jan Klee 15.11.17
%modified 1st July 2018

clear all

%% INPUT PARAMETERS
%DataFolder for Jan
%DataFolder= 'T:\arto\data\2016-0098-003\AD_Lev Experiment_final\analysis_Jan\';

%DataFolder for Arto
DataFolder='S:\rd-psy-TimeCell\Tompouce_8thJan18\data\2016-0098-003\AD_Lev Experiment_final\analysis_Jan';

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
wfsRaw=[];
%AutoCorrelograms=[]

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
    [wfsS,includeWFS,SpikeV2PS,spikeWidthS,DepthS,wfsRawS]=AD_waveFormAnalysis(spikes,0);
    wfs=cat(1,wfs,wfsS);
    includeWF=cat(2,includeWF,includeWFS);
    SpikeV2P=cat(1,SpikeV2P,SpikeV2PS);
    spikeWidth=cat(2,spikeWidth,spikeWidthS);
    wfsRaw=cat(1,wfsRaw,wfsRawS);

    % checks and saves if session was drug (1) or not (0)
    drug(cellcounter:cellcounter+length(FRmeansS)-1)= contains(sessions{s},'Lev');
    % checks and saves if session was from a transgenic (1) or not (0)
    genotype(cellcounter:cellcounter+length(FRmeansS)-1)=contains(sessions{s},transgenicAnimals);
    % keeps count of batch
    batch(cellcounter:cellcounter+length(FRmeansS)-1)=contains(sessions{s},batch1);
    % keeps count of animal ID
    animals(cellcounter:cellcounter+length(FRmeansS)-1)=str2num(sessions{s}(1:5));

     [AllCor]=AD_AutoCorrelograms;
     AutoCorrelograms(cellcounter:cellcounter+size(AllCor,1)-1,:)=AllCor;
      cellcounter=cellcounter+size(AllCor,1);

      area optioal
      area=(cellcounter:cellcounter+length(FRmeansS)-1)= contains(sessions{s},'M2');
    
   % change cell counter every iteration
    cellcounter=cellcounter+length(FRmeansS);
if exist('Power.mat')==0
load('ChannelMap128Unmapped.mat') 
FileName=['100_CH',num2str(chanMap(15)),'.continuous']  %% changes the filename on every iteration of the loop
[Data, Timestamps, Info] = load_open_ephys(FileName);
CutLow=500; %Low Cut off
Fn=30000; % Sampling Rate Raw Data
[b,a]   = butter(3,[CutLow]/(Fn/2),'low');
samplerate=2000;
dataRaw=filter(b,a,double(Data));
data=dataRaw(1:(30000/samplerate):end);
[Power,Phase] = TFA(Data',[1:100],5,2000); 
save ('Power.mat', 'Power');
save ('Phase.mat', 'Phase');
clearvars Power Phase
end
end

% Waveform Analysis on all extracted waveforms 
% determine NeuronType by cluster analysis of all cells
good=find(includeWF==1);
% X=horzcat(spikeWidth(good)',SpikeV2P(good));
% [idxInt,C] = kmeans(X,2);
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
% set manual threshold


%% Save
clearvars -except FRmeans FRall5m FRstd Animals drug batch genotype Depth includeWF neuronType SpikeV2P spikeWidth wfs AutoCorrelograms wfsRaw
%Data Folder for Jan !!!!
%cd('T:\arto\data\2016-0098-003\AD_Lev Experiment_final\analysis_Jan')

%Data Folder for Arto
cd('S:\rd-psy-TimeCell\Tompouce_8thJan18\data\2016-0098-003\AD_Lev Experiment_final\analysis_Jan')
save AD_FiringRates_004    

%%% analysis scatch
% stats 

clear all
cd('S:\rd-psy-TimeCell\Tompouce_8thJan18\data\2016-0098-003\AD_Lev Experiment_final\analysis_Jan')
load AD_FiringRates_004



%% plots
AbsRefra=sum(AutoCorrelograms(:,1:2),2);
PerRefra=(sum(AutoCorrelograms(:,1:2),2)./sum(AutoCorrelograms,2))*100;
Exclude=zeros(1,size(FRmeans,2));
Exclude(AbsRefra>50)=1;
Exclude(PerRefra>1)=1;
Exclude(includeWF==0)=1;
Exclude(neuronType>2)=1;
% I think 3-way anova (genotype,drug,neuronType) is the way to go here
% just ttest would not be enough because no correction for multiple
% comparisons
% if significant you can go on with individual ttests
[h,p]=kstest(FRmeans(Exclude==0))

%% twoway first
[p,tbl,stats] = anovan(FRmeans(Exclude==0&drug==0),{genotype(Exclude==0&drug==0);neuronType(Exclude==0&drug==0)},'model','interaction');


% Overall Effect of tg vs wt

names={'tg','wt'};
ylab={'Firing Rates Hz'};
ttl={'Firing Rates for tg / wt'};

data{1}=FRmeans(genotype==1&Exclude==0);
data{2}=FRmeans(genotype==0&Exclude==0);

AD_BeehivePlot(data,names,ylab,ttl)

% This effect is specific to tg pyramidal cells 

names={'Tg pyramidal','WT pyramidal'};
ylab={'Firing Rates Hz'};
ttl={'Pyramidal cells tg / wt'; 'in saline'};

data{1}=FRmeans(genotype==1&neuronType==2&drug==0&Exclude==0);
data{2}=FRmeans(genotype==0&neuronType==2&drug==0&Exclude==0);

AD_BeehivePlot(data,names,ylab,ttl)

% and there is no effect of genotype on interneurons 

names={'Tg Interneurons','WT Interneurons'};
ylab={'Firing Rates Hz'};
ttl={'Interneurons tg / wt'; 'in saline'};

data{1}=FRmeans(genotype==1&neuronType==1&drug==0&Exclude==0);
data{2}=FRmeans(genotype==0&neuronType==1&drug==0&Exclude==0);

AD_BeehivePlot(data,names,ylab,ttl)


% but levetiracetam can reverse Hypoactivity in Pyramdical cells in
% transgenic


names={'tg - pyr - leve','tg - pyr - sal'};
ylab={'Firing Rates Hz'};
ttl={'Firing Rates for tg - pyr - leve/sal'};

data{1}=FRmeans(genotype==1&neuronType==2&drug==1&Exclude==0);
data{2}=FRmeans(genotype==1&neuronType==2&drug==0&Exclude==0);

AD_BeehivePlot(data,names,ylab,ttl)