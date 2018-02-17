


% Function to analyse waveforms of neurons that where clustered with kilosort
% and phy. Uses a version of the getWaveforms function provided by the
% cortex lab. Input 1: number of Spikes that should be used to compute the
% mean waveforms. Input 2: 1 to plot result. Uncomment kmeans clustering
% section to determine neuron class automatically but this works better for
% bigger datasets posthoc. Jan Klee 15.11.17


function  [wfs,include,SpikeV2P,spikeWidth,Depth,wfMeanBest]=AD_waveFormAnalysis(SpikeNum,Fig)% spikev2p n1-P2
% SpikeNum=100;
%find good clusters
spike_clusters=readNPY('spike_clusters.npy');
spike_times=readNPY('spike_times.npy');
fileID = fopen('cluster_groups.csv','r');
delimiter = '\t';
startRow = 2;
formatSpec = '%f%s%[^\n\r]';
ClusterQual = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'HeaderLines' ,startRow-1, 'ReturnOnError', false);
fclose(fileID);
clearvars filename delimiter startRow formatSpec fileID ans;
GoodClusters=ClusterQual{1,1}(strmatch('good',ClusterQual{1,2}));


%use kilosort team script to extract waveforms&waveform mean
%use loop to avoid memory issues
load ChannelMap128Unmapped.mat
for ii=1:size(GoodClusters,1)
curpath=pwd;
gwfparams.dataDir = [curpath];    % KiloSort/Phy output folder
gwfparams.fileName = 'Data4Kilo.dat';         % .dat file containing the raw 
gwfparams.dataType = 'int16';            % Data type of .dat file (this should be BP filtered)
gwfparams.nCh = 128;                      % Number of channels that were streamed to disk in .dat file
gwfparams.wfWin = [-40 41];              % Number of samples before and after spiketime to include in waveform
gwfparams.nWf = SpikeNum;                    % Number of waveforms per unit to pull out
gwfparams.spikeTimes =    spike_times(find(ismember(spike_clusters,GoodClusters(ii)))); % Vector of cluster spike times (in samples) same length as .spikeClusters
gwfparams.spikeClusters = spike_clusters(find(ismember(spike_clusters,GoodClusters(ii)))); % Vector of cluster IDs (Phy nomenclature)   same length as .spikeTimes
wf = AD_getWaveForms(gwfparams)

wfMeans=wf.waveFormsMean;

%find channel with stongest signal for every cluster
[x,MinChannel]=min(squeeze(min(wfMeans,[],3)),[],2);
wfMeanBest(ii,:)=squeeze(wfMeans(1,MinChannel,:));
Depth(ii)=ycoords(MinChannel);
clearvars wf
ii
end



%subtract baseline
wfMeanBestBsCor=(wfMeanBest-mean(wfMeanBest(:,1:10),2));

%normalise
[Mins,MinT]=(min(wfMeanBestBsCor,[],2));
wfMeanBestNorm=wfMeanBestBsCor./abs(Mins);

%alinge Peaks
Offset=MinT-40;
pad=-1*MinT+max(MinT);
for i=1:size(GoodClusters,1)
 padpre=zeros(1,pad(i));
 padpost=zeros(1,max(pad)-pad(i));
 data=wfMeanBestNorm(i,:);
 wfspad(i,:)=horzcat(padpre,data,padpost);
end
wfs=wfspad(:,gwfparams.wfWin(1)+max(MinT):max(MinT)+gwfparams.wfWin(2));


%find strange waveforms
good=find(min(wfs(:,1:25),[],2)>-.1&max(wfs(:,1:25),[],2)<.1&max(wfs,[],2)<1);
include(1:size(wfs,1))=0;
include(good)=1;

%find max hyperpolarization within 1ms from peak and time to that
[hypPol,SpikeV2P]=max(wfs(:,40:end),[],2);

%find spike with at 30% of peak
for i=1:size(wfs,1)
    belowTresh=find(wfs(i,1:41)<-.33);
    if length(belowTresh>=1)
    firstTreshX(i)=belowTresh(1);
    else
    firstThresX(i)=41;
    end
    aboveTresh=find(wfs(i,40:end)>-.33);
    if aboveTresh>=1
    secTreshX(i)=aboveTresh(1)+40;
    else
    secTreshX(i)=41;
    end
end
spikeWidth=secTreshX-firstTreshX;


% %cluster with kmeans   
% X=horzcat(spikeWidth',SpikeV2P);
% [idxC,C] = kmeans(X,2);
% 
% neuronClass=zeros(1,size(GoodClusters,1));
% neuronClass(good(idxC==1))=1;
% neuronClass(good(idxC==2))=2;

if Fig==1
figure()
scatter(spikeWidth,SpikeV2P,'r' )
hold on
scatter(spikeWidth,SpikeV2P,'b' )

figure()

for i=1:size(good,1)

    plot(wfs(good(i),:),'r')
hold on
end
end
