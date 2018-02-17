%Function to extract firing rates from neurons that were clustered with
%kilosort and phy. Input: bins specifices over how many bins firing rates
%should be averaged. Jan Klee 15.11.17

function [FRperiods,FRmean,FRstd]=AD_FiringRates(bins)

%load good clusters
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

filename=['100_CH1.continuous']; 
[Data] = load_open_ephys(filename);

%determine number of possible 1min intervals
sets=floor(length(Data)/bins); 

for ii=1:length(GoodClusters)
    
TScell(1:length(Data))=0;
TScell(spike_times(find(spike_clusters==GoodClusters(ii))))=1;
Periods=reshape(TScell(1:sets*(bins)),sets,(bins));
FRperiods(ii,:)=(sum(Periods,2)/60)'; % matrix with different bins for all neurons 
FRmean(ii)=mean(FRperiods(ii,:)); % periods grand mean 
FRstd(ii)=std(FRperiods(ii,:)); % periods grand sd
end

