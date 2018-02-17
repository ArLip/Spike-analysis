
function [AllCor]=AD_AutoCorrelograms

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


%computes ratio of spikes during refractory period and shortly after
%refractory period
win=500;
AllCor=zeros(size(GoodClusters,1),win);
for i=1:size(GoodClusters,1)
    x=double(spike_times(find(spike_clusters==(GoodClusters(i)))));
    if size(x,1)>2
    InterSpikeIntervals=diff(x);
    [Xcor]=histcounts(InterSpikeIntervals,ceil(max(InterSpikeIntervals)/30));
    if length(Xcor)>win
    AllCor(i,:)=Xcor(1:win);
    else
    AllCor(i,1:length(Xcor))=Xcor;
    end
    meanRefra=mean(Xcor(1:2));
    if size(Xcor,2)>15
    meanAfterRefra=mean(Xcor(10:15));
    Refraclean(i)=meanRefra/meanAfterRefra;
    else
     Refraclean(i)=NaN;
    end
    else
    Refraclean(i)=NaN;
    end
end   


