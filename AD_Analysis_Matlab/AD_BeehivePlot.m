
% Function to create beehive plots of any kind of data in this case of firing 
% rates of AD neurons; Input data and names should be cells the same length
% Jan Klee 15.11.17


function []=AD_beehivePlot(data,names,ylab,ttl)
figure()
col=rand(size(data,2),3);

for ds=1:length(data)
    
    
    [n,edges]=histcounts(data{ds},length(data{ds}));
    sortedData=sort(data{ds});
    counter=1;
    for i=1:length(n)

    if n(i)>0
        if counter+n(i)<length(n)
        bindata=sortedData(counter:counter+n(i));
        else
        bindata=sortedData(counter:end);
        end
        xmin=((-n(i)/2)/10)+ds;
        xmax=((n(i)/2)/10)+ds;
        ni=n(i);
        rp=xmin+rand(1,ni)*(xmax-xmin);
        for ii=1:n(i)
            plot(rp(ii),bindata(ii),'color',col(ds,:),'marker','o','LineWidth',2);
            hold on
        end
    else
        continue
    end

 counter=counter+n(i);
 clearvars bindata rp
    end

%mean and error bars    
plot(ds,(mean(data{ds})),'color',[0 0 0],'marker','square','LineWidth',3);
errorbar(ds,mean(data{ds}),std(data{ds}/sqrt(length(data{ds}))),'color',[0 0 0],'LineWidth',3); 

if ds>1
%significants stars
[p,h]=ranksum(data{ds-1},data{ds})
if p>0.05
    txt=' n.s. ' ;
elseif p<0.001    
    txt=' *** ' ;
elseif p<0.01    
    txt=' ** ';
elseif p<0.05   
    txt=' * ';
end
    t=text(ds-.5,max(data{ds})/2,txt);
    s = t.FontSize;
    t.FontSize = 20;
end

%labels
xticks(1:length(data));
xticklabels(names);
ylabel(ylab);
title(ttl);
set(gca,'FontSize',20);
box off

end
