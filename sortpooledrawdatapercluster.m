sortedtraces=[];
astrocytefiguressortedcellsall3=[raw2314_2;raw2314_3;raw2314_4;raw2314_5;
    raw2331_1;raw2331_2;raw2331_3;raw2331_4;
    rawprimary_1;rawprimary_2;rawprimary_3;rawprimary_4];
for cluster=1:5
    thatclustersignal=astrocytefiguressortedcellsall3(astrocytefiguressortedcellsall3(:,end)==cluster,1:end-1);
    for i=1:size(thatclustersignal,1)
        if thatclustersignal(i,9)==1 && thatclustersignal(i,10)==1 
            sortedtraces=[sortedtraces;rawprimary_1(1,:)];
            rawprimary_1=rawprimary_1(2:end,:);
        end
        if thatclustersignal(i,9)==1 && thatclustersignal(i,10)==2 
            sortedtraces=[sortedtraces;rawprimary_2(1,:)];
            rawprimary_2=rawprimary_2(2:end,:);
        end        
        if thatclustersignal(i,9)==1 && thatclustersignal(i,10)==3 
            sortedtraces=[sortedtraces;rawprimary_3(1,:)];
            rawprimary_3=rawprimary_3(2:end,:);
        end        
        if thatclustersignal(i,9)==1 && thatclustersignal(i,10)==4 
            sortedtraces=[sortedtraces;rawprimary_4(1,:)];
            rawprimary_4=rawprimary_4(2:end,:);
        end        
        if thatclustersignal(i,9)==2331 && thatclustersignal(i,10)==1 
            sortedtraces=[sortedtraces;raw2331_1(1,:)];
            raw2331_1=raw2331_1(2:end,:);
        end        
        if thatclustersignal(i,9)==2331 && thatclustersignal(i,10)==2 
            sortedtraces=[sortedtraces;raw2331_2(1,:)];
            raw2331_2=raw2331_2(2:end,:);
        end        
        if thatclustersignal(i,9)==2331 && thatclustersignal(i,10)==3 
            sortedtraces=[sortedtraces;raw2331_3(1,:)];
            raw2331_3=raw2331_3(2:end,:);
        end        
        if thatclustersignal(i,9)==2331 && thatclustersignal(i,10)==4 
            sortedtraces=[sortedtraces;raw2331_4(1,:)];
            raw2331_4=raw2331_4(2:end,:);
        end        
        if thatclustersignal(i,9)==2314 && thatclustersignal(i,10)==2 
            sortedtraces=[sortedtraces;raw2314_2(1,:)];
            raw2314_2=raw2314_2(2:end,:);
        end        
        if thatclustersignal(i,9)==2314 && thatclustersignal(i,10)==3 
            sortedtraces=[sortedtraces;raw2314_3(1,:)];
            raw2314_3=raw2314_3(2:end,:);
        end        
        if thatclustersignal(i,9)==2314 && thatclustersignal(i,10)==4 
            sortedtraces=[sortedtraces;raw2314_4(1,:)];
            raw2314_4=raw2314_4(2:end,:);
        end        
        if thatclustersignal(i,9)==2314 && thatclustersignal(i,10)==5 
            sortedtraces=[sortedtraces;raw2314_5(1,:)];
            raw2314_5=raw2314_5(2:end,:);
        end        
    end
end
writematrix(sortedtraces,'C:\Users\lbinan\Desktop\astrocyte_figuresagain\bettersmoothing\sortedtracesall.csv')
