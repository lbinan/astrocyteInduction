astrocytefiguressortedcellsall3=[rawprimary_1;rawprimary_2;rawprimary_3;rawprimary_4;
    raw2331_1;raw2331_2;raw2331_3;raw2331_4;
    raw2314_2;raw2314_3;raw2314_4;raw2314_5;
    ];

sortedall=astrocytefiguressortedcellsall3(idx==1,:);
sortedall(:,end+1)=1;
for cluster=2:5
    addthis=astrocytefiguressortedcellsall3(idx==cluster,:);
    addthis(:,end+1)=cluster;
    sortedall=[sortedall;addthis];
end
smoothedsortedall=sortedall(:,1:end-1);
for i=1:size(sortedall,1)
    thistrace=sortedall(i,:);
newtrace=smooth(smooth(thistrace(1:end-1)-min(thistrace(1:end-1))),50,'sgolay');
smoothedsortedall(i,:)=newtrace;
end
features(:,13)=sortedall(:,end);
for mycluster=1:max(idx)
    %     subplot (max(idx),1,mycluster)
    cluster1cellsspikesdenoised=smoothedsortedall(features(:,13)==mycluster,1:(end-1));
    cluster1cellsspikesdenoised=cluster1cellsspikesdenoised(randperm(size(cluster1cellsspikesdenoised, 1)), :);
    cluster1cellsspikesdenoised=cluster1cellsspikesdenoised(1:min(60,size(cluster1cellsspikesdenoised,1)),:);
     maxvaluethis=max(cluster1cellsspikesdenoised,[],2);
     minvaluethis=min(cluster1cellsspikesdenoised,[],2);
     cluster1cellsspikesdenoised=cluster1cellsspikesdenoised-minvaluethis;
    maxvaluethis=1;
    cluster1cellsspikesdenoised=cluster1cellsspikesdenoised./maxvaluethis;
     maxvalue=1.8*prctile(cluster1cellsspikesdenoised,97,'all');

    if mycluster==1
        mytraces=uint8(255/maxvalue*cluster1cellsspikesdenoised);
    else
        mytraces=[mytraces;uint8(255/maxvalue*cluster1cellsspikesdenoised)];
       %imshow(uint8(255/maxvalue*cluster1cellsspikesdenoised))
    end
    imwrite(uint8(255/maxvalue*cluster1cellsspikesdenoised),strcat('C:\Users\lbinan\Desktop\astrocyte_figuresagain\bettersmoothing\tracecluster',num2str(mycluster),'_60.png'));
    if mycluster<max(idx)
        mytraces=[mytraces;200*ones(20, size(mytraces,2))];
    end
    time=[0:0.5:478.5*0.5];
    for i=1:60
        newtrace=cluster1cellsspikesdenoised(i,:);
        findpeaks(newtrace,time(1:end-1),'MinPeakHeight',10,'MinPeakProminence',10,'MinPeakDistance',6);
        hold on
    end
end
imshow(mytraces)
imwrite(uint8(255/maxvalue*cluster1cellsspikesdenoised),strcat('C:\Users\lbinan\Desktop\astrocyte_figuresagain\bettersmoothing\Finaltracecluster_60.png'));
