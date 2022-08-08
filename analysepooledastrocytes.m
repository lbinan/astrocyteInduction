features_2331=[features2331_1;features2331_2;features2331_3;features2331_4];
featuresprimary=[featuresprimary_1;featuresprimary_2;featuresprimary_3;featuresprimary_4];
activecellsprimary=[activecellsprimary_1;activecellsprimary_2;activecellsprimary_3;activecellsprimary_4];
activecells2331=[activecells2331_1;activecells2331_2;activecells2331_3;activecells2331_4];
% 
activecells=[activecellsprimary;activecells2331;activecells2314];
features=[featuresprimary;features_2331;features_2314];
activecells=[activecells(:,4:11),activecells(:,480:481)];
features=[features(:,4:11),features(:,480:481)];

clear activecellsprimary activecells2331 activecells2314 featuresprimary features_2331 features_2314
clear features2314_2 features2314_3 features2314_4 features2314_5 activecells2314_2 activecells2314_3 activecells2314_4 activecells2314_5
clear features2331_1 features2331_2 features2331_3 features2331_4
clear activecellsprimary_1 activecellsprimary_2 activecellsprimary_3 activecellsprimary_4
clear activecells2331_1 activecells2331_2 activecells2331_3 activecells2331_4 
clear featuresprimary_1 featuresprimary_2 featuresprimary_3 featuresprimary_4
kvalue=5;
opts = statset('Display','final');
[idx,C]=kmeans(activecells(:,1:8),kvalue,'Distance','correlation','Replicates',1000,'MaxIter',10000,'Options',opts);
sortedcellsall=activecells(idx==1,:);
sortedcellsall(:,end+1)=1;
for i=2:max(idx)
    addthis=activecells(idx==i,:);
    addthis(:,end+1)=i;
    sortedcellsall=[sortedcellsall;addthis];
end
writematrix(sortedcellsall,'C:\Users\lbinan\Desktop\astrocyte_figuresagain\bettersmoothing\sortedcellsall35clusters.csv')
shuflled=[];
for i=1:max(idx)
    temp=sortedcellsall(sortedcellsall(:,end)==i,:);
    shuflled=[shuflled;temp(randperm(size(temp, 1)), :)];
end
writematrix(sortedcellsall,'C:\Users\lbinan\Desktop\astrocyte_figuresagain\bettersmoothing\sortedcellsall35clustersshuflled.csv')

mydataall=sortedcellsall(:,1:8);
mycorrelationsall=corrcoef(mydataall');
evenhalf=mydataall(2:2:end,:);
oddhalf=mydataall(1:2:end,:);
mycorrelationsall=corrcoef(oddhalf');

figure,
imagesc(mycorrelationsall);
colormap('jet')
colorbar
ratios=zeros(12,kvalue);
globalratios=zeros(3,kvalue);
for cluster=1:5
    thiscluster=sortedcellsall(sortedcellsall(:,end)==cluster,:);
    primary=thiscluster(thiscluster(:,9)==1,:);
    induced=thiscluster(thiscluster(:,9)==2331,:);
    induced_2=thiscluster(thiscluster(:,9)==2314,:);
    globalratios(1,cluster)=size(primary,1);
    globalratios(2,cluster)=size(induced,1);
     globalratios(3,cluster)=size(induced_2,1);
    for dataset=1:4
        thisdatasetprimary=primary(primary(:,10)==dataset,:);
        thisdatasetinduced=induced(induced(:,10)==dataset,:);
        thisdatasetinduced_2=induced_2(induced_2(:,10)==dataset+1,:);
        ratios(dataset,cluster)=size(thisdatasetprimary,1);
        ratios(4+dataset,cluster)=size(thisdatasetinduced,1);
        ratios(8+dataset,cluster)=size(thisdatasetinduced_2,1);
    end
end

writematrix(ratios,strcat('C:\Users\lbinan\Desktop\astrocyte_figuresagain\bettersmoothing\ratiospersampleglobal',num2str(kvalue),'.csv'))
writematrix(globalratios,strcat('C:\Users\lbinan\Desktop\astrocyte_figuresagain\bettersmoothing\ratiospercelltypeglobal',num2str(kvalue),'.csv'))
clustersfeatures=zeros(max(idx),8);
clustersstd=zeros(max(idx),8);

for thiscluster=1:max(idx)
    clustersfeatures(thiscluster,:)=mean(features(idx==thiscluster,1:8),1);
    clustersstd(thiscluster,:)=std(features(idx==thiscluster,1:8),1);
end

writematrix(clustersfeatures,strcat('C:\Users\lbinan\Desktop\astrocyte_figuresagain\bettersmoothing\clustersfeaturesglobal',num2str(kvalue),'.csv'))
writematrix(clustersstd,strcat('C:\Users\lbinan\Desktop\astrocyte_figuresagain\bettersmoothing\clustersstdglobal',num2str(kvalue),'.csv'))

