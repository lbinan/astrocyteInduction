thosepaths=["/broad/clearylab/Users/Loic/7_19_2022astrocyte/20220719_121218_007_2314_1";
    "/broad/clearylab/Users/Loic/7_19_2022astrocyte/20220719_122834_968_2314_2";
    "/broad/clearylab/Users/Loic/7_19_2022astrocyte/20220719_131920_567_2314_3";
    "/broad/clearylab/Users/Loic/7_19_2022astrocyte/20220719_132416_439_2314_4";
    "/broad/clearylab/Users/Loic/7_19_2022astrocyte/20220719_133905_044_2314_5";
    "/broad/clearylab/Users/Loic/7_19_2022astrocyte/20220719_122223_184_2331_1";
    "/broad/clearylab/Users/Loic/7_19_2022astrocyte/20220719_130914_235_2331_2";
    "/broad/clearylab/Users/Loic/7_19_2022astrocyte/20220719_133406_783_2331_3";
    "/broad/clearylab/Users/Loic/7_19_2022astrocyte/20220719_134554_108_2331_4";
    "/broad/clearylab/Users/Loic/7_19_2022astrocyte/20220719_121728_223_primHuman_1";
    "/broad/clearylab/Users/Loic/7_19_2022astrocyte/20220719_131407_811_primHuman_2";
    "/broad/clearylab/Users/Loic/7_19_2022astrocyte/20220719_132912_256_primHuman_3";
    "/broad/clearylab/Users/Loic/7_19_2022astrocyte/20220719_135053_894_primHuman_4"];

%% loop on all datasets

for mypathnumber=1:size(thosepaths,1)
    if mypathnumber>9
        brightness=50;
        largethresh=800;
        dust=250;
    end
    if mypathnumber<10
        brightness=150;
        largethresh=1800;
        dust=500;
    end
    if mypathnumber<10
        brightness=120;
        largethresh=1800;
        dust=500;
    end
%load the data
myimage=imread(fullfile(mypath,'Time00000_Point0000_ChannelFITC-Penta_Seq0000.tiff'));
for i=2:9
    myimage(:,:,i)=imread(fullfile(mypath,strcat('Time0000',num2str(i-1),'_Point0000_ChannelFITC-Penta_Seq000',num2str(i-1),'.tiff')));
end
for i=10:10
    myimage(:,:,i)=imread(fullfile(mypath,strcat('Time00009_Point0000_ChannelFITC-Penta_Seq000',num2str(i-1),'.tiff')));
end
for i=11:99
    myimage(:,:,i)=imread(fullfile(mypath,strcat('Time000',num2str(i-1),'_Point0000_ChannelFITC-Penta_Seq00',num2str(i-1),'.tiff')));
end
for i=100:100
    myimage(:,:,i)=imread(fullfile(mypath,strcat('Time00099_Point0000_ChannelFITC-Penta_Seq00',num2str(i-1),'.tiff')));
end
for i=101:172
    myimage(:,:,i)=imread(fullfile(mypath,strcat('Time00',num2str(i-1),'_Point0000_ChannelFITC-Penta_Seq0',num2str(i-1),'.tiff')));
end
mypath=thosepaths(mypathnumber,:)
oldoutputs=strcat(mypath,'/outputsglobal/');
outputpath=strcat(mypath,'/newoutputsglobal/');
outputs=outputpath;
mkdir(outputpath)
%% create cell mask
tomakemask=imread(fullfile(mypath,'Time00000_Point0000_ChannelFITC-Penta_Seq0000.tiff'));
% imshow(myimage*60)
bw=imbinarize(tomakemask*brightness,'adaptive','sensitivity',0.57);%150 for 1st %50for priary
 figure, imshow(bw)
cleaned=bwareaopen(bw,50);
% figure, imshow(cleaned)
closed=imclose(cleaned,strel('disk',2));
% figure, imshow(closed)
opened=imopen(closed,strel('disk',3));
% figure, imshow(opened)
largestuff=bwareaopen(opened,largethresh);%800 for primary  *induced 1800
% figure, imshow(largestuff)
D = bwdist(~largestuff);
D = -D;
I3 = imhmin(D,0.9);
split = watershed(I3);
split(~largestuff)=0;
% imshow(3*split)
finalsplit=bwareaopen(logical(split),200);
% imshow(finalsplit)
BWcleanedfrom5=logical(opened-largestuff+finalsplit);
BWcleanedfrom5=bwareaopen(BWcleanedfrom5,150);%500 for induced %250 for primary
% figure, imshow(BWcleanedfrom5)
imwrite(BWcleanedfrom5,fullfile(outputs,'BWcleanedfrom5.tif'))
clear opened largestuff finalsplit I3 D cleaned bw tomakemask
BWcleanedfrom5=imread(fullfile(oldoutputs,'BWcleanedfrom5.tif'));
clear opened largestuff finalsplit I3 D cleaned bw tomakemask
%% find peaks, and extract features
stats=regionprops(BWcleanedfrom5,'Area','Centroid','PixelIdxList');
celltable=table2array(readtable(fullfile(oldoutputs,'cellmeanfluotable.csv')));
features=zeros(size(stats,1),size(celltable,2));
time=[0:0.5:478.5*0.5];
for i=1:size(stats,1)
firingtable(i,1)=stats(i).Area;
firingtable(i,2)=stats(i).Centroid(1);
firingtable(i,3)=stats(i).Centroid(2);
end

analysedtraceS=[];
denoisedtrace=[];
for thiscell=1:size(celltable,1)
    thistrace=celltable(thiscell,:);
    thismedian=median(thistrace);
    storemin=sign(min(thistrace))*min(thistrace);
    newtrace=smooth(smooth(thistrace-min(thistrace)),50,'sgolay');
    denoisedTrace=newtrace;
    denoisedtrace=[denoisedtrace;denoisedTrace];
    analysedtraceS=denoisedTrace;
    [pks, locs]=findpeaks(newtrace,time,'MinPeakHeight',10,'MinPeakProminence',10,'MinPeakDistance',6);
    savemedian=median(denoisedTrace);
    savenumberpfpeaks=size(pks,1);
    peakInterval = diff(locs);
    firinginterval=median(peakInterval);
    firingstrength=median(pks);
    duration=[];
    areaundercurve=[];
    risingtime=[];
    fallingtime=[];
    for thispeak=1:size(locs,2)
        derived=gradient(denoisedTrace(:)) ./ gradient(time(:));
        xd = (time(2:end)+time(1:(end-1)))/2;
        indexstart = find(derived(1:2*(locs(thispeak)-10))<0, 1, 'last');
        indexend = min(479,2*(locs(thispeak)+15)+find(derived(2*(locs(thispeak)+15): end)> 0, 1, 'first'));
        duration=[duration;(indexend-indexstart)*0.5];
        areaundercurve=[areaundercurve;sum(thistrace(indexstart:indexend))];
        risingtime=[risingtime;2*locs(thispeak)-indexstart];
        fallingtime=[fallingtime;indexend-2*locs(thispeak)];
    end
    features(thiscell,4)=savemedian;
    features(thiscell,5)=savenumberpfpeaks;
    features(thiscell,6)=firinginterval;
    features(thiscell,7)=firingstrength;
    features(thiscell,8)=median(duration(duration>8));
    features(thiscell,9)=median(areaundercurve(duration>8));
    features(thiscell,10)=median(risingtime(duration>8));
    features(thiscell,11)=median(fallingtime(duration>8));
end
features(isnan(features))=0;
writematrix(features,fullfile(outputs,'featuresnotfiltered.csv'))

%cluster
activecells=features;
activecells(:,4:11)=activecells(:,4:11)-mean(activecells(:,4:11),1);
activecells(:,4:11)=activecells(:,4:11)./std(activecells(:,4:11),1);
activecells(isnan(activecells))=0;

% [coeff,score,latent] = pca(activecells(:,4:11));
opts = statset('Display','final');

[idx,C]=kmeans(activecells(:,4:11),5,'Distance','correlation','Replicates',1000,'MaxIter',100000,'Options',opts);
sortedcells=activecells(idx==1,4:11);
sortedcells(:,end+1)=1;
for i=2:max(idx)
    addthis=activecells(idx==i,4:11);
    addthis(:,end+1)=1;
    sortedcells=[sortedcells;addthis];
end
writematrix(sortedcells,fullfile(outputs,'sortedcells.csv'))
mydata=sortedcells;
mycorrelations=corrcoef(mydata');

% figure,
% imagesc(mycorrelations);
% colormap('jet')
%% create image of sample traces

rgbImage = ind2rgb(mycorrelations, colormap('jet'));
imwrite(mycorrelations,fullfile(outputs,'mycorrelations.tif'))

writematrix(mycorrelations,fullfile(outputs,'mycorrelations.csv'))
clustersfeatures=zeros(max(idx),11);
clustersstd=zeros(max(idx),11);

for thiscluster=1:max(idx)
    clustersfeatures(thiscluster,:)=mean(features(idx==thiscluster,1:11),1);
    clustersstd(thiscluster,:)=std(features(idx==thiscluster,1:11),1);
end
writematrix(clustersfeatures,fullfile(outputs,'clustersfeatures.csv'))
writematrix(clustersstd,fullfile(outputs,'clustersstd.csv'))

% figure,
features(:,13)=idx;
for mycluster=1:max(idx)
    %     subplot (max(idx),1,mycluster)
    cluster1cellsspikesdenoised=denoisedtrace(features(:,13)==mycluster,:);
    maxvaluethis=max(cluster1cellsspikesdenoised,[],2);
    cluster1cellsspikesdenoised=cluster1cellsspikesdenoised./maxvaluethis;
    maxvalue=max(cluster1cellsspikesdenoised,[],'all');
    if mycluster==1
        mytraces=uint8(255/maxvalue*cluster1cellsspikesdenoised);
    else
        mytraces=[mytraces;uint8(255/maxvalue*cluster1cellsspikesdenoised)];
    % imshow(uint8(255/maxvalue*cluster1cellsspikesdenoised))
    end
    if mycluster<max(idx)
        mytraces=[mytraces;100*ones(20, size(mytraces,2))];
    end
end
imwrite(mytraces,fullfile(outputs,'mytraces.tif'))
clear firingtable features mytraces clustersfeatures clustersstd mycorrelations activecells sortedcells celltable cellmeanfluotable
end