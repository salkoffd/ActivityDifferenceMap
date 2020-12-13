function [adm] = ActivityDifferenceMap( responseHits, responseMisses, nIterations)
%Function by David Salkoff dovsalkoff@gmail.com
%ActivityDifferenceMap Calculate difference between two trial types on a
%pixel-by-pixel basis. Calculate significance of pixels and clusters with a
%permutation test (does not assume gaussian distribution of data). Clusters
%contain pixels with p<alpha (default 0.05) and which all have the same
%sign in the activity difference map i.e. clusters can either be negative
%or positive sign but cannot contain pixels of both sign.
%
%   Outputs
%   adm is a structure containing the following:
%   adm.map is the map of activity difference.
%   adm.pmap is the p-value estimate per pixel. (note: adjacent pixels with p<alpha could belong to different clusters because the adm value is opposite sign.)
%   adm.nPixelsPval is the p-value associated with the total number of pixels with significant difference (compared with null distribution)
%   adm.clusterSizePvals is the p-value associated with the size of each cluster
%   adm.clusterMassPvals is the p-value associated with the mass of each cluster
%   adm.clustersSigPos is a logical map showing where positive sign clusters are (Hit trial activity significantly larger than Miss trial activity)
%   adm.clustersSigPos is a logical map showing where negative sign clusters are (Hit trial activity significantly smaller than Miss trial activity)
%
%   Inputs
%   responseHits is a h x w x t matrix (3D) of hit trial responses (or
%   false-alarms). hxw is the height and width of image. t is the trial
%   number.
%   videosMisses is a h x w x t matrix (3D) of miss trial responses (or
%   correct rejections).
%   nIterations is the number of iterations for the permutations test (default 100)

%plot figures?
plot_figs = true;
%alpha level for significance
alpha = 0.05; %default 0.05

if nargin<3
    nIterations = 1000;
end

%calculate difference between hits and misses
diffMap = nanmean(responseHits,3) - nanmean(responseMisses,3);

%number of trials for each trial type
nTrials_hits = size(responseHits,3);
nTrials_misses = size(responseMisses,3);
%calculate statistical significance for each pixel
%pre-allocate pmap
pmap = nan(size(responseHits,1),size(responseHits,2)); %significance level that value is sig different from shuffled
% pre-allocate shuffled distributions for permutation test
%will be used to test each pixel value for significance of activity difference as well as measuring image statistics (cluster sizes, nPixels significant (uncorrected)).
ad_shuff_dist = nan(size(responseHits,1),size(responseHits,2),nIterations);
%chose nIterations sets of random numbers corresponding to random trials.
%For each iteration of the permutation, each per-pixel activity difference
%will be computed using the same set of trials.
labels_real = [true(1,nTrials_hits),false(1,nTrials_misses)]; %true trial labels to be permuted
labels_fake = false(nTrials_hits+nTrials_misses,nIterations);
for iPerm = 1:nIterations
    labels_fake(:,iPerm) = labels_real(randperm(nTrials_hits+nTrials_misses));
end
for i=1:size(responseHits,1)
    for ii=1:size(responseHits,2)
        if isnan(diffMap(i,ii))
            ad_shuff_dist(i,ii,:) = nan;
            pmap(i,ii) = 1;
        else
            vals_all = vertcat(squeeze(responseHits(i,ii,:)),squeeze(responseMisses(i,ii,:))); %all hit/miss values for this pixel
            for iPerm = 1:nIterations
                %generate shuffled dataset chosen from pool of all responses
                vals_hit_fake = vals_all(labels_fake(:,iPerm));
                vals_miss_fake = vals_all(~labels_fake(:,iPerm));
                ad_shuff_dist(i,ii,iPerm) = nanmean(vals_hit_fake) - nanmean(vals_miss_fake); %activity difference for shuffled values
            end
            %calculate p value of measured activity difference
            centile = nnz(squeeze(ad_shuff_dist(i,ii,:))<diffMap(i,ii))/nnz(~isnan(ad_shuff_dist(i,ii,:))); %percentile in decimal form. If no NaNs then denominator is equal to nIterations
            %p-value (two-tailed test for difference in means)
            if isnan(centile)
                pmap(i,ii) = 1;
            else
                pmap(i,ii) = 1-(2*abs(centile-.5));   %1-centile %two-tailed test
            end
        end
    end
end

%measure shuffled data image statistics using ad_shuff_dist
%for each value generated from the null distribution, calculate (uncorrected) significance
ad_shuff_sig = false(size(ad_shuff_dist)); %1 for significant, otherwise 0
for i=1:size(responseHits,1)
    for ii=1:size(responseHits,2)
        for iPerm = 1:nIterations
            %calculate p value of shuffled activity difference
            centile = nnz(squeeze(ad_shuff_dist(i,ii,:))<ad_shuff_dist(i,ii,iPerm))/nnz(~isnan(ad_shuff_dist(i,ii,:))); %percentile in decimal form. If no NaNs then denominator is equal to nIterations
            %p-value (two-tailed test for difference in means)
            if isnan(centile)
                ad_shuff_sig(i,ii,iPerm) = false;
            else
                ad_shuff_sig(i,ii,iPerm) = 1-(2*abs(centile-.5)) < alpha;   %1-centile %two-tailed test
            end
        end
    end
end
%measure total number of (uncorrected) significant pixels
nPixelsSig_dist = nan(1,nIterations);
for iPerm = 1:nIterations
    nPixelsSig_dist(iPerm) = nnz(ad_shuff_sig(:,:,iPerm));
end
%nPixelsSig_dist = sort(nPixelsSig_dist,'descend'); %sort for easy indexing
%nPixelsSigCritical = nPixelsSig_dist(floor(nIterations*alpha)+1); %critical value is c+1 largest member of permutation distribution where c=floor(alpha*N)
%calculate p-value for total number of significant (uncorrected) pixels
nPixelsPval = 1 - (nnz(nnz(pmap<alpha)>nPixelsSig_dist)/nIterations); %one-tailed test to see if number of pixels is larger than null hypothesis

%measure size and mass of clusters
clusterSizes_dist = nan(1,nIterations);   %sizes of largest permuted clusters
clusterMasses_dist = nan(1,nIterations);   %masses of most massive permuted clusters
for iPerm = 1:nIterations
    pic = ad_shuff_dist(:,:,iPerm);
    statsPositive = regionprops(ad_shuff_sig(:,:,iPerm)&pic>0,'Area','PixelIdxList'); %positive difference cluster stats
    statsNegative = regionprops(ad_shuff_sig(:,:,iPerm)&pic<0,'Area','PixelIdxList'); %negative difference cluster stats
    statsAreaAll = horzcat([statsPositive.Area],[statsNegative.Area]); %area of all positive and negative difference clusters
    if isempty(statsAreaAll)
        clusterSizes_dist(iPerm) = 0;
        clusterMasses_dist(iPerm) = 0;
    else
        clusterSizes_dist(iPerm) = max(statsAreaAll);
        %find heaviest cluster
        massesPositive = zeros(1,length([statsPositive.Area]));
        massesNegative = zeros(1,length([statsNegative.Area]));
        for iROI = 1:length([statsPositive.Area])
            massesPositive(iROI) = sum(pic(statsPositive(iROI).PixelIdxList));
        end
        for iROI = 1:length([statsNegative.Area])
            massesNegative(iROI) = sum(pic(statsNegative(iROI).PixelIdxList));
        end
        clusterMasses_dist(iPerm) = max([abs(massesPositive),abs(massesNegative)]); %absolute mass of heaviest cluster
    end
end
%cluster size critical value
clusterSizes_dist = sort(clusterSizes_dist,'descend'); %sort for easy indexing
clusterSizeSig = clusterSizes_dist(floor(nIterations*alpha)+1); %critical value is c+1 largest member of permutation distribution where c=floor(alpha*N)
%cluster mass critical value
clusterMasses_dist = sort(clusterMasses_dist,'descend'); %sort for easy indexing
clusterMassSig = clusterMasses_dist(floor(nIterations*alpha)+1);

%find size,mass of clusters in measured data
statsPositive = regionprops((pmap<alpha)&diffMap>0,'Area','PixelIdxList');
statsNegative = regionprops((pmap<alpha)&diffMap<0,'Area','PixelIdxList');
statsAreaAll =  horzcat([statsPositive.Area],[statsNegative.Area]); %area of all positive and negative difference clusters
clusterSizes = sort(statsAreaAll,'descend');
clusterMassesPositive = zeros(1,length([statsPositive.Area]));
clusterMassesNegative= zeros(1,length([statsNegative.Area]));
for iCluster = 1:length([statsPositive.Area])
    clusterMassesPositive(iCluster) = sum(diffMap(statsPositive(iCluster).PixelIdxList));
end
for iCluster = 1:length([statsNegative.Area])
    clusterMassesNegative(iCluster) = abs(sum(diffMap(statsNegative(iCluster).PixelIdxList))); %absolute mass
end
clusterMasses = sort([clusterMassesPositive,clusterMassesNegative],'descend');
% %delete list of cluster sizes, masses under critical value DON'T DO THIS BECAUSE THEN WE CAN'T FIGURE OUT WHICH CLUSTERS ARE SIGNIFICANT AFTER FUNCTION ENDS
% clusterSizes(clusterSizes<clusterSizeSig) = [];
% clusterMasses(clusterMasses<clusterMassSig) = [];
%calculate p-values of cluster sizes, masses
clusterSizePvals = zeros(1,length(clusterSizes)); %preallocate
for iCluster = 1:length(clusterSizes)
    clusterSizePvals(iCluster) = 1 - (sum(clusterSizes_dist<clusterSizes(iCluster))/nIterations);
end
clusterMassPvals = zeros(1,length(clusterMasses));
for iCluster = 1:length(clusterMasses)
    clusterMassPvals(iCluster) = 1 - (sum(clusterMasses_dist<clusterMasses(iCluster))/nIterations);
end

%find size of active clusters
%positive clusters
cc = bwconncomp((pmap<alpha)&(diffMap>0));
stats = regionprops(cc,'Area');
adm_clusterSizes = [stats.Area];   %size of positive clusters
%generate mask for only significantly large clusters
idx = find(adm_clusterSizes >= clusterSizeSig);
adm_sigClusters_positive = ismember(labelmatrix(cc), idx);
%negative clusters
cc = bwconncomp((pmap<alpha)&(diffMap<0));
stats = regionprops(cc,'Area');
adm_clusterSizes = [stats.Area];   %size of negative clusters
%generate mask for only significantly large clusters
idx = find(adm_clusterSizes >= clusterSizeSig);
adm_sigClusters_negative = ismember(labelmatrix(cc), idx);

%display significance
if nPixelsPval<alpha
    disp('Total number of active pixels significant');
else
    disp('Total number of active pixels NOT significant');
end
%how many significantly large, (and massive) clusters?
disp(['Found ' num2str(sum(clusterSizePvals<alpha)) ' significantly large clusters']);
disp(['Found ' num2str(sum(clusterMassPvals<alpha)) ' significantly massive clusters']);

%save data in structure
adm.map = diffMap;
adm.pmap = single(pmap);
adm.nPixelsPval = nPixelsPval;
adm.clusterSizeSig = clusterSizeSig;
adm.clusterMassSig = clusterMassSig;
adm.clusterSizePvals = single(clusterSizePvals);
adm.clusterMassPvals = single(clusterMassPvals);
adm.clustersSigPos = adm_sigClusters_positive;
adm.clustersSigNeg = adm_sigClusters_negative;
adm.nTrialA = nTrials_hits;
adm.nTrialB = nTrials_misses;
adm.nIterations = nIterations;

if plot_figs == true
    [B,~] = bwboundaries(adm_sigClusters_positive);
    [C,~] = bwboundaries(adm_sigClusters_negative);
    colorMax = ceil(max(diffMap(:))*100);
    colorMin = floor(min(diffMap(:))*100);
    figure
    imshow(diffMap*100,[colorMin colorMax],'InitialMagnification',1000);
    %set(h,'AlphaData',adm_sigClusters_positive+.4)
    hold on
    for k = 1:length(B)
        boundary = B{k};
        plot(boundary(:,2), boundary(:,1), 'k', 'LineWidth', 1)
    end
    for k = 1:length(C)
        boundary = C{k};
        plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 1)
    end
    title('Activity difference with outlined significantly large clusters')
    colormap(gca,jet(255))
    h2 = colorbar;
    ylabel(h2, 'dF/F (%)')
    hold off
end

end

