clear 

% load relevant data files 
load('Data/AllChirps.mat')

% go through all the syllables or "chirps" pulled from the virtual sound traces 
for ii = 1:size(myChirps, 2)
    chirp_sm = myChirps{ii}(1:6000);
    S = spectrogram(chirp_sm, 100); % get the spectrogram 
    S_abs = abs(S); 
    if ii == 1
        mySpecs = zeros(size(S_abs,1), size(S_abs,2), size(myChirps, 2));
    end
    mySpecs(:,:,ii) = flipud(S_abs); % save all the spectrograms in one 3D matrix 
end

% initialize a variable 
myCorrMat = zeros(size(myChirps, 2));

% populate the variable with the 2D correlation of each syllable or chirp spectrogram with every other one 
for ii = 1:size(myCorrMat, 1)
    for jj = 1:size(myCorrMat, 2)
        myCorrMat(ii,jj) = corr2(mySpecs(:,:,ii), mySpecs(:,:,jj));
    end
end

% use the correlation matrix for spectrograms to create a pdist vector 
iter = 0;
for ii = 1:(size(myCorrMat, 1) - 1)
    for jj = (ii+1):size(myCorrMat, 1)
        iter = iter + 1;
        myPDist(iter) = 1 - myCorrMat(ii,jj);
    end
end

% create a linkage with the pdist vector 
Z = linkage(myPDist);

% pre-define a maximum of number of syllable or chirp clusters to find (with > 90 chirps, 40 seems reasonable)
numOfClusters = 40;
meanClustDist = zeros(numOfClusters, 1);
stdClustDist = zeros(numOfClusters, 1);

% go through and perfrom hierarchical clustering with increasing number of clusters 
% and check mean and stdev cluster distance 
for ii = 1:numOfClusters
    T = cluster(Z, 'maxclust', ii);
    totalDistForEach = zeros(ii,1);
    for jj = 1:ii
        myIndices = find(T == jj);
        iter = 0;
        if size(myIndices, 1) > 1
            for kk = 1:size(myIndices, 1)
                for ll = kk+1:size(myIndices, 1)
                    iter = iter + 1; 
                    listOfDists(iter) = 1 - myCorrMat(myIndices(kk), myIndices(ll));
                end
            end
        else 
            listOfDists = 0;
        end
        totalDist = sum(listOfDists);
        totalDistForEach(jj,1) = totalDist;
    end
    meanClustDist(ii,1) = mean(totalDistForEach);
    stdClustDist(ii,1) = std(totalDistForEach);    
end

% plot mean +/- stdev cluster distance as a function of the number of clusters 
errorbar(meanClustDist, stdClustDist, 'o-')
title('Mean distance with number of clusters')
xlabel('Number of clusters')
ylabel('Mean within cluster distance')

% from the plot it's clear that more than 11 clusters show little to no change in mean cluster distance 
numOfClusters = 11;
T = cluster(Z, 'maxclust', numOfClusters);
chirpsByCluster = cell(numOfClusters, 1);
durationByCluster = cell(numOfClusters, 1);

% cluster syllables or chirps into 11 main clusters with each containing similar syllables
for ii = 1:numOfClusters 
    myIndices = find(T == ii);
    clustChirps = zeros(size(myIndices, 1), 6000);
    durationChirps = zeros(size(myIndices, 1), 1);
    for jj = 1:size(myIndices, 1)
        thisChirp = myChirps{myIndices(jj)};
        thisChirp = thisChirp(1:6000); 
        clustChirps(jj,:) = thisChirp;
        durationChirps(jj) = chirpDurations(myIndices(jj));
    end
    chirpsByCluster{ii} = clustChirps;
    durationByCluster{ii} = durationChirps;
end

% save clustered chirp data
save('Data/ChirpsByCluster.mat', 'chirpsByCluster', 'durationByCluster');        
        
        
        
