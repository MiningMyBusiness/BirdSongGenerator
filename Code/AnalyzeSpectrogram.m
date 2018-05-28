clear 

load('AllChirps.mat')

for ii = 1:size(myChirps, 2)
    chirp_sm = myChirps{ii}(1:6000);
    S = spectrogram(chirp_sm, 100);
    S_abs = abs(S); 
    thisFft = fft(chirp_sm);
    if ii == 1
        mySpecs = zeros(size(S_abs,1), size(S_abs,2), size(myChirps, 2));
        myFfts = zeros(size(myChirps, 2), size(chirp_sm, 1));
    end
    mySpecs(:,:,ii) = flipud(S_abs);
    myFfts(ii,:) = abs(thisFft)';
end

myCorrMat = zeros(size(myChirps, 2));

for ii = 1:size(myCorrMat, 1)
    for jj = 1:size(myCorrMat, 2)
        myCorrMat(ii,jj) = corr2(mySpecs(:,:,ii), mySpecs(:,:,jj));
    end
end


iter = 0;
for ii = 1:(size(myCorrMat, 1) - 1)
    for jj = (ii+1):size(myCorrMat, 1)
        iter = iter + 1;
        myPDist(iter) = 1 - myCorrMat(ii,jj);
    end
end
fftDist = pdist(myFfts);

Z = linkage(myPDist);
Z_fft = linkage(fftDist);

numOfClusters = 40;
meanClustDist = zeros(numOfClusters, 1);
stdClustDist = zeros(numOfClusters, 1);

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

errorbar(meanClustDist, stdClustDist, 'o-')
title('Mean distance with number of clusters')
xlabel('Number of clusters')
ylabel('Mean within cluster distance')

numOfClusters = 11;
T = cluster(Z, 'maxclust', numOfClusters);
chirpsByCluster = cell(numOfClusters, 1);
durationByCluster = cell(numOfClusters, 1);

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
        
save('ChirpsByCluster.mat', 'chirpsByCluster', 'durationByCluster');        
        
        
        