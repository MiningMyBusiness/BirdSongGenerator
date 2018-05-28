%% Pick out chirps in model generation

clear

dirName = 'BirdSongSims/';
birdSongFiles = dir([dirName '*.mat']);

Fs = 22050;  %sampling

numChirps = 0;
chirpDurations = 0;
for ii = 1:size(birdSongFiles, 1)
    load([dirName birdSongFiles(ii).name]);
    s1_sm = conv(abs(s1), ones(100,1), 'same')/100;
    myEnv = s1_sm > power(10,-4);
    myEnv_diff = myEnv(2:end) - myEnv(1:end-1);
    myEnds = find(myEnv_diff == -1);
    myStarts = find(myEnv_diff == 1);
    if isempty(myEnds) 
        myEnds = myStarts(1) + 6000;
    end
    if myEnds(1) < myStarts(1) 
        myEnds = myEnds(2:end);
    end
    if isempty(myEnds) 
        myEnds = myStarts(1) + 6000;
    end
    if myEnds(end) < myStarts(end) 
        myEnds = [myEnds; size(s1,1)];
    end
    for jj = 1:size(myStarts, 1)
        if (myEnds(jj) - myStarts(jj)) >= 1000 && (myEnds(jj) - myStarts(jj)) <= 6000
            numChirps = numChirps + 1;
            thisChirp = zeros(6000,1);
            chirpLength = myEnds(jj) - myStarts(jj) + 1;
            thisChirp(1:chirpLength) = s1(myStarts(jj):myEnds(jj)); 
            myChirps{numChirps} = thisChirp;  
            chirpDurations(numChirps) = chirpLength;
        end
    end       
end

save('AllChirps.mat', 'myChirps', 'chirpDurations')