% NAME:
% CreateHMMSongModel
% Creates a certain number of user-defined songs for a virtual bird whose 
%       songs are created by 
%
% REQUIRED INPUTS:
% numOfSongs            the number of songs the user wants generated 
% 
%
% OUTPUTS:
% myFullSong            a waveform which contains all generated songs. 
%
% numOfHiddenStates     the number of hidden states in the HMM model used
%                       to generate the syllable or chirp sequence 
%
% stateTransProb        the hidden state transition probability of the HMM
%
% stateTransCPD         the state transition CPD 
%
% obsStateProb          the observed state probability given a hidden state
%
% obsStateCPD           the observed state CPD 
%
% Fs                    the sampling frequency 
%
% sampMulti             multiplier slowing down playback frequency 
%            
% 
% NOTE:
% This function uses a data file call ChirpsByCluster.mat
%
% AUTHOR:
% Kiran D. Bhattacharyya (bhattacharyykiran12@gmail.com)
%
% License: MIT License 
%
% Example use: 
%           [myFullSong, numOfHiddenStates, stateTransProb, stateTransCPD,
%           obsStateProb, obsStateCPD, Fs, sampMulti] = CreateHMMSongModel(5);
%           sound(myFullSong, Fs/sampMulti); % play back generated songs 

function [myFullSong, numOfHiddenStates, stateTransProb, stateTransCPD, obsStateProb, obsStateCPD, Fs, sampMulti] = CreateHMMSongModel(numOfSongs)

load('Data/ChirpsByCluster.mat')

numOfClust = size(chirpsByCluster, 1);
myIncr = 1/numOfClust;

%% create random hidden markov model
numOfHiddenStates = round(13*rand(1) + 3); % generate number of hidden states 
stateTransProb = rand(numOfHiddenStates);
stateTransCPD = zeros(numOfHiddenStates);
for ii = 1:numOfHiddenStates % get state transition probability
    if ii == 1
        valMulti = 1/sum(stateTransProb(ii:end,ii));
        stateTransProb(ii:end,ii) = stateTransProb(ii:end,ii)*valMulti;
    else
        valMulti = (1 - sum(stateTransProb(1:ii-1,ii)))/sum(stateTransProb(ii:end,ii));
        stateTransProb(ii:end,ii) = stateTransProb(ii:end,ii)*valMulti;
    end
    valMulti = (1 - sum(stateTransProb(ii,1:ii)))/sum(stateTransProb(ii,ii+1:end));
    stateTransProb(ii,ii+1:end) = stateTransProb(ii,ii+1:end)*valMulti;
end
for ii = 1:numOfHiddenStates % get state transition CPD
    stateTransCPD(ii,:) = sum(stateTransProb(1:ii,:),1);
end
obsStateProb = zeros(numOfHiddenStates, numOfClust); % probability of a chirp given a state 
obsStateCPD = zeros(numOfHiddenStates, numOfClust); % cumulative probability at each chirp type 

for ii = 1:numOfHiddenStates % generate observed state probability 
    obsStateProb(ii,:) = rand(1, numOfClust);
    if sum(obsStateProb(ii,:)) > 1
        valMulti = 1/sum(obsStateProb(ii,:));
        obsStateProb(ii,:) = obsStateProb(ii,:)*valMulti;
    end
    for jj = 1:numOfClust
        obsStateCPD(ii,jj) = sum(obsStateProb(ii,1:jj));
    end
end

songEndProb = betarnd(3,3); % probability that a song will end after a syllable, sampled from beta distribution

%% create parameters controlling chirp duration and playback rate
Fs = 22050;
sampMulti = rand(1)*2 + 1.3; % change the fundamental frequency of the chirps by slowing playback rate  
timeDur = 0.009*randn(1) + 0.24; % time duration of a chirp 
chirpDur = round(timeDur/(1/Fs));
chirpDur_1 = round(chirpDur*0.95);
chirpDur_2 = chirpDur - chirpDur_1;
expDecay = exp(linspace(0, -6, chirpDur_2)); % add exponential decay to the end of chirps
myFullSong = 0; % initialize a variale to store chirps 

%% create the songs using the randomly generated hidden markov model and the parameters
%  equally likely to start in any hidden state 

for songNum = 1:numOfSongs % for each song asked by the user to generate 
    currState = ceil(rand(1)/(1/numOfHiddenStates)); % intialize state 

    listOfStates = [currState, 0];
    listOfChirps = [find(rand(1) < obsStateCPD(currState,:), 1), 0]; % intialize first chirp
    endState = 0;

    % enter while loop to create full song based on HMM
    iter = 1;
    while endState == 0 
        iter = iter + 1;
        currState = find(rand(1) < stateTransCPD(:,currState), 1);
        listOfStates(iter) = currState;
        listOfChirps(iter) = find(rand(1) < obsStateCPD(currState,:), 1);
        if rand(1) <= songEndProb
            endState = 1;
        end
    end

    % go through generated observed chirps and string them together 
    for ii = 1:size(listOfChirps, 2)
        thisChirpCluster = listOfChirps(ii); % get the chirp cluster 
        pickChirp = ceil(size(chirpsByCluster{thisChirpCluster}, 1)*rand(1)); % pick a chirp type from that cluster
        thisChirp = chirpsByCluster{thisChirpCluster}(pickChirp, 1:chirpDur); 
        thisDur = durationByCluster{thisChirpCluster}(pickChirp); % grab the duration of that chirp
        if thisDur >= chirpDur % apply exponential decay to avoid pops at the end of chirps 
            thisChirp(chirpDur_1:end) = thisChirp(chirpDur_1:end).*[expDecay 0];
        elseif thisDur < chirpDur 
            thisDur_1 = round(thisDur*0.95);
            thisDur_2 = thisDur - thisDur_1;
            thisDecay = exp(linspace(0,-6, thisDur_2));
            thisChirp(thisDur_1:thisDur) = thisChirp(thisDur_1:thisDur).*[thisDecay 0];
        end
        thisChirp = thisChirp*((randn(1)*0.2) + 1); % amplitude modulation for song to song variability 
        % pitch modulation
        pitchMod = (rand(1)*0.07) + 1; % multiple of pitch
        newLength = pitchMod*chirpDur; % pitch modulate by resampling the signal 
        newChirp = spline(1:1:chirpDur, thisChirp, linspace(1, chirpDur, newLength));
        if newLength < chirpDur
            newChirp = newChirp + zeros(1, chirpDur - size(newChirp, 2));
        elseif newLength > chirpDur
            newChirp = newChirp(1:chirpDur);
            newChirp(chirpDur_1:end) = newChirp(chirpDur_1:end).*[expDecay 0];
        end
        flipChirp = rand(1) > betarnd(7,7); % flip the chirp with some probability 
        if flipChirp
            thisChirp = fliplr(thisChirp);
        end
        myFullSong = [myFullSong; thisChirp'];
    end
    % add white noise the whole track to add "realness" of actual recording
    myFullSong = myFullSong + (wgn(size(myFullSong, 1), size(myFullSong, 2), 0.01)/2000);
    fillerLength = randn(1)*0.3 + 3; % define the silent space between songs 
    filler = wgn(round(fillerLength*chirpDur), 1, 0.01); % create white noise for silent space 
    filler = filler/2000;
    myFullSong = [myFullSong; filler];
end
    