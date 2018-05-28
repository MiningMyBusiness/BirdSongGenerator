%% Generate and save bird songs by sampling the parameter space of 
%  differential equations 

clear all

myRhos = -13:0.1:-5.5; % a close look at Figure 2a in Laje and Mindlin(2002), Diversity within a Birdsong, Physical Review
                       % will explain the choice of this range for the parameter

for ii = 1:size(myRhos, 2)
    tic;
    disp(['Performing bird syllable synthesis ' num2str(ii) ' of ' num2str(size(myRhos, 2)) '.'])
    [m1 s1]=singSyllable(myRhos(ii));
    thisRho = -myRhos(ii);
    save(['BirdSongSims/birdSim_' num2str(floor(thisRho)) '_' num2str(10*(thisRho - floor(thisRho))) '_paperParams.mat'], 's1')
    disp(['Synthesis took ' num2str(toc) ' seconds.'])
end
    
