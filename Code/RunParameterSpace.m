%% Generate and save bird songs by sampling the parameter space of 
%  differential equations 

clear all

myRhos = -13:0.1:-5.5;

for ii = 1:size(myRhos, 2)
    tic;
    disp(['Performing bird syllable synthesis ' num2str(ii) ' of ' num2str(size(myRhos, 2)) '.'])
    [m1 s1]=singSyllable_3(myRhos(ii));
    thisRho = -myRhos(ii);
    save(['BirdSongSims_3/birdSim_' num2str(floor(thisRho)) '_' num2str(10*(thisRho - floor(thisRho))) '_paperParams.mat'], 's1')
    disp(['Synthesis took ' num2str(toc) ' seconds.'])
end
    