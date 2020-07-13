% This script is to implement beam training strategies proposed in 
% Section III-C in [1] and generate similar simulation results in Fig. 6.
% The mixed beam search is only proposed for 1 RF chain without spatial
% multiplexing. Readers are advised to implement the 3GPP spatial 
% consistency Prodecure A in [2] to reproduce the simulation results in 
% [1]. The main purpose of this script is to demonstrate the 
% implementation of beam training process. The beam training process is 
% assumed to be noise-free. 
%
% [1] Narengerile,F.Alsaleem,J.Thomspon,T.Ratnarajah,"Low-Complexity Beam 
% Training for Tracking Spatially Consistent Millimeter Wave Channels",
% PIMRC, 2020.
% [2] 3GPP TR 38.901,"Study on channel model for frequencies from 0.5 to 
% 100 GHz," 2017.

%% Simulation parameters

clear;
close all;

% For repeatibility
stream = RandStream('mt19937ar','Seed',3);
RandStream.setGlobalStream(stream);

% System parameters
fc = 30e9; % carrier frequency, Hz
BW = 100e6; % bandwidth, Hz
speedUE = 1; % m/sec
updatePeriod = 0.5; % channel update period, second
timeEpoch = 30; % second
movingDis = speedUE*timeEpoch;
T = timeEpoch/updatePeriod; % number of locations
directionUE = -90; % horizontal
numSC = 64; % number of subcarriers 
numRFBS = 1; % number of BS RF chains
numRFUE = 1; % number of UE RF chains
nVar_dBm = -174+10*log10(BW)-10*log10(numSC); % noise power per subcarrier
nVar = 10^(nVar_dBm/10)*1e-3; % noise variance, watt
pT = 0.001; % transmit power, watt

disBSUE2D = 50; % BS-UE distance, meter
heightBS = 25; % meter
heightUE = 1.5; % meter
arraySizeBS = [8,8]; 
arraySizeUE = [8,8];
Nt = prod(arraySizeBS);
Nr = prod(arraySizeUE);
N = 20; % number of spatial clusters
velocityUE = speedUE*[cosd(directionUE);sind(directionUE);0];
velocityUE = repmat(velocityUE,1,T);
cellLayout = getCellLayout(disBSUE2D,heightBS,heightUE); % drop UE in the cell
locationBS = cellLayout.locationBS;
locationUE_t0 = cellLayout.locationUE;
locationUE_update = zeros(3,T);
locationUE_update(:,1) = locationUE_t0;
for t = 2:T
    locationUE_update(:,t) = velocityUE(:,t)*updatePeriod+locationUE_update(:,t-1);
end

% UE trajectory
figure();
dpUE = locationUE_update(1:2,end)-locationUE_update(1:2,1);
p1 = plot(locationBS(1),locationBS(2),'rpentagram','MarkerFaceColor','r','MarkerSize',20);hold on;
p2 = plot(locationUE_update(1,1),locationUE_update(2,1),'o','MarkerFaceColor','b','MarkerSize',12);hold on;
p3 = plot(locationUE_update(1,end),locationUE_update(2,end),'d','MarkerFaceColor','b','MarkerSize',12);hold on;
quiver(locationUE_update(1,1),locationUE_update(2,1),dpUE(1),dpUE(2),0,'MaxHeadSize',0.6);hold on;
dp1 = [0,20];dp2 = [20*cosd(30);-20*sind(30)];dp3 = [-20*cosd(30);-20*sind(30)];
quiver(0,0,dp1(1),dp1(2),0,'ShowArrowHead','off','LineStyle','--','Color','k');hold on;
quiver(0,0,dp2(1),dp2(2),0,'ShowArrowHead','off','LineStyle','--','Color','k');hold on;
quiver(0,0,dp3(1),dp3(2),0,'ShowArrowHead','off','LineStyle','--','Color','k');hold on;
xlim([-100,100]);ylim([-100,100]);grid on;

% DFT codebooks
[beamBS,beamAzBS,beamElBS,~,~] = getDFTCodebook3D(arraySizeBS(1),arraySizeBS(2),0.5,0.5);
[beamUE,beamAzUE,beamElUE,~,~] = getDFTCodebook3D(arraySizeUE(1),arraySizeUE(2),0.5,0.5);
numBeamBS = size(beamBS,2);
numBeamUE = size(beamUE,2);
numBeamAzBS = size(beamAzBS,2);
numBeamElBS = size(beamElBS,2);
numBeamAzUE = size(beamAzUE,2);
numBeamElUE = size(beamElUE,2);

% Blocker type
blockProb = 0.1;
blockDuration = 2; % number of samples blocked
blockCluster = zeros(N,T);
blockCluster(:,2:end) = blockage(T-1,N,blockDuration,blockProb);
blockCluster = logical(blockCluster);
% Without blockages
% blockCluster = false(N,T);

%% Channel matrices

% Load channel matrices for 1 Monte-Carlo run
load('chan.mat','chanMatrix');
chanMatrix_bl = chanMatrix; % applying blockages
for n = 1:N
    chanMatrix_bl(:,:,n,blockCluster(n,:)) = chanMatrix(:,:,n,blockCluster(n,:))*sqrt(0.01); % attenuated by 20 dB
end

% Frequency-domain channel
freqChan = zeros(Nr,Nt,numSC,T);
for t = 1:T
    H = chanMatrix_bl(:,:,:,t); % Nr x Nt x N
    for sc = 1:numSC
        chanPerSC = zeros(Nr,Nt);
        for n = 1:N
            chanPerSC = chanPerSC+1/sqrt(numSC)*H(:,:,n)*exp(-1j*2*pi*sc/numSC*(n-1)); % maintain the same power
        end
        freqChan(:,:,sc,t) = chanPerSC;
    end
end

% Calculate NLOS path loss
disBSUE3D_update = vecnorm(locationUE_update-repmat(locationBS,1,T),2); % 1 x T, meter
pathLoss = 32.4+20*log10(fc*1e-9)+30*log10(disBSUE3D_update); % dB, 1 x T
pR = pT./10.^(pathLoss./10); % average receive power, watt
SNR_linear = pR./nVar; % actual SNR per subcarrier

%% Beam training - mixed beam search

dr_Mixed = zeros(1,T);
drGap = 0.5;
beamTrain = zeros(1,T); % 0 - no beam training, 1 - Local Search 1, 2 - Local Search 2, 3 - Exhaustive Search

% 1st location at time t = 1
t_1 = 1;
[~,~,beam_Ex,~,~,~] = performBeamTraining(beamBS,1:numBeamBS,beamUE,1:numBeamUE,freqChan(:,:,:,t_1),pR(t_1)); 
dr_1 = computeDataRate(beam_Ex,freqChan(:,:,:,t_1),beamBS,beamUE,numRFUE,SNR_linear(t_1));
dr_Mixed(t_1) = dr_1;

% From time t = 2 to t = T
beam = repmat(beam_Ex,1,1,T); % assume always using the same beam pair and then replace with the selected beam pair
for t = 2:T
    dr = computeDataRate(beam(:,:,t),freqChan(:,:,:,t),beamBS,beamUE,numRFUE,SNR_linear(t));
    dr_Mixed(t) = dr;
    if dr_Mixed(t-1)-dr_Mixed(t) > drGap
        num = 1; % Local Search 1
        [~,~,beamTemp,~,~,~,~,~] = performBeamTraining_Next(num,beam(:,:,t-1),beamBS,beamAzBS,beamElBS,beamUE,beamAzUE,beamElUE,freqChan(:,:,:,t),pR(t));
        dr = computeDataRate(beamTemp,freqChan(:,:,:,t),beamBS,beamUE,numRFUE,SNR_linear(t));
        dr_Mixed(t) = dr;   
        beam(:,:,t:end) = repmat(beamTemp,1,1,(T-t)+1); % assume always using this beam pair and then replace with the selected beam pair
        beamTrain(t) = 1;      
        if dr_Mixed(t-1)-dr_Mixed(t) > drGap
            num = 2; % Local Search 2
            [~,~,beamTemp,~,~,~,~,~] = performBeamTraining_Next(num,beam(:,:,t-1),beamBS,beamAzBS,beamElBS,beamUE,beamAzUE,beamElUE,freqChan(:,:,:,t),pR(t));
            dr = computeDataRate(beamTemp,freqChan(:,:,:,t),beamBS,beamUE,numRFUE,SNR_linear(t));
            dr_Mixed(t) = dr;   
            beam(:,:,t:end) = repmat(beamTemp,1,1,(T-t)+1);
            beamTrain(t) = 2;      
            if dr_Mixed(t-1)-dr_Mixed(t) > drGap % Exhaustive Search
                [~,~,beamTemp,~,~,~] = performBeamTraining(beamBS,1:numBeamBS,beamUE,1:numBeamUE,freqChan(:,:,:,t),pR(t)); 
                dr = computeDataRate(beamTemp,freqChan(:,:,:,t),beamBS,beamUE,numRFUE,SNR_linear(t));
                dr_Mixed(t) = dr;   
                beam(:,:,t:end) = repmat(beamTemp,1,1,(T-t)+1);
                beamTrain(t) = 3;      
            end
        end
    end        
end

% For comparison - Exhaustive Search
[~,~,beam_compare,~,~,~] = performBeamTraining(beamBS,1:numBeamBS,beamUE,1:numBeamUE,freqChan,pR); 
dr_Ex = computeDataRate(beam_compare,freqChan,beamBS,beamUE,numRFUE,SNR_linear);

figure();
plot(dr_Mixed,'b-d');hold on;
plot(dr_Ex,'r-x');hold off;
legend('1 RF - Mixed Search','1 RF - Exhaustive Search');

%% Average number of beam measurements

beamTrainIdx = 0:3; % beam training method index
probSt = zeros(1,numel(beamTrainIdx));
for p = 1:numel(beamTrainIdx)
    probSt(p) = sum(logical(beamTrain==beamTrainIdx(p)))/T;       
end
numActivation = probSt.*T; % average number of activations per beam training method

numBM = zeros(1,T); % number of beam measurements required
numBM_perStage = [1,1+9^2,1+9^2+25^2,1+9^2+25^2+64^2];
for p = 1:numel(beamTrainIdx)   
    numBM(logical(beamTrain==beamTrainIdx(p))) = numBM_perStage(p);
end
aveNumBM = mean(numBM_perStage.*probSt); % average number of beam measurements per location
