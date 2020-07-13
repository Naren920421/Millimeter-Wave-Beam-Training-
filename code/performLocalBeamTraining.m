function [beamTable,rPower,beamPair,beamIdx_BS,beamIdx_UE,numTB,numRB,numBM] = performLocalBeamTraining(num,beamPrevious,beamBS,beamAzBS,beamElBS,beamUE,beamAzUE,beamElUE,freqChan,pR)

T = length(pR); % 1 x T
numTB = (2*num+1)^2; % 9 or 25
numRB = (2*num+1)^2; % 9 or 25
numBM = numTB*numRB; % 81 or 625

% Store beam selection per location
bestAz = zeros(1,2,T);
bestEl = zeros(1,2,T);

% Best beam at location 1
% BS
bestEl(1,1,1) = ceil(beamPrevious(1,1)./size(beamAzBS,2));
bestAz(1,1,1) = mod(beamPrevious(1,1),size(beamAzBS,2));
if bestAz(1,1,1) == 0
    bestAz(1,1,1) = size(beamAzBS,2);
end
% UE
bestEl(1,2,1) = ceil(beamPrevious(1,2)./size(beamAzUE,2));
bestAz(1,2,1) = mod(beamPrevious(1,2),size(beamAzUE,2));
if bestAz(1,2,1) == 0
    bestAz(1,2,1) = size(beamAzUE,2);
end

beamTable = zeros(numBM,2,T);
rPower = zeros(numTB,numRB,T);
beamPair = zeros(1,2,T);
beamIdx_BS = zeros(numTB,T);
beamIdx_UE = zeros(numRB,T);

for t = 2:T

    % Local codebook index
    % BS
    azBeamIdx_BS = getLocalBeamIdx(num,size(beamAzBS,2), bestAz(1,1,t-1));
    elBeamIdx_BS = getLocalBeamIdx(num,size(beamElBS,2), bestEl(1,1,t-1));
    beamIdx = (elBeamIdx_BS-1)*size(beamAzBS,2)+azBeamIdx_BS.';
    beamIdx_BS(:,t) = sort(beamIdx(:),'ascend'); % 9 or 25
    codebook_BS = beamBS(:,beamIdx_BS(:,t).');
    
    % UE
    azBeamIdx_UE = getLocalBeamIdx(num,size(beamAzUE,2), bestAz(1,2,t-1));
    elBeamIdx_UE = getLocalBeamIdx(num,size(beamElUE,2), bestEl(1,2,t-1));
    beamIdx = (elBeamIdx_UE-1)*size(beamAzBS,2)+azBeamIdx_UE.';
    beamIdx_UE(:,t) = sort(beamIdx(:),'ascend'); % 9 or 25
    codebook_UE = beamUE(:,beamIdx_UE(:,t).');
    
    [beamTable(:,:,t),rPower(:,:,t),beam,~,~,~] = performBeamTraining(codebook_BS,beamIdx_BS(:,t),codebook_UE,beamIdx_UE(:,t),freqChan(:,:,:,t),pR(t));

    % Azimuth and elevation beam index (for checking)
    elBS = ceil(beam(1)/size(beamAzBS,2));
    azBS = mod(beam(1),size(beamAzBS,2));
    azBS(logical(azBS)==0) = size(beamAzBS,2);
    
    elUE = ceil(beam(2)/size(beamAzUE,2));
    azUE = mod(beam(2),size(beamAzUE,2));
    azUE(logical(azUE)==0) = size(beamAzUE,2);

    bestAz(1,1,t) = azBS;
    bestAz(1,2,t) = azUE;
    bestEl(1,1,t) = elBS;
    bestEl(1,2,t) = elUE;
    
    beamPair(:,:,t) = beam;
    
end

beamTable(:,:,1) = [];
rPower(:,:,1) = [];
beamPair(:,:,1) = []; 
beamIdx_BS(:,1) = [];
beamIdx_UE(:,1) = [];

end

function localBeamIdx = getLocalBeamIdx(numLocalBeam,numTotalBeam,centerBeamIdx)

beamIdxVec = 1:numTotalBeam;
beamIdxVec = [beamIdxVec(numTotalBeam-numLocalBeam+1:end) beamIdxVec beamIdxVec(1:numLocalBeam)];
centerBeamIdx_new = centerBeamIdx+numLocalBeam;
localBeamIdx = beamIdxVec(centerBeamIdx_new-numLocalBeam:centerBeamIdx_new+numLocalBeam);
localBeamIdx = unique(localBeamIdx,'stable');

end