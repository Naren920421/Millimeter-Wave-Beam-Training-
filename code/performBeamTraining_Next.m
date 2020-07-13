function [beamTable,rPower,beamPair,beamIdx_BS,beamIdx_UE,numTB,numRB,numBM] = performBeamTraining_Next(num,beamPrevious,beamBS,beamAzBS,beamElBS,beamUE,beamAzUE,beamElUE,freqChan,pR)

numTB = (2*num+1)^2; % 9 or 25
numRB = (2*num+1)^2; % 9 or 25
numBM = numTB*numRB; % 81 or 625

% Store beam selection per location
bestAz = zeros(1,2);
bestEl = zeros(1,2);

% Store best beam at location 1
% BS
bestEl(1,1) = ceil(beamPrevious(1,1)./size(beamAzBS,2));
bestAz(1,1) = mod(beamPrevious(1,1),size(beamAzBS,2));
if bestAz(1,1) == 0
    bestAz(1,1) = size(beamAzBS,2);
end
% UE
bestEl(1,2) = ceil(beamPrevious(1,2)./size(beamAzUE,2));
bestAz(1,2) = mod(beamPrevious(1,2),size(beamAzUE,2));
if bestAz(1,2) == 0
    bestAz(1,2) = size(beamAzUE,2);
end

azBeamIdx_BS = getLocalBeamIdx(num,size(beamAzBS,2), bestAz(1,1));
elBeamIdx_BS = getLocalBeamIdx(num,size(beamElBS,2), bestEl(1,1));
beamIdx = (elBeamIdx_BS-1)*size(beamAzBS,2)+azBeamIdx_BS.';
beamIdx_BS = sort(beamIdx(:),'ascend'); % 9 or 25
codebook_BS = beamBS(:,beamIdx_BS.');

% UE
azBeamIdx_UE = getLocalBeamIdx(num,size(beamAzUE,2), bestAz(1,2));
elBeamIdx_UE = getLocalBeamIdx(num,size(beamElUE,2), bestEl(1,2));
beamIdx = (elBeamIdx_UE-1)*size(beamAzBS,2)+azBeamIdx_UE.';
beamIdx_UE = sort(beamIdx(:),'ascend'); % 9 or 25
codebook_UE = beamUE(:,beamIdx_UE.');

[beamTable,rPower,beamPair,~,~,~] = performBeamTraining(codebook_BS,beamIdx_BS,codebook_UE,beamIdx_UE,freqChan,pR);

end

function localBeamIdx = getLocalBeamIdx(numLocalBeam,numTotalBeam,centerBeamIdx)

beamIdxVec = 1:numTotalBeam;
beamIdxVec = [beamIdxVec(numTotalBeam-numLocalBeam+1:end) beamIdxVec beamIdxVec(1:numLocalBeam)];
centerBeamIdx_new = centerBeamIdx+numLocalBeam;
localBeamIdx = beamIdxVec(centerBeamIdx_new-numLocalBeam:centerBeamIdx_new+numLocalBeam);
localBeamIdx = unique(localBeamIdx,'stable');

end