function beamPair_SM = performLocalBeamTraining_SM(num,beamPrevious,beamBS,beamAzBS,beamElBS,beamUE,beamAzUE,beamElUE,numRFUE,freqChan,pR)

T = length(pR); % 1 x T
numTB = (2*num+1)^2; % 9 or 25
numRB = (2*num+1)^2; % 9 or 25

% Previous beampair information
elPrevious = zeros(numRFUE,2,T);
azPrevious = zeros(numRFUE,2,T);

elPrevious(:,1,1) = ceil(beamPrevious(:,1)./size(beamAzBS,2));
azPrevious(:,1,1) = mod(beamPrevious(:,1),size(beamAzBS,2));
azPrevious(logical(azPrevious(:,1,1)==0),1,1) = size(beamAzBS,2);

elPrevious(:,2,1) = ceil(beamPrevious(:,1)./size(beamAzUE,2));
azPrevious(:,2,1) = mod(beamPrevious(:,1),size(beamAzUE,2));
azPrevious(logical(azPrevious(:,2,1)==0),2,1) = size(beamAzUE,2);

beamPair_SM = zeros(numRFUE,2,T);
beamIdx_BS = zeros(numTB,T);
beamIdx_UE = zeros(numRB,T);
for t = 2:T

    % 1st beampair
    azBeamIdx_BS = getLocalBeamIdx(num,size(beamAzBS,2), azPrevious(1,1,t-1));
    elBeamIdx_BS = getLocalBeamIdx(num,size(beamElBS,2), elPrevious(1,1,t-1));
    beamIdx = (elBeamIdx_BS-1)*size(beamAzBS,2)+azBeamIdx_BS.';
    beamIdx_BS(:,t) = sort(beamIdx(:),'ascend'); % 9 or 25
    codebook_BS = beamBS(:,beamIdx_BS(:,t).');
    
    azBeamIdx_UE = getLocalBeamIdx(num,size(beamAzUE,2), azPrevious(1,2,t-1));
    elBeamIdx_UE = getLocalBeamIdx(num,size(beamElUE,2), elPrevious(1,2,t-1));
    beamIdx = (elBeamIdx_UE-1)*size(beamAzBS,2)+azBeamIdx_UE.';
    beamIdx_UE(:,t) = sort(beamIdx(:),'ascend'); % 9 or 25
    codebook_UE = beamUE(:,beamIdx_UE(:,t).');
    
    [~,~,beam1st,~,~,~] = performBeamTraining(codebook_BS,beamIdx_BS(:,t),codebook_UE,beamIdx_UE(:,t),freqChan(:,:,:,t),pR(t));
    beamPair_SM(1,:,t) = beam1st;
    
    % 2nd beampair
    azBeamIdx_BS = getLocalBeamIdx(num,size(beamAzBS,2), azPrevious(2,1,t-1));
    elBeamIdx_BS = getLocalBeamIdx(num,size(beamElBS,2), elPrevious(2,1,t-1));
    beamIdx = (elBeamIdx_BS-1)*size(beamAzBS,2)+azBeamIdx_BS.';
    beamIdx_BS(:,t) = sort(beamIdx(:),'ascend'); % 9 or 25
    codebook_BS = beamBS(:,beamIdx_BS(:,t).');
    
    azBeamIdx_UE = getLocalBeamIdx(num,size(beamAzUE,2), azPrevious(2,2,t-1));
    elBeamIdx_UE = getLocalBeamIdx(num,size(beamElUE,2), elPrevious(2,2,t-1));
    beamIdx = (elBeamIdx_UE-1)*size(beamAzBS,2)+azBeamIdx_UE.';
    beamIdx_UE(:,t) = sort(beamIdx(:),'ascend'); % 9 or 25
    codebook_UE = beamUE(:,beamIdx_UE(:,t).');
    
    [beamTable,~,beam2nd,~,~,~] = performBeamTraining(codebook_BS,beamIdx_BS(:,t),codebook_UE,beamIdx_UE(:,t),freqChan(:,:,:,t),pR(t));

    % Remove same beam
    if beam2nd(1) == beam1st(1)
        beamTable = beamTable(logical(beamTable(:,1)~=beam1st(1)),:);
    end
    if beam2nd(2) == beam1st(2)
        beamTable = beamTable(logical(beamTable(:,2)~=beam1st(2)),:);
    end
    beam2nd = beamTable(1,:); % 1 x 2
    beamPair_SM(2,:,t) = beam2nd;

    beam = [beam1st;beam2nd];
    elBS = ceil(beam(:,1)/size(beamAzBS,2));
    azBS = mod(beam(:,1),size(beamAzBS,2));
    azBS(logical(azBS)==0) = size(beamAzBS,2);
    
    elUE = ceil(beam(:,2)/size(beamAzUE,2));
    azUE = mod(beam(:,2),size(beamAzUE,2));
    azUE(logical(azUE)==0) = size(beamAzUE,2);

    azPrevious(:,:,t) = [azBS azUE];
    elPrevious(:,:,t) = [elBS elUE];
        
end

end

function localBeamIdx = getLocalBeamIdx(numLocalBeam,numTotalBeam,centerBeamIdx)

beamIdxVec = 1:numTotalBeam;
beamIdxVec = [beamIdxVec(numTotalBeam-numLocalBeam+1:end) beamIdxVec beamIdxVec(1:numLocalBeam)];
centerBeamIdx_new = centerBeamIdx+numLocalBeam;
localBeamIdx = beamIdxVec(centerBeamIdx_new-numLocalBeam:centerBeamIdx_new+numLocalBeam);
localBeamIdx = unique(localBeamIdx,'stable');

end