function dataRate = computeDataRate(beamSelection,freqChan,beamBS,beamUE,numRFRx,SNR_linear)

[~,~,numSC,T] = size(freqChan);

dataRateTemp = zeros(numSC,T);
for t = 1:T
    
    f = beamBS(:,beamSelection(:,1,t).');
    w = beamUE(:,beamSelection(:,2,t).');

    for sc = 1:numSC
        H = freqChan(:,:,sc,t);
        dataRateTemp(sc,t) = abs(log2(det(eye(numRFRx)+SNR_linear(t)/numRFRx*w'*H*f*f'*H'*w)));
    end
end

dataRate = mean(dataRateTemp,1); % averaged over subcarriers