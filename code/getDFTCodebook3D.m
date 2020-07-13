function [beam3D,beamAz,beamEl,beamAngleAz,beamAngleEl] = getDFTCodebook3D(numElementV,numElementH,eleSpacingV,eleSpacingH)

numAngle = numElementV*numElementH; % total number of resolvable angles

% Vertical - resolved elevation angles
virtAngleEl = (0:numElementV-1)./numElementV;
av = unique(real(asind(-virtAngleEl./eleSpacingV)),'stable');
bv = flip(-av(1:end-1));
beamAngleEl = unique([av bv],'stable');
beamEl = 1/sqrt(numElementV)*exp(1j*2*pi*((0:numElementV-1)-numElementV/2).'*virtAngleEl);

% Horizontal - resolved azimuth angles
virtAngleAz = (0:numElementH-1)./numElementH;
ah = unique(real(asind(-virtAngleAz./eleSpacingH)),'stable');
bh = flip(-ah(1:end-1));
beamAngleAz = unique([ah bh],'stable');
beamAz = 1/sqrt(numElementH)*exp(1j*2*pi*((0:numElementH-1)-numElementH/2).'*virtAngleAz); 

beam3D = zeros(numAngle,numAngle);
beamAngle = zeros(2,numAngle);
idx = 1;
for nv = 1:numel(beamAngleEl)
    for nh = 1:numel(beamAngleAz)
        beam2D = flip(beamEl(:,nv)).*beamAz(:,nh).'; % numV x numH
        beam3D(:,idx) = beam2D(:); % numElement x 1 [numV x 1;numV x 1;...]
        beamAngle(:,idx) = [beamAngleEl(nv);beamAngleAz(nh)];
        idx = idx+1;
    end
end
