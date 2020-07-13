function beam = performBeamTraining_SM(beamTable,numRFRx)

% The same beam cannot be selected twice.
% The first pair is selected for one RF chain by default.

beam = zeros(numRFRx,2);
beam(1,:) = beamTable(1,:);

tempBeamTable = beamTable;
for k = 2:numRFRx
    tempBeamTable = tempBeamTable(logical(tempBeamTable(:,1)~=beam(k-1,1)),:);
    tempBeamTable = tempBeamTable(logical(tempBeamTable(:,2)~=beam(k-1,2)),:);
    beam(k,:) = tempBeamTable(1,:);
end
