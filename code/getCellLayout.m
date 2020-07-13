function cellLayout = getCellLayout(disBSUE2D,heightBS,heightUE)

numUser = numel(disBSUE2D);
disBSUE3D = sqrt((heightBS-heightUE)^2+disBSUE2D.^2); 

% Direct path angles
% Zenith angle [0,180] degree
% Azimuth angle [-180,180] degree
ZOA = atand(disBSUE2D./(heightBS-heightUE));
ZOD = 180-ZOA;
AOD = (rand(1,numUser)-0.5)*360; 
AOA = zeros(1,numUser);
AOA(logical(AOD>0)) = AOD(logical(AOD>0))-180;
AOA(logical(AOD<0)) = AOD(logical(AOD<0))+180;

% User location
xUE = disBSUE2D.*cosd(AOD);
yUE = disBSUE2D.*sind(AOD); 
zUE = heightUE*ones(1,numUser); 
locationUE = [xUE;yUE;zUE];

% BS location
locationBS = [0;0;heightBS];

cellLayout = struct('disBSUE3D',disBSUE3D,...
    'locationBS',locationBS,...
    'locationUE',locationUE,...
    'ZOA',ZOA,...
    'ZOD',ZOD,...
    'AOD',AOD,...
    'AOA',AOA);
