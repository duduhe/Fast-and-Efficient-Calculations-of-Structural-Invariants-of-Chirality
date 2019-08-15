clear;clc;
close all;
%% Load Model data
disp('Model data loading...')
ptCloud_temp=load('01obj.mat');
ptCloud=ptCloud_temp.ptCloud;
clear ptCloud_temp;
surface.TRIV=ptCloud.f.v;
surface.X=ptCloud.v(:,1);
surface.Y=ptCloud.v(:,2);
surface.Z=ptCloud.v(:,3);
surface.I=ones(size(surface.Z,1),1)*(-1); 
disp('Finished');

%% Calculate chiral invariants
disp('Invariants calculating...')
try
    tempChiralInv=load('ChiralInv_Head.mat');
    ChiralInv=tempChiralInv.ChiralInv;
    clear tempChiralInv;
catch
    [ ChiralInv ] = ChiralInvariant( surface );
end

disp('Finished');

%% Visualization of S1-S5
disp('Visualization of S1-S5 : Processing...');
fViz_S(surface,ChiralInv);
disp('Finished');

%% symmetry plane detection (SPD)

%% SPD Method-1:Single Point Matching_Assigned Point
p_ISN=5942;% the point on ear in paper
disp(['SPD Method-1: Assigned Point ',num2str(p_ISN),' Processing...']);
fMatching_AssignedPoint(surface,ChiralInv,p_ISN);
disp('Finished');

%% SPD Method-1:Single Point Matching_All Points
disp('SPD Method-1: All Points Processing...');
fMatching_SinglePoint(surface,ChiralInv);
disp('Finished');

%% SPD Method-2:RandomSamplingMatching&PCA
samplingNum=1000;% the number of sampling points
disp(['SPD Method-2: ',num2str(samplingNum),' Points Processing...']);
fSPD_RSPCA(surface,ChiralInv,samplingNum);
disp('Finished');

%% PD Method-2:The Average of Rounds RandomSamplingMatching&PCA
samplingNum=200;% the number of sampling points
roundNum=10;% the num of round 
disp(['SPD Method-2-ave: ',num2str(samplingNum),' Points and ',num2str(roundNum),' Rounds Processing...']);
fSPD_RSPCA_Ave(surface,ChiralInv,samplingNum,roundNum);
disp('Finished');

%% SPD Method-2:RandomSamplingMatching&PCA
disp('SPD Method-2: Test Processing...');
StartNumofPoint=10; % n_start in paper
StepLength=10;
MaxNumofPoint=1000;  % n_end in paper
MaxRound=50;  % m in paper
fTest_RSPCA(surface,ChiralInv,StartNumofPoint,StepLength,MaxNumofPoint,MaxRound);

disp('THE END!');
