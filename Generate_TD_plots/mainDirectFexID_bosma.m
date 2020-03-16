% Default Main Code: The purpose of this script is to be the single
% executable script to run an example. 
% 
% As of right now the script has the mesh, nemoh, IRFs, and experimental
% data this script also sets up to simulate the FPS which encounters the
% wave first.

clear all
close all
clc

tic
rng(1)

indicies = [8,8];
for j = indicies(1):indicies(2)
%% Initialization and Data Generation
% generate path to needed folders
Paths.exampleDir = cd;
Paths.functionsDir = fullfile(Paths.exampleDir, '..', filesep,'..', filesep, 'Functions'); 
Paths.dataDir = fullfile(Paths.exampleDir, '..', filesep,'..', filesep, 'Data'); 
% Paths.functionsDir = fullfile(Paths.exampleDir, '..', filesep, 'Functions'); % original path for Dev/Test
Paths.outputDir = fullfile(Paths.exampleDir, filesep, 'Outputs'); % original path for Dev/Test
% Paths.dataDir = fullfile(Paths.exampleDir, '..', filesep, 'Data'); % original path for Dev/Test

% Generate GLOBAL VARIABLES
global Param Time Func Data SpectrumData Est SmoothCounter

% Simulation and data choices
name.body = {'FPS','OWC'};
bodyFunc = @generateBody_BosmaFPS; % this should be unique to the application
% bodyString: decide whether you want to run NEMOH
bodyString = 'check'; % options: skip, run, or check
% dataString decide wheather youy want to read bosmaData or load
dataString = 'check'; %  options: read, check
% measurements are the values evaluated by the EKF. These can be
% experiemental or numerical
savedDataName = 'RefSim_JP1.mat';
measString = 'loadExper'; % options: loadSim, loadExper, sim
% decide if you want to generate the functions
genFunString = 'gen'; % options: gen, check
bodyNumber = 1; % 1-> FPS

waveBool = [1]; % options: 1 -> generate wave, 0-> no wave, 2 -> regular wave with unit amplitude
mooringBool = [0]; % 0-> taught & 1-> catenary
RTSBool = [0]; % 0-> Do NOT implement RTS, 1-> run RTS Smoother
% this is used in the estimateVel function
Param.velEstimateVal = 1; % 0-> uses particleDynamics_spec 1-> uses method 1 2-> method 2...

% generate MATLAB Path
% addedPath = genpath('..\Functions');
addpath(Paths.functionsDir)
addpath(Paths.dataDir)
addpath('Outputs')
[status, msg, msgID] = mkdir('Outputs');


%% Enter Parameters
Param.numBodies = length(name.body);
Param.RTSBool = RTSBool;

% Choosing specific data
% Param.runNo = '176'; Param.wavePeriod = 2.69; Param.waveHs = 0.3; % regular wave with O-ring
% Param.runNo = '125'; Param.wavePeriod = 2.69; Param.waveHs = 0.3; %irregular wave TFI 30%
% Param.runNo = '113'; Param.wavePeriod = 3; Param.waveHs = 0.2; % regular wave TFI 30%
% Param.runNo = '099'; Param.wavePeriod = 2.06; Param.waveHs = 0.1; % regular wave TFI 30%
% Param.runNo = '117'; Param.wavePeriod = 4.74; Param.waveHs = 0.2; % regular wave TFI 30%
% Param.runNo = '162'; Param.wavePeriod = 3; Param.waveHs = 0.2; % regular wave O ring
% Param.runNo = '122'; Param.wavePeriod = 2.21; Param.waveHs = 0.3; % regular TFI 30%
% Param.runNo = '173'; Param.wavePeriod = 2.21; Param.waveHs = 0.3; % regular O-ring
% Param.runNo = '133'; Param.wavePeriod = 1.26; Param.waveHs = 0.1; % regular o ring
% oRingIndex = 11;
oRingIndex = j
waveRuns = [1,180,0.075,1.1,0.0081,0.0549;
            2,181,0.175,1.42,0,0;
            3,182,0.05,1.9,0,0;
            4,183,0.175,1.9,0,0;
            5,184,0.3,1.9,0,0;
            6,185,0.050,2.37,0,0;
            7,186,0.175,2.37,0,0;
            8,199,0.3,2.37,0,0;
            9,189,0.05,2.69,0,0;
            10,190,0.175,2.69,0,0;
            11,191,0.3,2.69,0.0171,0.0845;
            12,192,0.05,3.16,0,0;
            13,193,0.175,3.16,0,0;
            14,194,0.175,3.95,0,0;
            15,195,0.175,4.74,0,0;
            16,196,0.175,5.53,0,0;
            17,197,0.375,2.06,0,0;
            18,198,0.375,2.37,0,0];
Param.runNo = num2str(waveRuns(oRingIndex,2)); Param.wavePeriod = waveRuns(oRingIndex,4); Param.waveHs = waveRuns(oRingIndex,3);
load runStruct.mat

name.sensorFile = strcat('Run_', Param.runNo, '.csv');
name.waveFile = strcat('Run_', Param.runNo, '.txt');
name.motionFile = strcat('Run_', Param.runNo, '_6D.tsv'); 
name.runNo = Param.runNo;

Param.samplingFreq.motion = 200;
Param.samplingFreq.wave = 128;
Param.samplingFreq.loadCell = 128;
Param.peakFrequency = 1/Param.wavePeriod;
Param.waveLength = 9.81*Param.wavePeriod/(2*pi);

Time.subSampleNo = 1; % uses one out of every subSampleNo
Time.dT = Time.subSampleNo/Param.samplingFreq.motion;
Time.simStepsPerDT = 50;
% Time.simStepsPerDT = 10
Time.dtSim = Time.dT/Time.simStepsPerDT;
Time.start =50; % min 1 period into timespan
Time.end = 110;% max 1 period from the end of the timespan
Time.t = [Time.start:Time.dT:Time.end];

% Number of states and control inputs
Param.sizeRadSS = 3;
Param.numIDStates = 1;
Param.numStates = 2+Param.numIDStates+Param.sizeRadSS; % 2 for dynamics 
Param.numOutputs = 2; 
Param.numNonRad = Param.numStates-Param.sizeRadSS;

% initialize EstStates
% Both of these global variables are modified during the smoother_RTS funciton
Est.states = zeros(Param.numStates,length(Time.t));
Est.vel = zeros(1,length(Time.t)); % updated in the estimateVel function
SmoothCounter = 0;


% Param structure data Properties of the experiment
Param.linRad = 20;
Param.rho = 1000;
Param.g = 9.81;
Param.mass = 9.582; % 1/10 scale FPS
Param.upperPanelLength = 0.51/(sqrt(2)+1); % originally this was 0.50
Param.lowerPanelLength = 0.35/(sqrt(2)+1);
Param.avgPanelLength = (Param.upperPanelLength+Param.lowerPanelLength)/2;
Param.lowerArea = 2*(1+sqrt(2))*Param.lowerPanelLength^2;
Param.upperArea = 2*(1+sqrt(2))*Param.upperPanelLength^2;
Param.avgArea = 2*(1+sqrt(2))*Param.avgPanelLength^2;
Param.buoyDepth = 0.075; % specified as depth from surface
Param.waveNo = Param.peakFrequency*2*pi; % deep water assumption implicit
Param.attenFactor = exp(-Param.buoyDepth*Param.waveNo); % calculated on peak frequency
Param.attenFactor = 1; % this removes the attenfactor influence while keeping it for further implementation

Param.addedMass = 21.15;
Param.exciteTerm = 1800;
Param.weight = Param.mass*Param.g;
Param.totalInertia = Param.mass+Param.addedMass;
% Param.hydrostaticCoeff = Param.avgArea*Param.rho*Param.g;
Param.hydrostaticCoeff = Param.upperArea*Param.rho*Param.g;
% Param.heaveMoorStiff = 700; % N/m TFI 30%
% Param.heaveMoorStiff = 650; % N/m
Param.heaveMoorStiff = 200;
% Param.heaveMoorStiff = 100;

% ScaledHs = 10*Param.waveHs
% ScaledTs = sqrt(10)*Param.wavePeriod
Param.dispVol = 11322077e-9; % from solidworks for the FPS
Param.Delta = Param.dispVol*Param.rho;
Param.C_D = 0.4;
Param.dragFactor = Param.C_D*1/2*Param.rho*Param.upperArea;
Param.linearDrag = 60;



% Initial Condtions 
% Param.IC.est = [0.1, 0.1, 1,0,0,0,0];
% Param.IC.true = [0, 0, 0.35, 0,0,0,0]; % parameters used when simulating 'data'
% Param.IC.P = 0.1*eye(Param.numStates);
Param.IC.P = zeros(Param.numStates);
Param.IC.P(1,1) = 0.01; Param.IC.P(2,2) = 0.01; Param.IC.P(3,3) = 0.01;
% Param.IC.P = 1*eye(Param.numStates)


%% Data initializations ___________________________________________________
% Generating or loading body data
[bodyData] = check_generateBody(bodyString,name,bodyFunc,Paths,Param.sizeRadSS);
for k = 1:length(name.body)
    Param.body(k).A_r = bodyData.body(k).A;
%     Param.body(k).A_r = zeros(4,4);
    Param.body(k).B_r = bodyData.body(k).B;
%     Param.body(k).B_r = zeros(4,1);
    Param.body(k).C_r = bodyData.body(k).C;
%     Param.body(k).C_r = zeros(1,4);
%     Param.body(k).D_r = bodyData.body(k).D;
    Param.body(k).D_r = 0; % see Perez and Fossen body of work
end
% Generating or loading Bosma Data
[Data] = check_readBosmaData(dataString,name,Paths);
% Defining simulated spectrum data
[SpectrumData] = generateSpectrum(waveBool,bodyData);

% Determining Excitation Force Lumped Parameters
bodyIndex = 1; % this will call the FPS
vecTemp = SpectrumData.W - 2*pi*Param.peakFrequency;
[val,index] = min(abs(vecTemp));
Param.FexLumped = SpectrumData.body(bodyIndex).FeMag(index);

% ___________________________________________________________________
% EKF initializations
Param.options = odeset();
Param.R = (0.0006075)^2*eye(Param.numOutputs); % variance of the residuals
% Param.R = (0.001)*eye(Param.numOutputs);
% Param.R(1,1) = 0.0010147;
% Param.R(2,2) = 0.0001;
% Param.R = (0.0587)*eye(Param.numOutputs);
TF_mag = 5/2;
% Param.R(2,2) = Param.R(2,2);
% Param.R = zeros(Param.numOutputs)
% Param.R(1,1) = 0.003 % variance of the residuals
% Param.R(2,2) = TF_mag*Param.R(1,1);

% display('modified R')
% process noise covariance matrix
Param.GQGt = zeros(Param.numStates);

% determining the input wave variance
    for k = 1:length(Time.t)
        [pos,vel,acc] = particleDynamics_spec(Time.t(k));
        incidentWaveTemp(k) = pos;
    end
stdDevInputPos = std(incidentWaveTemp);

% Param.GQGt(1,1) = 0.0001;
% TF_mag = 5/2;
% Param.GQGt(2,2) = TF_mag*Param.GQGt(1,1);
% Param.GQGt(1,1) = 0.0001;
% Param.GQGt(2,2) = 0.0001;
% Param.GQGt(1,1) = 0.0139^2; % value from Determining Process noise scipt
% Param.GQGt(2,2) = 0.0828^2; % value from DeterminingProcessNoise script
% Param.GQGt(1,1) = 0.008^2; % value from DeterminingProcessNoise script
% Param.GQGt(2,2) = 0.05^2; % value from DeterminingProcessNoise script
% Param.GQGt(1,1) = 0.0666;
% Param.GQGt(2,2) = 0.0749;

% Param.GQGt(1,1) = 0.00052;
% Param.GQGt(2,2) = 0.0001679;

% Param.GQGt(1,1) = 0.0216^2; %run 11 info
% Param.GQGt(2,2) = 0.1005^2;
% Param.GQGt(3,3) = 600;

Param.GQGt(1,1) = runInfo(oRingIndex).posNoiseSTD^2;
Param.GQGt(2,2) = runInfo(oRingIndex).velNoiseSTD^2;
Param.GQGt(3,3) = (0.95*runInfo(oRingIndex).FexSTD)^2; % Used for first
% implementation
% Param.GQGt(3,3) = (0.75*runInfo(oRingIndex).FexSTD)^2; % Used for implementation with sub
% sample numbers
% Param.GQGt(3,3) = (0.9*runInfo(oRingIndex).FexSTD)^2; % Used for Estimation when noise is added.


% Param.GQGt(1,1) = 0.0001;
% Param.GQGt(2,2) = TF_mag*Param.GQGt(1,1);
% Param.GQGt(2,2) = 0.1;
% Param.GQGt(3,3) = 0.005;
% Param.GQGt(3,3) = 0.1;
% Param.GQGt(3,3) = 0.01;
% Param.GQGt(3,3) = 1.374*10^6/1000;
% Param.GQGt(3,3) = 500;
% Param.GQGt = Param.GQGt.*10
% Param.GQGt(3,3) = 2000;
% Param.GQGt(3,3) = (0.638338*100000);
% Param.GQGt(3,3) = 3000;
% Param.Q = Param.GQGt(1,1);
% Param.G = [1,TF_mag,Param.GQGt(3,3)/Param.Q,0,0,0,0]';
% Param.GQGt = Param.GQGt*0.0001; 
% Param.GQGt = Param.G*Param.Q*Param.G';
% Param.GQGt = Param.Q.*eye(Param.numStates);

% Artificial noise addition
NoiseFactor = 0;
Param.Ra = [NoiseFactor*Param.R(1,1),NoiseFactor*Param.R(2,2)];
Param.R = Param.R*2+NoiseFactor*Param.Ra;

Param.Ga = [0,0]; % artificial noise added as 'process noise'



%% Equations of motion {The Model}
% variables used in functions (These need to be changed if the order of SS realization is changed)
syms x1 x2 x3 x4 x5 x6 x7 U1 U2 U3 symTime; 
% vars = [x1,x2,x3,x4,x5,x6,x7]; % this needs to include all states used.
syms U1 U2 U3 symTime;
for k = 1:Param.numStates
    vars(k) = sym(strcat('x',num2str(k)));
end
bodyNumber = 1;
symVec = vars(1+Param.numNonRad:end);
% 
% dynamicsEq = [vars(2); 
%     1/Param.totalInertia*(excitationForce_spec(symTime) - (Param.heaveMoorStiff+Param.hydrostaticCoeff)*vars(1) - (Param.body(bodyNumber).C_r*symVec'+Param.body(bodyNumber).D_r*vars(2)));
%     0;
%     Param.body(bodyNumber).A_r*symVec'+Param.body(bodyNumber).B_r*vars(2)];
% dynamicsEq = [vars(2); % JP1 model, with relative velocity removed (Viscous drag is linear)
%     1/Param.totalInertia*(excitationForce_spec(symTime) - (Param.heaveMoorStiff+Param.hydrostaticCoeff)*vars(1) - Param.dragFactor*vars(3)*vars(2) - (Param.body(bodyNumber).C_r*symVec'+Param.body(bodyNumber).D_r*vars(2)));
%     0;
%     Param.body(bodyNumber).A_r*symVec'+Param.body(bodyNumber).B_r*vars(2)];
% dynamicsEq = [vars(2);
%     1/Param.totalInertia*(excitationForce_spec(symTime) - (Param.heaveMoorStiff+Param.hydrostaticCoeff)*vars(1) - 150*vars(2));
%     0];

% outputEq = [vars(1)];

% Generating functions
% [Functions] = check_genFunctions(genFunString,bodyNumber,vars,dynamicsEq,outputEq,Paths);
% Func = Functions;

% Loading data 
Param.IC.est = [0.1, 0.1, 0,zeros(1,Param.sizeRadSS)];
Param.IC.true = [0, 0, 0,zeros(1,Param.sizeRadSS)]; % parameters used when simulating 'data'

[savedData] = check_measurements(measString,savedDataName,Paths,[Time.start,Time.end]);
% measResidual = var(Param.residual)

% figure
% plot(Time.t,savedData(:,1),Time.t,savedData(:,2))
% legend('pos','vel')

% generating the functions that will be used to identify
clear dynamicsEq Functions
% dynamicsEq = [vars(2);
%     1/Param.totalInertia*(vars(3) - (Param.heaveMoorStiff+Param.hydrostaticCoeff)*vars(1) - 150*vars(2));
%     0];
% dynamicsEq = [vars(2); % includes quadradic drag with no relative
%     1/Param.totalInertia*(vars(3) - Param.dragFactor*Param.C_D*(vars(2))*abs(vars(2)) - (Param.heaveMoorStiff+Param.hydrostaticCoeff)*vars(1) - (Param.body(bodyNumber).C_r*symVec'+Param.body(bodyNumber).D_r*vars(2)));
%     0;
%     Param.body(bodyNumber).A_r*symVec'+Param.body(bodyNumber).B_r*vars(2)];
% dynamicsEq = [vars(2); % current best working
%     1/Param.totalInertia*(vars(3) - (Param.heaveMoorStiff+Param.hydrostaticCoeff)*vars(1) - Param.linearDrag*vars(2) - (Param.body(bodyNumber).C_r*symVec'+Param.body(bodyNumber).D_r*vars(2)));
%     0;
%     Param.body(bodyNumber).A_r*symVec'+Param.body(bodyNumber).B_r*vars(2)];
dynamicsEq = [vars(2); % nonlinear and Linear Drag
    1/Param.totalInertia*(vars(3) - (Param.heaveMoorStiff+Param.hydrostaticCoeff)*vars(1) - Param.linearDrag*vars(2) - Param.dragFactor*(vars(2)-U2)*abs(vars(2)-U2) - (Param.body(bodyNumber).C_r*symVec'+Param.body(bodyNumber).D_r*vars(2)));
    0;
    Param.body(bodyNumber).A_r*symVec'+Param.body(bodyNumber).B_r*vars(2)];

outputEq = [vars(1);vars(2)];


[Functions] = check_genFunctions(genFunString,bodyNumber,vars,dynamicsEq,outputEq,Paths);
Func = Functions;


% return
%% Run EKF to identify drag coefficient
Measurements = savedData;
[sim] = smoother_RTS(Func, Measurements(:,1:Param.numOutputs));

% storing EKF and RTS structures for plotting
EKF = sim.EKF; % excecute if RTS is active
if Param.RTSBool == 1
    RTS = sim.RTS;
end

% saveFileName = 'EKF4.mat';
% save(saveFileName,'EKF','Measurements')
% movefile(saveFileName,Paths.outputDir);


%% Plotting results
tempTime = sim.EKF.time.hat(1:Time.simStepsPerDT:end);
for k = 1:length(tempTime)
    
%     F_ex(k) = excitationForce_spec(tempTime(k)); % simulated
    F_ex(k) = excitationConv(tempTime(k),bodyData.body(1).FIRF,bodyData.body(1).irfTime,Data.wave.incidentFPS,Data.wave.time);
end
   
% calculating the residual sum of squares
num = sum((sim.EKF.x.hat(3,1:Time.simStepsPerDT:end)-F_ex).^2);
% den = sum((F_ex-mean(F_ex)).^2);
den = sum((F_ex).^2);

NMSE.Fex(j) = 1-(num/den)

% test for whiteness
stateError = (Measurements(1:end-1,1)-sim.EKF.x.hat(1,1:Time.simStepsPerDT:end)');
[m,n] = size(stateError);
point = 2;
% [rho,limit95]=whiteness_test(stateError(m-200:end),point)

if indicies(1) == indicies(2)
% plotting the autocorrelation
figure()
plot(tempTime,stateError)
axis([Time.start,Time.end,-0.02,0.02])
ylabel('position error')
title('plotting error')

figure()
autocorr(stateError)


figure()
hold on
plot(tempTime,F_ex,tempTime,sim.EKF.x.hat(3,1:Time.simStepsPerDT:end))
plot(tempTime,sim.EKF.x.hat(3,1:Time.simStepsPerDT:end)-sim.EKF.sig3(3,1:Time.simStepsPerDT:end),'m--')
plot(tempTime,sim.EKF.x.hat(3,1:Time.simStepsPerDT:end)+sim.EKF.sig3(3,1:Time.simStepsPerDT:end),'m--')
legend('F_ex true','F_ex est')
title('Fex comparison')
hold off

figure()
hold on
plot(tempTime,Measurements(1:end-1,1),tempTime,sim.EKF.x.hat(1,1:Time.simStepsPerDT:end))
plot(tempTime,sim.EKF.x.hat(1,1:Time.simStepsPerDT:end)-sim.EKF.sig3(1,1:Time.simStepsPerDT:end),'m--')
plot(tempTime,sim.EKF.x.hat(1,1:Time.simStepsPerDT:end)+sim.EKF.sig3(1,1:Time.simStepsPerDT:end),'m--')
title('Comparing position')
legend('meas','est')
hold off

figure()
hold on
plot(tempTime,Measurements(1:end-1,2),tempTime,sim.EKF.x.hat(2,1:Time.simStepsPerDT:end))
plot(tempTime,sim.EKF.x.hat(2,1:Time.simStepsPerDT:end)-sim.EKF.sig3(2,1:Time.simStepsPerDT:end),'m--')
plot(tempTime,sim.EKF.x.hat(2,1:Time.simStepsPerDT:end)+sim.EKF.sig3(2,1:Time.simStepsPerDT:end),'m--')
title('Comparing velocity')
legend('meas','est')
hold off


% % testing
%  figure()
% plot(tempTime,F_ex./1000,Data.wave.time,Data.wave.incidentFPS)
% legend('Fex','wave')
% axis([35,37.5,-0.1,0.1])

figure()
plot(Data.wave.time(1:end-1),Data.wave.incidentFPSVel,Time.t,Est.vel)
legend('meas','est')
title('Comparison of water particle velocity')
axis([Time.start,Time.end,-1,1])


Results.time = tempTime;
Results.measPos = Measurements(1:end-1);
Results.measVel = Measurements(2:end-1);
Results.FexTrue = F_ex;
Results.water.time = Data.wave.time(1:end-1);
Results.water.vel = Data.wave.incidentFPSVel;

Results.estFex = sim.EKF.x.hat(3,1:Time.simStepsPerDT:end);
Results.estPos = sim.EKF.x.hat(1,1:Time.simStepsPerDT:end);
Results.estVel = sim.EKF.x.hat(2,1:Time.simStepsPerDT:end);
Results.estWater.time = Time.t;
Results.estWater.vel = Est.vel

measWaveVel = interp1(Data.wave.time(1:end-1),Data.wave.incidentFPSVel,Time.t);
num = sum((Est.vel-measWaveVel).^2);
den = sum((measWaveVel).^2);
Results.NMSE.velEst = 1-(num/den);
Results.NMSE.Fex = NMSE.Fex;

save('DirectResults','Results')

end

measWaveVel = interp1(Data.wave.time(1:end-1),Data.wave.incidentFPSVel,Time.t);
num = sum((Est.vel-measWaveVel).^2);
% den = sum((F_ex-mean(F_ex)).^2);
den = sum((measWaveVel).^2);

NMSE.velEst(j) = 1-(num/den);
toc



end
