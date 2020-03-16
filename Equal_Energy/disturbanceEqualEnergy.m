% Default Main Code: The purpose of this script is to be the single
% executable script to run an example. 
% 
% As of right now the script has the mesh, nemoh, IRFs, and experimental
% data this script also sets up to simulate the FPS which encounters the
% wave first.

clear all
close all
clc


rng(1)

indices = [18,18];
for j = indices(1):indices(2)
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
Param.velEstimateVal = 2; % 0-> uses particleDynamics_spec 1-> uses the direct Fex Est approach  2-> disturbance model
Param.numFexComponents = 5;
Param.variableFexComponents = []; % enter which entries you want to be variable, [first second third]
Param.numVarComponents = length(Param.variableFexComponents); % specify which discretized frequencies you want variable, otherwise use [0]

% generate MATLAB Path
% addedPath = genpath('..\Functions');
addpath(Paths.functionsDir)
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
oRingIndex = j
waveRuns = [1,180,0.075,1.1;
            2,181,0.175,1.42;
            3,182,0.05,1.9;
            4,183,0.175,1.9;
            5,184,0.3,1.9;
            6,185,0.050,2.37;
            7,186,0.175,2.37;
            8,199,0.3,2.37;
            9,189,0.05,2.69;
            10,190,0.175,2.69;
            11,191,0.3,2.69;
            12,192,0.05,3.16;
            13,193,0.175,3.16;
            14,194,0.175,3.95;
            15,195,0.175,4.74;
            16,196,0.175,5.53;
            17,197,0.375,2.06;
            18,198,0.375,2.37];
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

Time.subSampleNo = 2; % uses one out of every subSampleNo
Time.dT = Time.subSampleNo/Param.samplingFreq.motion;
Time.simStepsPerDT = 50;
% Time.simStepsPerDT = 10
Time.dtSim = Time.dT/Time.simStepsPerDT;
Time.start = 50; % min 1 period into timespan
Time.end = 80;% max 1 period from the end of the timespan
Time.t = [Time.start:Time.dT:Time.end];

% Number of states and control inputs
Param.sizeRadSS = 4;
Param.numDynamics = 2;
Param.numIDStates = 2*Param.numFexComponents+Param.numVarComponents;
Param.numStates = Param.numDynamics+Param.numIDStates+Param.sizeRadSS; % 2 for dynamics 
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


Param.dispVol = 11322077e-9; % from solidworks for the FPS
Param.Delta = Param.dispVol*Param.rho;
Param.dragFactor = 1/2*Param.rho*Param.upperArea;
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
vecTemp = SpectrumData.W - 2*pi*Param.peakFrequency;
[val,index] = min(abs(vecTemp));
Param.FexLumped = SpectrumData.body(bodyNumber).FeMag(index);
% ___________________________________________________________________

% % EKF initializations
pwrThreshold = 0.003; % cutoff very low energy wave frequencies 
% [Param.DistVals] = chooseDistFreq(bodyNumber,pwrThreshold);

% for k =1:Param.numFexComponents
%     Param.W(k) = Param.DistVals.W(k);
%     Param.DF(k) = Param.DistVals.force(k);
% end

%% chooseDistFreq contents copied here, used for plotting
numComponents = Param.numFexComponents; % this will be the number of spectral divisions
[peakVal,peakIndex] = max(SpectrumData.pwr);
[val,firstThresIndex] = min(abs(SpectrumData.pwr(1:peakIndex) - pwrThreshold*peakVal));
[val,secondThreshIndex] = min(abs(SpectrumData.pwr(peakIndex:end) - pwrThreshold*peakVal));

% generating a new spectrum with only 'high' energy frequencies
num = 10000;
minOmega = SpectrumData.W(firstThresIndex);
maxOmega = SpectrumData.W(secondThreshIndex);
frequency = linspace(minOmega/(2*pi),maxOmega/(2*pi),num);
dW = (frequency(2)-frequency(1))*2*pi;

for j = 1:length(frequency) % loop to calculate specta
    A = frequency(j)/(Param.wavePeriod^(-1));
    S(j) = 5*(Param.waveHs)^2/(16*(Param.wavePeriod)^(-1))*(1/A^5)*exp(-5/4*A^(-4));
end

PowerSpec = S;
% Amplitude = sqrt(2*S.*dW); not applicable selecting few frequencies
Omega = frequency.*2*pi;

% using SpectrumData to determine total energy
totalArea = trapz(PowerSpec); % if evenly spaced, only the total area is necessary
componentArea = totalArea/numComponents;

% vec1 = SpecrtrumData.amp-maxAmp/numComponents;
% PeriodStart = periodRange(1);
% PeriodEnd = periodRange(2);

if numComponents == 1
    clear frequency Amplitude
    DistVals.frequency = Param.wavePeriod^(-1);
    DistVals.amplitude = Param.waveHs;   
    DistVals.force = DistVals.amplitude*Param.FexLumped;
else
%     
    divisionIndex = 1; % starting at the first division of the power
    divisionStart = 1; % staring the with the first index of spectrum Data
    divisionSave = divisionStart;
    for k = 1:length(PowerSpec)
        if trapz(PowerSpec(1:k)) >= divisionIndex*componentArea
            avgDivPwr = mean(PowerSpec(divisionStart:k));
            [val,index] = min(abs(PowerSpec(divisionStart:k)-avgDivPwr));
            spectrumIndices(divisionIndex) = index+divisionStart-1;
            divisionDf(divisionIndex) = (frequency(k)-frequency(divisionStart));
            
%             increments
            divisionStart = k+1;
            divisionIndex = divisionIndex + 1;
            divisionSave(divisionIndex) = divisionStart;
            
        end
        
    end
    
%     if statement in case there is numerical error and a last frequency is
%     not chosen
    if length(spectrumIndices) < numComponents
        avgDivPwr = mean(PowerSpec(divisionStart:k));
       [val,index] = min(abs(PowerSpec(divisionStart:k)-avgDivPwr));
        spectrumIndices(divisionIndex) = index+divisionStart-1;
       divisionDf(divisionIndex) = (frequency(k)-frequency(divisionStart));

    end
    


for k = 1:numComponents % cycle through the choices of frequency
    DistVals.power(k) = PowerSpec(spectrumIndices(k));
    DistVals.amplitude(k) = sqrt(2*DistVals.power(k)*divisionDf(k));
    DistVals.frequency(k) = frequency(spectrumIndices(k));
    FexGain(k) = interp1(SpectrumData.frequency,SpectrumData.body(bodyNumber).FeMag,frequency(spectrumIndices(k)));
    DistVals.force(k) = DistVals.amplitude(k)*FexGain(k);
end

end % terminates 1 frequency condition
DistVals.W = 2*pi.*DistVals.frequency;

%plotting the spectum

figure()
hold on
plot(frequency,S,'k',DistVals.frequency,DistVals.power,'ro')
for k = 1:Param.numFexComponents
plot([frequency(divisionSave(k)),frequency(divisionSave(k))],[S(divisionSave(k)),0],'b')
end

hold off
xlabel('Wave Frequency [Hz]')
ylabel('Wave Power [m^2/Hz]')
title('Wave Spectrum')
legend('Power Spectrum','Discrete Frequency','Area Boundaries')
axis([0.1,1.2,0,0.032])

end
