Sections \ref{sec:wecModel} and \ref{sec:experiment} 
describe the mathematical model and the experimental tests used in the development and implementation of the algorithms in this work, respectively.clear all
close all
clc

load runStruct.mat
for k = 2:18
HsVec(k-1) = runInfo(k).Hs;
TpVec(k-1) = runInfo(k).Ts;
sizeRatioVec(k-1) = runInfo(k).sizeRatio;
end

% load first runs 0.95*FexSTD
load NMSE_no4.mat; fullID(1) = NMSE; % DirectID
load NMSE_no5.mat; fullID(4) = NMSE; % 3rd order
load NMSE_no6.mat; fullID(3) = NMSE; % 2nd order
load NMSE_no7.mat; fullID(2) = NMSE; % 1st order
load NMSE_no8.mat; fullID(5) = NMSE; % 4th order
load NMSE_no9.mat; fullID(6) = NMSE; % 5th order

% boxplot of best case identification
figure()
boxplot([fullID(1).Fex(2:end)',fullID(2).Fex(2:end)',fullID(3).Fex(2:end)',fullID(4).Fex(2:end)',fullID(5).Fex(2:end)',fullID(6).Fex(2:end)'],'Label',{'Direct','Dist. O(1)','Dist. O(2)','Dist. O(3)','Dist. O(4)','Dist. O(5)'})
ylabel('GOF')
xlabel('Method')
title('Comparison of Estimation Methods')


%% plots relating wave parameters to NMSE
figure()
plot(HsVec,fullID(1).Fex(2:end),'kx',HsVec,fullID(2).Fex(2:end),'m+',HsVec,fullID(3).Fex(2:end),'bo',HsVec,fullID(4).Fex(2:end),'rd',HsVec,fullID(5).Fex(2:end),'gs',HsVec,fullID(6).Fex(2:end),'c^','markers',8)
title('GOF as a function of significant wave height')
legend('Direct','Dist. O(1)','Dist. O(2)','Dist. O(3)','Dist. O(4)','Dist. O(5)','Location','SouthEast')
ylabel('GOF')
xlabel('Significant Wave Height')
axis([0,0.4,0.5,1])

figure()
plot(TpVec,fullID(1).Fex(2:end),'kx',TpVec,fullID(2).Fex(2:end),'m+',TpVec,fullID(3).Fex(2:end),'bo',TpVec,fullID(4).Fex(2:end),'rd',TpVec,fullID(5).Fex(2:end),'gs',TpVec,fullID(6).Fex(2:end),'c^','markers',8)
title('GOF as a function of peak period')
ylabel('GOF')
xlabel('Peak Period [s]')
legend('Direct','Dist. O(1)','Dist. O(2)','Dist. O(3)','Dist. O(4)','Dist. O(5)','Location','SouthEast')

figure()
plot(sizeRatioVec,fullID(1).Fex(2:end),'kx',sizeRatioVec,fullID(2).Fex(2:end),'m+',sizeRatioVec,fullID(3).Fex(2:end),'bo',sizeRatioVec,fullID(4).Fex(2:end),'rs',sizeRatioVec,fullID(5).Fex(2:end),'gs',sizeRatioVec,fullID(6).Fex(2:end),'c^','markers',8)
title('GOF as a function of size ratio')
ylabel('GOF')
xlabel('Size ratio')
legend('Direct','Dist. O(1)','Dist. O(2)','Dist. O(3)','Dist. O(4)','Dist. O(5)','Location','SouthWest')

%Velocity comparison
figure()
boxplot([fullID(1).velEst(2:end)',fullID(2).velEst(2:end)',fullID(3).velEst(2:end)',fullID(4).velEst(2:end)',fullID(5).velEst(2:end)',fullID(6).velEst(2:end)'],'Label',{'Direct','Dist. O(1)','Dist. O(2)','Dist. O(3)','Dist. O(4)','Dist. O(5)'})
ylabel('GOF')
xlabel('Method')
title('Comparison of Estimated Velocity')


%% SubSampling comparison
load NMSE_no10.mat; subSamDist(1) = NMSE; % 200Hz
load NMSE_no11.mat; subSamDist(2) = NMSE; % 100Hz
load NMSE_no12.mat; subSamDist(3) = NMSE; % 40Hz
load NMSE_no13.mat; subSamDist(4) = NMSE; % 20Hz
load NMSE_no14.mat; subSamDist(5) = NMSE; % 10Hz
load NMSE_no15.mat; subSamDist(6) = NMSE; % 5Hz
load NMSE_no48.mat; subSamDist(7) = NMSE; % 2Hz


figure()
boxplot([subSamDist(1).Fex(2:end)',subSamDist(2).Fex(2:end)',subSamDist(3).Fex(2:end)',subSamDist(4).Fex(2:end)',subSamDist(5).Fex(2:end)',subSamDist(6).Fex(2:end)',subSamDist(7).Fex(2:end)'],'Label',{'200Hz','100Hz','40Hz','20Hz','10Hz','5Hz','2Hz'})
ylabel('GOF')
xlabel('Sampling Frequency')
title('Disturbance Method Sampling Rate Sensitivity')

load NMSE_no16.mat; subSamDir(1) = NMSE; % 200Hz
load NMSE_no17.mat; subSamDir(2) = NMSE; % 100Hz
load NMSE_no18.mat; subSamDir(3) = NMSE; % 40Hz
load NMSE_no19.mat; subSamDir(4) = NMSE; % 20Hz
load NMSE_no20.mat; subSamDir(5) = NMSE; % 10Hz
load NMSE_no21.mat; subSamDir(6) = NMSE; % 5Hz
load NMSE_no47.mat; subSamDir(7) = NMSE; % 2Hz

figure()
boxplot([subSamDir(1).Fex(2:end)',subSamDir(2).Fex(2:end)',subSamDir(3).Fex(2:end)',subSamDir(4).Fex(2:end)',subSamDir(5).Fex(2:end)',subSamDir(6).Fex(2:end)',subSamDir(7).Fex(2:end)'],'Label',{'200Hz','100Hz','40Hz','20Hz','10Hz','5Hz','2Hz'})
ylabel('GOF')
xlabel('Sampling Frequency')
title('Direct Method Sampling Rate Sensitivity')


%% Noise addition
load NMSE_no22.mat; noiseDir(1) = NMSE; % noiseFactor = 5
load NMSE_no23.mat; noiseDir(2) = NMSE; % noiseFactor = 10
load NMSE_no24.mat; noiseDir(3) = NMSE; % noiseFactor = 50
load NMSE_no29.mat; noiseDir(4) = NMSE; % noiseFactor = 75
load NMSE_no25.mat; noiseDir(5) = NMSE; % noiseFactor = 100
load NMSE_no28.mat; noiseDir(6) = NMSE; % noiseFactor = 300
load NMSE_no26.mat; noiseDir(7) = NMSE; % noiseFactor = 500

figure()
boxplot([noiseDir(1).Fex(2:end)',noiseDir(2).Fex(2:end)',noiseDir(3).Fex(2:end)',noiseDir(4).Fex(2:end)',noiseDir(5).Fex(2:end)',noiseDir(6).Fex(2:end)',noiseDir(7).Fex(2:end)'],'Label',{'5','10','50','75','100','300','500'})
ylabel('GOF')
xlabel('Noise Factor')
title('Direct Method Sampling Noise Sensitivity')


load NMSE_no30.mat; noiseDist(1) = NMSE; % noiseFactor = 5
load NMSE_no31.mat; noiseDist(2) = NMSE; % noiseFactor = 10
load NMSE_no32.mat; noiseDist(3) = NMSE; % noiseFactor = 50
load NMSE_no36.mat; noiseDist(4) = NMSE; % noiseFactor = 75
load NMSE_no33.mat; noiseDist(5) = NMSE; % noiseFactor = 100
load NMSE_no34.mat; noiseDist(6) = NMSE; % noiseFactor = 300
load NMSE_no35.mat; noiseDist(7) = NMSE; % noiseFactor = 500

figure()
boxplot([noiseDist(1).Fex(2:end)',noiseDist(2).Fex(2:end)',noiseDist(3).Fex(2:end)',noiseDist(4).Fex(2:end)',noiseDist(5).Fex(2:end)',noiseDist(6).Fex(2:end)',noiseDist(7).Fex(2:end)'],'Label',{'5','10','50','75','100','300','500'})
ylabel('GOF')
xlabel('Noise Factor')
title('Disturbance Method Sampling Noise Sensitivity')



%% Radiation Order 
load NMSE_no37.mat; radDir(1) = NMSE; % RadSS 1
load NMSE_no39.mat; radDir(2) = NMSE; % RadSS 2
load NMSE_no41.mat; radDir(3) = NMSE; % RadSS 3
load NMSE_no43.mat; radDir(4) = NMSE; % RadSS 4
load NMSE_no45.mat; radDir(5) = NMSE; % RadSS 5

figure()
boxplot([radDir(1).Fex(2:end)',radDir(2).Fex(2:end)',radDir(3).Fex(2:end)',radDir(4).Fex(2:end)',radDir(5).Fex(2:end)'],'Label',{'O(1)','O(2)','O(3)','O(4)','O(5)'})
ylabel('GOF')
xlabel('Radiation SS Realization Order')
title('Direct Method Radiation Order Sensitivity')


load NMSE_no38.mat; radDist(1) = NMSE; % RadSS 1
load NMSE_no40.mat; radDist(2) = NMSE; % RadSS 2
load NMSE_no42.mat; radDist(3) = NMSE; % RadSS 3
load NMSE_no44.mat; radDist(4) = NMSE; % RadSS 4
load NMSE_no46.mat; radDist(5) = NMSE; % RadSS 5

figure()
boxplot([radDist(1).Fex(2:end)',radDist(2).Fex(2:end)',radDist(3).Fex(2:end)',radDist(4).Fex(2:end)',radDist(5).Fex(2:end)'],'Label',{'O(1)','O(2)','O(3)','O(4)','O(5)'})
ylabel('GOF')
xlabel('Radiation SS Realization Order')
title('Disturbance Method Radiation Order Sensitivity')


