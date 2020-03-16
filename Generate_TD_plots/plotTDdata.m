clear all
close all
clc

load DirectResults.mat
DirectResults = Results; clear Results
load DistResults
DistResults = Results; clear Results

% plotting excitation force estimation
figure()
plot(DirectResults.time,DirectResults.FexTrue,'k--',DirectResults.time,DirectResults.estFex,'b',DistResults.time,DistResults.estFex,'r')
legend('Reference','Direct Est.','Disturb Est.')
ylabel('Excitation Force [N]')
xlabel('Time [s]')
title('Excitation Force Estimations')

% plotting the water velocity
timeStart = 50;
timeEnd = 110;
[val,IndexStart] = min(abs(DirectResults.water.time-timeStart));
[val,IndexEnd] = min(abs(DirectResults.water.time-timeEnd));


figure()
plot(DirectResults.water.time(IndexStart:IndexEnd),DirectResults.water.vel(IndexStart:IndexEnd),'k--',DirectResults.estWater.time,DirectResults.estWater.vel,'b',DistResults.estWater.time,DistResults.estWater.vel,'r')
legend('Reference','Direct Est.','Disturb Est.')
ylabel('Water Velocity [m/s]')
xlabel('Time [s]')
title('Incident Wave Estimations')