clear all
close all
clc

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

  Hs = waveRuns(:,3); % significant wave height
  Wp = 2*pi./waveRuns(:,4); % natural period

  M0 = (Hs/4).^2; % 0th moment of the spectrum
  M2 = 1.982.*M0.*Wp;% second moment of the spectrum
  
  Tz = sqrt(M0./M2);% zero crossing period
  
  %% Determining Bandwidth
num = 1000; % number of sine waves summed to generate irregular waves
minOmega = 0.15;
maxOmega = 30;
frequency = linspace(minOmega/(2*pi),maxOmega/(2*pi),num);
dW = (frequency(2)-frequency(1))*2*pi;
dF = frequency(2)-frequency(1);

index = 18;
for j = 1:length(frequency) % loop to calculate specta
    A = frequency(j)/(Wp(index)^(-1));
    S(j) = 5*(Hs(index))^2/(16*(Wp(index))^(-1))*(1/A^5)*exp(-5/4*A^(-4));
end

Smax = max(S)
plot(frequency,S-0.5*Smax)
axis([frequency(1),frequency(end),0,0.5*Smax])


%% Determining Resonant Period of WEC
% load FPS_ABFe_freq.mat
% 
% figure()
% B_BEM(1,:) = B(3,3,:);
% plot(omega,B_BEM(1,:))
% title('Radiation Damping')
% figure()
% plot(omega,squeeze(A(3,3,:)))
% title('Added Mass')
% figure()
% plot(omega,real(Fe(:,3)))
% title('real part of excitation force')
% figure()
% plot(omega,imag(Fe(:,3)))
% title('imag part of excitation force')
% figure()
% plot(omega,abs(Fe(:,3)))
% title('magnitude of Fex')
% figure()
% plot(omega,angle(Fe(:,3)))
% title('phase of Fex')
clear all

Rb = (.2706+0.25)/2; % intermediate big radius for FPS
Rs = (.1894+.175)/2; % small radius for FPS
z = [0.15,0.15,0,-0.07,-0.075,-0.075];
r = [0,Rb,Rb,Rs+(Rb-Rs)/3,Rs,0];
display('FPS: Enter the Following: 36, storage, 0.05, 500')
[Mass,Inertia,KH,XB,YB,ZB] = axiMesh(r,z,length(z));

num_omega = 50;
min_omega = 0.1;
max_omega = 10;
% rounding omega to the nearest radian/second
min_omega = round(min_omega,1); max_omega = round(max_omega,1); % Warning!!!
% watch these round functions, because the most common error I am getting
% is when the round function takes one of the input frequencies to be zero
omega = linspace(min_omega,max_omega,num_omega);

num_angle = 1;
min_angle = 0;
max_angle = 0;
angle_specified = linspace(min_angle,max_angle,num_angle);
depth = 0; % deep water waves
[A,B,Fe] = Nemoh(omega,angle_specified,depth);

%%
close all
figure()
subplot(2,1,1)
plot(omega,squeeze(B(3,3,:)))
title('Damping and Added Mass')
ylabel('Radiation Damping [kg/s]')
axis([0,10,0,105])
subplot(2,1,2)
plot(omega,squeeze(A(3,3,:)))
ylabel('Added Mass [kg]')
axis([0,10,0,45])
xlabel('Frequency [rad/s]')

figure()
subplot(2,1,1)
plot(omega,abs(Fe(:,3)))
axis([0,10,0,2200])
title('Excitation Force')
ylabel('Magnitude [N]')
subplot(2,1,2)
plot(omega,angle(Fe(:,3)))
axis([0,10,-2.3,0])
ylabel('Phase [rad]')
xlabel('Frequency [rad/s]')


%% Radiation Realization Recreating realizations
clear all
close all
clc
% 
% Rb = (.2706+0.25)/2; % intermediate big radius for FPS
% Rs = (.1894+.175)/2; % small radius for FPS
% z = [0.15,0.15,0,-0.07,-0.075,-0.075];
% r = [0,Rb,Rb,Rs+(Rb-Rs)/3,Rs,0];
% display('FPS: Enter the Following: 36, storage, 0.05, 500')
% [Mass,Inertia,KH,XB,YB,ZB] = axiMesh(r,z,length(z));
% 
% num_omega = 25;
% max_period =10;
% min_omega = 2*pi/max_period; 
% min_period = .15;
% max_omega = 2*pi/min_period;
% % rounding omega to the nearest radian/second
% min_omega = round(min_omega,1); max_omega = round(max_omega,1); % Warning!!!
% % watch these round functions, because the most common error I am getting
% % is when the round function takes one of the input frequencies to be zero
% omega = linspace(min_omega,max_omega,num_omega);
% 
% num_angle = 1;
% min_angle = 0;
% max_angle = 0;
% angle_specified = linspace(min_angle,max_angle,num_angle);
% depth = 0; % deep water waves
% [A,B,~] = Nemoh(omega,angle_specified,depth);

load('FPS_ABFe_freq')
b = squeeze(B(3,3,:));
plot(omega,b)

figure()
plot(omega,b,'bo')
hold on
dW = 0.01;
xx = 0:dW:max(omega); % desired Resolution
yy = spline(omega,b,xx); % creating cubic spline data to fill in gaps
plot(xx,yy)
hold off

% compute the impulse response function using trapazoidal integration
h = 0.001;
t = [0:h:2.5];

Sum = zeros(1,length(t));
for jj = 1:length(t)
    for j = 1:length(xx)-1
        
        Trap = 1/2*(yy(j)+yy(j)) * dW;
        Trap = Trap * cos(xx(j)*t(jj));
        Sum(jj) = Sum(jj)+Trap;
    end    
end
IRF = 2/pi*Sum;
figure
plot(t,IRF)

% O(1) model
load('FPS_ABCD_r1')
A1 = A_r; B1 = B_r; C1 = C_r; D1 = 0;
model1 = ss(A1,B1,C1,D1);
[IRF1,~] = impulse(model1,t);
fit1 = goodnessOfFit(IRF1, IRF', 'NMSE')
% O(2) model
load('FPS_ABCD_r2')
A2 = A_r; B2 = B_r; C2 = C_r; D2 = 0;
model2 = ss(A2,B2,C2,D2);
[IRF2,~] = impulse(model2,t);
fit2 = goodnessOfFit(IRF2, IRF', 'NMSE')
% O(3) model
load('FPS_ABCD_r3')
A3 = A_r; B3 = B_r; C3 = C_r; D3 = 0;
model3 = ss(A3,B3,C3,D3);
[IRF3,~] = impulse(model3,t);
fit3 = goodnessOfFit(IRF3, IRF', 'NMSE')
% O(4) model
load('FPS_ABCD_r4')
A4 = A_r; B4 = B_r; C4 = C_r; D4 = 0;
model4 = ss(A4,B4,C4,D4);
[IRF4,~] = impulse(model4,t);
fit4 = goodnessOfFit(IRF4, IRF', 'NMSE')
% O(5) model
load('FPS_ABCD_r5')
A5 = A_r; B5 = B_r; C5 = C_r; D5 = 0;
model5 = ss(A5,B5,C5,D5);
[IRF5,~] = impulse(model5,t);
fit5 = goodnessOfFit(IRF5, IRF', 'NMSE')

H5 = freqresp(model5,xx);
for j = 1:length(xx)
    Mag5(j) = abs(real(H5(1,1,j)));
end
figure()
plot(xx,yy,'b',xx,Mag5,'r.')
legend('Nemoh','realization')
legend('Nemoh','realization')
title('F.R.F. Reduced System')


figure()
plot(t,IRF,'k',t,IRF1,'g',t,IRF2,'r-.',t,IRF3,t,IRF4,'r',t,IRF5,'b:')
axis([0,2.5,-300,300])
legend('NEMOH IRF','O(1) Realization','O(2) Realization','O(3) Realization','O(4) Realization','O(5) Realization')



figure()
plot(t,IRF,'k',t,IRF1,'g',t,IRF2,'r-.',t,IRF3,t,IRF4,'r',t,IRF5,'b:')
legend('NEMOH IRF','O(1) Realization','O(2) Realization','O(3) Realization','O(4) Realization','O(5) Realization')
title('Radiation IRF')
ylabel('Radiation Force [N]')
xlabel('Time [s]')




% num_omega = 50;
% min_omega = 0.1;
% max_omega = 10;
% % rounding omega to the nearest radian/second
% min_omega = round(min_omega,1); max_omega = round(max_omega,1); % Warning!!!
% % watch these round functions, because the most common error I am getting
% % is when the round function takes one of the input frequencies to be zero
% omega = linspace(min_omega,max_omega,num_omega);
% 
% num_angle = 1;
% min_angle = 0;
% max_angle = 0;
% angle_specified = linspace(min_angle,max_angle,num_angle);
% depth = 0; % deep water waves
% [A,B,Fe] = Nemoh(omega,angle_specified,depth);
% % Creating SS realization of the Radiation Damping Coefficient
% b = squeeze(B(3,3,:));
% % frequency depended excitation force magnitude
% 
% figure()
% plot(omega,b,'bo')
% hold on
% dW = 0.01;
% xx = 0:dW:max(omega); % desired Resolution
% yy = spline(omega,b,xx); % creating cubic spline data to fill in gaps
% plot(xx,yy)
% hold off
% 
% % compute the impulse response function using trapazoidal integration
% % h = 0.002;
% h = 0.05
% t = [0:h:5];
% 
% Sum = zeros(1,length(t));
% for jj = 1:length(t)
%     for j = 1:length(xx)-1
%         
%         Trap = 1/2*(yy(j)+yy(j)) * dW;
%         Trap = Trap * cos(xx(j)*t(jj));
%         Sum(jj) = Sum(jj)+Trap;
%     end    
% end
% IRF = 2/pi*Sum;
% figure
% plot(t,IRF)
% 
% % O(1) model
% load('FPS_ABCD_r1')
% A1 = A_r; B1 = B_r; C1 = C_r; D1 = 0;
% model1 = ss(A1,B1,C1,D1);
% [IRF1,~] = impulse(model1,t);
% % O(2) model
% load('FPS_ABCD_r2')
% A2 = A_r; B2 = B_r; C2 = C_r; D2 = 0;
% model2 = ss(A2,B2,C2,D2);
% [IRF2,~] = impulse(model2,t);
% % O(3) model
% load('FPS_ABCD_r3')
% A3 = A_r; B3 = B_r; C3 = C_r; D3 = 0;
% model3 = ss(A3,B3,C3,D3);
% [IRF3,~] = impulse(model3,t);
% % O(4) model
% load('FPS_ABCD_r4')
% A4 = A_r; B4 = B_r; C4 = C_r; D4 = 0;
% model4 = ss(A4,B4,C4,D4);
% [IRF4,~] = impulse(model4,t);
% % O(5) model
% load('FPS_ABCD_r5')
% A5 = A_r; B5 = B_r; C5 = C_r; D5 = 0;
% model5 = ss(A5,B5,C5,D5);
% [IRF5,~] = impulse(model5,t);
% H5 = freqresp(model5,xx);
% for j = 1:length(xx)
%     Mag5(j) = abs(real(H5(1,1,j)));
% end
% figure()
% plot(xx,yy,'b',xx,Mag5,'r.')
% legend('Nemoh','realization')
% legend('Nemoh','realization')
% title('F.R.F. Reduced System')
% 
% 
% figure()
% plot(t,IRF,'k--',t,IRF1,t,IRF2,t,IRF3,t,IRF4,t,IRF5)
% axis([0,5,-300,450])
% legend('NEMOH IRF','O(1) Realization','O(2) Realization','O(3) Realization','O(4) Realization','O(5) Realization')







