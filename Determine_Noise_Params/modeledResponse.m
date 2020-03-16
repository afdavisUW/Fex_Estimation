function [qdot] = modeledResponse(time,q)
%This will be sued to calculate tyhe process noise covariance values for
%both states

global SpectrumData Param Data bodyData
FX = SpectrumData.body(1).FeMag;
FxPhase = SpectrumData.body(1).FePhase;
W = SpectrumData.W;
Amp = SpectrumData.amp;
Phi = SpectrumData.phase;
dW = SpectrumData.dW;

omega = SpectrumData.W;
Mag = SpectrumData.amp;
% Phi = SpectrumData.phase;
sum1 = 0; sum2 = 0; sum3 =0;

% [U1,U2,U3] = particleDynamics_spec(time);
Disp = Data.wave.incidentFPS;
DispTime = Data.wave.time;
temp = DispTime-time;
[val,index] = min(abs(temp));
U2 = Data.wave.incidentFPSVel(index);
F_exIRF = bodyData.body(1).FIRF;
IRFtime = bodyData.body(1).irfTime;
[exForce] = excitationConv(time,F_exIRF,IRFtime,Disp,DispTime);
bodyNumber = 1;


% qdot(1) = q(2); % Benchmark
% qdot(2) = 1/Param.totalInertia*(exForce - (Param.heaveMoorStiff+Param.hydrostaticCoeff)*q(1) - linearDrag*(q(2)-U2)*abs(q(2)-U2) - (Param.body(bodyNumber).C_r*q(3:6)+Param.body(bodyNumber).D_r*q(2)));
% qdot(3:6) = Param.body(bodyNumber).A_r*q(3:6)+Param.body(bodyNumber).B_r*q(2);

% qdot(1) = q(2); % linear coefficient
% qdot(2) = 1/Param.totalInertia*(exForce - (Param.heaveMoorStiff+Param.hydrostaticCoeff)*q(1) - Param.linearDrag*q(2) - 2*(Param.body(bodyNumber).C_r*q(3:6)+Param.body(bodyNumber).D_r*q(2)));
% qdot(3:6) = Param.body(bodyNumber).A_r*q(3:6)+Param.body(bodyNumber).B_r*q(2);

% qdot(1) = q(2); % includes quadradic drag
% qdot(2) = 1/Param.totalInertia*(exForce - (Param.heaveMoorStiff+Param.hydrostaticCoeff)*q(1) - Param.dragFactor*1*(q(2)-U2)*abs(q(2)-U2)- (Param.body(bodyNumber).C_r*q(3:6)+Param.body(bodyNumber).D_r*q(2)));
% qdot(3:6) = Param.body(bodyNumber).A_r*q(3:6)+Param.body(bodyNumber).B_r*q(2);

% qdot(1) = q(2); % all terms combined non relative quad drag
% qdot(2) = 1/Param.totalInertia*(exForce - (Param.heaveMoorStiff+Param.hydrostaticCoeff)*q(1) - Param.linearDrag*q(2) - Param.dragFactor*(q(2))*abs(q(2)) - 2*(Param.body(bodyNumber).C_r*q(3:6)+Param.body(bodyNumber).D_r*q(2)));
% qdot(3:6) = Param.body(bodyNumber).A_r*q(3:6)+Param.body(bodyNumber).B_r*q(2);

% qdot(1) = q(2); % all terms combined
% qdot(2) = 1/Param.totalInertia*(exForce - (Param.heaveMoorStiff+Param.hydrostaticCoeff)*q(1) - Param.linearDrag*q(2) - Param.dragFactor*(q(2)-U2)*abs(q(2)-U2) - 2*(Param.body(bodyNumber).C_r*q(3:6)+Param.body(bodyNumber).D_r*q(2)));
% qdot(3:6) = Param.body(bodyNumber).A_r*q(3:6)+Param.body(bodyNumber).B_r*q(2);

qdot(1) = q(2); % all terms combined added attenuation factor
qdot(2) = 1/Param.totalInertia*(exForce - (Param.heaveMoorStiff+Param.hydrostaticCoeff)*q(1) - Param.linearDrag*q(2) - Param.dragFactor*(q(2)-Param.attenFactor*U2)*abs(q(2)-Param.attenFactor*U2) - 2*(Param.body(bodyNumber).C_r*q(3:6)+Param.body(bodyNumber).D_r*q(2)));
qdot(3:6) = Param.body(bodyNumber).A_r*q(3:6)+Param.body(bodyNumber).B_r*q(2);


qdot = qdot';
end

