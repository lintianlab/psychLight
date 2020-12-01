%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Fiber photometry analysis for in vivo drug application
%Chunyang Dong
%Tian Lab, UC Davis
%9/4/2020
%
%Summary: Fit 405 to 465, calculate dff
%
%   Inputs
%       -digitone
%       -csv from doric system
%
%   Outputs
%       -dff, time, plots. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% doric photometry data basic script
% import data A2:E3403050 and select 'Numeric Matrix' from drop down menu.

data = animal3(:,:);

%%
%create separate matrices for 465nm and 405nm excitation channels
sig = data(:,1:3);
sig(:,2) = [];

iso = data(:,1:2);

%plot demodulated signal from 465nm and 405nm excitation channels
f1 = figure;
subplot(2,1,1);
plot(sig(:,1),sig(:,2));
title('Demodulated Green Sensor Signal')
xlabel('seconds')
ylabel('F')

subplot(2,1,2);
plot(iso(:,1),iso(:,2));
title('Demodulated Isosbestic Signal')
xlabel('seconds')
ylabel('F')

%% linearize & correct 465nmsignal using best fit line to 405nm data
temp_fit2 = fit(iso(:,1),iso(:,2),'exp2');
fit2 = fitlm(temp_fit2(sig(:,1)),sig(:,2));

dff(:,1) = sig(:,1);
dff(:,2) = (100*(sig(:,2)-(fit2.Fitted))./(fit2.Fitted));

% plotting linearized data
f2 = figure;
subplot(3,1,1);
plot(iso(:,1),iso(:,2));
hold on
plot(temp_fit2);
xlabel('sec')
ylabel('F)')
title('Fitting Isosbestic Signal with Biexponential Curve')

subplot(3,1,2);
plot(iso(:,1),(iso(:,2)-(fit2.Fitted))./(fit2.Fitted));
xlabel('sec')
ylabel('F')
title('Linearized Isosbestic Signal')
% 
subplot(3,1,3);
plot(dff(:,1),dff(:,2));
xlabel('sec')
ylabel('dF/F (%)')
title('Linearized Green Sensor Signal')

%% Downsample
dsdff= downsample(dff,1000);

%% Low pass filter
LPfilteredDFF = lowpass(dsdff,1/3);


