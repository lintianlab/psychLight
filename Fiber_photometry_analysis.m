%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Fiber photometry analysis
%Chunyang Dong
%Tian Lab, UC Davis
%09/20/2019
%
%Summary: Fit 405 to 465, calculate dff, chop entire trace into 15 same
%length trace with the help of digitone file. 
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
% make sure 'delimited' is checked in the upper right hand corner.
load('/Volumes/Drive/Tian Lab/Data Analysis/FP_analysis/matlab/digitone_mat/digitoneday3noshockpart.mat')
s2tone =digitoneday3noshockpart;

data = day33(:,:);
%data(:,2:3) = [];

% manually enter frame # of first '0' value in digi column, e.g.
%digi(:,1) = digi(:,1)-digi(9697,1);
digi = data(:,1:2);
digi(:,1) = digi(:,1)-digi(5004,1);
data(:,1) = digi(:,1);

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


%% chopping up the data
% change out the timestamp filename if need be-----------------------------------------------------------------------------------------------

s2bl = 15000; %number of datapoints to look at before tone
s2seg_dur = 120000; %number of datapoints to look at after tone

s2470 = dff;

ERF = []; %necessary if you re-run this part of the code changing baseline or seg_dur
for tone = 1:length(s2tone)
        temp = find(s2470(:,1) >= s2tone(tone)); %find all FP data indices that are greater than or equal to keystroke timestampe
        temp = temp(1); %take the first one
    for fiber = 1
        ERF(fiber).data(:,tone) = s2470(temp-s2bl:temp+s2seg_dur,fiber+1); %save in nice neat datastructure
        
        
        
        %note, these data were already corrected -- so, theoretically, the
        %differences in baseline activity are real. if you'd like to baseline
        %correct these data, you can subtract the mean of the bs period

        ERF(fiber).data(:,tone) = ERF(fiber).data(:,tone) - mean(ERF(fiber).data(1:s2bl,tone));
    end
    clear puff fiber temp
end

fps = 1/mean(diff(s2470(:,1))); %note, this way of calculating FPS fails if you started and stopped your recording (from the driver box)

% looking at system frame rate error (since we are taking mean fps)

f5 = figure;
plot(diff(s2470(:,1)));
xlabel('frame #')
ylabel('Period in ms')
title('Period vs. Frame #')

ERF_time = linspace(-s2bl/fps,s2seg_dur/fps,length(ERF(1).data));

%plot chopped up data (subplot1) and heatmap (subplot2)

% ERF.data(:, 3) = [];
% ERF.data(:, 5) = [];
% ERF.data(:, 8) = [];

f6 = figure;
for fiber = 1
    subplot(2,1,fiber)
    plot(ERF_time,ERF(fiber).data);
    xlabel('Time Relative to Tone (seconds)')
    ylabel('dF/F')
    hold on
    plot(ERF_time,mean(ERF(fiber).data,2),'k','LineWidth',2)
    %xlim([-100 1000])
    %pbaspect([0.8 1 1])
end

for fiber = 1
    subplot(2,1,fiber+1)
    imagesc(ERF_time,1:length(s2tone),ERF(fiber).data')
    xlabel('Time Relative to Tone (seconds)')
    ylabel('#Trails')
    %xlim([0 100])
    %pbaspect([0.8 1 1])
end

%% looking at digital input from TTL box, i.e. when shock is happening

digbl = 1000; 
digseg_dur = 120000; 
superdigitone = s2tone;

DIG = []; %necessary if you re-run this part of the code changing baseline or seg_dur-----------------------------------------------------------------------------------------------
for tone = 1:length(superdigitone)
        temp = find(digi(:,1) >= superdigitone(tone)); %find all FP data indices that are greater than or equal to keystroke timestampe-----------------------------------------------------------------------------------------------
        temp = temp(1); %take the first one
    for fiber = 1
        DIG(fiber).data(:,tone) = digi(temp-digbl:temp+digseg_dur,fiber+1); %save in nice neat datastructure
        


        DIG(fiber).data(:,tone) = DIG(fiber).data(:,tone) - mean(DIG(fiber).data(1:digbl,tone));
    end
    clear puff fiber temp
end

fps = 1/mean(diff(digi(:,1))); 
DIG_time = linspace(-digbl/fps,digseg_dur/fps,length(DIG(1).data));

f6 = figure;

for fiber = 1
    subplot(2,1,fiber)
    plot(DIG_time,DIG(fiber).data)
    xlabel('Time Relative to Tone (seconds)')
    ylabel('dF/F')
    hold on
    plot(DIG_time,mean(DIG(fiber).data,2),'k','LineWidth',1)
    xlim([0 100])
end

for fiber = 1
    subplot(2,1,fiber+1)
    imagesc(DIG_time,1:length(superdigitone),DIG(fiber).data')
    xlabel('Time Relative to Tone (seconds)')
    ylabel('dF/F')
end

%% plot digi input (shock) over chopped up data
f7 = figure;
for fiber = 1
    subplot(2,1,fiber)
    plot(ERF_time,ERF(fiber).data);
    xlabel('Time Relative to Tone (seconds)')
    ylabel('dF/F')
    hold on
    plot(ERF_time,mean(ERF(fiber).data,2),'k','LineWidth',2)
    xlim([-20 100])
end
hold on
for fiber = 1
    
    plot(DIG_time,DIG(fiber).data)
    xlabel('Time Relative to Tone (seconds)')
    ylabel('dF/F')
    hold on
    plot(DIG_time,mean(DIG(fiber).data,2),'r','LineWidth',2)
    hold off
end

%%plot 465 trace over 405 trace raw

f8 = figure;
 
plot(sig(:,1),sig(:,2));
hold on
plot(iso(:,1),iso(:,2));
hold off
xlabel('seconds')
ylabel('F')


