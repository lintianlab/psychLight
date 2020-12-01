function h = Blue_makeDFF_SNRheatmap(DFF,time,maxWindow,cbarOn)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%makeDFF_SNRheatmap
%Chunyang Dong
%Tian Lab, UC Davis
%08/21/2020
%
%Summary: This function makes a side-by-side heatmap of the input traces 
%for two DFF trace sets. The heatmap varies from blue to yellow with 0 set as
%green. Use this after getting DFF and time from FP analysis script. 
%
%   Inputs
%       -DFF - first set of DFF traces. All sets should be arranged with
%           time along dimension 2
%       -time - time vector corresponding to traces
%       -nTimePlot - number of time points to input as tick mark labels
%       -maxWindow - two element vector specifying times to average over if
%           applicable
%       -cbarOn - str or char entry for the suptitle
%
%   Outputs
%       -h - axis handle for the heatmap
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Handle input number
if nargin == 0%run sample data through
    load('DR-_zscore.mat')
elseif nargin<3
    maxWindow = [0,size(DFF,2)];
    cbarOn = 1;
end



Normalization = DFF ./ min(DFF(50000:52352));
DFF=-Normalization;

%dt = mean(diff(time)); 
%samplingrate = 1/dt;
LPfilteredDFF = lowpass(DFF,1/15);
DFF=LPfilteredDFF;


%Process the data
%Convert NaN values to zeros
DFF(isnan(DFF)) = 0;
%Get indices of maximum deflection in specified window in descending order
goodTime = time>maxWindow(1) & time<maxWindow(2);
maxDFF = mean(DFF(:,goodTime),2);
[~,DFF_SortIdx] = sort(maxDFF,'Descend');

% %Do color limits for RWB heatmap
% %color limits based on the higher of the mins and maxes from both matrices
% clims = [min(DFF(:)),max(DFF(:))];
% bluePercent = abs(clims(1))/(clims(2) + abs(clims(1)));
% redPercent = clims(2)/(clims(2) + abs(clims(1)));
% %Stretch the colorspace over the clims with zero as white
% redCol = [linspace(0,1,floor(bluePercent*256)),ones(1,ceil(256*redPercent))];
% greenCol = [linspace(0,1,floor(bluePercent*256)),linspace(1,0,ceil(256*redPercent))];
% blueCol = [ones(1,floor(bluePercent*256)),linspace(1,0,ceil(256*redPercent))];
% RWBcolormap = [redCol;greenCol;blueCol]';
% 
% %Do xtick labels
% ordRange = round(log10(max(time) - min(time))+0.2);%less than 2ex? count by 1ex-1 else count by 1ex
% xVals = unique(round(time,-(ordRange-1)));
% xIdxs = interp1(time,1:size(DFF,2),xVals);

%Make plot
h = imagesc(DFF(DFF_SortIdx,:));
caxis([-1 1]);
xlim([2884 111931])
% colormap(RWBcolormap);
xlabel('Time relative to tone /sec');
ylabel('SingleTrail');
xticks([2884 15002 39234 63466 87699 111931 135001]);
xticklabels({'-10','0','20','40','60','80','100'});
%set(gca,'xticklabels',compose('%d',xVals));
if cbarOn == 1
    c = colorbar;
    if nargin == 0
        c.Label.String = 'Normalized Z-Score';
        %return;
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%single trail analysis
avgnostim = mean(DFF(:,43000:48000)'); %average no stimulus region 0-27sec
avgstim = mean(DFF(:,49533:51350)'); %average shock region 27-28.5sec
%avgdec = mean(DFF(:,49533:51350)'); %average decrease after shock region 28.5-30sec
meanstim = mean(avgstim); %mean values of signal during stimulus of all trails
meannostim = mean(avgnostim); %mean values of signal during no stimulus of all trails

stdstim = std(avgstim); %std of signal during stimulus
stdnostim = std(avgnostim);
SEMnostim = stdnostim/sqrt(length(avgnostim));

diff = meanstim - meannostim; % difference between stim and no stim
SNR = diff^2/SEMnostim^2;
disp(SNR)

end 
