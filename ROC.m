%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ROC plot
%Chunyang Dong
%Tian Lab, UC Davis
%8/13/2020
%
%Summary: calculate frequency distribution of designated "Stimulated" region
%and "baseline" region. Apply threshold across histograms from right to
%left and plot true detection rate vs. false positive rate. Calculate d'
%with 2 different formulas, and calcuate area under the ROC curve. Use this
%script after FP analysis. 
%
%   Inputs
%       -DFF, time
%       
%
%   Outputs
%       -ROC curve, d', area under ROC curve.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%Inputs and set time phrame of stim, no stim,template%%%%%%%%%%
%load('/Volumes/Drive/Tian Lab/Data Analysis/FP_analysis/Single trail analysis and d prime/analsis/roc/BNST-.mat')
time=time; %CHANGE THIS, vector of times
DFF=DFF; %CHANGE THIS, where columns are timepoints, rows are trials

%Low pass filter cutoff freq=1/2 Hz
dt = mean(diff(time)); 
samplingrate = 1/dt;
LPfilteredDFF = lowpass(DFF,1/2,samplingrate);
DFF=LPfilteredDFF;

t_nostim=9000:11500; %CHANGE THIS, in the form of a range of column numbers where the 
                        %response is not occurring, i.e. baseline

t_stim= 50100:52600; %CHANGE THIS, in the form of  a range of column numbers where the 
                     %response is occurring

%%Look at Average DFF
averageDFF=mean(DFF, 1); %calculate average of DFF
plot(averageDFF) %plot index versus averageDFF to determine template indices
                 %you want
in_template= 50710:55810; %CHANGE THIS to pick a range of column numbers for the 
                         %template

templateav=averageDFF(in_template)-min(averageDFF(in_template)); %establishes template
nostimregion=averageDFF(t_nostim);
stimregion=averageDFF(t_stim);
                                   
figure %plots template overlay on averageDFF
subplot(3,1,1)
title('template')
hold on
plot(time, averageDFF)
plot(time(in_template), templateav, 'g')
hold off
subplot(3,1,2)
title('no stim')
hold on
plot(time, averageDFF)
plot(time(t_nostim), nostimregion, 'm')
hold off
subplot(3,1,3)
title('stim')
hold on
plot(time, averageDFF)
plot(time(t_stim), stimregion, 'r')
hold off

%%%%%%%%%%%%%%%%%%%ROC Analysis%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%


subplot(3,3, [1 2]) %plots raw traces and regions of stim or nostim
plot(time,DFF, '-k')
hold on
plot(time(t_nostim),(DFF(:,t_nostim)),'r.') %red traces out where there is no stim response
plot(time(t_stim),(DFF(:,t_stim)),'b.')
xlabel('Time (s)')
ylabel('dF/F')
hold off

subplot(3, 3, 3) %plots histogram
hold on
histogram(DFF(:,t_nostim),100, 'facecolor', 'r', 'edgecolor', 'r')
histogram(DFF(:,t_stim),100, 'facecolor','b', 'edgecolor', 'b')
xlabel('dF/F')
ylabel('occurence')
hold off
%plot([1 1]*mean(dFpeak),[0 80],'g--')
%xlim([-.1 0.1])
%SNR_dF = dFpeak/std(dF(t_nostim)) % this number will be close to 5

%This step convolves each row with the template
[n m]=size(DFF);
dF_filtmat=zeros(n, m);
for num=1:n
    dFtemp=DFF(num,:);
    dF_tempfilt= conv(dFtemp, templateav, 'same');
    dF_filtmat(num,:)= dF_tempfilt;
end
%dF_filtmat now holds the filtered traces


subplot(3,3, [4 5]) %plots filtered traces
hold on
plot(time, dF_filtmat')
plot(time(t_nostim),dF_filtmat(:,t_nostim),'r.')
plot(time(t_stim),dF_filtmat(:,t_stim), '-b')
xlabel('Time (s)')
ylabel('filtered dF/F')
hold off

subplot(3, 3, 6) %plots histogram of filtered traces
hold on
histogram(dF_filtmat(:,t_nostim),100, 'facecolor', 'r', 'edgecolor', 'r')
histogram(dF_filtmat(:,t_stim),100, 'facecolor','b', 'edgecolor', 'b')
xlabel('filtered dF/F')
ylabel('occurence')
hold off

 


%Create dF_filtmin, which holds locations of minima where response occurred
in_filtmin=zeros(n,1);
dF_minima=zeros(n,1);
A=dF_filtmat(:,t_stim);
for i=1:n
    currrow=A(i,:);
    actualrow=dF_filtmat(i,:);
    val = min(currrow);                                                        
    dF_minima(i)=val;
end
B=dF_filtmat(:,t_nostim);

dF_baseline=zeros(n,1);
for i=1:n
    currrow=B(i,:);
    actualrow=dF_filtmat(i,:);
    val = min(currrow);
    dF_baseline(i)=val;
end

SNR_dF_filt = abs(mean(dF_minima)/std(dF_baseline(:)))

subplot(3,3,[7 8]) %plots histograms of peak
hold on
h1=histogram(dF_baseline, 'facecolor', 'r', 'edgecolor', 'r');
h2=histogram(dF_minima, 'facecolor', 'g', 'edgecolor', 'g');
h1.Normalization = 'probability';
h1.BinWidth = 100;                                                              % number of bins, smaller number smaller the bar is 
h2.Normalization = 'probability';
h2.BinWidth = 100;
xlabel('filtered dF/F')
ylabel('probabiliy')
hold off


%%%%%
threshold = -50000:0.1:50000;                                                    %%range of the two histograms, and number of steps
dF_nostim = dF_baseline;
dF_stim = dF_minima;

nn = length(threshold);
false_pos = zeros(1,nn);
false_neg = zeros(1,nn);

%normalize=n.*length(t_nostim);
for i = 1:length(threshold)
    false_pos(i) = length(find(dF_nostim < threshold(i)))/n;        
    false_neg(i) = length(find(dF_stim > threshold(i)))/n;                   
end
subplot(3,3,9)
plot(false_pos,1-false_neg)
ylabel('True detection rate')
xlabel('False positive rate')
axis equal
axis([0 1 0 1])

Area=trapz(false_pos, 1-false_neg)                           
dPrimeFromArea = sqrt(2)*norminv(Area)

f1 = figure; 
plot(false_pos,1-false_neg)
ylabel('True detection rate')
xlabel('False positive rate')
axis equal
axis([-0.1 1.1 -0.1 1.1])
