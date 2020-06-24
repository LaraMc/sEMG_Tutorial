%%
% To run the examples, first place all supplementary EMG files in the same
% folder, and then change the current Matlab folder to this folder, using 
% the command "cd" with the address of the folder:
% e.g.
% cd('C:\Documents\MyFolder\')
% Then highlight the section of code you want to run, right-click and "Evaluate"

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Key Functions for Extracting Surface EMG Features
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   %%	Function 1: Surface EMG Power Spectral Density (Surface EMG in the Frequency Domain)
%
%   %%	Function 2: Low-pass Filtering Surface EMG
%
%   %%  Function 3: Surface EMG amplitude features - Calculate root-mean-square and average rectified value of EMG signal 
%
%   %%  Function 4: Surface EMG spectral features - Calculate the Mean Frequency and the Median Frequency of the surface EMG signal
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%	Function 1: Surface EMG Power Spectral Density (Surface EMG in the Frequency Domain)

clear all
pkg load signal                                         %Load Package containing some functions needed for the following lines

% load surface EMG data
load('Tutorial_EMG_sig.mat')
% or specify the full address of the directory the file is contained in  
% e.g.
% load('C:\Documents\MyFolder\Tutorial_EMG_sig.mat')
% or drag and drop file into the Matlab/Octave command window to load file

fs=2000;                                                %Sampling frequency
dt=1/fs;                                                %Time-step between samples
t=[0:length(EMG_sig_2k)-1].*dt;                         %Time vector of sampled signal

%Plot 5 seconds of the signal
T_period=5;                                             %Time period of the signal, 5 seconds

%Crop signal
EMG_sig_5s=EMG_sig_2k(find(t>=2,1,'first'):find(t<2+T_period,1,'last'),:);

figure;
subplot(1,2,1)
plot([0:length(EMG_sig_5s)-1].*dt,EMG_sig_5s)           %Plot time-domain signal

N=length(EMG_sig_5s);                                   %Number of samples in the signal, 10000
xlim([0 5])                                             %Show signal between 0 and 5 seconds
ylim([-1.5 1.5])
xlabel('Time (s)')
ylabel('mV')
title('Surface EMG - Time Domain')

%L is the length of the subsections of the total signal that will be
%analysed and averaged
L=round(0.5/dt);

% The number of points in the discrete Fourier Transform is chosen as the
% next power of 2 greater than the number of samples L
nfft = 2^nextpow2(L);                                    %i.e. 2^10 = 1024

%Note: this means that 24 zeros are appended to each signal epoch in the time
%domain to create a frequency bin spacing of fs/nfft = 0.98 Hz - approx 1Hz

%Overlap of successive signal epochs or segments
overlap = 0.5;                                           %This is equal to 50% overlap

%The number of samples in the overlap between successive epochs 
noverlap=round(overlap*nfft);

%The number of (overlapping) segments in the signal 
D=L-noverlap;

%Each successive segment starts D samples after the previous segment
N_avg=floor(((N-L)/D)+1);

%Calculate the power spectral density using a hamming window of length L
%(hamming window is used by default in pwelch), alternative windows are
%"hann" and "nuttallwin"
[Pxx,f]=pwelch(EMG_sig_5s-mean(EMG_sig_5s),hamming(L),overlap,nfft,fs);         % Jeremy Edit: overlap must be between 0 and 0.95
subplot(1,2,2)
plot(f,Pxx)
xlim([0 500])
ylim([0 3.75*1e-4])
title('Surface EMG - Frequency Domain')
xlabel('Frequency (Hz)')
ylabel('\muV^2/Hz')

%%	Function 2: Low-pass Filtering Surface EMG

clear all;
pkg load signal                                         %Load Package containing some functions needed for the following lines

% load surface EMG data
load('Tutorial_EMG_sig.mat')
fs=2000;                                                %Sampling frequency
dt=1/fs;                                                %Time-step between samples

t=[0:length(EMG_sig_2k)-1].*dt;                         %Time vector of sampled signal

% Low-pass filter example

%If the EMG signal was to be sampled at 500 Hz, it would first need to
%low-pass filtered at half that frequency (250 Hz) to remove
%the high-frequency components and prevent aliasing of the signal.

fc = 250;% Cutoff Frequency
FilterOrder = 4;   
ftype = 'low';

% Zero-Pole-Gain design, Butterworth filter
[b,a]=butter(FilterOrder,fc./(fs/2),ftype);        %Low-pass filter design

EMG_sig_2k_lowpass=filtfilt(b,a,EMG_sig_2k);

figure;
subplot(1,2,1)
plot(t,EMG_sig_2k)                                       %Plot original (unfiltered) signal
hold on
plot(t,EMG_sig_2k_lowpass,'-.')         %Plot low-pass filtered signal Jeremy Edit: Signal needs to be scaled
xlim([5 5.25])
ylim([-1.15 1.15])
legend({'Original Signal';'Low-Pass Filtered Signal'})
xlabel('Time (s)')
ylabel('mV')
title('Surface EMG - Time Domain')

% The low-pass filter shapes the power spectral density of the EMG signal

L=round(0.5/dt);                                        %Split the signal into subsections of 0.5 seconds
nfft = 2^nextpow2(L);                                   %Length of nfft window is the next power of 2 greater than L
overlap = 0.5;                                          %This is equal to 50% overlap
noverlap=round(overlap*nfft);                           %Overlap in number of samples

%Calculate the power spectral density of the original signal using Welch's method
[Pxx,f]=pwelch(EMG_sig_2k-mean(EMG_sig_2k),hamming(L),overlap,nfft,fs);
subplot(1,2,2)
plot(f,Pxx)
hold on;
%Calculate the power spectral density of the low-pass filtered signal using Welch's method
[Pxx_lowpass,f]=pwelch(EMG_sig_2k_lowpass-mean(EMG_sig_2k_lowpass),hamming(L),overlap,nfft,fs);
plot(f,Pxx_lowpass,'-.')
legend({'Original Signal';'Low-Pass Filtered Signal'})
xlim([0 500])
xlabel('Frequency (Hz)')
ylabel('\muV^2/Hz')
title('Surface EMG - Frequency Domain')

%% Function 3: Surface EMG amplitude features - Calculate root-mean-square and average rectified value of EMG signal

clear all;
pkg load signal                                         %Load Package containing some functions needed for the following lines

% load surface EMG data

load('Tutorial_AmpChange_EMG.mat')

fs=2000;                                                %Sampling frequency
dt=1/fs;                                                %Time-step between samples
t=[0:length(EMG_AmpChange_sig_2k)-1].*dt;               %Time vector of sampled signal

% Plot the surface EMG signal, then calculate the RMS and ARV of the signal
% for 3 different segments

figure;
plot(t,EMG_AmpChange_sig_2k)                            %Plot time-domain signal

ylim([-1.5 1.5])
xlim([0 65])
xlabel('Time (s)')
ylabel('Amplitude (mV)')
title('Raw sEMG Signal')

Segment1=[10/dt:((15/dt)-dt)];                          %Segment 1 is from 10 seconds to 15 seconds
Segment2=[30/dt:((35/dt)-dt)];                          %Segment 2 is from 30 seconds to 35 seconds
Segment3=[45/dt:((50/dt)-dt)];                          %Segment 3 is from 45 seconds to 50 seconds

%Calculate the root-mean-square (RMS) value for each segment and print

RMS1=sqrt((1/length(Segment1)).*sum(EMG_AmpChange_sig_2k(Segment1).^2));
%Alternatively, use built-in Matlab function: RMS1=rms(EMG_AmpChange_sig_2k(Segment1));
RMS2=sqrt((1/length(Segment2)).*sum(EMG_AmpChange_sig_2k(Segment2).^2));
RMS3=sqrt((1/length(Segment3)).*sum(EMG_AmpChange_sig_2k(Segment3).^2));

fprintf(1, 'RMS-EMG Segment 1: %f\nRMS-EMG Segment 2: %f\nRMS-EMG Segment 3: %f\n', [RMS1;RMS2;RMS3]);

%Calculate the average rectified value (ARV) value for each segment and print

ARV1=(1/length(Segment1)).*sum(abs(EMG_AmpChange_sig_2k(Segment1)));
ARV2=(1/length(Segment2)).*sum(abs(EMG_AmpChange_sig_2k(Segment2)));
ARV3=(1/length(Segment3)).*sum(abs(EMG_AmpChange_sig_2k(Segment3)));

fprintf(1, 'ARV-EMG Segment 1: %f\nARV-EMG Segment 2: %f\nARV-EMG Segment 3: %f\n', [ARV1;ARV2;ARV3]);

% To calculate the change in RMS or ARV of the signal over time, see
% Example (xv) in Tutorial Code

%% Function 4: Surface EMG spectral features - Calculate the Mean Frequency and the Median Frequency of the surface EMG signal

clear all;
pkg load signal                                         %Load Package containing some functions needed for the following lines

% EMG signals from the start and end of a fatiguing isometric contraction
load('Tutorial_Fatigue_EMG.mat')
fs=2000;
dt=1/fs;

%Subtract mean from EMG signal (DC voltage offset)
EMG_Start=EMG_Start-mean(EMG_Start);
EMG_End=EMG_End-mean(EMG_End);

%Calculate power spectra

L=round(0.25/dt);                                       %Split the signal into subsections of 0.25 seconds
nfft = 2^nextpow2(L);                                   %Length of nfft window is the next power of 2 greater than L
overlap = 0.5;                                          %This is equal to 50% overlap
noverlap=round(overlap*nfft);                           %Overlap in number of samples

%Calculate the power spectral density at the start of EMG signal using Welch's method
[P_Start,f1]=pwelch(EMG_Start,hamming(L),overlap,nfft,fs);
%Calculate the power spectral density at the end of EMG signal using Welch's method
[P_End,f2]=pwelch(EMG_End,hamming(L),overlap,nfft,fs);

% Calculate the MEAN FREQUENCY

% return the frequency widths of each frequency bin
width1 = mean(diff(f1));
width2 = mean(diff(f2));

% multiply the PSD by the width to get the power within each bin
P1 = bsxfun(@times,width1,P_Start);
P2 = bsxfun(@times,width2,P_End);

% find all elements within the specified range
idx1 = find(0<=f1 & f1<=fs/2);
idx2 = find(0<=f2 & f2<=fs/2);

% compute the total power within the range
pwr1 = sum(P_Start(idx1,:));
pwr2 = sum(P_End(idx2,:));

% compute the central moment within the range
cen_moment1 = sum(bsxfun(@times,P_Start(idx1,:),f1(idx1))) ./ pwr1;
cen_moment2 = sum(bsxfun(@times,P_End(idx2,:),f2(idx2))) ./ pwr2;

Mean_freq_Start=cen_moment1;
Mean_freq_End=cen_moment2;

%Calculate the mean frequency (MNF) value for each segment and print
fprintf(1, 'Mean Freq. (Start): %.2f\nMean Freq. (End): %.2f\n', [Mean_freq_Start;Mean_freq_End]);

% Calculate the MEDIAN FREQUENCY

% Find the area under the power spectral density
PSDarea_Start = cumtrapz(f1, P_Start);
PSDarea_End = cumtrapz(f2, P_End); %% Jeremy Edit:f1 and f2 are the same but it's confusing to use the same for both
%Find the frequency that divides the EMG power spectrum into two equal
%regions, each containing half of the total power within the signal epoch
MedFreq_Start = interp1(PSDarea_Start, f1, PSDarea_Start(end)/2);
MedFreq_End = interp1(PSDarea_End, f2, PSDarea_End(end)/2);

%Calculate the median frequency (MDF) value for each segment and print
fprintf(1, 'Median Freq. (Start): %.2f\nMedian Freq. (End): %.2f\n', [MedFreq_Start;MedFreq_End]);

figure;
h1=plot(f1,P_Start,'r');                                %Plot power spectral density at the start of EMG signal
hold on;
%Plot a line to indicate the location the mean frequency for the section at the start of signal
h2=plot([Mean_freq_Start Mean_freq_Start],[0 max(P_Start)],'r','linewidth',2.5,'linestyle','--');

[P_End,f2]=pwelch(EMG_End,hamming(L),overlap,nfft,fs); %Calculate power spectral density at the end of EMG signal
h3=plot(f2,P_End,'b');                                  %Plot power spectral density at the end of EMG signal
%Plot a line to indicate the location the mean frequency for the section at the end of signal
h4=plot([Mean_freq_End Mean_freq_End],[0 max(P_End)],'b','linewidth',2.5,'linestyle','--');

ylim([0 2*1e-3])
xlim([0 300])
xlabel('Frequency (Hz)')
ylabel('\muV^2/Hz')
legend([h1,h3,h2,h4],'Start of Contraction','End of Contraction','Mean Freq. (Start)','Mean Freq. (End)')
title('Surface EMG - Frequency Domain')
