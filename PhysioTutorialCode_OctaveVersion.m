%% Examples to illustrate concepts outlined in McManus et al., 2020
%
% In the contents the relevant section of the paper that corresponds to the
% examples is noted, though the tutorial does not exactly follow the order 
% of the paper sections.
% where examples correspond to a figure, it is noted
%
%%
% To run the examples, first place all supplementary EMG files in the same
% folder, and then change the current Matlab folder to this folder, using 
% the command "cd" with the address of the folder:
% e.g.
% cd('C:\Documents\MyFolder\')
% Then highlight the section of code you want to run, right-click and "Evaluate"

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Contents
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%	%%	Section 3.1: Description of a signal (corresponding Section in paper)
%                   
%       Example (i) Plot a sine wave, see Figure 3(c)
%       Example (ii) Plot sine waves of different amplitudes, frequencies and phases, see Figure 4 
%
%   %%  Section 3.2: Frequency domain analysis
%
%       Example (iii) Plot sine/square wave in Time and Frequency Domain, see Figure 5
%
%   %%  Section 4.4: Sampling
%
%       Example (iv) Show effect of different sampling rates on the frequency spectrum, see Figure 8
%
%   %%  Section A.2: Zero-padding
%
%       Example (v) Show effect of zero-padding on the spacing between frequency bins, see Figure A1
%       Example (vi) Illustrate difference between frequency resolution and the spacing between frequency bins
%
%   %%	Section 5.1 (see also Section A.1): Power Spectral Density Estimation (Surface EMG in the Frequency Domain)
%
%       Example (vii) Calculate power spectrum for different epoch lengths and overlaps, see Figure 9
%
%       Example (viii) Effect of multiplying signal by different window functions, see Figure 9 (c)
%
%   %%	Section 5.2: Filtering
%
%       Example (ix) Illustrate the attenuation/roll-off characteristics of different filters, see Figure 10 (a)
%       Example (x) Examples of low-pass, high-pass, band-pass and notch filters, see Figure 10 (a)
%
%   %%	Section 4.2: Choice of Electrode: 
%
%       Example (xi) Demonstrate effect of Inter-Electrode Distance (IED), see Figure 6
%
%   %%	Section 4.3: Choice of Amplifier
%
%       Example (xii) Demonstrate signal amplification with differential filter, see Figure 7 (c)
%       Example (xiii) Calculate CMRR
%
%   %%  Section 5.3: Surface EMG amplitude features
%
%       Example (xiv) Plot gait data, see Figure 11
%       Example (xv) Calculate the Root-Mean Square value and the Average Rectified value of an EMG signal
%       Example (xvi) Calculate the average of the surface EMG signal for 3 separate segments
%
%   %% Section 5.4: Surface EMG spectral features
%
%       Example (xvii) Calculate the Mean Frequency and the Median Frequency of the surface EMG signal, see Figure 12
%       Example (xviii) Calculating the coherence between two EMG signals, Figure A2
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Section 3: Key concepts in signal processing

%% Section 3.1: Description of a signal

%% Example (i) Plot a sine wave (see Figure 3(c))
clear all;
A=1;                                                    %Amplitude of the sine wave
phi=0;                                                  %Phase of the sine wave
f=10;                                                   %Frequency of the sine wave

t=[0:0.001:0.35];                                       %Time vector of samples
y=A.*sin((2*pi*f*t)+phi);                               %Value of signal at each time-point           
figure;
plot(t,y)
xlabel('Time (s)')
ylabel('y(t)')

%%  Example (ii) Plot sine waves of different amplitudes, frequencies and phases, see Figure 4

clear all;
% Change amplitude
A1=1; %Amplitude of the sine wave 1
A2=0.5; %Amplitude of the sine wave 2
phi=0; %Phase of the sine wave
f=10; %Frequency of the sine wave
t=[0:0.001:0.35]; %Time vector of samples

figure;
plot(t,A1.*sin((2*pi*f*t)+phi));
hold on;
plot(t,A2.*sin((2*pi*f*t)+phi));
xlabel('Time (s)')
ylabel('y(t)')
legend({'A = 1'; 'A = 0.5'})

%Change frequency
A1=1;                                                    %Amplitude of the sine wave
phi=0;                                                   %Phase of the sine wave
f1=10;                                                   %Frequency of the sine wave 1
f2=60;                                                   %Frequency of the sine wave 2

figure;
plot(t,A1.*sin((2*pi*f1*t)+phi));
hold on;
plot(t,A1.*sin((2*pi*f2*t)+phi));
xlabel('Time (s)')
ylabel('y(t)')
legend({'f = 10'; 'f = 60'})

%Change phase
A1=1;                                                    %Amplitude of the sine wave
phi1=0;                                                  %Phase of the sine wave 1

%Delay signal by 90 degrees
phi2=-pi/2;                                              %Phase of the sine wave 2 in radians - pi = 180 degrees; 2*pi = 360 degrees

%Shift signal forward in time by 180 degrees
phi3=pi;                                                 %Phase of the sine wave 2 - pi = 180 degrees; 2*pi = 360 degrees
f=10;                                                    %Frequency of the sine wave

figure;
plot(t,A1.*sin((2*pi*f*t)+phi1));
hold on;
plot(t,A1.*sin((2*pi*f*t)+phi2));
plot(t,A1.*sin((2*pi*f*t)+phi3));
xlabel('Time (s)')
ylabel('y(t)')
legend({'\phi = 0'; '\phi = \pi/2 = 90^{\circ}'; '\phi = \pi = 180^{\circ}'})

%% Section 3.2: Frequency domain analysis

%%  Example (iii) Plot sine/square wave in Time and Frequency Domain, see Figure 5

%Time Domain
clear all
A=1;                                                    %Amplitude of the sine/triangle wave
phi=0;                                                  %Phase of the sine wave
f=10;                                                   %Frequency of the sine/triangle wave
phi2=pi/2;                                              %Phase of the triangle wave
T_period=2;                                             %Time period of the signal, 2 seconds
dt=0.001;                                               %Sampling period of the signal, samples taken every 0.001 seconds
t=[0:dt:T_period-dt];                                   %Time vector of samples
fs=1/dt;                                                %Sampling frequency, equal to 1/Sampling Period.

wav_sine=A.*sin((2*pi*f*t)+phi);
wav_triangle=sawtooth(2*pi*f*t+phi2,0.5);

figure;
subplot(1,2,1)
plot(t,wav_sine);
xlabel('Time (s)')
ylabel('y(t)')
legend({'Sine'})
title('Time Domain - Sine Wave')
xlim([0 0.35])                                           %Only show 0.35 seconds

subplot(1,2,2)
plot(t,wav_triangle,'r');
xlabel('Time (s)')
ylabel('y(t)')
legend({'Triangle'})
title('Time Domain - Triangle Wave')
xlim([0 0.35])

%Frequency Domain
Y1 = fft(wav_sine);                                       %computes the discrete Fourier transform using a fast Fourier transform (FFT) algorithm.
Y2 = fft(wav_triangle);

L1=length(wav_sine);                                      %Length of the sine wave in samples
L2=length(wav_triangle);                                  %Length of the triangle wave in samples

%Generate the two-sided amplitude frequency spectrum
P1_2side = abs(Y1/L1);
P2_2side = abs(Y2/L2);

%Generate the one-sided amplitude frequency spectrum
P1_1side = P1_2side(1:L1/2+1);
P1_1side(2:end-1) = 2*P1_1side(2:end-1);

P2_1side = P2_2side(1:L2/2+1);
P2_1side(2:end-1) = 2*P2_1side(2:end-1);

%Create the frequency vector for plotting
f1 = fs*(0:(L1/2))/L1;
f2 = fs*(0:(L2/2))/L2;

%Generate the amplitude spectrum
figure;
subplot(1,2,1)
plot(f1,P1_1side) 
 ylim([0 1.05])
 xlim([0 100])
xlabel('Frequency (Hz)')
ylabel('|Y(f)|')
title('Frequency Domain - Amplitude Spectrum Sine Wave')

subplot(1,2,2)
plot(f2,P2_1side,'r') 
 xlim([0 100])
 ylim([0 1.05])
xlabel('Frequency (Hz)')
ylabel('|Y(f)|')
title('Frequency Domain - Amplitude Spectrum Triangle Wave')

%Generate the Power spectrum
figure;
subplot(1,2,1)
plot(f1,P1_1side.^2) 
 ylim([0 1.05])
 xlim([0 100])
xlabel('Frequency (Hz)')
ylabel('|Y(f)|^{2}')
title('Frequency Domain - Power Spectrum Sine Wave')

subplot(1,2,2)
plot(f2,P2_1side.^2,'r') 
 xlim([0 100])
 ylim([0 1.05])
xlabel('Frequency (Hz)')
ylabel('|Y(f)|^{2}')
title('Frequency Domain - Power Spectrum Triangle Wave')

%Note: see Example (vii) for time and frequency domain representations of 
%real EMG signals.

%% Section 4.4: Sampling

%%  Example (iv) Show effect of different sampling rates on the frequency spectrum, see Figure 8

%Time Domain
clear all
A=1;                                                    %Amplitude of the sine wave
phi=-pi/2;                                              %Phase of the sine wave
f=15;                                                   %Frequency of the sine wave
T_period=2;                                             %Time period of the signal, 2 seconds

%Approximate continuous signal, very high sampling rate 10000Hz
dt=0.0001;                                              %Sampling period of the signal, samples taken every 0.001 seconds
t=[0:dt:T_period-dt];                                   %Time vector of samples
fs=1/dt;                                                %Sampling frequency, equal to 1/Sampling Period.
sine_orig=A.*sin((2*pi*f*t)+phi);                       %Original (continuous signal)

%Highest frequency in the signal is 15Hz so must be sampled at 30 samples/s
%or more

fs_good=40;                                             %Adequate sampling frequency (40 Hz)
dt_good=1/fs_good;
fs_bad=25;                                              %Low Sampling frequency (25 Hz) (distorts signal)
dt_bad=1/fs_bad;

t_good=[0:dt_good:T_period-dt];                         %Time vector of samples
sine_good=A.*sin((2*pi*f*t_good)+phi);

%Frequency Domain
figure;
subplot(1,2,1)
plot(t,sine_orig)
hold on
stem(t_good,sine_good)
xlim([0 0.5])
ylim([-1.5 1.5])
legend({'Original Signal';'Signal Sampled at 40Hz'})
xlabel('Time (s)')
ylabel('y(t)')

%Frequency Domain
Y1 = fft(sine_good);                                    %computes the discrete Fourier transform using a fast Fourier transform (FFT) algorithm.
L1=length(sine_good);                                   %Length of signal in samples

P1_2side = abs(Y1/L1);                                  %Compute the two-sided amplitude spectrum

P1_1side = P1_2side(1:L1/2+1);                          %Compute the single-sided amplitude spectrum
P1_1side(2:end-1) = 2*P1_1side(2:end-1);

f1 = fs_good*(0:(L1/2))/L1;                             %Vector of all the frequencies contained in the spectrum

subplot(1,2,2)
plot(f1,P1_1side)                                       %Plot the amplitude spectrum
text(.05 ,0.9,['Correct: 15 Hz'],'Units','normalized','Interpreter','none')
xlabel('Frequency (Hz)')
ylabel('|Y(f)|')
ylim([0 1.25])
xlim([0 20])

%Time Domain
t_bad=[0:dt_bad:T_period-dt];                            %Time vector of samples
sine_bad=A.*sin((2*pi*f*t_bad)+phi);                     %Plot signal at each time-point      
figure;
subplot(1,2,1)
plot(t,sine_orig)
hold on
stem(t_bad,sine_bad)
legend({'Original Signal';'Signal Sampled at 25Hz'})
xlim([0 0.5])
ylim([-1.5 1.5])
xlabel('Time (s)')
ylabel('y(t)')

%Frequency Domain

Y2 = fft(sine_bad);                                       %computes the discrete Fourier transform using a fast Fourier transform (FFT) algorithm.
L2=length(sine_bad);                                      %Length of signal in samples
P2_2side = abs(Y2/L2);                                    %Compute the two-sided amplitude spectrum
P2_1side = P2_2side(1:L2/2+1);                            %Compute the single-sided amplitude spectrum
P2_1side(2:end-1) = 2*P2_1side(2:end-1);
f2 = fs_bad*(0:(L2/2))/L2;                                %Vector of all the frequencies contained in the spectrum

subplot(1,2,2)
plot(f2,P2_1side,'r') 
xlabel('Frequency (Hz)')
ylabel('|Y(f)|')
ylim([0 1.25])
xlim([0 20])
text(.05 ,0.9,['Incorrect: 10 Hz'],'Units','normalized','Interpreter','none')
xlabel('Frequency (Hz)')
ylabel('|Y(f)|')

%% Section A.2: Zero-padding

%%  Example (v) Show effect of zero-padding on the spacing between frequency bins, see Figure A1

%Time Domain
clear all
A=1;                                                    %Amplitude of the sine wave
phi=0;                                                  %Phase of the sine wave
f1=8;                                                   %Frequency of the sine wave
f2=30;                                                  %Frequency of the sine wave
T_period=0.25;                                          %Time period of the signal, 0.25 seconds

dt=0.001;                                               %Sampling period of the signal, samples taken every 0.001 seconds
t=[0:dt:T_period-dt];                                   %Time vector of samples
fs=1/dt;                                                %Sampling frequency, equal to 1/Sampling Period.

%A sinusoid that consists of a 8Hz sine wave and a 30 Hz sine wave
sine_orig=A.*sin((2*pi*f1*t)+phi)+A.*sin((2*pi*f2*t)+phi);

%Frequency Resolution is 1/T_period = 4 Hz

%Frequency Domain
figure;
subplot(1,2,1)
plot(t,sine_orig)
hold on
plot(t,sine_orig,'o','MarkerSize',4)
xlim([0 0.25])
ylim([-2.5 2.5])
xlabel('Time (s)')
ylabel('y(t)')
title('Sinusoid comprised of 2 sine waves, both A = 1')

%Frequency Domain
nfft=length(sine_orig);
Y1 = fft(sine_orig,nfft);                               %computes the discrete Fourier transform using a fast Fourier transform (FFT) algorithm.
L1=length(sine_orig);                                   %Length of signal in samples
P1_2side = abs(Y1/L1);                                  %Compute the two-sided spectrum
P1_1side = P1_2side(1:L1/2+1);                          %Compute the single-sided spectrum
P1_1side(2:end-1) = 2*P1_1side(2:end-1);
f1 = fs*(0:(L1/2))/L1;                                  %Vector of all the frequencies contained in the spectrum
subplot(1,2,2)
plot(f1,P1_1side) 
hold on;
plot(f1,P1_1side,'g.','MarkerSize',20)

%Correct amplitude level is 1
plot([0 100],[1 1],'k--')
xlabel('Frequency (Hz)')
ylabel('|Y(f)|')
ylim([0 1.25])
xlim([-1 50])
title('Incorrectly Estimated Amplitude, \DeltaR_{NFFT} = 4Hz')

%Increase the length of the FFT window used by adding zeros to the signal
%Frequency Domain
sine_orig_zeropad=[zeros(1,round(length(sine_orig)/2)) sine_orig zeros(1,round(length(sine_orig)/2))];
t_zeropad=[0:length(sine_orig_zeropad)-1].*dt;          %zeropad time vector

figure;
subplot(1,2,1)
plot(t_zeropad,sine_orig_zeropad)
hold on
plot(t_zeropad,sine_orig_zeropad,'o','MarkerSize',4)
xlim([0 0.5])
ylim([-2.5 2.5])
xlabel('Time (s)')
ylabel('y(t)')

%Frequency Domain
nfft=length(sine_orig)*2;
%Setting nfft to a length greater than the length of the original signal
%will pad the original signal with trailing zeros until the total number of samples in the
%signal reaches length nfft. This is equivalent to calculating the fft of 'sine_orig_zeropad'

Y1 = fft(sine_orig,nfft);                               %computes the discrete Fourier transform using a fast Fourier transform (FFT) algorithm. 
L1=length(sine_orig);                                   %Length of signal in samples
P1_2side = abs(Y1/L1);                                  %Compute the two-sided spectrum
P1_1side = P1_2side(1:nfft/2+1);                        %Compute single-sided spectrum
P1_1side(2:end-1) = 2*P1_1side(2:end-1);
f1 = fs*(0:(nfft/2))/nfft;                              %Vector of all the frequencies contained in the spectrum
subplot(1,2,2)
plot(f1,P1_1side) 
hold on;
plot(f1,P1_1side,'g.','MarkerSize',20)

%Correct amplitude level is 1
plot([0 100],[1 1],'k--')
xlabel('Frequency (Hz)')
ylabel('|Y(f)|')
ylim([0 1.25])
xlim([-1 50])
title('Correctly Estimated Amplitude, \DeltaR_{NFFT} = 2Hz')

%% Example (vi) Illustrate difference between frequency resolution and the spacing between frequency bins (Advanced Topic)

clear all
A=1;                                                    %Amplitude of the sine wave
phi=0;                                                  %Phase of the sine wave
f1=30;                                                  %Frequency of the sine wave
f2=32;                                                  %Frequency of the sine wave
T_period=0.25;                                          %Time period of the signal, 0.25 seconds

dt=0.001;                                               %Sampling period of the signal, samples taken every 0.001 seconds
t=[0:dt:T_period-dt];                                   %Time vector of samples
fs=1/dt;                                                %Sampling frequency, equal to 1/Sampling Period.

%A sinusoid that consists of a 30Hz sine wave and a 32 Hz sine wave
sine_orig=A.*sin((2*pi*f1*t)+phi)+A.*sin((2*pi*f2*t)+phi);

%Frequency Resolution is 1/T_period = 4 Hz

%Frequency Domain - cannot distinguish between 30Hz and 32Hz sine waves
figure;
nfft=length(sine_orig);                                 %Length of the signal that the FFT will be applied to
Y1 = fft(sine_orig,nfft);                               %computes the discrete Fourier transform using a fast Fourier transform (FFT) algorithm.    
L1=length(sine_orig);                                   %Length of signal in samples
P1_2side = abs(Y1/L1);                                  %Compute the two-sided spectrum
P1_1side = P1_2side(1:L1/2+1);                          %Compute single-sided spectrum
P1_1side(2:end-1) = 2*P1_1side(2:end-1);
freq = fs*(0:(L1/2))/L1;                                %Vector of all the frequencies contained in the spectrum
subplot(2,2,1)
plot(freq,P1_1side) 
hold on;
plot(freq,P1_1side,'g.','MarkerSize',20)
xlabel('Frequency (Hz)')
ylabel('|Y(f)|')
ylim([0 1.5])
xlim([-1 50])
text(.05 ,0.85,['\DeltaR = 4Hz'],'Units','normalized')  %Frequency resolution
text(.05 ,0.75,['L = 0.25 s'],'Units','normalized')     %Signal length, i.e. T_period
%Before Zero-Padding, frequency resolution is equal to the spacing between frequency bins
title('Before Zero-Padding, \DeltaR = \DeltaR_{NFFT} = 4Hz')

%Correct amplitude level is 1
plot([0 100],[1 1],'k--')

%Zero-pad signal - still cannot distinguish between 30Hz and 32Hz sine waves

%Length of the signal that the FFT will be applied to, i.e. twice the
%length of the original signal, so the original signal will be zero-padded
nfft=length(sine_orig)*2;

Y1 = fft(sine_orig,nfft);                               %computes the discrete Fourier transform using a fast Fourier transform (FFT) algorithm.     
L1=length(sine_orig);                                   %Length of original signal in samples
P1_2side = abs(Y1/L1);                                  %Compute the two-sided spectrum
P1_1side = P1_2side(1:nfft/2+1);                        %Compute single-sided spectrum
P1_1side(2:end-1) = 2*P1_1side(2:end-1);
freq = fs*(0:(nfft/2))/nfft;                            %Vector of all the frequencies contained in the spectrum
subplot(2,2,2)
plot(freq,P1_1side) 
hold on;
plot(freq,P1_1side,'g.','MarkerSize',20)
xlabel('Frequency (Hz)')
ylabel('|Y(f)|')
ylim([0 1.5])
xlim([-1 50])
text(.05 ,0.85,['\DeltaR = 4Hz'],'Units','normalized')
text(.05 ,0.75,['L = 0.25 s'],'Units','normalized')     %Signal length, i.e. T_period
%After Zero-Padding, frequency resolution is still 4Hz but the spacing between frequency bins is 2 Hz
title('After Zero-Padding, \DeltaR_{NFFT} = 2 Hz')

%Correct amplitude level is 1
plot([0 100],[1 1],'k--')

%Increase the length of the signal to 0.5 s so that the frequency
%resolution is 2Hz

T_period=0.5;                                           %Time period of the signal, 0.5 seconds
t_2HzRes=[0:dt:T_period-dt];                            %Time vector of samples
fs=1/dt;                                                %Sampling frequency, equal to 1/Sampling Period.

%A sinusoid that consists of a 30Hz sine wave and a 32 Hz sine wave
sine_2HzRes=A.*sin((2*pi*f1*t_2HzRes)+phi)+A.*sin((2*pi*f2*t_2HzRes)+phi);

nfft=length(sine_2HzRes);                               %Length of the signal that the FFT will be applied to
Y1 = fft(sine_2HzRes,nfft);                             %computes the discrete Fourier transform using a fast Fourier transform (FFT) algorithm. 
L1=length(sine_2HzRes);                                 %Length of original signal in samples
P1_2side = abs(Y1/L1);                                  %Compute the two-sided spectrum
P1_1side = P1_2side(1:L1/2+1);                          %Compute single-sided spectrum
P1_1side(2:end-1) = 2*P1_1side(2:end-1);
freq = fs*(0:(L1/2))/L1;                                %Vector of all the frequencies contained in the spectrum
subplot(2,2,3)
plot(freq,P1_1side) 
hold on;
plot(freq,P1_1side,'g.','MarkerSize',20)
xlabel('Frequency (Hz)')
ylabel('|Y(f)|')
ylim([0 1.5])
xlim([-1 50])
text(.05 ,0.85,['\DeltaR = 2Hz'],'Units','normalized')
text(.05 ,0.75,['L = 0.5 s'],'Units','normalized')      %Signal length, i.e. new T_period

%Increasing the length of the signal results in better frequency resolution (2 Hz)
title('Increase Signal Length, Decrease \DeltaR')

%Correct amplitude level is 1
plot([0 100],[1 1],'k--')

%Increase the length of the signal to 1 s so that the frequency
%resolution is 1Hz

T_period=1;                                             %Time period of the signal, 1 second
t_2HzRes=[0:dt:T_period-dt];                            %Time vector of samples
fs=1/dt;                                                %Sampling frequency, equal to 1/Sampling Period.

%A sinusoid that consists of a 30Hz sine wave and a 32 Hz sine wave
sine_1HzRes=A.*sin((2*pi*f1*t_2HzRes)+phi)+A.*sin((2*pi*f2*t_2HzRes)+phi);

nfft=length(sine_1HzRes);                               %Length of the signal that the FFT will be applied to
Y1 = fft(sine_1HzRes,nfft);                             %computes the discrete Fourier transform using a fast Fourier transform (FFT) algorithm. 
L1=length(sine_1HzRes);                                 %Length of original signal in samples
P1_2side = abs(Y1/L1);                                  %Compute the two-sided spectrum
P1_1side = P1_2side(1:L1/2+1);
P1_1side(2:end-1) = 2*P1_1side(2:end-1);
freq = fs*(0:(L1/2))/L1;                                %Vector of all the frequencies contained in the spectrum
subplot(2,2,4)
plot(freq,P1_1side) 
hold on;
plot(freq,P1_1side,'g.','MarkerSize',20)
xlabel('Frequency (Hz)')
ylabel('|Y(f)|')
ylim([0 1.5])
xlim([-1 50])
text(.05 ,0.85,['\DeltaR = 1Hz'],'Units','normalized')
text(.05 ,0.75,['L = 1 s'],'Units','normalized')        %Signal length, i.e. new T_period
%Increasing the signal length gives a better frequency resolution, 1Hz
title('Increase Signal Length Further, Decrease \DeltaR Further')

%Correct amplitude level is 1
plot([0 100],[1 1],'k--')

%% Section 5.1 (see also Section A.1): Power Spectral Density Estimation 

%% Example (vii) Calculate power spectrum for different epoch lengths and overlaps, see Figure 9

clear all
pkg load signal    
% Power spectrum is generated by averaging spectra of smaller signal epochs

load('Tutorial_EMG_sig.mat')
% or specify the full address of the directory the file is contained in  
% e.g.
% load('C:\Documents\MyFolder\Tutorial_EMG_sig.mat')
% or drag and drop file into the Matlab/Octave command window to load file

fs=2000;                                                %Sampling frequency
dt=1/fs;
t=[0:length(EMG_sig_2k)-1].*dt;
T_period=0.5;                                           %Time period of the signal, 0.5 seconds

% figure;plot(t,EMG_sig_2k)                             %Plot the whole EMG signal

%Plot 0.5 seconds of the signal
EMG_sig_short=EMG_sig_2k(find(t>=2,1,'first'):find(t<2+T_period,1,'last'),:);
figure;
subplot(1,2,1)
plot([0:length(EMG_sig_short)-1].*dt,EMG_sig_short)
xlim([0 0.5])
ylim([-1.5 1.5])
title('Time Domain')
xlabel('Time (s)')
ylabel('mV')

%Frequency Domain

%Typically the mean of the signal is subtracted from the total signal as this 
%voltage offset is generated by electrical noise
Y1 = fft(EMG_sig_short-mean(EMG_sig_short));            %Calculate the discrete Fourier Transform of the signal
L1=length(EMG_sig_short);                               %Find the length of the signal in samples

P1_2side = abs(Y1/L1);                                  %Compute the two-sided amplitude spectrum

P1_1side = P1_2side(1:L1/2+1);                          %Compute the single-sided spectrum   
P1_1side(2:end-1) = 2*P1_1side(2:end-1);

f1 = int32((1/dt)*(0:(L1/2))/L1);                       %Vector of all the frequencies contained in the spectrum

subplot(1,2,2)
plot(f1,((P1_1side.^2)./(1/(L1*dt))))                   %Compute the power spectral density, in this case measured in Volts-squared per Hertz          

hold on;
xlabel('Frequency (Hz)')
ylabel('\muV^2/Hz')
title('Frequency Domain')
xlim([0 500])

%Plot 5 seconds of the signal
T_period=5;                                             %Time period of the signal, 5 seconds
EMG_sig_5s=EMG_sig_2k(find(t>=2,1,'first'):find(t<2+T_period,1,'last'),:);
figure;
subplot(1,3,1)
plot([0:length(EMG_sig_5s)-1].*dt,EMG_sig_5s)

N=length(EMG_sig_5s);                                   %Number of samples in the signal, 10000

% %Uncomment to plot the segments
% %Show 0.5s segments overlapping by 50%
% overlap = 0.5;L=round(0.5/dt);
% hold on
% clear x_data;clear y_data;
% L=round(0.5/dt);
% L_ovlp=floor((N/L-1)/(1-overlap))+1;
% x_data(1:2,[1:L_ovlp]')=[(round(L*(1-overlap))*dt)*([1:L_ovlp]'-1) (round(L*(1-overlap))*dt*([1:L_ovlp]'-1))+L*dt]';
% y_data=ones(size(x_data,1),size(x_data,2)).*1.1;
% y_data(:,1:2:length(x_data))=ones.*1.15;
% plot(x_data,y_data)
% 
% %Show 0.2s segments overlapping by 25%
% overlap = 0.25;L=round(0.2/dt);
% hold on
% clear x_data;clear y_data;
% L_ovlp=floor((N/L-1)/(1-overlap))+1;
% x_data(1:2,[1:L_ovlp]')=[(round(L*(1-overlap))*dt)*([1:L_ovlp]'-1) (round(L*(1-overlap))*dt*([1:L_ovlp]'-1))+L*dt]';
% y_data=ones(size(x_data,1),size(x_data,2)).*-1.1;
% y_data(:,1:2:length(x_data))=ones.*-1.15;
% plot(x_data,y_data)
 
xlim([0 5])
ylim([-1.5 1.5])
xlabel('Time (s)')
ylabel('mV')

%L is the length of the subsections of the total signal that will be
%analysed and averaged
L=round(0.5/dt);

% The number of points in the discrete Fourier Transform is chosen as the
% next power of 2 greater than the number of samples L
nfft = 2^nextpow2(L); %i.e. 2^10 = 1024

%Note: this means that 24 zeros are appended to each signal epoch in the time
%domain to create a frequency bin spacing of fs/nfft = 0.98 Hz - approx 1Hz

%Overlap of successive signal epochs or segments
overlap = 0.5; %This is equal to 50% overlap

%The number of samples in the overlap between successive epochs 
noverlap=round(overlap*nfft);

%The number of (overlapping) segments in the signal 
D=L-noverlap;

%Each successive segment starts D samples after the previous segment
N_avg=floor(((N-L)/D)+1);

%Calculate the power spectral density using a hamming window of length L
%(hamming window is used by default in pwelch), alternative windows are
%"hann" and "nuttallwin"
[Pxx,f]=pwelch(EMG_sig_5s-mean(EMG_sig_5s),hamming(L),overlap,nfft,fs);
subplot(1,3,2)
plot(f,Pxx)
xlim([0 500])
ylim([0 3.75*1e-4])
title(['N_avg = ',num2str(N_avg)],'interpreter','none')
xlabel('Frequency (Hz)')
ylabel('\muV^2/Hz')

%Decrease L to create a smoother spectrum, lower frequency spacing
L2=round(0.2/dt);
nfft2 = 2^nextpow2(L2); %i.e. 2^9 = 512
overlap2 = 0.25;                                          %This is equal to 25% overlap
noverlap2=round(overlap2*L2);                             %Overlap in number of samples

%The number of (overlapping) segments in the signal 
D2=L2-noverlap2;

%The number of (overlapping) segments in the signal 
N_avg=floor(((N-L2)/D2)+1);

[Pxx2,f2]=pwelch(EMG_sig_5s-mean(EMG_sig_5s),hamming(L2),overlap2,nfft2,fs);
subplot(1,3,3)
plot(f2,Pxx2)
xlim([0 500])
ylim([0 3.75*1e-4])
title(['N_avg = ',num2str(N_avg)],'interpreter','none')
xlabel('Frequency (Hz)')
ylabel('\muV^2/Hz')


%% Example (viii) Effect of multiplying signal by different window functions, see Figure 9 (c)
clear all;
load('Tutorial_EMG_sig.mat')
fs=2000;%Sampling frequency
dt=1/fs;
t=[0:length(EMG_sig_2k)-1].*dt;
T_period=0.5;                                           %Time period of the signal, 0.5 seconds
EMG_sig_short=EMG_sig_2k(find(t>=2,1,'first'):find(t<2+T_period,1,'last'),:);

figure;
subplot(3,3,1)
plot([0:length(EMG_sig_short)-1].*dt,EMG_sig_short)
ylim([-1.15 1.15])
xlim([0 0.5])
xlabel('Time (s)')
ylabel('mV')
subplot(3,3,2)
plot(hamming(length(EMG_sig_short)))                    %Plot a Hamming window
ylim([0 1.15])
xlabel('Samples')
ylabel('Amplitude')
subplot(3,3,3)
%Multiply signal segment by a Hamming window
plot([0:length(EMG_sig_short)-1].*dt,hamming(length(EMG_sig_short)).*EMG_sig_short)
xlim([0 0.501])
ylim([-1.15 1.15])
xlabel('Time (s)')
ylabel('mV')
subplot(3,3,4)
plot([0:length(EMG_sig_short)-1].*dt,EMG_sig_short)
ylim([-1.15 1.15])
xlim([0 0.5])
xlabel('Time (s)')
ylabel('mV')
subplot(3,3,5)
plot(hann(length(EMG_sig_short)))                      %Plot a Hanning window
ylim([0 1.15])
xlabel('Samples')
ylabel('Amplitude')
subplot(3,3,6)
%Multiply signal segment by a Hanning window
plot([0:length(EMG_sig_short)-1].*dt,hann(length(EMG_sig_short)).*EMG_sig_short)
xlim([0 0.501])
ylim([-1.15 1.15])
xlabel('Time (s)')
ylabel('mV')
subplot(3,3,7)
plot([0:length(EMG_sig_short)-1].*dt,EMG_sig_short)
ylim([-1.15 1.15])
xlim([0 0.5])
xlabel('Time (s)')
ylabel('mV')
subplot(3,3,8)
plot(nuttallwin(length(EMG_sig_short)))                %Plot a Nuttall window
ylim([0 1.15])
xlabel('Samples')
ylabel('Amplitude')
subplot(3,3,9)
%Multiply signal segment by a Nuttall window
plot([0:length(EMG_sig_short)-1].*dt,nuttallwin(length(EMG_sig_short)).*EMG_sig_short)
xlim([0 0.501])
ylim([-1.15 1.15])
xlabel('Time (s)')
ylabel('mV')

%% Section 5.2: Filtering, Figure 10

%% Example (ix) Illustrate the attenuation/roll-off characteristics of different filters, see Figure 10 (a)
clear all;
% Filter examples: Butterworth and Chebyshev Type I
% Two filter examples with a different order and roll-off

fc = 10;                                                %Cutoff Frequency
FilterOrder1 = 4;                                       %Higher order filters will have steeper roll-off            
FilterOrder2 = 6;         
ftype = 'low';                                          %Demonstrate filter frequency response plot for a low-pass filter

fs=2000;                                                %Sampling frequency
dt=1/fs;

% Zero-Pole-Gain design, Butterworth filter
[b_butter,a_butter] =  butter(FilterOrder1,fc./(fs/2),ftype);

% Chebyshev Type I filter
passbandripple=0.5;
[b_cheby1,a_cheby1] = cheby1(FilterOrder2,passbandripple,fc./(fs/2),ftype);

FreqPts=1000;
[h1,f1] = freqz(b_butter,a_butter,FreqPts,fs);                       %Frequency response of filter
[h2,f2] = freqz(b_cheby1,a_cheby1,FreqPts,fs);

%The following plot shows how much the filter attenuates each frequency
%component of the signal. For a low-pass filter there will be more
%attenuation at the higher frequencies (more negative) but the lower
%frequencies will pass through, i.e. 0 dB (no attenuation)
figure;
subplot(1,2,1)
f1_log=log10(f1);                                       %Plot Frequency on a log scale
dB=mag2db(abs(h1));                                     %Convert magnitude to decibels (dB)
plot(f1_log,dB)
ylim([-100 2])
xlim([0 2.5])
grid on;hold on;
set(gca,'XTickLabel',{'10^0','10^1','10^2','10^3'})
set(gca,'XTick',[0:1:3])
text(.55 ,0.85,[{'~80 dB Roll-off';'i.e. Slope'}],'Units','normalized')

%Plot a line between the first and second decades (a decade is a unit for measuring frequency ratios on a log scale)
line([1 1 2],[dB(find(f1_log>=1,1,'first')) dB(find(f1_log>=2,1,'first')) dB(find(f1_log>=2,1,'first'))],'color','k','linestyle','--')
xlabel('Frequency (Hz, log-scale)')
ylabel('Magnitude (dB)')
title([num2str(FilterOrder1),'th Order Butterworth Filter'])

subplot(1,2,2)
plot(log10(f2),mag2db(abs(h2)))
ylim([-100 2])
xlim([0 2.5])
grid on;hold on;
set(gca,'XTickLabel',{'10^0','10^1','10^2','10^3'})
set(gca,'XTick',[0:1:3])
text(.55 ,0.85,[{'~160 dB Roll-off';'Steeper Slope'}],'Units','normalized')
title([num2str(FilterOrder2),'th Order Chebyshev Type I Filter'])

%% Example (x) Examples of low-pass, high-pass, band-pass and notch filters, see Figure 10 (a)

clear all;
load('Tutorial_EMG_sig.mat')
fs=2000;                                                %Sampling frequency
dt=1/fs;

t=[0:length(EMG_sig_2k)-1].*dt;

% Low-pass filter example

%If the EMG signal was to be sampled at 500 Hz, it would first need to
%low-pass filtered at half that frequency (250 Hz) to remove
%the high-frequency components and prevent aliasing of the signal.
fc = 250;                                               %Cutoff Frequency
FilterOrder = 4;   
ftype = 'low';

% Zero-Pole-Gain design, Butterworth filter
[b,a] =  butter(FilterOrder,fc./(fs/2),ftype);

EMG_sig_2k_lowpass = filtfilt(b,a,EMG_sig_2k);        %Filter the signal
figure;
subplot(1,2,1)
plot(t,EMG_sig_2k)
hold on
plot(t,EMG_sig_2k_lowpass,'-.')
xlim([5 5.25])
ylim([-1.15 1.15])
legend({'Original Signal';'Low-Pass Filtered Signal'})
xlabel('Time (s)')
ylabel('mV')

% The low-pass filter shapes the power spectral density of the EMG signal

L=round(0.5/dt);                                         %Split the signal into subsections of 0.5 seconds
nfft = 2^nextpow2(L);                                    %Length of nfft window is the next power of 2 greater than L
overlap = 0.5;                                           %This is equal to 50% overlap

%Calculate the power spectral density of original signal using Welch's method
[Pxx,f]=pwelch(EMG_sig_2k-mean(EMG_sig_2k),hamming(L),overlap,nfft,fs);
subplot(1,2,2)
plot(f,Pxx)
hold on;
%Calculate the power spectral density of low-pass filtered signal using Welch's method
[Pxx_lowpass,f]=pwelch(EMG_sig_2k_lowpass-mean(EMG_sig_2k_lowpass),hamming(L),overlap,nfft,fs);
plot(f,Pxx_lowpass,'-.')
legend({'Original Signal';'Low-Pass Filtered Signal'})
xlim([0 500])
xlabel('Frequency (Hz)')
ylabel('\muV^2/Hz')

% High-pass filter example

%High-pass filter the EMG signal to remove noise components below 20 Hz

fc = 20;                                                 %Cutoff Frequency
FilterOrder = 4;   
ftype = 'high';

clear b;clear a;
% Zero-Pole-Gain design, Butterworth filter
[b,a] =  butter(FilterOrder,fc./(fs/2),ftype);

EMG_sig_2k_highpass = filtfilt(b,a,EMG_sig_2k);
figure;
subplot(1,2,1)
plot(t,EMG_sig_2k)
hold on
plot(t,EMG_sig_2k_highpass,'-.')
xlim([5 5.25])
ylim([-1.15 1.15])
legend({'Original Signal';'High-Pass Filtered Signal'})
xlabel('Time (s)')
ylabel('mV')

% The high-pass filter shapes the power spectral density of the EMG signal

L=round(0.5/dt);                                         %Split the signal into subsections of 0.5 seconds
nfft = 2^nextpow2(L);                                    %Length of nfft window is the next power of 2 greater than L
overlap = 0.5;                                           %This is equal to 50% overlap

%Calculate the power spectral density of original signal using Welch's method
[Pxx,f]=pwelch(EMG_sig_2k-mean(EMG_sig_2k),hamming(L),overlap,nfft,fs);
subplot(1,2,2)
plot(f,Pxx)
hold on;
%Calculate the power spectral density of high-pass filtered signal using Welch's method
[Pxx_highpass,f]=pwelch(EMG_sig_2k_highpass-mean(EMG_sig_2k_highpass),hamming(L),overlap,nfft,fs);
plot(f,Pxx_highpass,'-.')
legend({'Original Signal';'High-Pass Filtered Signal'})
xlim([0 500])
xlabel('Frequency (Hz)')
ylabel('\muV^2/Hz')

% Band-pass filter example

%Remove noise components below 20Hz and above 500 Hz in the EMG signal
fc = [20 500];                                            %Band-pass filter the EMG signal between 20 and 500 Hz
FilterOrder = 4;   
ftype = 'bandpass';

clear b;clear a;

% Zero-Pole-Gain design, Butterworth filter
[b,a] =  butter(FilterOrder,fc./(fs/2),ftype);

EMG_sig_2k_bandpass = filtfilt(b,a,EMG_sig_2k);
figure;
subplot(1,2,1)
plot(t,EMG_sig_2k)
hold on
plot(t,EMG_sig_2k_bandpass,'-.')
xlim([5 5.25])
ylim([-1.15 1.15])
legend({'Original Signal';'Band-Pass Filtered Signal'})
xlabel('Time (s)')
ylabel('mV')

% The band-pass filter shapes the power spectral density of the EMG signal

L=round(0.5/dt);                                         %Split the signal into subsections of 0.5 seconds
nfft = 2^nextpow2(L);                                    %Length of nfft window is the next power of 2 greater than L
overlap = 0.5;                                           %This is equal to 50% overlap

%Calculate the power spectral density of original signal using Welch's method
[Pxx,f]=pwelch(EMG_sig_2k-mean(EMG_sig_2k),hamming(L),overlap,nfft,fs);
subplot(1,2,2)
plot(f,Pxx)
hold on;
%Calculate the power spectral density of band-pass filtered signal using Welch's method
[Pxx_bandpass,f]=pwelch(EMG_sig_2k_bandpass-mean(EMG_sig_2k_bandpass),hamming(L),overlap,nfft,fs);
plot(f,Pxx_bandpass,'-.')
legend({'Original Signal';'Band-Pass Filtered Signal'})
xlim([0 500])
xlabel('Frequency (Hz)')
ylabel('\muV^2/Hz')

% Notch filter example

%First add noise (50 Hz sine wave) to the signal to simulate power line interference
%Note that in real signals power line interference will also be present at
%higher harmonics (100 Hz, 150 Hz etc.)
f_noise=50;
PowerLine=0.025.*sin((2*pi*f_noise*t));                     %In mV
EMG_sig_2k_noise=EMG_sig_2k+PowerLine';
fc = [48 52];                                               %Notch filter the EMG signal between 20 and 500 Hz
FilterOrder = 4;   
ftype = 'stop';

clear b;clear a;

% Zero-Pole-Gain design, Butterworth filter
[b,a] =  butter(FilterOrder,fc./(fs/2),ftype);

EMG_sig_2k_notch = filtfilt(b,a,EMG_sig_2k_noise);
figure;
subplot(1,2,1)
plot(t,EMG_sig_2k_noise)
hold on
plot(t,EMG_sig_2k_notch,'-.')
xlim([5 5.25])
ylim([-1.15 1.15])
legend({'Original Signal';'Band-Pass Filtered Signal'})
xlabel('Time (s)')
ylabel('mV')

% The notch filter shapes the power spectral density of the EMG signal

L=round(0.5/dt);                                         %Split the signal into subsections of 0.5 seconds
nfft = 2^nextpow2(L);                                    %Length of nfft window is the next power of 2 greater than L
overlap = 0.5;                                           %This is equal to 50% overlap

%Calculate the power spectral density of signal with Power Line Interference using Welch's method
[Pxx,f]=pwelch(EMG_sig_2k_noise-mean(EMG_sig_2k_noise),hamming(L),overlap,nfft,fs);
subplot(1,2,2)
plot(f,Pxx)
hold on;
%Calculate the power spectral density of the notch filtered signal using Welch's method
[Pxx_notch,f]=pwelch(EMG_sig_2k_notch-mean(EMG_sig_2k_notch),hamming(L),overlap,nfft,fs);
plot(f,Pxx_notch,'-.')
legend({'Original Signal';'Notch Filtered Signal'})
xlim([0 500])
xlabel('Frequency (Hz)')
ylabel('\muV^2/Hz')

%% Section 4.2: Choice of Electrode: 

%% Example (xi) Demonstrate effect of Inter-Electrode Distance (IED), see Figure 6

%Two data-sets with different IEDs
clear all;
fs=2000;                                                    %Sampling frequency

dt=1/fs;
%IED = 1cm
load('Tutorial_Gait_EMG.mat')
t1=[0:length(EMG_Gait_sig_2k)-1].*dt;

%IED = 0.5cm
load('Tutorial_AmpChange_EMG.mat')
t2=[0:length(EMG_AmpChange_sig_2k)-1].*dt;

IED_1=EMG_Gait_sig_2k(t1>=7.45&t1<7.7);                     %Crop signal lengths
IED_0_5=EMG_AmpChange_sig_2k(t2>=48.5&t2<48.75);
L1=length(IED_1);
L2=length(IED_0_5);

Y_1 = fft(IED_1-mean(IED_1));                               %computes the discrete Fourier transform using a fast Fourier transform (FFT) algorithm.
Y_0_5 = fft(IED_0_5-mean(IED_0_5));

%Generate the two-sided amplitude frequency spectra
P1_2side = abs(Y_1/L1);
P2_2side = abs(Y_0_5/L2);

%Generate the one-sided amplitude frequency spectra
P1_1side = P1_2side(1:L1/2+1);
P1_1side(2:end-1) = 2*P1_1side(2:end-1);

P2_1side = P2_2side(1:L2/2+1);
P2_1side(2:end-1) = 2*P2_1side(2:end-1);

%Generate vectors of frequencies contained within the spectra
f1 = fs*(0:(L1/2))/L1;
f2 = fs*(0:(L2/2))/L2;

figure;
subplot(2,2,1)
plot([0:length(IED_1)-1].*dt,IED_1)
xlabel('Time (s)')
ylabel('Amplitude (mV)')
title('IED = 1 cm')
ylim([-1 1])
xlim([0 0.25])
subplot(2,2,2)
plot(f1,(P1_1side.^2)./(1/(L1*dt)))
xlabel('Frequency (Hz)')
ylabel('\muV^2/Hz')
title('Power Spectrum - Signal 1')
xlim([0 500])
subplot(2,2,3)
plot([0:length(IED_0_5)-1].*dt,IED_0_5,'r')
xlabel('Time (s)')
ylabel('Amplitude (mV)')
title('IED = 0.5 cm')
ylim([-0.75 0.75])
xlim([0 0.25])
subplot(2,2,4)
plot(f2,(P2_1side.^2)./(1/(L2*dt)),'r') 
xlabel('Frequency (Hz)')
ylabel('\muV^2/Hz')
title('Power Spectrum - Signal 2')
xlim([0 500])

% Note that this comparison for illustrative purposes only and that the power spectra of the EMG signals
% cannot be directly compared, as they were recorded in different muscles, under different conditions, 
% with different electrodes

%% Section 4.3: Choice of Amplifier

%% Example (xii) Demonstrate signal amplification with differential filter, see Figure 7 (c)
clear all;
%Generate simple sine waves to demonstrate signal amplification using a
%differential amplifier

T_period=1;                                                 %Time period of the signal, 2 seconds
dt=0.001;                                                   %Sampling period of the signal, samples taken every 0.001 seconds
t=[0:dt:T_period-dt];                                       %Time vector of samples

%Create two sine wave signals, to represent what is detected at each input
%to the amplifier, i.e. at V+ and V-

V_neg=0.75.*sin((2*pi*10*t));
V_pos=sin((2*pi*15*t)+pi/2);

%The difference between these two inputs, Vd, is what is amplified by the
%amplifier

Vd=V_neg-V_pos;

figure;
subplot(2,3,1)
plot(t,V_neg)
xlim([0 0.5])
ylim([-2.5 2.5])
title('V_{-}')
xlabel('Time (s)')
ylabel('Amplitude')
subplot(2,3,4)
plot(t,V_pos)
xlim([0 0.5])
ylim([-2.5 2.5])
title('V_{+}')
xlabel('Time (s)')
ylabel('Amplitude')
subplot(3,3,5);
plot(t,Vd)
xlim([0 0.5])
ylim([-2.25 2.25])
title('V_{d}')
xlabel('Time (s)')
ylabel('Amplitude')
%Set the amplifier gain, Ad, and the voltage output by the amplifier is
%equal to Vd multiplied by the gain

Ad=2;
Vout=Vd*Ad;

subplot(1,3,3);
plot(t,Vout)
xlim([0 0.5])
ylim([-4 4])
title('V_{out}')
xlabel('Time (s)')
ylabel('Amplitude')

 %% Example (xiii) Calculate CMRR
 
% Calculate the common-mode rejection ratio (CMRR) of an amplifier given
% the differential gain, Ad, and the common-mode gain, Acm

Ad=40;
Acm=2;
CMRR=20*log(abs(Ad/Acm));


%% Section 5: EMG Signal Analysis Techniques

%% Section 5.3: Surface EMG amplitude features

 %% Example (xiv) Plot gait data, see Figure 11

% Plot EMG gait data from the soleus muscle during walking, with heel
% strike and toe-off information.

clear all;
load('Tutorial_Gait_EMG.mat')
fs=2000;                                                    %Sampling frequency
dt=1/fs;

t=[0:length(EMG_Gait_sig_2k)-1].*dt;
figure;
% see Figure 11 (a)
subplot(4,1,1)
h1=plot(t,EMG_Gait_sig_2k);                                 %Plot raw EMG data during walking
hold on;
h2=plot(HS'*[1 1],[-1 1],'k--');                            %Plot heel strike information
h3=plot(TO'*[1 1],[-1 1],'r--');                            %Plot toe-off information
xlabel('Time (s)')
ylabel('Amplitude (mV)')
title('Raw sEMG Signal')
legend([h1,h2(1),h3(1)],'EMG Signal','Heel Strike','Toe off')

%Plot the absolute value of the EMG gait data
% see Figure 11 (b)
subplot(4,1,2)
plot(t,abs(EMG_Gait_sig_2k-mean(EMG_Gait_sig_2k)))          %Plot rectified EMG data during walking
xlabel('Time (s)')
ylabel('Amplitude (mV)')
title('Rectified sEMG Signal')
ylim([-0.1 1])

%Apply the Teager Kaiser Energy Operator to the EMG signal, a signal
%transformation that helps to identify when the muscle is "active"

%First remove frequency components below 30Hz and above 300 Hz in the EMG signal
fc = [30 300];% Band-pass filter the EMG signal between 30 and 300 Hz
FilterOrder = 4;   
ftype = 'bandpass';


% Zero-Pole-Gain design, Butterworth filter
[b,a] =  butter(FilterOrder,fc./(fs/2),ftype);

%Filter signal to remove frequency components below 30Hz and above 300 Hz in the EMG signal
bp_sig = filtfilt(b,a,EMG_Gait_sig_2k-mean(EMG_Gait_sig_2k));

% Applying the TKEO operator to the filtered signal
tk=[];
for idx = 2:length(bp_sig)-1
    tk(idx-1) = bp_sig(idx)^2-bp_sig(idx-1)*bp_sig(idx+1);
end

%Use a 50Hz low pass filter to extract the outline or "shape" of the EMG signal

fc = 50;                                                    %Cutoff Frequency
FilterOrder = 4;   
ftype = 'low';
% Zero-Pole-Gain design, Butterworth filter
clear b;clear a;
[b,a] =  butter(FilterOrder,fc./(fs/2),ftype);

tk_env = filtfilt(b,a,abs(tk));                           %low-pass filtering to get the shape/outline of the teager kaiser signal

tk_env_norm = tk_env/max(abs(tk_env));                      %Normalising the teager kaiser envelope

t2 = 1/fs:1/fs:length(tk_env_norm)/fs;                      %time vector to be used in the final plot

%Plot the outline or "shape" of the EMG signal
% see Figure 11 (c)

subplot(4,1,3)
plot(t2,tk_env_norm)                                        %Plot teager kaiser envelope
xlabel('Time (s)')
ylabel('Amplitude (mV)')
title('sEMG Signal Envelope')
ylim([-0.1 1.1])

%Calculate when the muscle is active, by setting a threshold for the
%activation. Whenever the outline of the EMG signal is above this
%threshold the muscle is active

n = 1.5; %number of standard deviations greater than the mean for the threshold, this can be adjust according to data

thresh = mean(tk_env_norm(round(fs*0.5):round(fs*7.5)))+n*std(tk_env_norm(round(fs*0.5):round(fs*7.5))); % threshold value

activation = zeros(1,length(tk_env));                        %preallocating a variable to determine the activation times

% if the EMG signal exceeds the threshold for more than 10 samples, the muscle is
% considered active

for idx = 10:length(tk_env_norm)-(0.5*fs)-1
    if all(tk_env_norm(idx-9:idx) > thresh)                  % Checking for when 10 samples in a row are greater than the threshold
        activation(idx-9:idx+(0.1*fs)) = 1;                  % If there are 10 samples in a row over the threshold, muscle is active for next 0.1 seconds
    end
end

hold on;

%Plot when the muscle is "OFF" (activation = 0) and "ON" (activation = 1) 
plot(t2,activation)


%Find the moving average amplitude of the EMG signal on the rectified or absolute
%value of the EMG signal, using a window 0.2 seconds in length with an overlap of 25%

absEMG=abs(EMG_Gait_sig_2k-mean(EMG_Gait_sig_2k));           %absolute value of the EMG gait data
WinLen=0.2;                                                  %in seconds, length of the moving average window
ovlp=0.25;                                                   % a 25% overlap
TimeStep=WinLen*ovlp;                                        %in seconds, denotes the degree of overlap between each successive window,
            
StartWin=[1:(TimeStep*fs):size(absEMG,1)-(WinLen*fs)];       %Index for the start of each window
rmsfunc=[];
for index=StartWin
    rmsfunc((index-1)/(TimeStep*fs)+1,:)=mean(absEMG(index:(index+WinLen*fs),:));
end

%Low-pass filtering the EMG signal is an alternative way of calculating a "moving average"
% see Figure 11 (d)
subplot(4,1,4)
h1=plot((StartWin/fs)+(WinLen/2),rmsfunc);
hold on;

%Lower cut-off frequencies will result in a smoother outline or shape.
fc = 5;                                                      %Cutoff Frequency
ftype = 'low';
clear b;clear a;
[b,a] =  butter(FilterOrder,fc./(fs/2),ftype);

%To calculate the moving average using a low-pass filter, the function 
%filtfilt is used, as this applies the filter twice to the EMG data in both 
%the forward and reverse directions (i.e. zero-phase digital filtering).
 
EMG_av = filtfilt(b,a,absEMG);                             %low-pass filtering the rectified EMG to obtain a moving average
h2=plot(t,EMG_av);

%If the filter is not applied twice,i.e. using sosfilt, there will be a delay in the moving
%average window

EMG_av_delay = filter(b,a,absEMG);  

h3=plot(t,EMG_av_delay,'linestyle','--');
legend([h1,h2,h3],'Moving Average','Low-Pass Filter','Low-Pass Filter (with delay)');
xlabel('Time (s)')
ylabel('Amplitude (mV)')
title('Moving Average Value sEMG Signal')
ylim([-0.05 0.25])

%% Example (xv) Calculate the Root-Mean Square value and the Average Rectified value of an EMG signal

%Sample surface EMG signal shows an increase in amplitude

clear all;
load('Tutorial_AmpChange_EMG.mat')
fs=2000;
dt=1/fs;
t=[0:length(EMG_AmpChange_sig_2k)-1].*dt;

%Plot Raw sEMG signal
figure;
subplot(3,1,1)
plot(t,EMG_AmpChange_sig_2k) ;
ylim([-1.5 1.5])
xlim([0 65])
xlabel('Time (s)')
ylabel('Amplitude (mV)')
title('Raw sEMG Signal')
%Plot Rectified sEMG signal
subplot(3,1,2)
plot(t,abs(EMG_AmpChange_sig_2k))
ylim([-0.1 1.5])
xlim([0 65])
xlabel('Time (s)')
ylabel('Amplitude (mV)')
title('Rectified sEMG Signal')

%Calculate moving average of the rectified surface EMG signal, with a window length
%of 2 seconds and an overlap of 50% between successive windows
absEMG=abs(EMG_AmpChange_sig_2k-mean(EMG_AmpChange_sig_2k)); %absolute value of the EMG force data (rectified value)
WinLen=2;                                                    %in seconds, length of the moving average window
ovlp=0.5;                                                    % a 50% overlap
TimeStep=WinLen*ovlp;                                        %in seconds, denotes the degree of overlap between each successive window
            
StartWin=[1:(TimeStep*fs):size(absEMG,1)-(WinLen*fs)];       %Index for the start of each window

%Calculate moving average 
rmsfunc=[];
for index=StartWin
    rmsfunc((index-1)/(TimeStep*fs)+1,:)=mean(absEMG(index:(index+WinLen*fs),:));
end

subplot(3,1,3)
plot(StartWin/fs,rmsfunc)
ylim([-0.025 0.2])
xlim([0 65])
xlabel('Time (s)')
ylabel('Amplitude (mV)')
title('Moving Average Value sEMG Signal')

%% Example (xvi) Calculate the average of the surface EMG signal for 3 separate segments
clear all
load('Tutorial_AmpChange_EMG.mat')
fs=2000;
dt=1/fs;
t=[0:length(EMG_AmpChange_sig_2k)-1].*dt;

Segment1=[10/dt:((15/dt)-dt)];                               %Segment 1 is from 10 seconds to 15 seconds
Segment2=[30/dt:((35/dt)-dt)];                               %Segment 2 is from 30 seconds to 35 seconds
Segment3=[45/dt:((50/dt)-dt)];                               %Segment 3 is from 45 seconds to 50 seconds

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

%% Section 5.4: Surface EMG spectral features

%% Example (xvii) Calculate the Mean Frequency and the Median Frequency of the surface EMG signal, see Figure 12

clear all;
% EMG signals from the start and end of a fatiguing isometric contraction
load('Tutorial_Fatigue_EMG.mat')
fs=2000;
dt=1/fs;

%Subtract mean from EMG signal (DC voltage offset)
EMG_Start=EMG_Start-mean(EMG_Start);
EMG_End=EMG_End-mean(EMG_End);

figure;
subplot(2,2,1)
plot([0:length(EMG_Start)-1].*dt,EMG_Start)

% Calculate the RMS and ARV amplitude of the sEMG signal at the start of the contraction
RMS_Start=sqrt((1/length(EMG_Start)).*sum(EMG_Start.^2));
ARV_Start=(1/length(EMG_Start)).*sum(abs(EMG_Start));

text(.55 ,0.85,{['RMS = ',num2str(RMS_Start,2)];['ARV = ',num2str(ARV_Start,2)]},'Units','normalized')

xlabel('Time (s)')
ylabel('Amplitude (mV)')
xlim([0 5])
ylim([-2 2])

% Change limits to zoom in
% ylim([-1 1])
% xlim([2.5 2.75])

subplot(2,2,3)
plot([0:length(EMG_End)-1].*dt,EMG_End)


% Calculate the RMS and ARV amplitude of the sEMG signal at the end of the contraction
RMS_End=sqrt((1/length(EMG_End)).*sum(EMG_End.^2));
ARV_End=(1/length(EMG_End)).*sum(abs(EMG_End));

text(.55 ,0.85,[{['RMS = ',num2str(RMS_End,2)];['ARV = ',num2str(ARV_End,2)]}],'Units','normalized')

xlabel('Time (s)')
ylabel('Amplitude (mV)')
xlim([0 5])
ylim([-2 2])

% Change limits to zoom in
% ylim([-1 1])
% xlim([2.5 2.75])

L=round(0.25/dt);                                            %Split the signal into subsections of 0.25 seconds
nfft = 2^nextpow2(L);                                        %Length of nfft window is the next power of 2 greater than L
overlap = 0.5;                                               %This is equal to 50% overlap

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

subplot(2,2,[2,4])
h1=plot(f1,P_Start,'r');
hold on;
%Plot a line to indicate the Mean Frequency at the start of the contraction on the power spectral density plot
h2=plot([Mean_freq_Start Mean_freq_Start],[0 max(P_Start)],'r','linewidth',2.5,'linestyle','--');

%Calculate the power spectral density at the end of the signal using Welch's method
[P_End,f2]=pwelch(EMG_End,hamming(L),overlap,nfft,fs);
h3=plot(f2,P_End,'b');

%Plot a line to indicate the Mean Frequency at the end of the contraction on the power spectral density plot
h4=plot([Mean_freq_End Mean_freq_End],[0 max(P_End)],'b','linewidth',2.5,'linestyle','--');

ylim([0 2*1e-3])
xlim([0 300])
xlabel('Frequency (Hz)')
ylabel('\muV^2/Hz')
legend([h1,h3,h2,h4],'Start of Contraction','End of Contraction','Mean Freq. (Start)','Mean Freq. (End)')

%Calculate the mean frequency (MNF) value for each segment and print
fprintf(1, 'Mean Freq. (Start): %.2f\nMean Freq. (End): %.2f\n', [Mean_freq_Start;Mean_freq_End]);

% Calculate the median frequency

% Find the area under the power spectral density
PSDarea_Start = cumtrapz(f1, P_Start);
PSDarea_End = cumtrapz(f2, P_End); %% Jeremy Edit:f1 and f2 are the same but it's confusing to use the same for both
%Find the frequency that divides the EMG power spectrum into two equal
%regions, each containing half of the total power within the signal epoch
MedFreq_Start = interp1(PSDarea_Start, f1, PSDarea_Start(end)/2);
MedFreq_End = interp1(PSDarea_End, f2, PSDarea_End(end)/2);ithin the signal epoch
MedFreq = interp1(PSDarea, f1, PSDarea(end)/2);

%Calculate the median frequency (MDF) value for each segment and print
fprintf(1, 'Median Freq. (Start): %.2f\nMedian Freq. (End): %.2f\n', [Median_freq_Start;Median_freq_End]);

%% Example (xviii) Calculating the coherence between two EMG signals, Figure A2
clear all;
% EMG signals recorded from the biceps brachii and the brachioradialis during a
% sustained isometric contraction at 50% MVC
load('Tutorial_Coherence_EMG.mat')
fs=1000;
dt=1/fs;

figure;
% Plot the time domain signals
subplot(2,3,1)
plot([0:length(EMG_BB)-1].*dt,EMG_BB)
ylim([-1 1])
xlim([2.5 42.5])
xlabel('Time (s)')
ylabel('Amplitude (mV)')

subplot(2,3,4)
plot([0:length(EMG_BR)-1].*dt,EMG_BR)
ylim([-1 1])
xlim([2.5 42.5])
xlabel('Time (s)')
ylabel('Amplitude (mV)')


%Crop the signals to analyze the steady-state EMG in the middle of the
%contraction where there is minimal recruitment and de-recruitment of motor
%units.
EMG_BB_crop=EMG_BB(int32(17.5/dt):int32(((40/dt)-dt)));
EMG_BR_crop=EMG_BR(int32(17.5/dt):int32(((40/dt)-dt)));

% Plot the power spectral densities of each EMG signal

L=round(0.5/dt);                                             %Split the signal into subsections of 0.5 seconds
nfft = 2^nextpow2(L);                                        %Length of nfft window is the next power of 2 greater than L
overlap = 0.75;                                              %This is equal to 75% overlap
noverlap=round(overlap*nfft);                                %Overlap in samples

%Calculate the power spectral density for signals recorded from BB and BR using Welch's method
[P_BB,f1]=pwelch(EMG_BB_crop,hamming(L),overlap,nfft,fs);
[P_BR,f2]=pwelch(EMG_BR_crop,hamming(L),overlap,nfft,fs);

subplot(2,3,2)
plot(f1,P_BB)
ylim([0 4*1e-4])
xlim([0 300])
xlabel('Frequency (Hz)')
ylabel('\muV^2/Hz')

subplot(2,3,5)
plot(f1,P_BR)
ylim([0 3*1e-4])
xlim([0 300])
xlabel('Frequency (Hz)')
ylabel('\muV^2/Hz')

% Calculate the magnitude-squared coherence
L=round(1/dt);                                               %Split the signal into subsections of 1 second
nfft = 2^nextpow2(L);                                        %Length of nfft window is the next power of 2 greater than L
overlap = 0.75;                                              %This is equal to 75% overlap
noverlap=round(overlap*nfft);                                %Overlap in samples

%Calculate the magnitude-squared coherence between signals recorded from BB and BR
[c, f] = mscohere(EMG_BB_crop,EMG_BR_crop,hamming(nfft),overlap,nfft,fs);

subplot(2,3,[3 6])
plot(f,c)
xlim([0 80])
ylim([0 0.5])
xlabel('Frequency (Hz)')
ylabel('Magnitude-Squared Coherence')
