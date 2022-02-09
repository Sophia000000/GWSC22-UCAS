close all; clear; clc;
%% Plot the AM signal
% Signal parameters
f0 = 100;
phi0 = 0;
f1 = 2;
A = 10;
% Base frequency after 1 sec is 
maxFreq = f0+2*f1;
samplFreq = 5*maxFreq;
samplIntrvl = 1/samplFreq;

% Time samples
timeVec = 0:samplIntrvl:1.0;
% Number of samples
nSamples = length(timeVec);

% Generate the signal
sigVec = crcbgenamsig(timeVec,A,[f1, f0, phi0]);

%Plot the signal 
figure;
plot(timeVec,sigVec,'Marker','.','MarkerSize',12);

%Plot the periodogram
%--------------
%Length of data 
dataLen = timeVec(end)-timeVec(1);
%DFT sample corresponding to Nyquist frequency
kNyq = floor(nSamples/2)+1;
% Positive Fourier frequencies
posFreq = (0:(kNyq-1))*(1/dataLen);
% FFT of signal
fftSig = fft(sigVec);
% Discard negative frequencies
fftSig = fftSig(1:kNyq);

%Plot periodogram
figure;
plot(posFreq,abs(fftSig));

%Plot a spectrogram
%----------------
winLen = 0.05;%sec
ovrlp = 0.025;%sec
%Convert to integer number of samples 
winLenSmpls = floor(winLen*samplFreq);
ovrlpSmpls = floor(ovrlp*samplFreq);
[S,F,T]=spectrogram(sigVec,winLenSmpls,ovrlpSmpls,[],samplFreq);
figure;
imagesc(T,F,abs(S)); axis xy;
xlabel('Time (sec)');
ylabel('Frequency (Hz)');

%% %% Remove lower side band
% Design high pass filter
filtOrdr = 30;
b = fir1(filtOrdr,(maxFreq)/(samplFreq/2),'high');
% Apply filter
filtSig = fftfilt(b,sigVec);

% Plots
figure;
hold on;
plot(timeVec,sigVec);
plot(timeVec,filtSig);
%% Generating AM signal
function sigVec = crcbgenamsig(dataX,snr,amCoefs)
% Generate a AM signal
% S = crcbgenamsig(X,SNR,C)
% Generates a AM signal S. X is the vector of
% time stamps at which the samples of the signal are to be computed. SNR is
% the matched filtering signal-to-noise ratio of S and C is the vector of
% three coefficients [f1, f0, phi0] that parametrize the amplitude and phase of the signal:
% cos(2*pi*f1*t)*sin(2*pi*(f0*t+phi0)). 

%Soumya D. Mohanty, May 2018
%DianWei Wang modified, Feb 2022

ampVec = cos(2*pi*amCoefs(1)*dataX);
phaseVec = amCoefs(2)*dataX + amCoefs(3);
sigVec = ampVec.*sin(2*pi*phaseVec);
sigVec = snr*sigVec/norm(sigVec);
end