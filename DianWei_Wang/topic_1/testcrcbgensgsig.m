close all; clear; clc;
%% Plot the sin-Gaussin signal
% Signal parameters
t0 = 0.5;
sigma = 0.2;
phi0 = 0;
f0 = 25;
f1 = 3;
A = 10;
% Instantaneous frequency at t0 sec is 
maxFreq = f0;
samplFreq = 5*maxFreq;
samplIntrvl = 1/samplFreq;

% Time samples
timeVec = 0:samplIntrvl:1.0;
% Number of samples
nSamples = length(timeVec);

% Generate the signal
sigVec = crcbgensgsig(timeVec,A,[t0, sigma, phi0, f0]);

%Plot the signal 
figure;
plot(timeVec,sigVec,'Marker','.','MarkerSize',24);

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
winLen = 0.1;%sec
ovrlp = 0.05;%sec
%Convert to integer number of samples 
winLenSmpls = floor(winLen*samplFreq);
ovrlpSmpls = floor(ovrlp*samplFreq);
[S,F,T]=spectrogram(sigVec,winLenSmpls,ovrlpSmpls,[],samplFreq);
figure;
imagesc(T,F,abs(S)); axis xy;
xlabel('Time (sec)');
ylabel('Frequency (Hz)');

%% %% Retain frequencies below a threshold frequency
% Design low pass filter
filtOrdr = 30;
b1 = fir1(filtOrdr,(maxFreq - 0/(sigma))/(samplFreq/2));
b2 = fir1(filtOrdr,(maxFreq - 0.5/(sigma))/(samplFreq/2));
b3 = fir1(filtOrdr,(maxFreq - 2/(sigma))/(samplFreq/2));
% Apply filter
filtSig1 = fftfilt(b1,sigVec);
filtSig2 = fftfilt(b2,sigVec);
filtSig3 = fftfilt(b3,sigVec);

% Plots
figure;
hold on;
plot(timeVec,sigVec);
plot(timeVec,filtSig1);
plot(timeVec,filtSig2);
plot(timeVec,filtSig3);
%% Generating sin-Gaussian signal
function sigVec = crcbgensgsig(dataX,snr,sgCoefs)
% Generate a sin-Gaussian signal
% S = crcbgensgsig(X,SNR,C)
% Generates a sin-Gaussian signal S. X is the vector of
% time stamps at which the samples of the signal are to be computed. SNR is
% the matched filtering signal-to-noise ratio of S and C is the vector of
% four coefficients [t0, sigma, phi0, f0] that parametrize the amplitude and phase of the signal:
% exp(-((t - t0)^2) / (2*sigma^2)) * sin(2*pi*(f0*t + phi0)). 

%Soumya D. Mohanty, May 2018
%DianWei Wang modified, Feb 2022

ampVec = exp(- (dataX-sgCoefs(1)).^2 / (2 * sgCoefs(2)^2));
phaseVec = sgCoefs(3) + sgCoefs(4)*dataX;
sigVec = ampVec.*sin(2*pi*phaseVec);
sigVec = snr*sigVec/norm(sigVec);
end