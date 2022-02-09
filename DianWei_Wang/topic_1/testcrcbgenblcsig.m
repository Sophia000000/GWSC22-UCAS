close all; clear; clc;
%% Plot the box-linear chirp signal
% Signal parameters
t_a = 2;
L = 3;
phi0=0;
f0=10;
f1=5;
A = 10;
% Instantaneous frequency after (t_a + L) sec is 
maxFreq = f0+2*f1*L;
samplFreq = 5*maxFreq;
samplIntrvl = 1/samplFreq;

% Time samples
timeVec = 0:samplIntrvl:(t_a + L + 0);
% Number of samples
nSamples = length(timeVec);

% Generate the signal
sigVec = crcbgenlcsig(timeVec,A,[t_a, L, phi0, f0, f1]);

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
winLen = 0.2;%sec
ovrlp = 0.1;%sec
%Convert to integer number of samples 
winLenSmpls = floor(winLen*samplFreq);
ovrlpSmpls = floor(ovrlp*samplFreq);
[S,F,T]=spectrogram(sigVec,winLenSmpls,ovrlpSmpls,[],samplFreq);
figure;
imagesc(T,F,abs(S)); axis xy;
xlabel('Time (sec)');
ylabel('Frequency (Hz)');

%% %% Remove frequencies above half the central frequency
%{
% Design low pass filter
filtOrdr = 30;
b = fir1(filtOrdr,(f0+f1 * L)/(samplFreq/2));
% Apply filter
filtSig = fftfilt(b,sigVec);

% Plots
figure;
hold on;
plot(timeVec,sigVec);
plot(timeVec,filtSig);
%}
%% Generating box-linear chirp signal
function sigVec = crcbgenlcsig(dataX,snr,blcCoefs)
% Generate a box-linear chirp signal
% S = crcbgenlcsig(X,SNR,C)
% Generates a linear chirp signal S. X is the vector of
% time stamps at which the samples of the signal are to be computed. SNR is
% the matched filtering signal-to-noise ratio of S and C is the vector of
% three coefficients [t_a, L, phi0, f0, f1] that parametrize the phase of the signal:
% phi0+f0*t+f1*t^2. 

%Soumya D. Mohanty, May 2018
%DianWei Wang modified, Feb 2022

ampVec = rectpulse(blcCoefs(1), blcCoefs(1) + blcCoefs(2), dataX);
phaseVec = blcCoefs(3) + blcCoefs(4)*dataX + blcCoefs(5)*dataX.^2;
sigVec = ampVec.*sin(2*pi*phaseVec);
sigVec = snr*sigVec/norm(sigVec);
end
function y = rectpulse(a, b, x)
y = (b > x & x > a) + 0.0;
end