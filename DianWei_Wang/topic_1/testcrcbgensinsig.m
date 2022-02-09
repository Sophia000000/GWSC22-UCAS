close all; clear; clc;
%% Plot the sinusoidal signal
% Signal parameters
phi0 = 0;
f0 = 25;
A = 10;
% Instantaneous frequency after 1 sec is 
maxFreq = f0;
samplFreq = 5*maxFreq;
samplIntrvl = 1/samplFreq;

% Time samples
timeVec = 0:samplIntrvl:1.0;
% Number of samples
nSamples = length(timeVec);

% Generate the signal
sigVec = crcbgensinsig(timeVec,A,[phi0,f0]);

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

%% Generating sinusoidal signal
function sigVec = crcbgensinsig(dataX,snr,Coefs)
% Generate a sinusoidal signal
% S = crcbgenlcsig(X,SNR,C)
% Generates a sinusoidal signal S. X is the vector of
% time stamps at which the samples of the signal are to be computed. SNR is
% the matched filtering signal-to-noise ratio of S and C is the vector of
% two coefficients [phi0, f0] that parametrize the phase of the signal:
% phi0+f0*t. 

%Soumya D. Mohanty, May 2018
%DianWei Wang modified, Feb 2022

phaseVec = Coefs(1) + Coefs(2)*dataX;
sigVec = sin(2*pi*phaseVec);
sigVec = snr*sigVec/norm(sigVec);
end