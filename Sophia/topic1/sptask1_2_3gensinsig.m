%% task1 ��task2�� task3  ���������źŲ���ͼ���ԱȲ�ͬ����Ƶ�ʵĽ��������Ƶ��ͼ
%% Signal parameters
f0=50;
phi0=pi;
A = 10;%�����
% Instantaneous frequency after 1 sec is 
maxFreq = f0;
samplFreq = 5*maxFreq;%����Ƶ��������maxFreq������
samplFreq2 = 0.5*maxFreq;%�������Ƶ�ʲ���
samplIntrvl = 1/samplFreq;
samplIntrvl2 = 1/samplFreq2;
% Time samples
timeVec = 0:samplIntrvl:1.0;
timeVec2 = 0:samplIntrvl2:1.0;
% Number of samples
nSamples = length(timeVec);
nSamples2 = length(timeVec2);

% Generate the signal
sigVec = gensinsig(timeVec,A,[f0,phi0]);
sigVec2 = gensinsig(timeVec2,A,[f0,phi0]);

%Plot the signal 
figure;
plot(timeVec,sigVec,'Marker','.','MarkerSize',24);
hold;
plot(timeVec2,sigVec2,'Marker','*','MarkerSize',8);

%% task3 Plot the periodogram��magnitude of the FFT��
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