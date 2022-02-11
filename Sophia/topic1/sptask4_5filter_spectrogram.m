%% Task 4��task 5 ��ͨ����ͨ����ͨ�˲���ʱƵͼ %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ���������������źŵ��Ӷ��ɵ��ź�
sampFreq = 1024;
nSamples = 2048;
timeVec = (0:(nSamples-1))/sampFreq;

% Signal parameters
A1=10;
f1=100;
phi1=0;

A2=5;
f2=200;
phi2=pi/6;

A3=2.5;
f3=300;
phi3=pi/4;

% Signal length
sigLen = (nSamples-1)/sampFreq;
%Maximum frequency �ο�˹��Ƶ�� 
maxFreq1= f1;
maxFreq2= f2;
maxFreq3= f3;
maxFreq=max([maxFreq1,maxFreq2,maxFreq3]);
disp(['The maximum frequency of the discrete time sinusoid  is ', num2str(maxFreq)]);

% Generate signal
sigVec1 = gensinsig(timeVec,A1,[f1,phi1]);
sigVec2 = gensinsig(timeVec,A2,[f2,phi2]);
sigVec3 = gensinsig(timeVec,A3,[f3,phi3]);
sigVec = sigVec1 + sigVec2 + sigVec3;

%% �˲�
% Use Matlab's `fir1` function to design 3 different filters such that 
% filter #i alows only sigVeci to pass through.

% Remove frequencies such that sigVec1/sigVec2/sigVec3 can pass through
% Design low/band/high pass filters
%Length of data 

filt0rdr = 30;%�˲�������
b1 = fir1( filt0rdr, ((maxFreq1+maxFreq2)/2)/(sampFreq/2), 'low' );%bΪ�˲�������
b2 = fir1( filt0rdr, [((maxFreq1+maxFreq2)/2)/(sampFreq/2),((maxFreq2+maxFreq3)/2)/(sampFreq/2)], 'bandpass' );
b3 = fir1( filt0rdr, ((maxFreq2+maxFreq3)/2)/(sampFreq/2), 'high' );
% Apply filters �˲�
filtSig1 = fftfilt(b1, sigVec);
filtSig2 = fftfilt(b2, sigVec);
filtSig3 = fftfilt(b3, sigVec);

%% �˲�����ź�������Ҷ�任
dataLen = timeVec(end)-timeVec(1);
%DFT sample corresponding to Nyquist frequency
kNyq = floor(nSamples/2)+1;
% Positive Fourier frequencies
posFreq = (0:(kNyq-1))*(1/dataLen);

% FFT of signals  �˲����źŸ���Ҷ�任
fftSig1 = fft(filtSig1);
fftSig2 = fft(filtSig2);
fftSig3 = fft(filtSig3);
% Discard negative frequencies
fftSig1 = fftSig1(1:kNyq);
fftSig2 = fftSig2(1:kNyq);
fftSig3 = fftSig3(1:kNyq);
%% Plots
% figure;
% hold on;
% plot(timeVec,sigVec);
% plot(timeVec,filtSig1);
% plot(timeVec,filtSig2);
% plot(timeVec,filtSig3);
% title('ԭʼ�ź����˲����ź�')

figure;
subplot(2,3,1);
plot(timeVec,sigVec,timeVec,filtSig1);
title('��ͨ�˲�')
subplot(2,3,2);
plot(timeVec,sigVec,timeVec,filtSig2);
title('��ͨ�˲�')
subplot(2,3,3);
plot(timeVec,sigVec,timeVec,filtSig3);
title('��ͨ�˲�')

subplot(2,3,4);
plot(posFreq,abs(fftSig1));
title('��ͨ�˲�Ƶ��')
xlabel('Frequency (Hz)');

subplot(2,3,5);
plot(posFreq,abs(fftSig2));
title('��ͨ�˲�Ƶ��')
xlabel('Frequency (Hz)');

subplot(2,3,6);
plot(posFreq,abs(fftSig3));
title('��ͨ�˲�Ƶ��')
xlabel('Frequency (Hz)');

%% task5 ʱƵͼ
%% ԭʼ�ź���ͼ
figure;
subplot(3,1,1);
plot(timeVec,sigVec);
title('ԭʼ�ź�');xlabel('Time (sec)');

% FFT of signals  ԭʼ�źŸ���Ҷ�任
fftSig = fft(sigVec);
% Discard negative frequencies
fftSig = fftSig(1:kNyq);

subplot(3,1,2);
plot(posFreq,abs(fftSig));
title('ԭʼ�ź�Ƶ��');xlabel('Frequency (Hz)');

subplot(3,1,3);
spectrogram(sigVec, 256,250,[],sampFreq);
title('ԭʼ�ź�ʱƵͼ');

