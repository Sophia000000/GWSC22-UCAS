%% Plot the quadratic chirp signal
% Signal parameters
a1=10;
a2=3;
A = 10;
% Instantaneous frequency after 1 sec is 
maxFreq = a1+2*a2;
samplFreq = 5*maxFreq;
samplIntrvl = 1/samplFreq;

% Time samples
timeVec = 0:samplIntrvl:1.0;
% Number of samples
nSamples = length(timeVec);

% Generate the signal
[sigVec,sigVec1,sigVec2] = signal1(timeVec,A,[a1,a2]);

%Plot the signal 
figure;
plot(timeVec,sigVec,'Marker','.','MarkerSize',24);
figure;
hold on;
plot(timeVec,sigVec1);
plot(timeVec,sin(sigVec2));
plot(timeVec,a1+2*a2*timeVec);
legend("a(t)","sin(\phi(t))","d\phi(t)/dt")
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

% Time samples
the_fre=5*kNyq;
the_Fre=0.5*kNyq;
samplIntrvl=1/the_fre;
samplIntrv2=1/the_Fre;
timeVec = 0:samplIntrvl:1.0;
timeVec2 = 0:samplIntrv2:1.0;
% Number of samples
nSamples = length(timeVec);
nSamples2 = length(timeVec2);

% Generate the signal
[sigVec5,sigVec15,sigVec25] = signal1(timeVec,A,[a1,a2]);
[sigVecA2,sigVec12,sigVec22] = signal1(timeVec2,A,[a1,a2]);
figure;
hold on;
plot(timeVec,sigVec5,'Marker','.','MarkerSize',24);
plot(timeVec2,sigVecA2,'Marker','.','MarkerSize',24);
legend('5*Nyquist','0.5*Nyquist')
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
