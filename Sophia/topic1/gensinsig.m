function sigVec = gensinisg(dataX,snr,qcCoefs)
% Generate a sinusoidal signal
% S = CRCBGENQSIG(X,SNR,C)
% Generates a sinusoidal signal S.
% X is the vector of time stamps.
% SNR is the matched filtering signal-to-noise ratio of S. 
% C is the vector of three coefficients [f0, phi0] 
% phi(t)=f0*t

%…Ú∆º, 2022.2.8

phaseVec = qcCoefs(1)*dataX;
sigVec = sin(2*pi*phaseVec + qcCoefs(2));
sigVec = snr*sigVec/norm(sigVec);