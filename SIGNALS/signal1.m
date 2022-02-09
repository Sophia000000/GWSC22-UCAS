function [sigVec, sigVec1, sigVec2] = signal1(dataX,snr,qcCoefs)
% Generate a signal
% S = signal1(X,SNR,C)
% Generates a quadratic chirp signal S. X is the vector of
% time stamps at which the samples of the signal are to be computed. SNR is
% the matched filtering signal-to-noise ratio of S and C is the vector of
% three coefficients [a1, a2, a3] that parametrize the phase of the signal:
% cos(t)*sin(2*pi*(a*t+b*t.^2)). 

%Xiaotong Wei, Fre 2022
sigVec1=cos(dataX);
sigVec2=qcCoefs(1)*dataX + qcCoefs(2)*dataX.^2;
sigVec = cos(dataX).*sin(2*pi*(qcCoefs(1)*dataX + qcCoefs(2)*dataX.^2));

sigVec = snr*sigVec/norm(sigVec);


