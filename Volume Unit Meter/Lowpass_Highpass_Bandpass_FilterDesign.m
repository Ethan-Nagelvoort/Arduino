clearvars;
close all;
%highpass
Fc = 8000;  % in hz, highpass 3 dB frequency in Hz.change for each filter, this is cutoff 8000
Fs = 44100; % Sample rate

FcNorm = Fc/(Fs/2); % Normalized

[z,p,k] = butter(4,FcNorm, 'high');%returns the zero,poles, and gain. 4 order filter
[sos,g] = zp2sos(z,p,k);%takes zeros poles, creates finds a second order matrix(coefficients) with gain
h = fvtool(sos, 'Analysis', 'Freq');% 
h.Name = 'Highpass Filter Mag and Phase';
h;
isstable(round(32768*sos));
fprintf("highpass coeff");
sos
b01 = sos(1,1); % Section 1, b0 
b11 = sos(1,2); % Section 1, b1
b21 = sos(1,3); % Section 1, b2
a01 = sos(1,4); % Section 1, a0
a11 = sos(1,5); % Section 1, a1
a21 = sos(1,6); % Section 1, a2


b02 = sos(2,1); % Section 2, b0
b12 = sos(2,2); % Section 2, b1
b22 = sos(2,3); % Section 2, b2
a02 = sos(2,4); % Section 2, a0
a12 = sos(2,5); % Section 2, a1
a22 = sos(2,6); % Section 2, a2

% Quantize to 16 bits
b11Q = round(32768 * b11);
b11Qa = round(b11Q/2);
b11Qb = b11Q - b11Qa;

a11Q = round(32768 * a11);
a11Qa = round(a11Q/2);
a11Qb = a11Q - a11Qa;

b12Q = round(32768 * b12);
b12Qa = round(b12Q/2);
b12Qb = b11Q - b12Qa;

a12Q = round(32768 * a12);
a12Qa = round(a12Q/2);
a12Qb = a12Q - a12Qa;

a21Q = round(32768 * a21);
a22Q = round(32768 * a22);

TestFreq = 10000;%change for all filters
Ampl = 0.95;
S1Reg0 = 0; S1Reg1 = 0; S2Reg0 = 0; S2Reg1 = 0;

NumPts = 100000;
S1Reg0Array = zeros(NumPts, 1); S1Reg1Array = zeros(NumPts, 1);
S2Reg0Array = zeros(NumPts, 1); S2Reg1Array = zeros(NumPts, 1);
Sec1FilterOutArray = zeros(NumPts, 1);
Sec2FilterOutArray = zeros(NumPts, 1);

for i = 1 : NumPts
   Inp = 2^11 * Ampl * sin(2*pi*TestFreq * i/Fs) + 0;
%     Inp = 2e9;
    
    % Section 1 - Denominator
    NextReg0 = Inp - (S1Reg0*a11Qa + S1Reg0*a11Qb + S1Reg1*a21Q)/(2^15);         
 
    % Section 1 - Numerator
	Sec1FilterOutput = NextReg0 + S1Reg0 + S1Reg0 + S1Reg1;

    Sec1FilterOutArray(i) = Sec1FilterOutput;
    
    % Shift Section 1
    S1Reg1 = S1Reg0; S1Reg0 = NextReg0;
    S1Reg0Array(i,1) = S1Reg0;  S1Reg1Array(i,1) = S1Reg1;
    
	% Section 2 - Denominator   
	NextReg0 = Sec1FilterOutput - (S2Reg0 * a12Qa + S2Reg0 * a12Qb + S2Reg1 * a22Q)/(2^15);         
    
    % Section 2 - Numerator
    Sec2FilterOutput = NextReg0 + S2Reg0 + S2Reg0 + S2Reg1;         
    Sec2FilterOutArray(i) = Sec2FilterOutput;
    
    % Shift Section 2
    S2Reg1 = S2Reg0; S2Reg0 = NextReg0;


    
    S2Reg0Array(i,1) = S2Reg0;  S2Reg1Array(i,1) = S2Reg1;
end
%     
% Plot bit growth.    
% plot(1:NumPts, log(S1Reg0Array(1:NumPts, 1))./log(2)); 
% plot(1:NumPts, log(S1Reg1Array(1:NumPts, 1))./log(2)); 
% plot(1:NumPts, log(S2Reg0Array(1:NumPts, 1))./log(2)); 
% plot(1:NumPts, log(abs(S2Reg0Array(1:NumPts, 1)))./log(2)); 

% Plot register values
% plot(1:NumPts, S1Reg0Array(1:NumPts, 1)); 
% plot(1:NumPts, S2Reg0Array(1:NumPts, 1)); 



% Plot output
% plot(1:NumPts, Sec1FilterOutArray(1:NumPts, 1)); 

plot(1:NumPts, Sec2FilterOutArray(1:NumPts, 1)*g); 
title("Test for Highpass");  
%  
%lowpass
Fc = 8000;  % in hz, Lowpass 3 dB frequency in Hz.change for each filter, this is cutoff 8000
Fs = 44100; % Sample rate

FcNorm = Fc/(Fs/2); % Normalized

[z,p,k] = butter(4,FcNorm);%returns the zero,poles, and gain. 4 order filter
[sos,g] = zp2sos(z,p,k);%takes zeros poles, creates finds a second order matrix(coefficients) with gain
h = fvtool(sos, 'Analysis', 'Freq');% 
h.Name = 'Lowpass Filter Mag and Phase';
h
isstable(round(32768*sos));
fprintf("lowpass coeff");
sos
b01 = sos(1,1); % Section 1, b0 
b11 = sos(1,2); % Section 1, b1
b21 = sos(1,3); % Section 1, b2
a01 = sos(1,4); % Section 1, a0
a11 = sos(1,5); % Section 1, a1
a21 = sos(1,6); % Section 1, a2


b02 = sos(2,1); % Section 2, b0
b12 = sos(2,2); % Section 2, b1
b22 = sos(2,3); % Section 2, b2
a02 = sos(2,4); % Section 2, a0
a12 = sos(2,5); % Section 2, a1
a22 = sos(2,6); % Section 2, a2

% Quantize to 16 bits
b11Q = round(32768 * b11);
b11Qa = round(b11Q/2);
b11Qb = b11Q - b11Qa;

a11Q = round(32768 * a11);
a11Qa = round(a11Q/2);
a11Qb = a11Q - a11Qa;

b12Q = round(32768 * b12);
b12Qa = round(b12Q/2);
b12Qb = b11Q - b12Qa;

a12Q = round(32768 * a12);
a12Qa = round(a12Q/2);
a12Qb = a12Q - a12Qa;

a21Q = round(32768 * a21);
a22Q = round(32768 * a22);

TestFreq = 10000;%change for all filters
Ampl = 0.95;
S1Reg0 = 0; S1Reg1 = 0; S2Reg0 = 0; S2Reg1 = 0;

NumPts = 100000;
S1Reg0Array = zeros(NumPts, 1); S1Reg1Array = zeros(NumPts, 1);
S2Reg0Array = zeros(NumPts, 1); S2Reg1Array = zeros(NumPts, 1);
Sec1FilterOutArray = zeros(NumPts, 1);
Sec2FilterOutArray = zeros(NumPts, 1);

for i = 1 : NumPts
   Inp = 2^11 * Ampl * sin(2*pi*TestFreq * i/Fs) + 0;
%     Inp = 2e9;
    
    % Section 1 - Denominator
    NextReg0 = Inp - (S1Reg0*a11Qa + S1Reg0*a11Qb + S1Reg1*a21Q)/(2^15);         
 
    % Section 1 - Numerator
	Sec1FilterOutput = NextReg0 + S1Reg0 + S1Reg0 + S1Reg1;

    Sec1FilterOutArray(i) = Sec1FilterOutput;
    
    % Shift Section 1
    S1Reg1 = S1Reg0; S1Reg0 = NextReg0;
    S1Reg0Array(i,1) = S1Reg0;  S1Reg1Array(i,1) = S1Reg1;
    
	% Section 2 - Denominator   
	NextReg0 = Sec1FilterOutput - (S2Reg0 * a12Qa + S2Reg0 * a12Qb + S2Reg1 * a22Q)/(2^15);         
    
    % Section 2 - Numerator
    Sec2FilterOutput = NextReg0 + S2Reg0 + S2Reg0 + S2Reg1;         
    Sec2FilterOutArray(i) = Sec2FilterOutput;
    
    % Shift Section 2
    S2Reg1 = S2Reg0; S2Reg0 = NextReg0;


    
    S2Reg0Array(i,1) = S2Reg0;  S2Reg1Array(i,1) = S2Reg1;
end

plot(1:NumPts, Sec2FilterOutArray(1:NumPts, 1)*g);
title("Test for Lowpass");
% fvtool(sos, 'Fs', 44100);

%bandpass

%Fc = 8000;  % in hz, Lowpass 3 dB frequency in Hz.change for each filter, this is cutoff
Fs = 44100; % Sample rate

FcNorm = [4000 6000]/(Fs/2); % Normalized

[z,p,k] = butter(2,FcNorm, 'bandpass');%returns the zero,poles, and gain. 4 order filter
[sos,g] = zp2sos(z,p,k);%takes zeros poles, creates finds a second order matrix(coefficients) with gain
h = fvtool(sos, 'Analysis', 'Freq');% 
h.Name = 'Bandpass Filter Mag and Phase';
h
isstable(round(32768*sos));
fprintf("bandpass coeff");
sos
b01 = sos(1,1); % Section 1, b0 
b11 = sos(1,2); % Section 1, b1
b21 = sos(1,3); % Section 1, b2
a01 = sos(1,4); % Section 1, a0
a11 = sos(1,5); % Section 1, a1
a21 = sos(1,6); % Section 1, a2


b02 = sos(2,1); % Section 2, b0
b12 = sos(2,2); % Section 2, b1
b22 = sos(2,3); % Section 2, b2
a02 = sos(2,4); % Section 2, a0
a12 = sos(2,5); % Section 2, a1
a22 = sos(2,6); % Section 2, a2

% Quantize to 16 bits
b11Q = round(32768 * b11);
b11Qa = round(b11Q/2);
b11Qb = b11Q - b11Qa;

a11Q = round(32768 * a11);
a11Qa = round(a11Q/2);
a11Qb = a11Q - a11Qa;

b12Q = round(32768 * b12);
b12Qa = round(b12Q/2);
b12Qb = b11Q - b12Qa;

a12Q = round(32768 * a12);
a12Qa = round(a12Q/2);
a12Qb = a12Q - a12Qa;

a21Q = round(32768 * a21);
a22Q = round(32768 * a22);

TestFreq = 5000;%change for all filters
Ampl = 0.95;
S1Reg0 = 0; S1Reg1 = 0; S2Reg0 = 0; S2Reg1 = 0;

NumPts = 100000;
S1Reg0Array = zeros(NumPts, 1); S1Reg1Array = zeros(NumPts, 1);
S2Reg0Array = zeros(NumPts, 1); S2Reg1Array = zeros(NumPts, 1);
Sec1FilterOutArray = zeros(NumPts, 1);
Sec2FilterOutArray = zeros(NumPts, 1);

for i = 1 : NumPts
   Inp = 2^11 * Ampl * sin(2*pi*TestFreq * i/Fs) + 0;
%     Inp = 2e9;
    
    % Section 1 - Denominator
    NextReg0 = Inp - (S1Reg0*a11Qa + S1Reg0*a11Qb + S1Reg1*a21Q)/(2^15);         
 
    % Section 1 - Numerator
	Sec1FilterOutput = NextReg0 + S1Reg0 + S1Reg0 + S1Reg1;

    Sec1FilterOutArray(i) = Sec1FilterOutput;
    
    % Shift Section 1
    S1Reg1 = S1Reg0; S1Reg0 = NextReg0;
    S1Reg0Array(i,1) = S1Reg0;  S1Reg1Array(i,1) = S1Reg1;
    
	% Section 2 - Denominator   
	NextReg0 = Sec1FilterOutput - (S2Reg0 * a12Qa + S2Reg0 * a12Qb + S2Reg1 * a22Q)/(2^15);         
    
    % Section 2 - Numerator
    Sec2FilterOutput = NextReg0 + S2Reg0 + S2Reg0 + S2Reg1;         
    Sec2FilterOutArray(i) = Sec2FilterOutput;
    
    % Shift Section 2
    S2Reg1 = S2Reg0; S2Reg0 = NextReg0;


    
    S2Reg0Array(i,1) = S2Reg0;  S2Reg1Array(i,1) = S2Reg1;
end
plot(1:NumPts, Sec2FilterOutArray(1:NumPts, 1)*g);
title("Test for Bandpass");

% % fvtool(sos, 'Fs', 44100);