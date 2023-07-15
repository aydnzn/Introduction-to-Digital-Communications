clear;
%AYDIN UZUN
%2015401210
%EE 477 HW#1
%Please install Communications Toolbox to run this code. Because this code
%has some toolbox specific functions and classes.
%% BPSK
%a)
%set frequency
f_bpsk=5;
%set the number of distinct symbols in BPSK
M_bpsk=2;
% % comm.PSKModulator creates a modulator System object,
% % MODULATOR. This object modulates the input signal using the M-ary phase
% % shift keying (M-PSK) method. Use binary symbolmapping. 0 = PhaseOffset
modulator_obj_bpsk = comm.PSKModulator(M_bpsk, 0,'SymbolMapping', 'binary');
% constellation(OBJ) generates a constellation plot.
constellation(modulator_obj_bpsk);
%b) plot symbols
% Y = step(H,X) modulates input data, X, with the PSK modulator System object, 
% H. It returns the baseband modulated output, Y. Depending on the value of
% the BitInput property, input X can be an integer or bit valued column 
%vector with numeric, logical, or fixed-point data types. 
% [0;1] is my symbols array, because I used binary symbolmapping.
% Baseband modulated output is created.
baseband_mod_out_symbol_bpsk = step(modulator_obj_bpsk,[0;1]);
%create a time space for plot
time_space_for_symbols_bpsk = linspace(0,2,2000);
% rectangular pulse shaped version of modulated output
Rect_symbol_bpsk = rectpulse(baseband_mod_out_symbol_bpsk,1000);
%This is my symbol on IQ plane. 
IQ_symbol_bpsk = real(Rect_symbol_bpsk).*cos(2*pi*f_bpsk*time_space_for_symbols_bpsk)'+imag(Rect_symbol_bpsk).*sin(2*pi*f_bpsk*time_space_for_symbols_bpsk)';

figure(2);
subplot(2,1,1);
plot(time_space_for_symbols_bpsk,IQ_symbol_bpsk);
xlabel('time(s)');
ylabel('amplitude(V)');
title('waveforms to represent the symbols in the constellation diagram');

%c) plot stream
%set s
s = [1;0;1;1;1;1;0;0];
% create timespace to plot the stream
time_space_for_stream_bpsk = linspace(0,8,8000);
% baseband modulated output of stream
baseband_mod_out_stream_bpsk = step(modulator_obj_bpsk,s);
%Rectangular shaped version
Rect_stream_bpsk = rectpulse(baseband_mod_out_stream_bpsk,1000);
%the stream on IQ plane
IQ_stream_bpsk = real(Rect_stream_bpsk).*cos(2*pi*f_bpsk*time_space_for_stream_bpsk)'+imag(Rect_stream_bpsk).*sin(2*pi*f_bpsk*time_space_for_stream_bpsk)';
figure(2);
subplot(2,1,2);
plot(time_space_for_stream_bpsk,IQ_stream_bpsk);
xlabel('time(s)');
ylabel('amplitude(V)');
title('waveforms to represent modulated pulse stream');

clear;
%% QPSK
clear;
% The algorithm to plot signals and diagrams is the same for QPSK, 4PAM and
% 16QAM
%set frequency
f_qpsk=5;
%set the number of distinct symbols in QPSK
M_qpsk=4;
% % comm.PSKModulator creates a modulator System object,
% % MODULATOR. This object modulates the input signal using the M-ary phase
% % shift keying (M-PSK) method. Use binary symbolmapping. 0 = PhaseOffset
modulator_obj_qpsk = comm.PSKModulator(M_qpsk, 0,'SymbolMapping', 'binary');
% constellation(OBJ) generates a constellation plot.
constellation(modulator_obj_qpsk);
%b) plot stmbols
% Y = step(H,X) modulates input data, X, with the PSK modulator System object, 
% H. It returns the baseband modulated output, Y. Depending on the value of
% the BitInput property, input X can be an integer or bit valued column 
%vector with numeric, logical, or fixed-point data types.

% Baseband modulated output is created.
baseband_mod_out_symbol_qpsk = step(modulator_obj_qpsk,[0;1;2;3]);
%create a time space for plot
time_space_for_symbols_qpsk = linspace(0,4,4000);
% rectangular pulse shaped version of modulated output
Rect_symbol_qpsk = rectpulse(baseband_mod_out_symbol_qpsk,1000);
%the stream on IQ plane
IQ_symbol_qpsk = real(Rect_symbol_qpsk).*cos(2*pi*f_qpsk*time_space_for_symbols_qpsk)'+imag(Rect_symbol_qpsk).*sin(2*pi*f_qpsk*time_space_for_symbols_qpsk)';
figure(4);
subplot(2,1,1);
plot(time_space_for_symbols_qpsk,IQ_symbol_qpsk);
xlabel('time(s)');
ylabel('amplitude(V)');
title('waveforms to represent the symbols in the constellation diagram');

%c) plot stream
% The stream was [1;0;1;1;1;1;0;0]
% correspondence between new stream and old stream 
% 1;0 = 2, 1;1 = 3 , 0;0 =0,  0;1 = 1
s = [2;3;3;0];
% create timespace to plot the stream
time_space_for_stream_qpsk = linspace(0,4,4000);
% baseband modulated output of stream
baseband_mod_out_stream_qpsk = step(modulator_obj_qpsk,s);
%Rectangular shaped version
Rect_stream_qpsk = rectpulse(baseband_mod_out_stream_qpsk,1000);
% the stream on IQ plane
IQ_stream_qpsk = real(Rect_stream_qpsk).*cos(2*pi*f_qpsk*time_space_for_stream_qpsk)'+imag(Rect_stream_qpsk).*sin(2*pi*f_qpsk*time_space_for_stream_qpsk)';
figure(4);
subplot(2,1,2);
plot(time_space_for_stream_qpsk,IQ_stream_qpsk);
xlabel('time(s)');
ylabel('amplitude(V)');
title('waveforms to represent modulated pulse stream');

clear;
%% 4_PAM
clear;
%set frequency
f_4pam=5;
%set the number of distinct symbols in 4PAM
M_4pam=4;
%The PAMModulator object modulates using M-ary pulse amplitude modulation. 
%The output is a baseband representation of the modulated signal. 
%The M-ary number parameter, M, represents the number of points in the 
%signal constellation and requires an even integer.
% binary symbopl mapping is used and the average power is normalized to
% unity.
modulator_obj_4pam = comm.PAMModulator(M_4pam,'SymbolMapping', 'binary', 'NormalizationMethod','Average Power');
constellation(modulator_obj_4pam);
% constellation(OBJ) generates a constellation plot.

% b) plot symbols
% Y = step(H,X) modulates input data, X, with the PSK modulator System object, 
% H. It returns the baseband modulated output, Y. Depending on the value of
% the BitInput property, input X can be an integer or bit valued column 
%vector with numeric, logical, or fixed-point data types.

% Baseband modulated output is created.
baseband_mod_out_symbol_4pam = step(modulator_obj_4pam,[0;1;2;3]);
%create a time space for plot
time_space_for_symbols_4pam = linspace(0,4,4000);
% rectangular pulse shaped version of modulated output
Rect_symbol_4pam = rectpulse(baseband_mod_out_symbol_4pam,1000);
%the stream on IQ plane, stream consists of each symbol. 
IQ_symbol_4pam = real(Rect_symbol_4pam).*cos(2*pi*f_4pam*time_space_for_symbols_4pam)'+imag(Rect_symbol_4pam).*sin(2*pi*f_4pam*time_space_for_symbols_4pam)';
figure(6);
subplot(2,1,1);
plot(time_space_for_symbols_4pam,IQ_symbol_4pam);
xlabel('time(s)');
ylabel('amplitude(V)');
title('waveforms to represent the symbols in the constellation diagram');


%c) plot stream
% The stream was [1;0;1;1;1;1;0;0]
% correspondence between new stream and old stream 
% 1;0 = 2, 1;1 = 3 , 0;0 =0,  0;1 = 1
s = [2;3;3;0];
% create timespace to plot the stream
time_space_for_stream_4pam = linspace(0,4,4000);
% baseband modulated output of stream
baseband_mod_out_stream_4pam = step(modulator_obj_4pam,s);
%Rectangular shaped version
Rect_stream_4pam = rectpulse(baseband_mod_out_stream_4pam,1000);
% the stream on IQ plane
IQ_stream_4pam = real(Rect_stream_4pam).*cos(2*pi*f_4pam*time_space_for_stream_4pam)'+imag(Rect_stream_4pam).*sin(2*pi*f_4pam*time_space_for_stream_4pam)';
figure(6);
subplot(2,1,2);
plot(time_space_for_stream_4pam,IQ_stream_4pam);
xlabel('time(s)');
ylabel('amplitude(V)');
title('waveforms to represent modulated pulse stream');
clear;

%% 16-QAM
%set frequency
f_16qam=5;
%set the number of distinct symbols in 16-QAM
M_16qam=16;
%The RectangularQAMModulator object modulates using M-ary quadrature amplitude
%modulation with a constellation on a rectangular lattice. The output is a 
%baseband representation of the modulated signal.
% binary symbol mapping is used and the average power is normalized to
% unity.
modulator_obj_16qam = comm.RectangularQAMModulator(M_16qam,'SymbolMapping', 'binary', 'NormalizationMethod','Average Power');
constellation(modulator_obj_16qam);
% constellation(OBJ) generates a constellation plot.
%b) plot symbols
% Y = step(H,X) modulates input data, X, with the PSK modulator System object, 
% H. It returns the baseband modulated output, Y. Depending on the value of
% the BitInput property, input X can be an integer or bit valued column 
%vector with numeric, logical, or fixed-point data types.

% Baseband modulated output is created. [0:15]' is used, because there are
% 16 distinct symbols
baseband_mod_out_symbol_16qam = step(modulator_obj_16qam,[0:15]');
%create a time space for plot
time_space_for_symbols_16qam = linspace(0,16,16000);
%Rectangular shaped version
Rect_symbol_16qam = rectpulse(baseband_mod_out_symbol_16qam,1000);
%the stream on IQ plane, stream consists of each symbol. 
IQ_symbol_16qam = real(Rect_symbol_16qam).*cos(2*pi*f_16qam*time_space_for_symbols_16qam)'+imag(Rect_symbol_16qam).*sin(2*pi*f_16qam*time_space_for_symbols_16qam)';
figure(8);
subplot(2,1,1);
plot(time_space_for_symbols_16qam,IQ_symbol_16qam);
xlabel('time(s)');
ylabel('amplitude(V)');
title('waveforms to represent the symbols in the constellation diagram');

% c) plot stream (modulated signal)
% The stream was [1;0;1;1;1;1;0;0]
% correspondence between new stream and old stream 
% 1;0;1;1 = 11, 1;1;0:0= 12
s = [11;12];
% create timespace to plot the stream
time_space_for_stream_16qam = linspace(0,2,2000);
% baseband modulated output of stream
baseband_mod_out_stream_16qam = step(modulator_obj_16qam,s);
%Rectangular shaped version
Rect_stream_16qam = rectpulse(baseband_mod_out_stream_16qam,1000);
% the stream on IQ plane
IQ_stream_16qam = real(Rect_stream_16qam).*cos(2*pi*f_16qam*time_space_for_stream_16qam)'+imag(Rect_stream_16qam).*sin(2*pi*f_16qam*time_space_for_stream_16qam)';
figure(8);
subplot(2,1,2);
plot(time_space_for_stream_16qam,IQ_stream_16qam);
xlabel('time(s)');
ylabel('amplitude(V)');
title('waveforms to represent modulated pulse stream');
clear;

%% BFSK
clear;
% set frequency
f_bfsk=5;
%set the number of distinct symbols in BPSK
M_bfsk=2;
% The FSKModulator object modulates using the M-ary frequency shift keying method. 
% The output is a baseband representation of the modulated signal.
% Binary Symbol mapping is used.
modulator_obj_bfsk = comm.FSKModulator(M_bfsk,50,'SymbolMapping', 'binary');
% But one cannot plot with constellation function. Therefore I used a
% method from official mathworks website.
% https://www.mathworks.com/matlabcentral/answers/405733-how-do-i-plot-the-constellation-diagram-for-a-binary-fsk-modulator
 modulator_obj_bfsk.SamplesPerSymbol = 1;
 data = [1;0];
 modSignal = modulator_obj_bfsk(data);  
 scope = comm.ConstellationDiagram('ShowReferenceConstellation', false);
 scope(modSignal);
 %b) plot symbols
 % time space for one symbol
time_space_for_one_symbol_bfsk = linspace(0,1,1000);
% time space two plor two symbols side by side.
time_space_for_symbols_bfsk = linspace(0,2,2000);

% fo represents waveform for binary '0'
% f1 represents waveform for binary '1'
% the difference between frequencies is 1, because symbolling period is 1s.
f0 = cos(2*pi*f_bfsk*time_space_for_one_symbol_bfsk);
f1 = cos(2*pi*(f_bfsk+1)*time_space_for_one_symbol_bfsk);
% concatanated symbol representations
IQ_symbols_bfsk = [f0,f1];
figure(10);
subplot(2,1,1);
plot( time_space_for_symbols_bfsk, IQ_symbols_bfsk );
xlabel('time(s)');
ylabel('amplitude(V)');
title('waveforms to represent the symbols in the constellation diagram');

% c) plot stream
% The stream was [1;0;1;1;1;1;0;0]
% correspondence between new stream and old stream 
% 0 = f0, 1= f1
s = [f1,f0,f1,f1,f1,f1,f0,f0];
%create timespace to plot stream 
time_space_for_stream_bfsk = linspace(0,8,8000);
figure(10);
subplot(2,1,2);
plot(time_space_for_stream_bfsk,s);
xlabel('time(s)');
ylabel('amplitude(V)');
title('waveforms to represent modulated pulse stream');
clear;
%% PART2
%a)
p1 = [1;0;0;0];
p2 = [0;1;0;0];
p3 = [0;0;1;0];
p4 = [0;0;0;1];
%forms an orthonormal basis because
NORM = [norm(p1); norm(p2) ;norm(p3) ;norm(p4) ];
% is ones array
%and
DOT = [p1.*p2 ; p1.*p3; p1.*p4; p2.*p3; p2.*p4; p3.*p4 ];
% is zeros array.
s1 = 2*p1 - p2 - p3 - p4;
s2 = -2*p1 + p2 +p3;
s3 = p1 -p2 +p3 -p4;
s4 = p1 - p2 -2*p3 +p4;
%b)
%Gram-Schmidt Method
v1 = s1;
v2 = s2 - (dot(s2,v1)/dot(v1,v1))*v1;
v3 = s3 - (dot(s3,v1)/dot(v1,v1))*v1 - (dot(s3,v2)/dot(v2,v2))*v2;
v4 = s4 - (dot(s4,v1)/dot(v1,v1))*v1 - (dot(s4,v2)/dot(v2,v2))*v2 - (dot(s4,v3)/dot(v3,v3))*v3;
%Normalizing
v1_normalized = v1/norm(v1); 
v2_normalized = v2/norm(v2); 
v3_normalized = v3/norm(v3); 
v4_normalized = v4/norm(v4);
% The calculated and normalized vectors build an orthonormal basis.

Basis_Matrix = [v1_normalized , v2_normalized, v3_normalized, v4_normalized];
% Si = Basis_Matrix * Coefficients
% Coefficients_i = inv(Basis_Matrix) * Si
Inv_Basis_Matrix = inv(Basis_Matrix);
s1_coef = Inv_Basis_Matrix * s1;
s2_coef = Inv_Basis_Matrix * s2;
s3_coef = Inv_Basis_Matrix * s3;
s4_coef = Inv_Basis_Matrix * s4;
