% Mother & Fetus hearbeat detection
%DB : https://physionet.org/content/fpcgdb/1.0.0/
clear all; close all; clc;

Fs = 333;
Time = 40;
NumSamp = Time * Fs;
%wfdb2mat('fetal_PCG_p01_GW_36'); %convert .dat to .mat
load fetal_PCG_p01_GW_36m.mat;
load Hd; 


%% generate synthetic Mom's heartbeat

x1 = 3.5*ecg(2700).'; % gen synth ECG signal
y1 = sgolayfilt(kron(ones(1,ceil(NumSamp/2700)+1),x1),0,21); % repeat for NumSamp length and smooth
n = 1:Time*Fs';
del = round(2700*rand(1)); % pick a random offset
mhb = y1(n + del)'; %construct the ecg signal from some offset
t = 1/Fs:1/Fs:Time';


%% generate synthetic Fetus heartbeat

x2 = 0.25*ecg(1725);
y2 = sgolayfilt(kron(ones(1,ceil(NumSamp/1725)+1),x2),0,17);
del = round(1725*rand(1));
fhb = y2(n + del)';


%% The measured signal

Wopt = [0 1.0 -0.5 -0.8 1.0  -0.1 0.2 -0.3 0.6 0.1];
%Wopt = rand(1,10);

data = val(1, 1:13320) ;  
d = filter(Wopt,1, data); 


%% Measured Mom's heartbeat (for reference)

x = mhb + 0.02*randn(size(mhb));

subplot(3,2,1); plot(t,x);
axis([0 20 -4 4]);
grid;
xlabel('Time [sec]');
ylabel('Voltage [mV]');
title('Reference Signal');


%% Applying the filter

h = dsp.LMSFilter('Length',9,'Method','Normalized LMS', ...
              'StepSize',0.001);
d = d';
[y,e] = step(h, d, x);


subplot(3,2,2); plot(t,d,'c');
axis([0 20.0 -4 400]);
grid;
xlabel('Time [sec]');
ylabel('Voltage [mV]');
title('Original Signal ');
legend('Signal');

%% Extract the mother heartbeat 
% Nous avons essayé de transformer cette partie de code 
% pour qu'elle fonctionne avec nos données mais le résultat ne peut qu'être 
% peu convaincant car la variable 'x', contenant le signal mecgsyn, est utilisée.

subplot(3,2,3); plot(t,e,'r');
axis([0 20 -0.5 0.5]);
grid on;
xlabel('Time [sec]');
ylabel('Voltage [mV]');
title('Mother Heartbeat');
legend('signal of mother');


%% Counting the peaks

[num,den] = fir1(100,100/2000);
filt_e = filter(Hd,e);
subplot(3,2,4); plot(t,fhb,'r'); hold on; plot(t,filt_e,'b');
axis([0 20 -2 2]);
grid on;
xlabel('Time [sec]');
ylabel('Voltage [mV]');
title('Filtered signal');
legend('Fetus reference','Fetus filtered');


thresh = 4*mean(abs(filt_e))*ones(size(filt_e));
peak_e = (filt_e >= thresh);
edge_e = (diff([0; peak_e]) >0);
subplot(3,2,5); plot(t,filt_e,'c'); hold on; plot(t,thresh,'r'); plot(t,peak_e,'b');
xlabel('Time [sec]');
ylabel('Voltage [mV]');
title('Peak detection');
legend('Filtered fetus','Threshold','Peak marker', 'Location','Southwest');
axis([0 20 -2 2]);


subplot(3,2,6); plot(t,filt_e,'r'); hold on; plot(t,edge_e,'b'); plot(0,0,'w');
xlabel('Time [sec]');
ylabel('Voltage [mV]');
title('Reconstructed fetus signal');
legend('Fetus Sig','Edge marker', 'Location','SouthEast');
axis([0 20 -2 2]);