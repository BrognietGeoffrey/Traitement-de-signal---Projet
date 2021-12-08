%%
clc; close all; clear all;
% wfdb2mat('fetal_PCG_p01_GW_36'); %conversion

load('fetal_PCG_p01_GW_36m.mat'); %https://physionet.org/content/fpcgdb/1.0.0/#files-panel

x = 0:0.003:1199.999; %-> 1 echantillons ttes les 0.003sec
%executer ->[tm,signal,Fs,siginfo] = rdmat('fetal_PCG_p01_GW_36m');
%(après avoir load) et ca correspond à la variable tm

%%
h=val(1,1:2000);
Fs=333;
t=(0:length(h)-1)/Fs;

%% Déterminer le bruit

Fn = Fs/2;
ecg_raw = h;
L = numel(ecg_raw);
time = linspace(0, L, L)/Fs;
N = 2^nextpow2(L);
FTs = fft(ecg_raw,N)/L;
Fv = linspace(0, 1, N/2+1)*Fn;
Iv = 1:numel(Fv);
figure
plot(Fv, abs(FTs(Iv))*2)
grid

%% test denoising with sgolayfilt
    
rd = 3;
fl = 77;
smtlb = sgolayfilt(h,rd,fl);

figure;
subplot(2,1,1)
plot(t,h)

title('Original in time response')
grid

subplot(2,1,2)
plot(t,smtlb)

title('Filtered using  Savitzky-Golay in time response')
grid

%% test denoising with butter

% [b1,a1]=butter(1,[0.3 0.35],'stop'); 
[b1,a1]=butter(1,0.06,'low');

butter_filter1=filter(b1,a1,h);
figure;
subplot(2,1,1)
plot(t,h)
grid;
subplot(2,1,2)
plot(t,butter_filter1)
grid;
title('filtering using butterworth in time domain')

%% autre test

d1 = designfilt('lowpassiir','FilterOrder',12, ...
    'HalfPowerFrequency',0.15,'DesignMethod','butter');
y = filtfilt(d1,h);

figure
subplot(1,1,1)
plot(h)
hold on
plot(y,'LineWidth',2)
legend('Noisy ECG','Zero-Phase Filtering')

%% Séparation fecg-mecg
