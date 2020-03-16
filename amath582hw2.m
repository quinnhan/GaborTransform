%%%%%%%%%%%%%%%%%%%%%% HW2: Gabor transform %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Han Song
% AMATH 582

% Part 1: Handel
clear all; close all; clc
load handel
v = y'/2;
figure(1)
plot((1:length(v))/Fs,v);
xlabel('Time [sec]');
ylabel('Amplitude');
title('Handel');
% This is in time domain

% p8 = audioplayer(v,Fs);
% playblocking(p8);

L = 9;
n = length(v);
t2 = (1:n+1)/Fs;
t = t2(1:n); 
k = (2*pi/L)*[0:n/2-1 -n/2:-1];
ks = fftshift(k);

vt = fft(v(1:end-1));
figure(2)
subplot(2,1,1)
plot(t,v)
title('Time Domain')
xlabel('time[sec]')
ylabel('Amplitude')
subplot(2,1,2)
plot(ks,fftshift(vt))
title('Frequency Domain')
xlabel('frequency [Hz]')
ylabel('n')
% Gabor

tslide = 0:0.1:9;
s = zeros(length(tslide),length(t)-1);
for j = 1:length(tslide)
    %g = exp(-(t-tslide(j)).^10); 
    %g = exp(-0.01*(t-tslide(j)).^10); %bigger window 
    %g = exp(-100*(t-tslide(j)).^10); %smaller window: localizing in time
    g = (1-(t-tslide(j)).^2).*exp(-(t-tslide(j)).^2); %Mexican hat
    %g = (1-0.1*(t-tslide(j)).^2).*exp(-0.1*(t-tslide(j)).^2); %Mexican hat big
    %g = (1-10*(t-tslide(j)).^2).*exp(-10*(t-tslide(j)).^2); %Mexican hat small
    %g = exp(-(t-tslide(j)).^2); %Gaussian Wavelet
    %g = exp(-0.1*(t-tslide(j)).^2); %Gaussian Big
    %g = exp(-10*(t-tslide(j)).^2); %Gaussian Small
    %g = abs(t-tslide(j)) <= 0.5; %Shannon function
    %g = 0.1*abs(t-tslide(j)) <= 0.5; %Shannon function big 
    %g = 10*abs(t-tslide(j)) <= 0.5; %Shannon function small
    
    vf = g.*v; % signal filter
    vft = fft(vf(1:end-1));
    s(j,:) = abs(fftshift(vft)); %save fourier transform of the signal
    
    figure(3)
    subplot(3,1,1)
    plot(t,v,'k-',t,g,'r-')
    title('Signal and Filter')
    xlabel('time')
    ylabel('Amplitude')
    legend('signal','filter','Location','North')
    subplot(3,1,2)
    plot(t,vf,'k-')
    xlabel('time')
    ylabel('Filtered Signal Amplitude')
    axis([0 9 -0.5 0.5])
    subplot(3,1,3)
    plot(fftshift(k),abs(fftshift(vft))/max(abs(fftshift(vft))),'b')
    title('FFT of Filtered Signal')
    xlabel('Frequency [Hz]')
    ylabel('n')
    axis([-200 200 0 0.2])
    pause(0.2)
end

figure(4)
pcolor(tslide,fftshift(k),s.'),shading interp, colormap(hot)
axis([0 7 -8000 8000])
title('Spectrogram using Mexican Hat Filter')


%% Part 2: Piano

L = 16;
n = 701440;
t2 = linspace(0,L,n+1);
t = t2(1:n); % make n+1 points and throw out the last point
k = (2*pi/L)*[0:n/2-1 -n/2:-1];

t_p=16; % record time in seconds
y=audioread('music1.wav'); 
Fs=length(y)/t_p;
y = y.';
figure(5)
plot((1:length(y))/Fs,y);
xlabel('Time [sec]'); ylabel('Amplitude');
title('Mary had a little lamb (piano)'); drawnow

figure(6)
tslide = 0:0.1:16; 
Max = zeros(length(tslide),1);
s = zeros(length(tslide), length(t));
for j = 1:length(tslide)
    g = exp(-20*(t-tslide(j)).^2); %simple Gaussian
    Sf = g.*y; % signal filter
    Sft = fft(Sf);
    [f_max,ind] = max(Sft);
    Max(j) = k(ind);
    s(j,:) = abs(fftshift(Sft));
 
    subplot(3,1,1)
    plot(t,y,'k-',t,g,'r-')
    title('Signal and Filter')
    xlabel('time')
    ylabel('Amplitude')
    legend('signal','filter','Location','North')
    subplot(3,1,2)
    plot(t,Sf,'k-')
    xlabel('time')
    ylabel('Filtered Signal Amplitude')
    axis([0 16 -1 1])
    subplot(3,1,3)
    plot(fftshift(k),abs(fftshift(Sft))/max(abs(fftshift(Sft))),'b')
    title('FFT of Filtered Signal')
    xlabel('Frequency [Hz]')
    ylabel('n')
    axis([-100 100 0 0.1])
    pause(0.2)
end
figure(7)
plot(tslide,abs(Max./(2*pi)))
xlabel('Time [sec]'); ylabel('Hz');
title('Mary had a little lamb (piano)')
figure(8)
pcolor(tslide,fftshift(k),s.'),shading interp, colormap(hot)
axis([0 7 -10000 10000])
title('Spectrogram of Mary had a little lamb (piano)')

% Record
L = 14;
n = 627712;
t2 = linspace(0,L,n+1);
t = t2(1:n); % make n+1 points and throw out the last point
k = (2*pi/L)*[0:n/2-1 -n/2:-1];

figure(9)
tr_rec=14; % record time in seconds
y2=audioread('music2.wav'); Fs2=length(y2)/tr_rec;
y2 = y2.';
plot((1:length(y2))/Fs2,y2);
xlabel('Time [sec]'); ylabel('Amplitude');
title('Mary had a little lamb (recorder)');

tslide = 0:0.1:14; 
Max = zeros(length(tslide),1);
s = zeros(length(tslide), length(t));
for j = 1:length(tslide)
    g = exp(-20*(t-tslide(j)).^2); % simple Gaussian
    Sf2 = g.*y2; % signal filter
    Sft2 = fft(Sf2);
    [f_max2,ind] = max(Sft2);
    Max(j) = k(ind);
    s(j,:) = abs(fftshift(Sft2));
end
figure(10)
plot(tslide,abs(Max./(2*pi)))
xlabel('Time [sec]'); ylabel('Hz');
title('Mary had a little lamb (recorder)')
figure(11)
pcolor(tslide,fftshift(k),s.'),shading interp, colormap(hot)
axis([0 7 -10000 10000])
title('Spectrogram of Mary had a little lamb (recorder)')