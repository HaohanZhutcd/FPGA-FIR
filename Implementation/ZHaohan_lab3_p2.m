clear all;
close all;


filename = 'audio.wav';
% filename = 'speech_20.wav';
[x,Fs] = audioread(filename); % Read in audio file

nfft = 2^10;
X = fft(x, nfft);
fstep = Fs/nfft; 
fvec = fstep*(0: nfft/2-1);
fresp = 2*abs(X(1:nfft/2));

%=====================================================
% fine frequency part work area begin
%=====================================================
% I is index which point to the maximum position 
% [I] = argmax(fvec, fresp);
% f_max = max(fresp);
% fstep is interval
% fvec are the fft samples duration
% frequency_I = fstep * I;

% disp(f_max);
% disp(I);
% disp(frequency_I);

%=====================================================
% fine frequency part work area over
%=====================================================

% Hd = FIR_Hp;
% Hd = FIR_Speech;
Hd = FIR_Audio;
y = double(filter(Hd,x));
% sound(y,Fs);
Hd.Arithmetic = 'fixed';
coewrite(Hd);
figure(1);
plot(fvec,fresp)
title('Single-Sided Amplitude Spectrum of x(t)')
xlabel('Frequency (Hz)')
ylabel('|X(f)|')


figure(2);
Y = fft(y, nfft);
fresp_Y = 2*abs(Y(1:nfft/2));
plot(fvec,fresp_Y)
title(' Spectrum of Y(t)')
xlabel('Frequency (Hz)')
ylabel('|Y(f)|')
grid on;
% sound(y, Fs)
audiowrite("clean.wav", y, Fs);

% %=====================================================
% % fine frequency part work area begin
% %=====================================================
% function [index] = argmax(fvec, fresp)
%     % initial index
%     index = 0;
%     % find maximum in the fresquency
%     fmax = max(fresp);
%     % when detection meet the maximum
%     % there is the index of noise
%     for i = 1 : length(fvec)
%         f_index = fresp(i);
%         if f_index == fmax
%             index = i;
%         end
%     end
% end
% %=====================================================
% % fine frequency part work area begin
% %=====================================================

