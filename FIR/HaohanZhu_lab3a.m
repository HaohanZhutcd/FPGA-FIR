[x,Fs] = audioread('speech_20.wav'); % Read in audio file

% Plot to find out which frequency to remove from signal

% Transform length, specified as [] or a nonnegative integer scalar. 
% Specifying a positive integer scalar for the transform length can 
% increase the performance of fft. The length is typically specified 
% as a power of 2 or a value that can be factored into a product of 
% small prime numbers. If n is less than the length of the signal, 
% then fft ignores the remaining signal values past the nth entry and 
% returns the truncated result. If n is 0, then fft returns an empty matrix.
nfft = 2^10;
X = fft(x, nfft);
fstep = Fs/nfft; 
fvec = fstep*(0: nfft/2-1);
fresp = 2*abs(X(1:nfft/2));

%=====================================================
% fine frequency part work area begin
%=====================================================
% I is index which point to the maximum position 
[I] = argmax(fvec, fresp);
f_max = max(fresp);
% fstep is interval
% fvec are the fft samples duration
frequency_I = fstep * I;

disp(f_max);
disp(I);
disp(frequency_I);

% F_max = fresp(I)
% fmax = max(fresp)
% [fmax, I] = max(fresp);
%=====================================================
% fine frequency part work area over
%=====================================================

Hd = filtershow;
y = filter(Hd,x);
% sound(y,Fs);

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

audiowrite("clean.wav", y, Fs);

Hd_percision = filtershow1;
z = eval(filter(Hd,x));

figure(3);
Z = fft(z, nfft);
fresp_Z = 2*abs(Z(1:nfft/2));
plot(fvec,fresp_Z)
title(' Spectrum of Y_p(t)')
xlabel('Frequency (Hz)')
ylabel('|Y_p(f)|')
audiowrite("clean_P.wav", y, Fs);
% sound(Z,Fs);

%=====================================================
% fine frequency part work area begin
%=====================================================
function [index] = argmax(fvec, fresp)
    % initial index
    index = 0;
    % find maximum in the fresquency
    fmax = max(fresp);
    % when detection meet the maximum
    % there is the index of noise
    for i = 1 : length(fvec)
        f_index = fresp(i);
        if f_index == fmax
            index = i;
        end
    end
end
%=====================================================
% fine frequency part work area begin
%=====================================================

