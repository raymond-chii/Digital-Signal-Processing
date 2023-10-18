clc
clear
close all

%%%% Lei(Raymond) Chi DSP ps01

%%question 2

Fs = 20 * 10^6;
N = 512; 
collected_samples = 500;
t = linspace(0, (collected_samples-1)/Fs, collected_samples); 
x = sin(2*pi*(6*10^6)*t);
x_padded = [x, zeros(1, N - collected_samples)];
dft = fft(x_padded, N);

mag = abs(dft)
k = 0:N-1;

% a 
bin_freq = k*Fs/N; 
bin_spacing = Fs/N; 

% b 
figure;  
plot(k,mag);
title("mag of dft vs. k indices");
hold on;
grid on; 
% read on the graph 155 and 359


% c

rectangle = rectwin(collected_samples);
scale_rectangle = rectangle / sqrt(collected_samples); 
x_rectangle = x .* scale_rectangle'; 
dft2 = fft(x_rectangle, N);
fe = abs(dft2);

Cheb = chebwin(collected_samples, 30);
ener_Cheb = sum (Cheb.^ 2); 
scale_Cheb = Cheb / sqrt(ener_Cheb); 
x_Cheb = x .* scale_Cheb'; 
dft3 = fft(x_Cheb,N);
ef = abs(dft3);

figure;
subplot(2,1,1); 
plot(bin_freq, fe); 
title("rectangle window of 512 points dft");
xlabel("Frequency(Hz)");
ylabel("Mag of dft"); 
grid on; 
subplot(2,1,2); 
plot(bin_freq, ef);
title("Chebyshev window"); 
xlabel("Frequency(Hz)");
ylabel("Mag of dft"); 
grid on; 

 

% d

Ps = rms(x)^2; 
Pn = Ps/100; 
noise_power = sqrt(Pn) * randn(1, 500);
noise_x = noise_power + x;
rect_x = fft(noise_x/2);
wind_x = noise_x' .* Cheb; 
X_Cheb = fft(wind_x);

x_rect_shift = fftshift(20*log10(abs(rect_x)))
x_Cheb_shift = fftshift(20*log10(abs(x_Cheb)));

k_in = (-collected_samples/2):(collected_samples/2-1)

figure; 
subplot(2,1,1);
plot(k_in, x_rect_shift);
title("rect window of dft");
xlabel("Frequency(Hz)");
ylabel("Mag of dft"); 
grid on; 
subplot(2,1,2)
plot(k_in, x_Cheb_shift);
title("Chebyshev window of dft");
xlabel("Frequency(Hz)");
ylabel("Mag of dft"); 
grid on; 

figure; 
plot(k_in, x_rect_shift);
hold;
plot(k_in, x_Cheb_shift);
title("superimposed graph")
xlabel("Frequency(Hz)");
ylabel("Mag of dft"); 
legend("Rect Wind","Chebyshev Wind")
grid on; 
xlim([145,165])


%%question 3

% a
Fs = 44.1*10^3; 
N0 = round(Fs/2); 
G4 =393; 
A4 = 440; 
D5 = 587.33;
S = N0 *50 - 1; 
t = (0:S)/Fs; 
x_1 = cos(2*pi*t*G4); 
x_2 = cos(2*pi*t*A4);
x_3 = cos(2*pi*t*D5);
x_M = [x_1;x_2;x_3]; 
bitmask = [0 1 1; 1 0 1; 1 1 0];
iBM = zeros(3,S+1);
r = randi([1,3], 1, S+1);


for i = 1: S
    iBM(:, i ) = bitmask(:,r(i)); 



end


iBM = iBM .* x_M; 
clean_n = sum(iBM); 
Ps = rms(clean_n) ^2; 
Pn = Ps/10^4; 
noise_power = rand(1, S+1) * sqrt(Pn) ;
q_notes = noise_power + clean_n


% b


wind = hamming(N0);
n = 2^nextpow2(N0);
m = N0/2
[periodogram, g] = pwelch(q_notes, wind, m, N, Fs)


% c

figure; 
subplot(2,1,1)
plot(g, 10*log10(periodogram));
title("periodogram(decibel)");
xlabel("Frequency(Hz)");
ylabel("power/Frequency(dB)"); 
grid on; 
subplot(2,1,2)
plot(g, 10*log10(periodogram));
xlim([300,600])
title("periodogram(decibel)(zoomed in)");
xlabel("Frequency(Hz)");
ylabel("power/Frequency(dB)"); 
grid on; 


% d

[X, Y, Z] = spectrogram(q_notes, wind, m, n, Fs);
Xdb = 10*log10(abs(X));
[time, frequency] = meshgrid(Z, Y);
figure; 
surf(frequency , time, Xdb);
colormap winter;

