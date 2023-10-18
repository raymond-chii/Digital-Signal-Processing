clc 
clear 
close all
%%%% Lei (Raymond) Chi ps02



f_pass1 = 9*10^6;
f_pass2 = 12.5*10^6;
f_stop1 = 9.5*10^6;
f_stop2 = 12*10^6;
Pv = 1.5;
Sa = 40; 
f_sample = 40 * 10^6; 

norm_fpass = [f_pass1, f_pass2] / (f_sample/2);
norm_fstop = [f_stop1, f_stop2] / (f_sample/2);


freq = linspace(0, 20*10^6, 1000); 
freqMega = freq/10^6;
w = 2*pi*freq;
w_pass = 2 * pi * norm_fpass * (f_sample/2);
w_stop = 2 * pi * norm_fstop * (f_sample/2);

% butter

[n_butter_dig, Wn_butter_dig] = buttord(norm_fpass, norm_fstop, Pv, Sa);
[b_butter_dig, a_butter_dig, ~] = butter(n_butter_dig, Wn_butter_dig, 'stop');
[h_butter_dig, ~] = freqz(poly(b_butter_dig), poly(a_butter_dig), freq, f_sample);
dB_butter_dig = 20 * log10(abs(h_butter_dig));
phase_butter_dig = unwrap(angle(h_butter_dig)) * 180 / pi;

[n_butter_ana, Wn_butter_ana] = buttord(w_pass, w_stop, Pv, Sa, "s");
[b_butter_ana, a_butter_ana, ~] = butter(n_butter_ana, Wn_butter_ana, 'stop', "s");
[h_butter_ana, ~]= freqs(poly(b_butter_ana), poly(a_butter_ana), w);
dB_butter_ana = 20 * log10(abs(h_butter_ana));
phase_butter_ana = unwrap(angle(h_butter_ana)) * 180 / pi;

% cheby1

[n_cheby1_dig, Wn_cheby1_dig] = cheb1ord(norm_fpass, norm_fstop, Pv, Sa);
[b_cheby1_dig, a_cheby1_dig, ~] = cheby1(n_cheby1_dig,Pv , Wn_cheby1_dig, 'stop');
[h_cheby1_dig, ~]= freqz(poly(b_cheby1_dig), poly(a_cheby1_dig), freq,f_sample);
dB_cheby1_dig = 20 * log10(abs(h_cheby1_dig));
phase_cheby1_dig = unwrap(angle(h_cheby1_dig)) * 180 / pi;

[n_cheby1_ana, Wn_cheby1_ana] = cheb1ord(w_pass, w_stop, Pv, Sa, "s");
[b_cheby1_ana, a_cheby1_ana, ~] = cheby1(n_cheby1_ana, Pv, Wn_cheby1_ana, 'stop', 's');
[h_cheby1_ana, ~]= freqs(poly(b_cheby1_ana), poly(a_cheby1_ana), w);
dB_cheby1_ana = 20 * log10(abs(h_cheby1_ana));
phase_cheby1_ana = unwrap(angle(h_cheby1_ana)) * 180 / pi;

% cheby2

[n_cheby2_dig, Wn_cheby2_dig] = cheb2ord(norm_fpass, norm_fstop, Pv, Sa);
[b_cheby2_dig, a_cheby2_dig, ~] = cheby2(n_cheby2_dig, Sa, Wn_cheby2_dig, 'stop');
[h_cheby2_dig, ~]= freqz(poly(b_cheby2_dig), poly(a_cheby2_dig), freq, f_sample);
dB_cheby2_dig = 20 * log10(abs(h_cheby2_dig));
phase_cheby2_dig = unwrap(angle(h_cheby2_dig)) * 180 / pi;

[n_cheby2_ana, Wn_cheby2_ana] = cheb2ord(norm_fpass, norm_fstop, Pv, Sa, 's');
[b_cheby2_ana, a_cheby2_ana, ~] = cheby2(n_cheby2_ana, Sa, Wn_cheby2_ana, 'stop', 's');
[h_cheby2_ana, ~]= freqs(poly(b_cheby2_ana), poly(a_cheby2_ana), w);
dB_cheby2_ana = 20 * log10(abs(h_cheby2_ana));
phase_cheby2_ana = unwrap(angle(h_cheby2_ana)) * 180 / pi;

% ellip

[n_ellip_dig, Wn_ellip_dig] = ellipord(norm_fpass, norm_fstop, Pv, Sa);
[b_ellip_dig, a_ellip_dig, ~] = ellip(n_ellip_dig, Pv, Sa, norm_fpass, 'stop');
[h_ellip_dig, ~]= freqz(poly(b_ellip_dig), poly(a_ellip_dig), freq, f_sample);
dB_ellip_dig = 20 * log10(abs(h_ellip_dig));
phase_ellip_dig = unwrap(angle(h_ellip_dig)) * 180 / pi;

[n_ellip_analog, Wn_ellip_analog] = ellipord(w_pass, w_stop, Pv, Sa, 's');
[b_ellip_ana, a_ellip_ana, ~] = ellip(n_ellip_analog, Pv, Sa, w_pass, 'stop', 's');
[h_ellip_ana, ~]= freqs(poly(b_ellip_ana), poly(a_ellip_ana), w);
dB_ellip_ana = 20 * log10(abs(h_ellip_ana));
phase_ellip_ana = unwrap(angle(h_ellip_ana)) * 180 / pi;

%part a

Digital_butterworth_filter_order = n_butter_dig * 2
Analog_butterworth_filter_order = n_butter_ana * 2

Digital_chebyshev1_filter_order = n_cheby1_dig * 2
Analog_chebyshev1_filter_order = n_cheby1_ana * 2

Digital_chebyshev2_filter_order = n_cheby2_dig * 2
Analog_chebyshev2_filter_order = n_cheby2_ana * 2

Digital_ellipitic_filter_order = n_ellip_dig * 2
Analog_ellipitic_filter_order = n_ellip_analog * 2

%part b

figure;
zplane(b_butter_dig, a_butter_dig);
title('Digital Butterworth pole-zero plot');
grid on;

figure;
zplane(b_butter_ana, a_butter_ana);
title('Analog Butterworth pole-zero plot');
grid on; 

figure;
zplane(b_cheby1_dig, a_cheby1_dig);
title('Digital Chebyshev 1 pole-zero plot');
grid on;

figure;
zplane(b_cheby1_ana, a_cheby1_ana);
title('Analog Chebyshev 1 pole-zero plot');
grid on; 

figure;
zplane(b_cheby2_dig, a_cheby2_dig);
title('Digital Chebyshev 2 pole-zero plot');
grid on;

figure;
zplane(b_cheby2_ana, a_cheby2_ana);
title('Analog Chebyshev 2 pole-zero plot');
grid on; 

figure;
zplane(b_ellip_dig, a_ellip_dig);
title('Digital Ellipitic pole-zero plot');
grid on;

figure;
zplane(b_ellip_ana, a_ellip_ana);
title('Analog Ellipitic pole-zero plot');
grid on; 

% part c
figure; 
subplot(2, 1, 1);
plot(freqMega, dB_butter_dig, 'b');
title('Digital Butterworth Magnitude Response');
xlabel('Frequency (MHz)');
ylabel('Magnitude (dB)');
ylim([-50, 2]);

subplot(2, 1, 2);
plot(freqMega, phase_butter_dig, 'r');
title('Phase Response');
xlabel('Frequency (MHz)');
ylabel('Phase (degrees)');

figure;
subplot(2, 1, 1);
plot(freqMega, dB_butter_ana, 'b');
title('Analog Butterworth Magnitude Response');
xlabel('Frequency (MHz)');
ylabel('Magnitude (dB)');
ylim([-50, 2]);

subplot(2, 1, 2);
plot(freqMega, phase_butter_ana, 'r');
title('Phase Response');
xlabel('Frequency (MHz)');
ylabel('Phase (degrees)');

figure;
subplot(2, 1, 1);
plot(freqMega, dB_cheby1_dig, 'b');
title('Digital Cheby1 Magnitude Response');
xlabel('Frequency (MHz)');
ylabel('Magnitude (dB)');
ylim([-50, 2]);

subplot(2, 1, 2);
plot(freqMega, phase_cheby1_dig, 'r');
title('Phase Response');
xlabel('Frequency (MHz)');
ylabel('Phase (degrees)');

figure;
subplot(2, 1, 1);
plot(freqMega, dB_cheby1_ana, 'b');
title('Analog Cheby1 Magnitude Response');
xlabel('Frequency (MHz)');
ylabel('Magnitude (dB)');
ylim([-50, 2]);

subplot(2, 1, 2);
plot(freqMega, phase_cheby1_ana, 'r');
title('Phase Response');
xlabel('Frequency (MHz)');
ylabel('Phase (degrees)');

figure;
subplot(2, 1, 1);
plot(freqMega, dB_cheby2_dig, 'b');
title('Digital cheby2 Magnitude Response');
xlabel('Frequency (MHz)');
ylabel('Magnitude (dB)');
ylim([-50, 2]);

subplot(2, 1, 2);
plot(freqMega, phase_cheby2_dig, 'r');
title('Phase Response');
xlabel('Frequency (MHz)');
ylabel('Phase (degrees)');

figure;
subplot(2, 1, 1);
plot(freqMega, dB_cheby2_ana, 'b');
title('Analog cheby2 Magnitude Response');
xlabel('Frequency (MHz)');
ylabel('Magnitude (dB)');
ylim([-50, 2]);

subplot(2, 1, 2);
plot(freqMega, phase_cheby2_ana, 'r');
title('Phase Response');
xlabel('Frequency (MHz)');
ylabel('Phase (degrees)');

figure;
subplot(2, 1, 1);
plot(freqMega, dB_ellip_dig, 'b');
title('Digital ellip Magnitude Response');
xlabel('Frequency (MHz)');
ylabel('Magnitude (dB)');
ylim([-50, 2]);

subplot(2, 1, 2);
plot(freqMega, phase_ellip_dig, 'r');
title('Phase Response');
xlabel('Frequency (MHz)');
ylabel('Phase (degrees)');

figure;
subplot(2, 1, 1);
plot(freqMega, dB_ellip_ana, 'b');
title('Analog ellip Magnitude Response');
xlabel('Frequency (MHz)');
ylabel('Magnitude (dB)');
ylim([-50, 2]);

subplot(2, 1, 2);
plot(freqMega, phase_ellip_ana, 'r');
title('Phase Response');
xlabel('Frequency (MHz)');
ylabel('Phase (degrees)');

% part d


gain_at_fpass1 = f_pass1/ ((20 * 10^6)/1000);
gain_at_fpass2 = f_pass2/ ((20 * 10^6)/1000);
gain_at_fstop1 = f_stop1/ ((20 * 10^6)/1000);
gain_at_fstop2 = f_stop2/ ((20 * 10^6)/1000);

dig_butter_gain_at_9MHz = dB_butter_dig(gain_at_fpass1)
dig_butter_gain_at_9_5MHz = dB_butter_dig(gain_at_fstop1)
dig_butter_gain_at_12MHz = dB_butter_dig(gain_at_fstop2)
dig_butter_gain_at_12_5MHz = dB_butter_dig(gain_at_fpass2)

ana_butter_gain_at_9MHz = dB_butter_ana(gain_at_fpass1)
ana_butter_gain_at_9_5MHz = dB_butter_ana(gain_at_fstop1)
ana_butter_gain_at_12MHz = dB_butter_ana(gain_at_fstop2)
ana_butter_gain_at_12_5MHz = dB_butter_ana(gain_at_fpass2)

dig_cheby1_gain_at_9MHz = dB_cheby1_dig(gain_at_fpass1)
dig_cheby1_gain_at_9_5MHz = dB_cheby1_dig(gain_at_fstop1)
dig_cheby1_gain_at_12MHz = dB_cheby1_dig(gain_at_fstop2)
dig_cheby1_gain_at_12_5MHz = dB_cheby1_dig(gain_at_fpass2)

ana_cheby1_gain_at_9MHz = dB_cheby1_ana(gain_at_fpass1)
ana_cheby1_gain_at_9_5MHz = dB_cheby1_ana(gain_at_fstop1)
ana_cheby1_gain_at_12MHz = dB_cheby1_ana(gain_at_fstop2)
ana_cheby1_gain_at_12_5MHz = dB_cheby1_ana(gain_at_fpass2)

dig_cheby2_gain_at_9MHz = dB_cheby2_dig(gain_at_fpass1)
dig_cheby2_gain_at_9_5MHz = dB_cheby2_dig(gain_at_fstop1)
dig_cheby2_gain_at_12MHz = dB_cheby2_dig(gain_at_fstop2)
dig_cheby2_gain_at_12_5MHz = dB_cheby2_dig(gain_at_fpass2)

ana_cheby2_gain_at_9MHz = dB_cheby2_ana(gain_at_fpass1)
ana_cheby2_gain_at_9_5MHz = dB_cheby2_ana(gain_at_fstop1)
ana_cheby2_gain_at_12MHz = dB_cheby2_ana(gain_at_fstop2)
ana_cheby2_gain_at_12_5MHz = dB_cheby2_ana(gain_at_fpass2)

dig_ellip_gain_at_9MHz = dB_ellip_dig(gain_at_fpass1)
dig_ellip_gain_at_9_5MHz = dB_ellip_dig(gain_at_fstop1)
dig_ellip_gain_at_12MHz = dB_ellip_dig(gain_at_fstop2)
dig_ellip_gain_at_12_5MHz = dB_ellip_dig(gain_at_fpass2)

ana_ellip_gain_at_9MHz = dB_ellip_ana(gain_at_fpass1)
ana_ellip_gain_at_9_5MHz = dB_ellip_ana(gain_at_fstop1)
ana_ellip_gain_at_12MHz = dB_ellip_ana(gain_at_fstop2)
ana_ellip_gain_at_12_5MHz = dB_ellip_ana(gain_at_fpass2)

% part a 
f_cuts = [f_pass1 f_stop1 f_stop2 f_pass2];

dev_pass = (10^(Pv/20)-1)/(10^(Pv/20)+1);

dev_stop = 10 ^ (-Sa/20); 
devs = [dev_pass dev_stop dev_pass];
magn = [1 0 1]; 

[n, Wn, beta, ftype] = kaiserord(f_cuts, magn, devs, f_sample);

[n_eq, fo, ao, w] = firpmord(f_cuts, magn, devs, f_sample);

% part b
n = n + rem(n,2); 
kaiser_filter_length = n + 1
b = fir1(n,Wn,ftype, kaiser(kaiser_filter_length,beta));

a = firpm(n_eq,fo,ao,w);
Equirpple_filter_length = n_eq + 1 


figure; 
zplane(b, 1);
title('Kaiser window design')
grid on; 

figure;
zplane(a, 1);
title('Equiripple FIR filter design')
grid on; 


[h,f] = freqz(b,1,1024,f_sample); 
figure; 
plot(f, 20*log10(abs(h)));
title('Kaiser window design')
grid on; 

[h_eq, f_eq] =freqz(a, 1, 1024, f_sample);
figure;
plot(f_eq, 20*log10(abs(h_eq)));
title('Equiripple FIR filter design')
grid on;

figure;
stem(b);
title('Stem plot for the Kaiser filter coefficients');
xlim([0 n]);
grid on;

figure;
stem(a);
title('Stem plot of Equiripple FIR Filter coefficients');
xlim([0 n_eq]);
grid on;

% part c


% w(1)
% w(2)
weight_ratio = w(1)/w(2)
deviations_ratio = dev_stop/dev_pass

difference_in_ratio = weight_ratio - deviations_ratio % basically eqaul
% part d

[k_pb_var, k_ps_gain] = max_min_gain(b, f_cuts, f_sample);

Kaiser_passband_variation_dB = k_pb_var
Kaiser_peakstopband_gain_dB = k_ps_gain

[eq_pb_var, eq_ps_gain] = max_min_gain(a, f_cuts, f_sample);

equiripple_passband_variation_dB = k_pb_var
equiripple_peakstopband_gain_dB = k_ps_gain

function [i, j] = max_min_gain(c, f_cuts, f_sample)
    
    f_pass1 = f_cuts(1);
    f_pass2 = f_cuts(4);
    f_stop1 = f_cuts(2);
    f_stop2 = f_cuts(3);
    freq1 = linspace(0, f_pass1, 1000);
    freq2 = linspace(f_pass2, 20*10^6, 1000);
    freq3 = linspace(f_stop1, f_stop2, 1000);

    p_gain = 0;
    min_gain = inf;

    for k = freq1
        z = exp(2 * pi * k * 1i / f_sample);
        h = polyval(c, z);
        acc_gain = abs(h);
        
        if acc_gain < min_gain
            min_gain = acc_gain;
        end
        if acc_gain > p_gain
            p_gain = acc_gain;
        end
    end

    for k = freq2
        z = exp(2 * pi * k * 1i / f_sample); 
        h2 = polyval(c, z);
        acc_gain = abs(h2);

        if acc_gain < min_gain
            min_gain = acc_gain;
        end
        if acc_gain > p_gain
            p_gain = acc_gain;
        end
        
    end
    i = 20*log10(p_gain)-20*log10(min_gain);

    p_gain = 0;
    for k = freq3
        z = exp(2 * pi * k * 1i / f_sample); 
        h3 = polyval(c, z);
        acc_gain = abs(h3);

        if acc_gain > p_gain
            p_gain = acc_gain;
        end
    end

    j = 20 * log10(p_gain);

end



