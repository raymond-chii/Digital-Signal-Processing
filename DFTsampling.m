clc
clear
close all

%%%% Lei(Raymond) Chi
% a
hlap = (1/6) * [1 4 1; 4 -20 4; 1 4 1];
hx = (1/4) * [1 0 -1; 2 0 -2; 1 0 -1];
hy = (1/4) * [-1 -2 -1; 0 0 0; 1 2 1];
hlapsob = -filter2(hx, hx, 'full') - filter2(hy, hy, 'full')
% b

[Hlap,x_lap_frequency,y_lap_frequency] = freqz2(hlap);
[Hlapsob,~] = freqz2(hlapsob);

% c
figure;
colormap(jet);
contour(x_lap_frequency, y_lap_frequency, real(Hlap), 'Fill', 'on');
title('Contour Plot of H_{Lap}');
xlabel('Frequency');
ylabel('Frequency');

figure;
colormap(jet);
contour(x_lap_frequency, y_lap_frequency, real(Hlapsob), 'Fill', 'on');
title('Contour Plot of H_{LapSob}');
xlabel('Frequency');
ylabel('Frequency');

figure;
surf(x_lap_frequency, y_lap_frequency, real(Hlap), 'EdgeColor', 'none');
title('Surface Plot of H_{Lap}');
xlabel('Frequency');
ylabel('Frequency');
zlabel('Magnitude');

figure;
surf(x_lap_frequency, y_lap_frequency, abs(Hlapsob), 'EdgeColor', 'none');
title('Surface Plot of H_{LapSob}');
xlabel('Frequency');
ylabel('Frequency');
zlabel('Magnitude');

% Hlap => bandstop
% Hlapsob => Bandpass

% They are isotrpoic due to that both have rotational symmetry

% d

load("Rodanimg.mat");
load("LilyImg.mat");

lily = Lilyx;
rodan = Rodanx;

lap_filt_lily = filter2(hlap, lily);
lap_filt_rodan = filter2(hlap, rodan);

lapsob_filt_lily = filter2(hlapsob, lily);
lapsob_filt_rodan = filter2(hlapsob, rodan);

figure;
colormap('gray');
title('Original');
subplot(1,2,1);
image(lily);
subplot(1,2,2);
image(rodan);

figure;
colormap('gray');
title('H_{Lap}');
subplot(1,2,1);
image(lap_filt_lily);
subplot(1,2,2);
image(lap_filt_rodan);



figure;
colormap('gray');
title('H_{LapSob}');
subplot(1,2,1);
image(lapsob_filt_lily);
subplot(1,2,2);
image(lapsob_filt_rodan);

% b
figure;
subplot(1,2,1);
image(Lilyx);
title('Original Lily');
colormap('gray');
subplot(1,2,2);
image(rodan);
title('Original Rodan');
colormap('gray');

lilyUpsample = upsample(Lilyx);
rodanUpsample = upsample(Rodanx);

figure;
subplot(1,2,1);
image(lilyUpsample);
title('Upsampled Lily');
colormap('gray');
subplot(1,2,2);
image(rodanUpsample);
title('Upsampled Rodan');
colormap('gray');

% c

lily_fft = abs(fftshift(fft2(Lilyx)));
rodan_fft = abs(fftshift(fft2(Rodanx)));

lilyUpsample_fft = fftshift(fft2(lilyUpsample));
rodanUpsample_fft = fftshift(fft2(rodanUpsample));

figure;
image(abs(lily_fft));
colormap('gray');
title('DFT lily');


figure;
image(abs(lilyUpsample_fft));
colormap('gray');
title('DFT Upsampled lily');

figure;
image(abs(rodan_fft));
colormap('gray');
title('DFT rodan');

figure;
image(abs(rodanUpsample_fft));
colormap('gray');
title('DFT Upsampled rodan');

%Phase distortion


function M = upsample(m)
    n = length(m);
    M = zeros(2*n);
    M(1:2:end, 1:2:end) = m;
end
