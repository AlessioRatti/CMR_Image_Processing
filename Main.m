%% Biomedical Image Processing Assignment
clear all
clc

%% Task 1
% Load and display the original CMR image.

% load
cmr_str = load('CMR_slice.mat');
cmr = cmr_str.I_noise;
fs = 1/0.7; % [mm^(-1)]
fn = fs/2*1000; % [m^(-1)]

% plot
ax_sd = imref2d(size(cmr)); % defining the axes object for spatial domain
ax_sd.XWorldLimits = [0 (size(cmr,1))/fs]; % [mm]
ax_sd.YWorldLimits = [0 (size(cmr,2))/fs]; % [mm]
figure(1)
imshow(cmr,ax_sd,[]); title('Imported CMR'); xlabel('x [mm]');ylabel('y [mm]'); set(gca, 'XAxisLocation', 'top');

%% Task 2
% Illustrate magnitude and angle of 2D Discrete Fourier Transform (DFT) of
% the image in a unique figure and using jet colormap. Shift the
% zero-frequency component to the center of the output. When plotting 2D
% DFT magnitude and phase in the spatial frequency domain, axes must always
% be expressed in spatial frequency units from now on. Specify axis marks
% and label.

% rescaling the image in the range [0 1]
%cmr_res = (cmr - min(cmr(:))) / (max(cmr(:) - min(cmr(:)))); % maximum in the whole matrix
cmr_res = mat2gray(cmr, [min(cmr(:)) max(cmr(:))]);

% frequency domain
cmr_fft = fft2(cmr_res);    % 2D fast fourier transform

% magnitude
cmr_mag = abs(cmr_fft);
%cmr_mag_res = (cmr_mag - min(min(cmr_mag))) / (max(cmr_mag(:))-min(cmr_mag(:))); % new range [0 1]
%cmr_mag_res = rescale(10*log10(cmr_mag));
cmr_mag_res = 10*log10(cmr_mag);
% without rescaling it, the magnitude plot is a uniform dark red square

% phase
% returns the phase angle in the interval [-pi,pi]
cmr_phase = angle(cmr_fft); % the original phase range is [-pi pi]
% cmr_phase_res = (cmr_phase - min(min(cmr_phase))) / (max(cmr_phase(:)) - min(cmr_phase(:))); % new range [0 1]
% rescale(cmr_phase); % rescale between [0 1]

% plotting
ax_fd = imref2d([size(cmr,1) size(cmr,2)]); % defining the axes object for sp. frequency domain
ax_fd.XWorldLimits = [-fn fn];
ax_fd.YWorldLimits = [-fn fn];

figure(2)

subplot(1,2,1);
imshow(fftshift(cmr_mag_res),ax_fd,[],'Colormap',jet);
title('CMR - 2D FFT magnitude'); xlabel('Spatial frequency [cycles/m]'); ylabel('Spatial frequency [cycles/m]');
c = colorbar;
c.Label.String = 'Magnitude [dB]';

subplot(1,2,2);
imshow(fftshift(cmr_phase),ax_fd,[-pi pi],'Colormap',jet); 
title('CMR - 2D FFT phase'); xlabel('Spatial frequency [cycles/m]'); ylabel('Spatial frequency [cycles/m]');
c = colorbar;
c.Ticks = [-pi 0 pi];
c.TickLabels ={'-\pi ','0','+\pi '};
c.Label.String = 'Phase [rad]';

%% Task 3
% The image is corrupted by noise, characterize the noise and provide a
% detailed description of its properties in both spatial and spatial
% frequency domains.

% The image is affected by the superimposition of a periodic noise (from 
% Image Lab #3: 'repetitive intensity pattern added to the original image 
% and spectral peaks in the frequency domain').
%
% Strong sinusoidal noise: ~387 cycles/m in the x-axis direction
% ~ 276 cycles/m in the y-axis direction
% in fact we can actually sse how the angle of the sinusoid with respect to
% the x-axis is slightly greater than 45 degrees, exaplining how the
% frequency on the x-axis is higher than on the y-axis
%
% The two peaks are at [140 126] and [42 56], the center of the DFT is [91
% 91]. Values of the peaks are in samples (convert to spatial frequencies.

%% Task 4
% Reduce the noise through denoising operations in the frequency domain
% (correction step #1). Show the CMR image before and after correction #1
% and the removed noise in both spatial and frequency domain (2D DFT
% magnitude only). Display all images in a unique figure. Axes must be
% expressed in spatial or spatial frequency units. Describe the denoising
% approach used.

% FILTERING with notch
% find the location of the noise spikes.
threshold = 30;
spikes = abs(fftshift(cmr_fft)) > threshold; % binary image

% exclude the central DC spike
spikes(60:120,50:135) = 0;

% plotting
figure(3)
imshow(spikes,ax_fd,[]), title('Notch filter "peaks"');
xlabel('Spatial frequency [cycles/m]'); ylabel('Spatial frequency [cycles/m]');
% check if the central peak has been removed

% filter/mask the spectrum
cmr_filt_fft = fftshift(cmr_fft); % saving the original fft before applying the notch
cmr_filt_fft(spikes) = 0;
cmr_filt = ifft2(ifftshift(cmr_filt_fft));

% extract the noise as mathematical difference
noise = cmr_res - cmr_filt;
noise_fft = fft2(noise);

% plot
figure(4)
subplot(2,3,1);
imshow(cmr_res,ax_sd,[]), title('Original CMR'), xlabel('x [mm]'), ylabel('y [mm]');
set(gca, 'XAxisLocation', 'top');
subplot(2,3,2);
imshow(cmr_filt,ax_sd,[]), title(strcat('CMR filtered with 2D notch filter')), xlabel('x [mm]'), ylabel('y [mm]');
set(gca, 'XAxisLocation', 'top');
subplot(2,3,3);
imshow(noise,ax_sd,[]), title('Removed noise'), xlabel('x [mm]'), ylabel('y [mm]');
set(gca, 'XAxisLocation', 'top');
subplot(2,3,4);
imshow(fftshift(cmr_mag_res),ax_fd,[min(cmr_mag_res(:)) max(cmr_mag_res(:))],'Colormap',jet);
title('CMR - 2D FFT mag'); xlabel('Spatial frequency [cycles/m]'); ylabel('Spatial frequency [cycles/m]');
c = colorbar;
c.Label.String = 'Magnitude [dB]';
subplot(2,3,5);
imshow(10*log10(abs(cmr_filt_fft)),ax_fd,[min(cmr_mag_res(:)) max(cmr_mag_res(:))],'Colormap',jet);
title('CMR denoised - 2D FFT mag'); xlabel('Spatial frequency [cycles/m]'); ylabel('Spatial frequency [cycles/m]');
c = colorbar;
c.Label.String = 'Magnitude [dB]';
subplot(2,3,6);
imshow(10*log10(fftshift(abs(noise_fft))),ax_fd,[min(cmr_mag_res(:)) max(cmr_mag_res(:))],'Colormap',jet);
title('Removed noise - 2D FFT mag'); xlabel('Spatial frequency [cycles/m]'); ylabel('Spatial frequency [cycles/m]');
c = colorbar;
c.Label.String = 'Magnitude [dB]';

%% Task 5
% Reduce any residual noise through denoising operations in the spatial 
% domain (correction step #2). In the same figure, show the CMR image 
% before and after correction #2 and the removed noise. Show the CMR image 
% before and after correction #2 and the removed noise in both spatial and 
% spatial frequency domain (2D DFT magnitude only). Display all images in a
% unique figure. Axes must always be expressed in spatial or spatial 
% frequency units.

% gaussian filter
sigma = 1;
cmr_filt_h2 = imgaussfilt(cmr_filt, sigma);

% frequency domain
cmr_filt_h2_fft = fft2(cmr_filt_h2); 
cmr_filt_h2_fft_mag = abs(cmr_filt_h2_fft);

% extract the noise as mathematical difference
noise = cmr_filt - cmr_filt_h2;
noise_fft_mag = abs(fft2(noise));

% plot
figure(6)
subplot(2,3,1);
imshow(cmr_filt,ax_sd,[]), title('CMR filtered in frequency only'), xlabel('x [mm]'), ylabel('y [mm]');
set(gca, 'XAxisLocation', 'top');

subplot(2,3,2);
imshow(cmr_filt_h2,ax_sd,[]), title('CMR w/ notch and Gaussian filtering'), xlabel('x [mm]'), ylabel('y [mm]');
set(gca, 'XAxisLocation', 'top');

subplot(2,3,3);
imshow(noise,ax_sd,[]), title('Removed noise after step #2'), xlabel('x [mm]'), ylabel('y [mm]');
set(gca, 'XAxisLocation', 'top');

subplot(2,3,4);
imshow(10*log10(abs(cmr_filt_fft)),ax_fd,[min(cmr_mag_res(:)) max(cmr_mag_res(:))],'Colormap',jet);
title('After step #1 - 2D FFT mag'); xlabel('Spatial frequency [cycles/m]'); ylabel('Spatial frequency [cycles/m]');
c = colorbar;
c.Label.String = 'Magnitude [dB]';

subplot(2,3,5);
imshow(10*log10(fftshift(cmr_filt_h2_fft_mag)),ax_fd,[min(cmr_mag_res(:)) max(cmr_mag_res(:))],'Colormap',jet);
title('After step #1 and #2 - 2D FFT mag'); xlabel('Spatial frequency [cycles/m]'); ylabel('Spatial frequency [cycles/m]');
c = colorbar;
c.Label.String = 'Magnitude [dB]';

subplot(2,3,6);
imshow(10*log10(fftshift(noise_fft_mag)),ax_fd,[min(cmr_mag_res(:)) max(cmr_mag_res(:))],'Colormap',jet);
title('Noise only - 2D FFT mag'); xlabel('Spatial frequency [cycles/m]'); ylabel('Spatial frequency [cycles/m]');
c = colorbar;
c.Label.String = 'Magnitude [dB]';

%% Task 6
% Maximize the image contrast by processing the histogram of the denoised
% image (after correction step #2) in order to maximize the difference
% between myocardium and blood. In the same figure, show the CMR image
% before and after the operation and the effect of the contrast enhancement
% and their corresponding histograms.

% Moreover, describe the type of histogram transformation you performed and
% calculate the CNR (expressed in the appropriate intensity units) before
% and after the transformation.

nbins = 64; % Matlab default value = 64
% when visualizing it, bins number does not seem to affect the result

cmr_eq = histeq(cmr_filt_h2,nbins);
equalization = cmr_filt_h2 - cmr_eq; % effect of the contrast enhancement

side = 5; % pixels

xc1 = 100;
yc1 = 85;

xc2 = 110;
yc2 = 107;

x1 = (xc1) - (side-1)/2; % pixels
y1 = (yc1) - (side-1)/2; % pixels
x2 = (xc2) + (side-1)/2; % pixels
y2 = (yc2) + (side-1)/2; % pixels

figure(7)
subplot(2,3,1); imshow(cmr_filt,ax_sd,[]);
title('Filterd image'); xlabel('x [mm]');ylabel('y [mm]');set(gca, 'XAxisLocation', 'top');
hold on
rectangle('Position',[x1/fs y1/fs side/fs side/fs], 'EdgeColor','r');
hold on
rectangle('Position',[x2/fs y2/fs side/fs side/fs], 'EdgeColor','g');
subplot(2,3,2); imshow(cmr_eq,ax_sd,[]);
hold on
rectangle('Position',[x1/fs y1/fs side/fs side/fs], 'EdgeColor','r');
hold on
rectangle('Position',[x2/fs y2/fs side/fs side/fs], 'EdgeColor','g');
title('Equalized image'); xlabel('x [mm]');ylabel('y [mm]');set(gca, 'XAxisLocation', 'top');
subplot(2,3,3); imshow(equalization,ax_sd,[]);
title('Effect of contrast enhancement'); xlabel('x [mm]');ylabel('y [mm]');set(gca, 'XAxisLocation', 'top');
subplot(2,3,4); imhist(cmr_filt,nbins);
xlabel({'','','Intensity'}); ylabel('# pixels'); title('Filterd image histogram');
subplot(2,3,5); imhist(cmr_eq,nbins);
xlabel({'','','Intensity'}); ylabel('# pixels'); title('Equalized image histogram');
subplot(2,3,6); imhist(equalization,nbins);
xlabel({'','','Intensity'}); ylabel('# pixels'); title('Contrast enhancement histogram');

% calculating the CNR value
sub_blood_filt = imcrop(cmr_filt_h2,[x1 y1 side side]);
sub_myo_filt = imcrop(cmr_filt_h2,[x2 y2 side side]);

sub_blood_eq = imcrop(cmr_eq,[x1 y1 side side]);
sub_myo_eq = imcrop(cmr_eq,[x2 y2 side side]);

mu_blood_filt = mean(sub_blood_filt, 'all')
mu_myo_filt = mean(sub_myo_filt, 'all')
sigma_myo_filt = std(sub_myo_filt, [], 'all')

mu_blood_eq = mean(sub_blood_eq, 'all')
mu_myo_eq = mean(sub_myo_eq, 'all')
sigma_myo_eq = std(sub_myo_eq, [], 'all')

CNR_filt = (mu_blood_filt - mu_myo_filt)/sigma_myo_filt
CNR_eq = (mu_blood_eq - mu_myo_eq)/sigma_myo_eq

%% Task 7
% Select a subregion in the restored image (61 x 61 matrix) including the
% left ventricle section. The coordinates of region of interest center are:
% xc = 74.2 mm, yc = 65.1 mm. Display the cropped image in one figure and
% its boundaries in the original image in a second figure.

xc = 74.2; % [mm]
yc = 65.1; % [mm]

side = 61; % pixels

x1 = (xc*fs) - (side-1)/2; % pixels
y1 = (yc*fs) - (side-1)/2; % pixels
x2 = (xc*fs) + (side-1)/2; % pixels
y2 = (yc*fs) + (side-1)/2; % pixels

subregion = imcrop(cmr_eq,[x1 y1 side side]);

ax_sd_c = imref2d(size(subregion)); % defining the axes object for spatial domain
ax_sd_c.XWorldLimits = [x1/fs-1/(2*fs) x2/fs+1/(2*fs)]; % [mm]
ax_sd_c.YWorldLimits = [y1/fs-1/(2*fs) y2/fs+1/(2*fs)]; % [mm]

figure(8)
subplot(1,2,1); imshow(subregion,ax_sd_c,[]);
title('Cropped region'); xlabel('x [mm]');ylabel('y [mm]');set(gca, 'XAxisLocation', 'top');
subplot(1,2,2); imshow(cmr_eq,ax_sd,[]);
hold on
rectangle('Position',[x1/fs-1/(2*fs) y1/fs-1/(2*fs) side/fs side/fs], 'EdgeColor','r'); % [mm]
title('Cropped region boundaries'); xlabel('x [mm]');ylabel('y [mm]');set(gca, 'XAxisLocation', 'top');

%% Task 8
% Using automatic techniques, segment the blood in left ventricle section
% and display the obtained binary image. Provide a detailed description of
% the used technique. Calculate the blood section area in mm2.

%maximum image value
I_max = max(subregion(:))

%binarization
Th = I_max/2;
Th = graythresh(subregion);
subregion_bin = imbinarize(subregion,Th);

stats = regionprops('table',subregion_bin,'Centroid','MajorAxisLength','MinorAxisLength', 'Area');
centers = stats.Centroid;

% adding red color to the blood
L = bwlabel(subregion_bin);
L = (L==3);
subregion_blood = subregion_bin+L+1;
map = [0 0 0; 1 1 1; 1 0 0];

figure(9)
imshow(subregion_blood, ax_sd_c, map)
hold on
plot(centers(:,1)/fs+ax_sd_c.XWorldLimits(1)-1/(fs*2),centers(:,2)/fs+ax_sd_c.YWorldLimits(1)-1/(fs*2),'b*')
hold off
title('Image segmentation');xlabel('x [mm]');ylabel('y [mm]');set(gca, 'XAxisLocation', 'top');

area = sum(subregion_blood == 3, 'all')*(1/fs)^2

%% Task 8
% Using automatic techniques, identify the edge of the segmented region and
% show it superimposed to the other structures using the red color. Provide
% a detailed description of the used edge detection approach

Th = 0.007;
edges = edge(subregion,'log', Th);

map = (0:255)'./255;
map = [map map map];
map = [map; 1 0 0];

subregion_edges = subregion;
subregion_edges(subregion_edges == 1) = 0.999;
subregion_edges = subregion_edges+edges;

figure(10)
imshow(subregion_edges*length(map), ax_sd_c, map)
title('Detected edges with LoG operator');xlabel('x [mm]');ylabel('y [mm]');set(gca, 'XAxisLocation', 'top');
