%% Anna Fellin
% Using RAFEEF GARBI - ELEC 421 - DIGITAL SIGNAL PROCESSING - 2018 Matlab
% as base
% Hough/Line detection: https://www.mathworks.com/help/images/hough-transform.html
close all; clear all; clc;

addpath 'C:\Users\Anna\Documents\UBC Classes\Personal Projects\Motion Photos';

%% Mathematical Model of the Motion Blur
mn=0; st=0.0;% 0.0 0.001 0.1  %noise level

USE_PSEUDO=0; % use psuedo filter
USE_WIENER=1; % use wiener filter
NSR = 0.1;
    
%parameter T (observation time) and motion rate
%guess for now: todo later, adjust

%reading Image
figure; imagesc(imread('blur.jpg')); axis image;
RGB=im2double(imread('blur.jpg'));
I = rgb2gray(RGB);

% Use Hough Function to detect lines
BW = edge(I);
[H,theta,rho] = hough(BW);
P = houghpeaks(H,5,'threshold',ceil(0.3*max(H(:))));
lines = houghlines(BW,theta,rho,P,'FillGap',15,'MinLength',7);

figure, imshow(I), hold on
max_len = 0;
for k = 1:length(lines)
   xy = [lines(k).point1; lines(k).point2];
   plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');

   % Plot beginnings and ends of lines
   plot(xy(1,1),xy(1,2),'x','LineWidth',2,'Color','yellow');
   plot(xy(2,1),xy(2,2),'x','LineWidth',2,'Color','red');

   % Determine the endpoints of the longest line segment
   len = norm(lines(k).point1 - lines(k).point2);
   if ( len > max_len)
      max_len = len;
      xy_long = xy;
   end
end

% take longest line
dX = xy_long(2,1) - xy_long(1,1);
dY = xy_long(2,2) - xy_long(1,2);

T=0.5;
ax = -dX; ay = dY;

%%
%generating frequencies for the blurring model
u=linspace(-0.5,0.5,size(I,2));
v=linspace(-0.5,0.5,size(I,1));

[U,V]=meshgrid(u,v);

H=(T./(pi*(U*ax + V*ay))).*sin(pi*(U*ax + V*ay)).*exp(-1i*pi*(U*ax + V*ay)); %Blurring in X and Y direction
figure; imagesc(abs(H)); title('Frequency Spectrum of the Filter');

%% Applying Inverse Filter
% filtering in f domain
I_f=fft2(I);
I_motion_f=fftshift(I_f);

% figure; imagesc(real(I_motion_f)); axis image; title('Freq of Image');

%% Reconstruction
InvFilt = 1./H;

if USE_WIENER
    InvFilt = abs(H).^2./(H.*(abs(H).^2 + NSR));
end

if USE_PSEUDO
    threshold = 0.025;
    InvFilt(abs(H)<threshold)=0;
end

I_recon_fn=I_motion_f.*InvFilt;
I_recon=ifft2(ifftshift(I_recon_fn));

figure;
subplot(121), imagesc(I); colormap(gray); title('Original Image');
subplot(122); imagesc(abs(I_recon)); colormap(gray); title('Reconstruction')

%% Lets sharpen it back up with laplacian
mask = [0 -1 0; -1 5 -1; 0 -1 0];
figure;
imagesc(conv2(abs(I_recon),mask)); title('Sharper'); colormap(gray); axis image;