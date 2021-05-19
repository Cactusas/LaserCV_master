clear; close all; clc;
main();

function main
%imshow(log(abs(img)),[]); <--- I got this expression from the link below
%https://www.youtube.com/watch?v=nWGirrTlpnk&t=879s

percents = zeros(200,3);
noise = zeros(200,3);

%CONST
%Image dimensions

img = imread('Picture 16.jpg');
height = size(img,1);
width = size(img,2);

green_img = imread('Picture 16_1.png');
%green_img = imresize(green_img, [height width]);
green_img = green_line(green_img);

% img = img(:,:,1);
% img = imresize(img, [height width]);
% img = double(img);
% img = FFT(img);
% img = switch_quadrants(img);

img = img(:,:,1);
img = fft2(img);
img = fftshift(img);

for i = 1:200
    
    
    %Filtering related constants
    filter_radius = i;
    butter_order = 15;
    
    %Frequency filtering
    filter_ideal = Ideal_HP_filter(height, width, filter_radius);
    filter_gauss = Gaussian_HP_filter(height, width, filter_radius);
    filter_butter = Butterworth_HP_filter(height, width, filter_radius, butter_order);
    
    img_ideal = img.*filter_ideal; %Apply filter
    img_gauss = img.*filter_gauss; %Apply filter
    img_butter = img.*filter_butter; %Apply filter
    
    img_ideal = fftshift(img_ideal); %Switch back quadrants
    img_gauss = fftshift(img_gauss); %Switch back quadrants
    img_butter = fftshift(img_butter); %Switch back quadrants
    
    img_ideal = ifft2(img_ideal); %Inverse FFT
    img_ideal = uint8(real(img_ideal));
    img_gauss = ifft2(img_gauss); %Inverse FFT
    img_gauss = uint8(real(img_gauss));
    img_butter = ifft2(img_butter); %Inverse FFT
    img_butter = uint8(real(img_butter));
    
    img_ideal = threshold(img_ideal, 0);
    img_gauss = threshold(img_gauss, 0);
    img_butter = threshold(img_butter, 0);
%            figure(1);
%            imshow(img_ideal);
%            title('Ideal');
%            figure(2);
%            imshow(img_gauss);
%            title('Gauss');
%            figure(3);
%            imshow(img_butter);
%            title('Butter');
    
    
    percents(i,1) = filtered_laser(green_img, img_ideal);
    percents(i,2) = filtered_laser(green_img, img_gauss);
    percents(i,3) = filtered_laser(green_img, img_butter);
    
    
    noise(i,1) = calc_noise(green_img, img_ideal);
    noise(i,2) = calc_noise(green_img, img_gauss);
    noise(i,3) = calc_noise(green_img, img_butter);
end


%We are looking for a red line, so use R channel
%       img = img(:,:,1);
%
%       %Cooley-Tukey algorithm limitation, power of 2
%       img = imresize(img, [height width]);
%
%
%       %FFT
%       img = double(img);
%       img = FFT(img);
%
%       img = switch_quadrants(img);
%
%       %Frequency filtering
%       filter_ideal = Ideal_HP_filter(height, width, 20);
%       filter_gauss = Gaussian_HP_filter(height, width, 20);
%       filter_butter = Butterworth_HP_filter(height, width, 20, 3);
%
%
%       img_ideal = img.*filter_ideal; %Apply filter
%       img_gauss = img.*filter_gauss; %Apply filter
%       img_butter = img.*filter_butter; %Apply filter
%
%       %Back to spatial domain
%       img_ideal = switch_quadrants(img_ideal); %Switch back quadrants
%       img_gauss = switch_quadrants(img_gauss); %Switch back quadrants
%       img_butter = switch_quadrants(img_butter); %Switch back quadrants
%
%       img_ideal = ifft2(img_ideal); %Inverse FFT
%       img_ideal = uint8(real(img_ideal));
%       img_gauss = ifft2(img_gauss); %Inverse FFT
%       img_gauss = uint8(real(img_gauss));
%       img_butter = ifft2(img_butter); %Inverse FFT
%       img_butter = uint8(real(img_butter));
%
%       figure(11);
%       imshow(img_ideal);
%       title("Ideal HP filtered");
%       figure(12);
%       imshow(img_gauss);
%       title("Gaussian HP filtered");
%       figure(13);
%       imshow(img_butter);
%       title("Butterworth HP filtered");
%
%       img = threshold(img, threshold_val);

%Detect lines, it will draw them on figure(2)
%Hough_transform(img);

%You can also pass your custom threshold value, for example:
%Hough_transform(img, 50);
end

function ret = green_line(img)
height = size(img,1);
width = size(img,2);
ret = zeros(height, width);

for i = 1:height
    for j = 1:width
        R = img(i,j,1);
        G = img(i,j,2);
        B = img(i,j,3);
        
        if (R == 0 && B == 0 && G == 255)
            ret(i,j) = 1;
        end
    end
end
end

function ret = filtered_laser(green_img, img)
[rows, cols] = find(green_img == 1);
len = length(rows);

filtered = 0;

for i = 1:len
    val = img(rows(i), cols(i));
    if (val == 0)
        filtered = filtered + 1;
    end
end
percent = filtered * 100 / len;
ret = percent;
end

function ret = calc_noise(green_img, img)
[rows, cols] = find(green_img == 1);
len = length(rows);

tmp = img;
for i = 1:len
    tmp(rows(i), cols(i)) = 0;
end
[rows2, cols2] = find(tmp == 1);
ret = length(rows2);
end


%Switches quadrants of an image
%img - image to switch qaudrants
%ret - image with switched quadrants
function ret = switch_quadrants(img)
%Image dimensions
height = size(img,1);
width = size(img,2);

%Simply get all 4 quadrants
Q1 = img(1:height/2, 1:width/2);
Q2 = img(1:height/2, width/2+1:width);
Q3 = img(height/2+1:height, 1:width/2);
Q4 = img(height/2+1:height, width/2+1:width);

%Now from a new image with quadrants switched
ret = [Q4 Q3; Q2 Q1];
end

%Creates High Pass filter for spectrum
%height, width - spectrum dimensions
%rad - filter radius
%ret - created new filter
function ret = Ideal_HP_filter(height, width, rad)
center = [height/2 width/2]; %Center of spectrum
ret = ones(height, width);

%Basically just draw a circle in the center
for i = center(1)-rad:center(1)+rad
    for j = center(2)-rad:center(2)+rad
        if (euclidean_distance([j,i], center) < rad)
            ret(j,i) = 0;
        end
    end
end
end

function ret = Gaussian_HP_filter(height, width, rad)
center = [height/2 width/2]; %Center of spectrum
ret = ones(height, width);

for i = 1:width
    for j = 1:height
        dist = euclidean_distance([j,i], center);
        ret(j,i) = 1 - exp(-dist^2/(2*rad^2));
    end
end
end

function ret = Butterworth_HP_filter(height, width, rad, n)
center = [height/2 width/2]; %Center of spectrum
ret = ones(height, width);

for i = 1:width
    for j = 1:height
        dist = euclidean_distance([j,i], center);
        ret(j,i) = 1-1/(1+(dist/rad)^(2*n));
    end
end
end

%Calculates distance between two points
%point1, point2 - [x y] coordinates
%dist - calculated distance between two points
function dist = euclidean_distance(point1, point2)
dist = sqrt((point1(1) - point2(1))^2 + (point1(2) - point2(2))^2);
end

%Performs 1D FFT
%singal - signal to perfrom FFT
%ret - FFT coefficients of signal
function ret = FFT__(signal)
len = length(signal);

%Last division of signal, return what we got
if len == 1
    ret(1) = signal(1);
    return;
end

%Don't mix up, matlab indexing starts from 1, FFT formula summation
%starts from 0
Ek = FFT__(signal(1:2:len)); %FFT of even indexes
Ok = FFT__(signal(2:2:len)); %FFT of odd indexes
len_h = len*0.5; %Half length of signal

for k = 1:len_h
    ret(k) = Ek(k) + exp(-2i*pi*(k-1)/len)*Ok(k);
    ret(k + len_h) = Ek(k) - exp(-2i*pi*(k-1)/len)*Ok(k);
end

end

%Performs FFT on the image
%img - image to perform FFT
%ret - image in frequency domain
function ret = FFT(img)
height = size(img,1);
width = size(img,2);
tmp = zeros(height, width);
for i = 1:height
    tmp(:,i) = FFT__(img(1:height,i));
end
for i = 1:width
    tmp(i,:) = FFT__(tmp(i,1:width));
end
ret = tmp;
end

%Thresholds a picutre with (max - thr)
%img - picture to threshold
%thr - threshold value (max - thr)
%ret - thresholded binary picture
function ret = threshold(img, thr)
%Image dimensions
height = size(img,1);
width = size(img,2);

%Get maximum value in an image
max_col = max(img);
max_total = max(max_col);

ret = zeros(height, width);
%thr = max_total - thr;

%We don't want negative values
%     if (thr <= 0)
%        error("Threshold is too big. Max = %d", max_total-1);
%     end

%Simply iterate through pixels and threshold them
for i = 1:height
    for j = 1:width
        val = img(i,j);
        if (val > thr)
            ret(i,j) = 1;
        end
    end
end
end

%Performs Hough transform and draw lines. It will draw the lines on
%figure(1)
%img - binary image to detect lines from
%thr - how many points a line must consist of
%Ref.: https://nl.mathworks.com/matlabcentral/fileexchange/4983-line-detection-via-standard-hough-transform
%Might be useful: https://en.wikipedia.org/wiki/Hough_transform
%Also might be useful: https://www.youtube.com/watch?v=4zHbI-fFIlI (6 mins)
function Hough_transform(img, thr)
%Image dimensions
height = size(img,1);
width = size(img,2);

%Hough space resolution
p_step = 1;
teta_step = 1;

%Hough space size
p_size = sqrt(height^2 + width^2);
teta_size = 180; %Degrees

%Form our empty Hough space
p = -p_size:p_step:p_size;
teta = 0:teta_step:teta_size;
hough_space = zeros(length(p),length(teta)); %Accumulator

%Get 1's points from binary image
[y_bin, x_bin] = find(img);
pts_cnt = size(x_bin);

%Loop through the points
for i = 1:pts_cnt
    ind_teta = 0; %Teta index in our Hough Space
    for i_teta = teta %Loop through our tetas in Hough Space
        ind_teta = ind_teta+1; %Increment teta index
        %Line Normal Form, calculate length of support line
        p_val = x_bin(i)*cosd(i_teta)+y_bin(i)*sind(i_teta);
        %Get distances from our calculated length to all lengths in Hough Space
        dists = abs(p_val - p);
        %Pick the shortest distance (the closest match to what we just calculated)
        min_dist = min(dists);
        %Find the index of that shortest distance
        p_ind = find(dists == min_dist);
        %Increment value in Hough Space at particular angle and calculated length
        hough_space(p_ind,ind_teta) = hough_space(p_ind,ind_teta)+1;
    end
end
%To make it simple, we just iterate through all the tetas and increment
%value by 1 at calculated length (P) of support line in Hough Space.
%Since a line in Hough Space can be defined as a POINT, more times we get
%the same length of support line, the brighter local maxima will be.

%If threshold was not provided, pick the best one
thresh = max(hough_space);
thresh = max(thresh);
if (nargin == 2)
    thresh = thr;
end

regional_max = imregionalmax(hough_space); %Find regional maximums in Hough Space
[max_p, max_teta] = find(regional_max == 1); %Get the indexes of those maximums
hough_thr = hough_space - thresh; %Threshold our Hough Space
p_detect = [];
teta_detect = [];
for i = 1:length(max_p) %Loop through possible p
    if (hough_thr(max_p(i),max_teta(i)) >= 0) %Check if it meets our threshold
        %Append new lines
        p_detect = [p_detect; p(max_p(i))];
        teta_detect = [teta_detect; teta(max_teta(i))];
    end
end

%How many lines have we detected
lines_cnt = numel(p_detect);

x0 = 1;
xend = height;

%Draw the lines
for i = 1 : lines_cnt
    r = p_detect(i);
    th = teta_detect(i);
    
    %If theta is 0, just draw a vertical line
    if (th == 0)
        %figure(1);
        %line([r r], [1 height], 'Color', 'red');
    else
        %Starting y coordinate (derived from Normal Line Form)
        y0 = (-sind(th)/cosd(th))*x0 + (r / cosd(th));
        
        %Ending y coordinate (derived from Normal Line Form)
        yend = (-sind(th)/cosd(th))*xend + (r / cosd(th));
        
        %Draw the line
        %figure(1);
        %line( [y0, yend], [x0, xend], 'Color', 'red', 'LineWidth', 2);
    end
end

%Visualise calculated Hough Space
% figure('NumberTitle', 'off', 'Name', '4. Calculated Hough Space');
% tmp = uint8(hough_space);
% tmp = imresize(tmp, [512 512]);
% tmp = rescale(tmp, 0,255); %Rescale values to make image brighter
% imshow(uint8(tmp));
% xlabel(['0 \leq \Theta \leq ', num2str(teta_size), ' (teta in degrees)']);
% ylabel([num2str(-p_size), ' \leq P \leq ', num2str(p_size)]);
end





















