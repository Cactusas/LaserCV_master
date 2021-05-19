clear; close all; clc;
main();

function main
%We had good results on 640x360 and 1280x1024 photos
%with these absolute values. Just saving.
% HP_filter_size = 50;
% kernel_size = 11;
% diag_kernel_size = 31;
% fill_gap = 20;
% min_length = 100;
% median_filter_size = 4;

%Constants percentage to image size
const_perc.HP_filter_size = 10;
const_perc.kernel_size = 1;
const_perc.diag_kernel_size = 0.7;
const_perc.fill_gap = 5;
const_perc.min_length = 20;
const_perc.fix_size = 1;

% Action what our code will process:
% 0 - photo
% 1 - video
% 2 - experiment with photo
action = 0;
switch action
    case 0 %Photo
        img_org = imread('C:\Users\ivano\Desktop\Magistrinis\LaserCV_master\320x180\Picture 33.jpg');
        
        [const, kernels, kernels_diag, HP_filter] = init_detection(img_org, const_perc);
        line = perform_detection(img_org, const, kernels, kernels_diag, HP_filter);
        figure(1); imshow(img_org); hold on;
        plot_line(line);
        
    case 1 %Video
        reader = VideoReader('C-s_L-d.webm');
        start_frame = 110;
        img_org = readFrame(reader);
        [const, kernels, kernels_diag, HP_filter] = init_detection(img_org, const_perc);
        %Just read some starting frames
        for i = 1:start_frame
            readFrame(reader);
        end
        
        while hasFrame(reader)
            img_org = readFrame(reader);
            line = perform_detection(img_org, const, kernels, kernels_diag, HP_filter);
            figure(1); imshow(img_org); hold on;
            plot_line(line);
        end
        
    case 2 %Experiment
        img_org = imread('Picture HD.jpg');
        figure(1); imshow(img_org); hold on;
        xi = zeros(1,2); yi = zeros(1,2);
        for i=1:2
            [xi(i), yi(i)] = ginput(1);
            plot(xi(i),yi(i),'x','LineWidth',2,'Color','red');
        end
        plot([xi(1) xi(2)],[yi(1) yi(2)],'LineWidth',2,'Color','red');
        
        [const, kernels, kernels_diag, HP_filter] = init_detection(img_org, const_perc);
        line = perform_detection(img_org, const, kernels, kernels_diag, HP_filter);
        plot_line(line);
        
        %Calculate angle error
        user_angle = points_to_angle(xi, yi);
        detected_angle = points_to_angle([line.point1(1) line.point2(1)]...
            ,[line.point1(2) line.point2(2)]);
        angle_error = abs(user_angle - detected_angle);
        
        %Calculate detected points error
        %Calculate required parameters for line equation
        m = (yi(2)-yi(1))/(xi(2)-xi(1));
        b = yi(1) - m*xi(1);

        %Populate line equation with x's and y's
        x = 1:size(img_org,2);
        y = m*x + b;
        
        %Calculate distances between one of detected points and
        %points that user marked
        distances = zeros(1,size(img_org,2));
        for i=1:length(y)
            distances(i) = euclidean_distance(line.point1, [x(i) y(i)]);
        end
        point_error = min(distances);
        
        fprintf('Angle error: %f\n', angle_error);
        fprintf('Point error: %f\n', point_error);
        
    otherwise
        error('Unknown action %d', action);
end
end

% Perform all the laser line detection algorithm.
% img_org - RGB image where to detect laser line.
% const - structure of const values.
% kernels, kernels_diag - arrays of rotated kernels.
% HP_filter - Gaussian High Pass filter.
% ret - detected laser line.
function ret = perform_detection(img_org, const, kernels, kernels_diag, HP_filter)
%ret = best_line([2 2 2 2 2 2 2 2 2]);%ASDDASDASDASDDDASDASDAS
fprintf('Detecting...\n');
%Apply Fourier Transform to red channel of the image
img_red = img_org(:,:,1);
img_fft = fft2(img_red);
img_fft = fftshift(img_fft);

%Apply High Pass Filter to image frequency spectrum
fprintf('Filtering frequencies... ');
img_filtered = img_fft.*HP_filter;
img_filtered = fftshift(img_filtered);
img_filtered = ifft2(img_filtered);
img_filtered = uint8(real(img_filtered));
fprintf('Done\n');

%Get all the lines based on rotating kernel
fprintf('Detecting lines on kernel rotation... ');
lines = get_lines(img_filtered, kernels, const.fill_gap, const.min_length);
fprintf('Done\n');

%Get half of the lines with biggest point count
intensities = [lines(1:end).intensity];
[~,ind] = maxk(intensities, ceil(length(lines)*0.5));
intense_lines = lines(ind);

%Get pixel values of each line
fprintf('Acquiring pixels of lines for fixing points... ');
intense_lines = add_lines_pixels(intense_lines, kernels_diag, img_org, const.fix_size);
fprintf('Done\n');
%Fix line points and clean pixels then again add pixels
for i=1:length(intense_lines)
    intense_lines(i) = fix_line_points(intense_lines(i), const.fix_size);
    intense_lines(i).Rvals = [];
    intense_lines(i).Gvals = [];
    intense_lines(i).Bvals = [];
    intense_lines(i).Xvals = [];
    intense_lines(i).Yvals = [];
end
fprintf('Reacquiring pixels of lines... ');
intense_lines = add_lines_pixels(intense_lines, kernels_diag, img_org, 0);
fprintf('Done\n');
fprintf('Possible lines found: %d\n', length(intense_lines));

%Choose the best line based on particular rules
fprintf('Searching for laser line... ');
line = best_line(intense_lines);
fprintf('Done\n');

%Place line points in the center of line
line = fix_line_points(line, const.fix_size);
fprintf('Line detected at points: [%d %d] [%d %d] ', ...
    line.point1(1), line.point1(2), line.point2(1), line.point2(2));
fprintf('with angle %f\n', points_to_angle([line.point1(1), line.point2(1)],...
    [line.point1(2), line.point2(2)]))
ret = line;
fprintf('-------------------------------------------------------\n');
end

% Perform initialization of laser line detection
% img - image to initialize detection to.
% const_percents - structure of constants percentage.
% const - structure of const values.
% kernels, kernels_diag - arrays of rotated kernels.
% HP_filter - Gaussian High Pass filter.
function [const, kernels, kernels_diag, HP_filter] = init_detection(img, const_percents)
fprintf('Starting detector initialization...\n');
const = init_constants(const_percents, img);

kernel = generate_kernel(const.kernel_size);
kernels = kernel_rotate(kernel, false);
kernels_cnt = size(kernels,3);
fprintf('Rotated kernels count %d ', kernels_cnt);
if (kernels_cnt > 45)
    n = round(kernels_cnt/45);
    kernels = kernels(:,:,n:n:end);
    fprintf('was reduced to %d', size(kernels,3));
end
fprintf('\n');

kernel_diag = generate_kernel(const.diag_kernel_size);
kernels_diag = kernel_rotate(kernel_diag, false);

height = size(img,1);
width = size(img,2);
HP_filter = Gaussian_HP_filter(height, width, const.HP_filter_size);
fprintf('Initialization done\n');
fprintf('-------------------------------------------------------\n');
end

% Calculate constant values based on image size.
% const - structure constant percentages by image size.
% img - image to calculate constants from.
function ret = init_constants(const, img)

height = size(img,1);
width = size(img,2);
frame_size = sqrt(height^2 + width^2);
frame_size = real(frame_size);

fprintf('Frame given with size %dx%d\n', width, height);

fprintf('Initializing detector constants with sizes: \n');
%Simply calculate value from percent given and round up to required number
ret.HP_filter_size = round(const.HP_filter_size*frame_size/100);
ret.kernel_size = round_odd(const.kernel_size*frame_size/100);
ret.diag_kernel_size = round_odd(const.diag_kernel_size*frame_size/100);
ret.fill_gap = round(const.fill_gap*frame_size/100);
ret.min_length = round(const.min_length*frame_size/100);
ret.fix_size = round(const.fix_size*frame_size/100);

fprintf('High Pass filter size: %d\n', ret.HP_filter_size);
fprintf('Kernel size: %d\n', ret.kernel_size);
fprintf('Diagonal kernel size: %d\n', ret.diag_kernel_size);
fprintf('Hough Transform fill gap: %d\n', ret.fill_gap);
fprintf('Hough Transform minimum line length: %d\n', ret.min_length);
fprintf('Fixing line points area size: %d\n', ret.fix_size);
end

% Round to nearest odd integer.
% S - number to round.
function S = round_odd(S)
idx = mod(S,2)<1;
S = floor(S);
S(idx) = S(idx)+1;
end

% Move line points to brightest spot (line center)
% line - line to move.
% fix_size - line length to average
% ret - line with fixed points
function ret = fix_line_points(line, fix_size)
Rpt1 = line.Rvals(1:fix_size,:);
Rpt2 = line.Rvals(end-fix_size+1:end,:);

% a1 = Rpt1;
% a1(:,:,2) = line.Gvals(1:fix_size,:);
% a1(:,:,3) = line.Bvals(1:fix_size,:);
% 
% a2 = Rpt2;
% a2(:,:,2) = line.Gvals(end-fix_size+1:end,:);
% a2(:,:,3) = line.Bvals(end-fix_size+1:end,:);
% 
% 
% figure(1);
% imshow(a1,'InitialMagnification','fit');
% figure(2);
% imshow(a2,'InitialMagnification','fit');
Rpt1 = mean(Rpt1);
Rpt2 = mean(Rpt2);

[~, ind1] = max(Rpt1(:));
[~, ind2] = max(Rpt2(:));
line.point1(1) = line.Xvals(1,ind1);
line.point1(2) = line.Yvals(1,ind1);
line.point2(1) = line.Xvals(end,ind2);
line.point2(2) = line.Yvals(end,ind2);

ret = line;
end

% Choose the best line to match red laser line from array of lines.
% lines - array of lines.
% ret - chosen best line.
function ret = best_line(lines)
max_coef = 0;
best_line_ind = 0;
for i = 1:length(lines)
    Rvals = lines(i).Rvals;
    Gvals = lines(i).Gvals;
    Bvals = lines(i).Bvals;


%          ab=[];
%          ab(:,:,1)= uint8(Rvals);
%          ab(:,:,2)= uint8(Gvals);
%          ab(:,:,3)= uint8(Bvals);
%           imshow(uint8(ab));
    
    
    %Divide line into 3 parts - rise/center/fall
    %     edge_len = fix(size(Rvals,2)/3);
    %     center_len = mod(size(Rvals,2),3);
    %
    %     risingR = mean(Rvals(:,1:edge_len));
    %     centerR = mean(Rvals(:,edge_len+1:edge_len+edge_len+center_len));
    %     fallingR = mean(Rvals(:,edge_len+edge_len+center_len+1:end));
    %
    %     risingG = mean(Gvals(1:edge_len));
    %     centerG = mean(Gvals(edge_len+1:edge_len+edge_len+center_len));
    %     fallingG = mean(Gvals(edge_len+edge_len+center_len+1:end));
    %
    %     risingB = mean(Bvals(1:edge_len));
    %     centerB = mean(Bvals(edge_len+1:edge_len+edge_len+center_len));
    %     fallingB = mean(Bvals(edge_len+edge_len+center_len+1:end));
    %Divide line into 3 parts - rise/center/fall
    
%     Rvals = 1:11;
%     Gvals = 1:11;
%     Bvals = 1:11;
    
    width = size(Rvals,2);
    sides = round(width/3);
    center = width - sides*2;
    if (center > sides)
       center = center - 2 ;
       sides = sides + 1;
    end
    
    Cvals = [];
    Cvals(:,:,1)= uint8(Rvals);
    Cvals(:,:,2)= uint8(Gvals);
    Cvals(:,:,3)= uint8(Bvals);
    Cvals= rgb2gray(Cvals);
    
    %risingC = mean(diff(Cvals(1:sides),1,2));
    %centerC = mean(Cvals(sides+1:sides+center));
    %fallingC = mean(diff(Cvals(sides+center+1:end),1,2));
   
    risingR = mean(Rvals(1:sides));
    centerR = mean(Rvals(sides+1:sides+center));
    fallingR = mean(Rvals(sides+center+1:end));
    
    risingG = mean(Gvals(1:sides));
    centerG = mean(Gvals(sides+1:sides+center));
    fallingG = mean(Gvals(sides+center+1:end));
    
    risingB = mean(Bvals(1:sides));
    centerB = mean(Bvals(sides+1:sides+center));
    fallingB = mean(Bvals(sides+center+1:end));
    
    risingC = (risingR + risingG + risingB)/3;
    fallingC = (fallingR + fallingG + risingB)/3;
    
%      i
    %if (risingC > 0 && fallingC < 0)
        if (risingG + risingB == 0)
            rise_avg = risingR;
            fall_avg = fallingR;
        else
            rise_avg = (risingR-((risingG + risingB)/2))*3000;
            fall_avg = (fallingR-((fallingG + fallingB)/2))*3000;
        end
        
        center_avg = (centerR + centerG + centerB)/3;
        coef = (rise_avg + fall_avg + center_avg);
        
        %coef = centerR;
        if (coef > max_coef)
            max_coef = coef;
            best_line_ind = i;
        end
        
    %end
    
end


%abs(rise_avg-fall_avg)
%max_coef
ret = lines(best_line_ind);
end

% Add pixel values of lines to line structure array.
% lines - lines structure array.
% kernels_diag - array of another (bigger) rotated kernels, used for
%                extracting pixel values diagonally to line.
% img - RGB image to extract pixel values from.
% ret - structure of lines with pixel values.
function ret = add_lines_pixels(lines, kernels_diag, img, fix)
height = size(img,1);
width = size(img,2);
for i = 1:length(lines)
    %lines(i) = fix_line_points(lines(i), 10);
    x = [lines(i).point1(1) lines(i).point2(1)];
    y = [lines(i).point1(2) lines(i).point2(2)];
    
    %Get diagonal kernel
    %diagonal_degrees = rotation_to_degrees(kernels, lines(i).rotation) + 90;
    %diagonal = degrees_to_kernel(kernels_diag, diagonal_degrees);
    
    diagonal_degrees = points_to_angle(x,y) + 90;
    diagonal = degrees_to_kernel(kernels_diag, diagonal_degrees);
    
    %Get line that will be used to extract pixel values.
    test_line = kernel_to_line(diagonal);
    center_line = ceil(size(test_line, 1)/2);
    
    [cx, cy, ~] = improfile(img, x, y);
    cx = round(cx);
    cy = round(cy);
    
    if (fix ~= 0)
        %Acquire only pixels required for line point fix
        k = [1:fix length(cx)-fix:length(cx)];
    else
        k = 1:length(cx);
    end
    
    for j=k
        cnt = 1;
        for ci = 1:size(test_line, 1)
            for cj = 1:size(test_line, 2)
                if (test_line(ci, cj) == 0)
                    continue;
                end
                y=cy(j)-center_line+ci;
                x=cx(j)-center_line+cj;
                if (x < 1 || x > width || y < 1 || y > height)
                    continue;
                end
                %Extract pixel values.
                lines(i).Rvals(j,cnt) = img(cy(j)-center_line+ci, cx(j)-center_line+cj,1);
                lines(i).Gvals(j,cnt) = img(cy(j)-center_line+ci, cx(j)-center_line+cj,2);
                lines(i).Bvals(j,cnt) = img(cy(j)-center_line+ci, cx(j)-center_line+cj,3);
                lines(i).Xvals(j,cnt) = cx(j)-center_line+cj;
                lines(i).Yvals(j,cnt) = cy(j)-center_line+ci;
                %Just for visualization
                %img(cy(j)-center_line+ci, cx(j)-center_line+cj) = 255;
                cnt = cnt + 1;
            end
        end
    end
    %imshow(img);
end
ret = lines;
end

% Get all the lines from image based on rotating kernel.
% img - image to search lines from.
% kernels - rotated kernels array.
% fill_gap, min_length - Hough Transform corresponding parameters.
function ret = get_lines(img, kernels, fill_gap, min_length)
%Structure initialization
lines(1).point1 = [0 0];
lines(1).point2 = [0 0];
lines(1).rotation = 0;
line_saved = false;
frame_size = 10; %To remove points from the very edges of binary image.
height = size(img,1);
width = size(img,2);

% imshow(img);
% k = degrees_to_kernel(kernels, 105);
% filtered = imfilter(img, k);
% imshow(filtered);
% filtered = imadjust(filtered); 
% imshow(filtered);
% img_bin  = filtered > 127;
% imshow(img_bin);
% img_bin = medfilt2(img_bin,[2 2]);
% imshow(img_bin);
% 
% k = degrees_to_kernel(kernels, 90);
% SE = kernel_to_strel(k);
% SE = strel('line',17, 45);
% img_bin = imopen(img_bin, SE);
% imshow(img_bin);
%imshow(mat2gray(k),'InitialMagnification','fit');

for i=1:size(kernels,3)
    kernel = kernels(:,:,i);
    
    %Apply kernel
    filtered = imfilter(img, kernel);
    
    %filtered = imgaussfilt(filtered,[2 2]);
    %imshow(filtered);
    %Top-hat
    SE = kernel_to_strel(kernel);
    %degrees = rotation_to_degrees(kernels, i);
    %SE = strel('line',size(kernels,1), degrees+90);
    
    filtered = imadjust(filtered);   
    %imshow(top_hat);

    %Threshold and median filter
    img_bin = filtered > 127;
    %imshow(img_bin);
    img_bin = medfilt2(img_bin,[2 2]);
%     figure(3);
%     imshow(img_bin);
    img_bin = imopen(img_bin, SE);
    %imshow(img_bin);
    
    %Remove white points from edges of the frame.
    img_bin(:,1:frame_size)=0;
    img_bin(:,width-frame_size:end)=0;
    img_bin(1:frame_size,:)=0;
    img_bin(height-frame_size:end,:)=0;
    
%     figure(1);
%     imshow(img_bin);
%     figure(2);
%     imshow(mat2gray(kernel),'InitialMagnification','fit');
    
    %---------Save lines and kernel rotation---------
    hough_lines = perform_hough(img_bin, fill_gap, min_length);
    for j=1:length(hough_lines)
        if line_saved == false
            lines(1).point1 = hough_lines(j).point1;
            lines(1).point2 = hough_lines(j).point2;
            lines(1).intensity = hough_lines(j).intensity;
            lines(1).rotation = i;
            line_saved = true;
        else
            lines(end+1).point1 = hough_lines(j).point1;
            lines(end).point2 = hough_lines(j).point2;
            lines(end).intensity = hough_lines(j).intensity;
            lines(end).rotation = i;
        end
    end
end
ret = lines;
end

% Perform Hough Transformation on binary image with provided parameters.
% img_bin - binary image to serach lines.
% fill_gap - minimum distance between two line segments to merge that
%            line into one line
% min_length - minimum length of line to consider it to be a line.
% ret  - structures of lines found with their parameters.
function ret = perform_hough(img_bin, fill_gap, min_length)
[H,T,R] = hough(img_bin,'RhoResolution',1,'Theta',-90:1:89);
P  = houghpeaks(H,30);
lines = houghlines(img_bin,T,R,P,'FillGap',fill_gap,'MinLength',min_length);

for i = 1:length(lines)
    theta = lines(i).theta;
    rho = lines(i).rho;
    T_theta = T==theta;
    R_rho = R==rho;
    hough_intensity = H(R_rho, T_theta);
    lines(i).intensity = hough_intensity;
end

ret=lines;
end

% Rotate kernel in all possible rotations.
% kernel - kernel in horizontal position.
% visualize - if set to true, visualization of rotation will be shown.
% ret - all the rotated kernels in one array.
function ret = kernel_rotate(kernel, visualize)
% Check if kernel dimensions are the same. 
if (size(kernel,1) ~= size(kernel,2))
   error('Kernel dimensions must be the same.'); 
end

% Check if kernel size is an odd number.
dim = size(kernel,1);
if (mod(dim,2) == 0)
    error('Kernel dimensions must be an odd number.'); 
end

center = round(dim/2);
rot = center;

kernel_org = kernel;

%Draw frame of zeros on kernel so we would have enough space for rotation.
kernel=conv2(kernel,[0,0,0;0,1,0;0,0,0]);
hor = kernel; %Save initial horizontal kernel.
hor_rot = [];
k=1;
%45 degrees clockwise from horizontal
for i=1:rot
    for j=rot:-1:i+1
        kernel(center+i:center+i+2, center+j) = kernel(center+i-1:center+i+1, center+j);
        kernel(center+i-1,center+j) = 0;
        
        kernel(center-i:center-i+2, center-j+2) = kernel(center-i+1:center-i+3, center-j+2);
        kernel(center-i+3,center-j+2) = 0;
        
        hor_rot(:,:,k) = kernel;
        k = k+1;
    end
end

kernel=conv2(kernel_org,[0,0,0;0,1,0;0,0,0]);
hor_inv = [];
k=1;
%45 degrees counter-clockwise from horizontal
for i=1:rot
    for j=rot:-1:i+1
        kernel(center+i:center+i+2, center-j+2) = kernel(center+i-1:center+i+1, center-j+2);
        kernel(center+i-1,center-j+2) = 0;
        
        kernel(center-i:center-i+2, center+j) = kernel(center-i+1:center-i+3, center+j);
        kernel(center-i+3,center+j) = 0;
        
        hor_inv(:,:,k) = kernel;
        k=k+1;
    end
end


%---------------------------Make kernel vertical---------------------------
kernel=rot90(kernel_org);
kernel=conv2(kernel,[0,0,0;0,1,0;0,0,0]);

ver = kernel;
ver_inv = [];
k=1;
%45 degrees counter-clockwise from vertical
for i=1:rot
    for j=rot:-1:i+1
        kernel(center+j, center+i:center+i+2) = kernel(center+j, center+i-1:center+i+1);
        kernel(center+j, center+i-1) = 0;
        
        kernel(center-j+2, center-i:center-i+2) = kernel(center-j+2, center-i+1:center-i+3);
        kernel(center-j+2, center-i+3) = 0;
        
        ver_inv(:,:,k) = kernel;
        k=k+1;
    end
end

kernel=rot90(kernel_org);
kernel=conv2(kernel,[0,0,0;0,1,0;0,0,0]);

ver_rot = [];
k=1;
%45 degrees clockwise from vertical
for i=1:rot
    for j=rot:-1:i+1
        kernel(center-j+2, center+i:center+i+2) = kernel(center-j+2, center+i-1:center+i+1);
        kernel(center-j+2, center+i-1) = 0;
        
        kernel(center+j, center-i:center-i+2) = kernel(center+j, center-i+1:center-i+3);
        kernel(center+j, center-i+3) = 0;
        
        ver_rot(:,:,k) = kernel;
        k=k+1;
    end
end

%Now put all the kernels in order in a single array.
kernels = hor;
kernels(:,:,end+1:end+size(hor_rot,3)) = hor_rot;
kernels(:,:,end+1:end+size(ver_inv,3)) = flip(ver_inv,3);
kernels(:,:,end+1:end+size(ver,3)) = ver;
kernels(:,:,end+1:end+size(ver_rot,3)) = ver_rot;
kernels(:,:,end+1:end+size(hor_inv,3)) = flip(hor_inv,3);
ret = kernels;

%---------------VIZUALISATION---------------
if (visualize == true)
    figure(99);
    viz_kernels = kernels;
    for i=1:size(kernels,3)
        mat=kernels(:,:,i);
        mat(mat==2)=255;
        mat(mat==-1)=127;
        viz_kernels(:,:,i)=mat;
    end
    
    for i=1:size(viz_kernels,3)
        imshow(mat2gray(viz_kernels(:,:,i)),'InitialMagnification','fit');
        pause(1);
    end
end
end

% Generates horizontal kernel with provided size.
% kernel_size - kernel size, has to be an odd number.
% kernel - generated kernel.
function kernel = generate_kernel(kernel_size)
%Check if size provided is an odd number.
if (mod(kernel_size, 2) == 0)
   error('Kernel size must be odd number') ;
end

kernel = zeros(kernel_size,kernel_size);
center = ceil(kernel_size/2);
kernel(center,:) = 2;
kernel(center - 1,:) = -1;
kernel(center + 1,:) = -1;
end

% Converts kernel to Structuring Element for Top-Hat Filter.
% kernel - kernel to convert.
% se - Structuring Element for Top-Hat Filter.
function se = kernel_to_strel(kernel)
se = kernel;
se(1,:) = [];
se(end,:) = [];
se(:,1) = [];
se(:,end) = [];
se(se == -1) = 0;

se( ~any(se,2), : ) = [];  %rows
se( :, ~any(se,1) ) = [];  %columns
se = strel(se);
end

% Converts kernel used for convolution to binary matrix.
% kernel - kernel to convert.
% ret - binary matrix converted from kernel.
function ret = kernel_to_line(kernel)
ret = kernel;
%We only need to know where 2's were in kernel.
ret(ret == -1) = 0;
ret(ret == 2) = 1;
%Remove the frame from matrix which was addded during kernel generation.
ret(1,:) = [];
ret(end,:) = [];
ret(:,1) = [];
ret(:,end) = [];
end

% Generates Gaussian High Pass filter with provided radius.
% height, width - dimensions of image to filter.
% rad - radius of Gaussian HP filter.
% ret - Gaussian High Pass filter.
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

% Simple Euclidean distance calculates distance between 2 points.
% point1, point2 - [x y] coordinates of points.
% dist - calulated distance between 2 points.
function dist = euclidean_distance(point1, point2)
dist = sqrt((point1(1) - point2(1))^2 + (point1(2) - point2(2))^2);
end

% ------------------------UNUSED------------------------ %
% Returns diagonal kernel based on rotation.
% kernels - array of rotated kernels.
% rotation - rotation number in kernels array.
% ret - diagonal kernel.
% function ret = diagonal_kernel(kernels, rotation)
% kernels_cnt = size(kernels,3);
% diagonal = kernels_cnt/2 + rotation;
% if (diagonal > kernels_cnt)
%     diagonal = diagonal - kernels_cnt;
% end
% ret = kernels(:,:,diagonal);
% end
% ------------------------UNUSED------------------------ %

% Returns degrees based on kernel rotation number.
% kernels - array of rotated kernels.
% rotation - rotation number in kernels array.
% degrees - numerical value of kernel rotation in degrees.
function degrees = rotation_to_degrees(kernels, rotation)
kernels_cnt = size(kernels,3);
rot_res = 180/kernels_cnt; %How many degrees per rotation
degrees_arr = 0:rot_res:180-rot_res;
degrees = degrees_arr(rotation);
end

% Returns kernel based on degrees it is rotated.
% kernels - array of rotated kernels.
% degrees - kernel rotation in degrees.
% kernel - kernel rotated in provided degrees.
function kernel = degrees_to_kernel(kernels, degrees)
%We take how many degrees kernel rotated form horizontal position
if (degrees > 180)
    degrees = degrees - 180;
end

kernels_cnt = size(kernels,3);
rot_res = 180/kernels_cnt; %How many degrees per rotation
degrees_arr = 0:rot_res:180-rot_res;
[~,degrees_ind] = min(abs(degrees_arr - degrees)); %Nearest kernel
kernel = kernels(:,:,degrees_ind);
end

% Plot line on an image.
% img - image where to plot the line.
% line - line to plot.
function plot_line(line)
xy = [line.point1; line.point2];
plot(xy(:,1),xy(:,2),'LineWidth',4,'Color','blue');
%Plot beginnings and ends of lines
plot(xy(1,1),xy(1,2),'x','LineWidth',2,'Color','blue');
plot(xy(2,1),xy(2,2),'x','LineWidth',2,'Color','blue');
end

% Calculate line rotation angle from 2 points.
% Angle show how many degrees line is rotated from horizontal position.
% x[2], y[2] - point coordinates.
% angle - calculated angle.
function angle = points_to_angle(x, y)
x1=x(1);
y1=y(1);
x2=x(2);
y2=y(2);
slope = (y2 - y1) ./ (x2 - x1);
angle = atand(slope);
if (angle < 0)
    angle = 180 + angle;
end
end
