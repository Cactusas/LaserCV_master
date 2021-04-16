clear; close all; clc;
main();

function main

%VIDEO%
% file = 'C-d_L-s.webm';
% reader = VideoReader(file);
% for i = 1:50
%     readFrame(reader);
% end
% for z = 1:200
% img_org = readFrame(reader);
%VIDEO%

img_org = imread('Picture 18.jpg');
height = size(img_org,1);
width = size(img_org,2);

img_red = img_org(:,:,1);
img_fft = fft2(img_red);
img_fft = fftshift(img_fft);

gauss_filter = Gaussian_HP_filter(height, width, 5);
img_filtered = img_fft.*gauss_filter;
img_filtered = fftshift(img_filtered);
img_filtered = ifft2(img_filtered);
img_filtered = uint8(real(img_filtered));

%img_fft = log(abs(img_fft));
%img_fft = uint8(255 * mat2gray(img_fft));

kernel = generate_kernel(11);
kernels = kernel_rotate(kernel, false);
kernels_cnt = size(kernels,3);

kernel_big = generate_kernel(21);
kernels_big = kernel_rotate(kernel_big, false);

viz = false;
lines(1).point1 = [1 1];
lines(1).point2 = [1 1];
lines(1).rotation = 1;
line_saved = false;
for i=1:size(kernels,3)
    kernel = kernels(:,:,i);
    
    %Apply kernel
    filtered = imfilter(img_filtered, kernel);
    
    %Top-hat
    SE = kernel_to_strel(kernel);
    top_hat = imtophat(filtered, SE);
    top_hat = imadjust(top_hat);
    
    %Threshold and median filter
    img_bin = top_hat > 127;
    img_bin = medfilt2(img_bin,[4 4]);
    
    %Remove white points from edges of the frame
    %Stupid workaround?????????
    img_bin(:,1:10)=0;
    img_bin(:,width-10:end)=0;
    img_bin(1:10,:)=0;
    img_bin(height-10:end,:)=0;
    
    %---------Save lines and kernel rotation---------
    hough_lines = test_hough(img_bin, img_org, false);
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
    
    if (viz)
        figure(1);
        imshow(mat2gray(kernel),'InitialMagnification','fit');
        figure(2);
        imshow(img_bin, 'InitialMagnification', 500);
        %axis on, axis normal, hold on;
        %plot(y,x,'s','color','white');
        pause(1);
    end
end

intensities = [lines(1:end).intensity];
[dummy,ind] = maxk(intensities, floor(length(lines)*0.5));
intense_lines = lines(ind);

%figure(45);
%imshow(img_org); hold on;
%Analyze line with diagonal kernel
for i = 1:length(intense_lines)
    x = [intense_lines(i).point1(1) intense_lines(i).point2(1)];
    y = [intense_lines(i).point1(2) intense_lines(i).point2(2)];
    
    kernel = kernels(:,:,intense_lines(i).rotation);
    %diagonal_1 = diagonal_kernel(kernels, intense_lines(i).rotation);
    diagonal_degrees = rotation_to_degrees(kernels, intense_lines(i).rotation) + 90;
    diagonal = degrees_to_kernel(kernels_big, diagonal_degrees);
    
    %figure(1);
    %imshow(mat2gray(kernel),'InitialMagnification','fit');
    %figure(2);
    %imshow(mat2gray(diagonal),'InitialMagnification','fit');
    
    
    test_line = kernel_to_line(diagonal);
    center_line = ceil(size(test_line, 1)/2);
    intense_lines(i).vals = [];
    
    [cx, cy, c] = improfile(img_org, x, y);
    cx = round(cx);
    cy = round(cy);
    for j = 1:length(cx)
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
                intense_lines(i).Rvals(j,cnt) = img_org(cy(j)-center_line+ci, cx(j)-center_line+cj,1);
                intense_lines(i).Gvals(j,cnt) = img_org(cy(j)-center_line+ci, cx(j)-center_line+cj,2);
                intense_lines(i).Bvals(j,cnt) = img_org(cy(j)-center_line+ci, cx(j)-center_line+cj,3);
                %img_red(cy(j)-center_line+ci, cx(j)-center_line+cj) = 255;
                cnt = cnt + 1;
            end
        end
        %figure(475);
        %imshow(img_red);
    end
    
    %     xy = [x(1) y(1); x(2) y(2)];
    %     plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');
    %     % Plot beginnings and ends of lines
    %     plot(xy(1,1),xy(1,2),'x','LineWidth',2,'Color','yellow');
    %     plot(xy(2,1),xy(2,2),'x','LineWidth',2,'Color','red');
    
end
%        figure(475);
%imshow(img_red);

% x = [pt1(1) pt2(1)];
% y = [pt1(2) pt2(2)];
% [cx, cy, c] = improfile(img_red, x, y);
% cx = round(cx);
% cy = round(cy);
% test_line = kernel_to_line(diagonal);
% for i = 1:length(cx)
%     x = cx(i);
%     y = cy(i);
%     center_line = ceil(size(test_line, 1)/2);
%     
%     for ci = 1:size(test_line, 1)
%         for cj = 1:size(test_line, 2)
%             if (test_line(ci, cj) == 0)
%                 continue;
%             end
%             img_red(y-center_line+ci, x-center_line+cj) = 0;
%         end
%     end
% end

figure(99);
imshow(img_org); hold on;
max_bright = 0;
max_intensity = 0;
best_kernel = 0;
pt1 = []; pt2 = [];
rf_diff = 0;
for i = 1:length(intense_lines)
    x = [intense_lines(i).point1(1) intense_lines(i).point2(1)];
    y = [intense_lines(i).point1(2) intense_lines(i).point2(2)];
    
    Rvals = intense_lines(i).Rvals;
    Gvals = intense_lines(i).Gvals;
    Bvals = intense_lines(i).Bvals;
    Rvals_avg = mean(Rvals,1);
    Gvals_avg = mean(Gvals,1);
    Bvals_avg = mean(Bvals,1);
    
    %ac = Rvals+(Gvals+Bvals);
    ac = Rvals;
    ac = mean(ac,1);
    
    
    edge_len = fix(length(Rvals_avg)/3);
    center_len = mod(length(Rvals_avg),3);
    risingR = mean(Rvals(1:edge_len));
    centerR = mean(Rvals(edge_len+1:edge_len+edge_len+center_len));
    fallingR = mean(Rvals(edge_len+edge_len+center_len+1:end));
    
    risingG = mean(Gvals(1:edge_len));
    centerG = mean(Gvals(edge_len+1:edge_len+edge_len+center_len));
    fallingG = mean(Gvals(edge_len+edge_len+center_len+1:end));
    
    risingB = mean(Bvals(1:edge_len));
    centerB = mean(Bvals(edge_len+1:edge_len+edge_len+center_len));
    fallingB = mean(Bvals(edge_len+edge_len+center_len+1:end));
    
    a = diff(Rvals_avg);
    %xn = [min(Rvals_avg):(max(Rvals_avg)-min(Rvals_avg))/(length(a)-1):max(Rvals_avg)];
    xn = [-length(a)/2:length(a)/2];
    yn = normpdf(xn,0,1);
    %figure(1); hold on;
    %plot(xn,yn);
    
    b = corrcoef(ac, yn);
    cof = abs(b(2,1));
    %plot(a)
    
    %if (center < falling && center > rising)
    %    pt1 = intense_lines(i).point1;
    %    pt2 = intense_lines(i).point2;
    %end
    
    %avg_vals = mean(vals(:));
    pixels = improfile(img_org, x, y);
    
    R = mean(pixels(:,:,1));
    G = mean(pixels(:,:,2));
    B = mean(pixels(:,:,3));
    avg = R/(G+B);
    %avg = mean(pixels);
    
    rise_avg = 0;
    fall_avg = 0;
    if (risingG + risingB == 0)
        rise_avg = risingR;
        fall_avg = fallingR;
    else
        rise_avg = risingR-((risingG + risingB)/2);
        fall_avg = fallingR-((fallingG + fallingB)/2);
    end
    center_avg = (centerR + centerG + centerB)/3;
    coef = rise_avg + fall_avg + center_avg;
    rise_avg
    fall_avg
    center_avg
    %if (coef == inf)
     %  continue;
    %end
    
%     %if (rising < center && falling < center)
%     if (cof > 0.7)
%         pt1 = intense_lines(i).point1;
%         pt2 = intense_lines(i).point2;
%         
%         
%         xy = [pt1; pt2];
%         plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');
%         % Plot beginnings and ends of lines
%         plot(xy(1,1),xy(1,2),'x','LineWidth',2,'Color','yellow');
%         plot(xy(2,1),xy(2,2),'x','LineWidth',2,'Color','red');
%     end
    
    if (coef > max_bright)
        max_bright = coef;
        best_kernel = intense_lines(i).rotation;
        pt1 = intense_lines(i).point1;
        pt2 = intense_lines(i).point2;
        rf_diff = rise_avg - fall_avg;
    end
end
max_bright
rf_diff
lines=[];

%kernel = kernels(:,:,best_kernel);
%figure(1);
%imshow(mat2gray(kernel),'InitialMagnification','fit');
%filtered = imfilter(img_filtered, kernel);
% figure(99);
% imshow(img_org); hold on;
% 
 xy = [pt1; pt2];
 plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','magenta');
 % Plot beginnings and ends of lines
 plot(xy(1,1),xy(1,2),'x','LineWidth',2,'Color','yellow');
 plot(xy(2,1),xy(2,2),'x','LineWidth',2,'Color','red');

%Calculating detected line angle
x1=pt1(1);
y1=pt1(2);
x2=pt2(1);
y2=pt2(2);
slope = (y2 - y1) ./ (x2 - x1);
angle = atand(slope);
if (angle < 0)
   angle = 180 + angle; 
end

%Calculate kernel rotation angle resolution
rot_res = 180/kernels_cnt;
tmp = 0;
corrected_kernel = 0;
for i=1:kernels_cnt
    tmp = tmp + rot_res;
    if (tmp > angle)
        corrected_kernel = i;
        break;
    end
end

kernel = kernels(:,:,corrected_kernel);
SE = kernel_to_strel(kernel);
SE_neigh = uint8(SE.Neighborhood);
SE_neigh(SE_neigh == 0) = 254;
%img_red(101:100+size(SE_neigh,1), 101:100+size(SE_neigh,1)) = SE_neigh;

%figure(80);
%imshow(mat2gray(SE.Neighborhood),'InitialMagnification','fit');
%figure(88);
%imshow(mat2gray(diagonal),'InitialMagnification','fit');

diagonal = diagonal_kernel(kernels, corrected_kernel);
SE = kernel_to_strel(diagonal);
SE_neigh = uint8(SE.Neighborhood);
SE_neigh(SE_neigh == 0) = 255;

%VIDEO%
%end
%VIDEO%
end

function ret = test_hough(img_bin, img_org, viz)
[H,T,R] = hough(img_bin,'RhoResolution',1,'Theta',-90:180/62:89);
P  = houghpeaks(H,30);
lines = houghlines(img_bin,T,R,P,'FillGap',20,'MinLength',100);

for i = 1:length(lines)
    theta = lines(i).theta;
    rho = lines(i).rho;
    T_theta = T==theta;
    R_rho = R==rho;
    hough_intensity = H(R_rho, T_theta);
    lines(i).intensity = hough_intensity;
end

if viz == true
    figure(30);
    imshow(img_org); hold on;
    for k = 1:length(lines)
        xy = [lines(k).point1; lines(k).point2];
        plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');
        
        % Plot beginnings and ends of lines
        plot(xy(1,1),xy(1,2),'x','LineWidth',2,'Color','yellow');
        plot(xy(2,1),xy(2,2),'x','LineWidth',2,'Color','red');
    end
    figure(10);
    imshow(H,[],'XData',T,'YData',R,'InitialMagnification','fit');
    axis on, axis normal, hold on;
    plot(T(P(:,2)),R(P(:,1)),'s','color','white');
end

ret=lines;
end


function ret = kernel_rotate(kernel, visualize)
dim = size(kernel,1);
center = round(dim/2);
rot = center;

kernel_org = kernel;

%Draw frame of zeros
kernel=conv2(kernel,[0,0,0;0,1,0;0,0,0]);
hor = kernel;
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

hor_inv = [];
%45 degrees counter-clockwise from horizontal
kernel=conv2(kernel_org,[0,0,0;0,1,0;0,0,0]);
k=1;
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
%45 degrees clockwise from vertical
k=1;
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

kernels = hor;
kernels(:,:,end+1:end+size(hor_rot,3)) = hor_rot;
kernels(:,:,end+1:end+size(ver_inv,3)) = flip(ver_inv,3);
kernels(:,:,end+1:end+size(ver,3)) = ver;
kernels(:,:,end+1:end+size(ver_rot,3)) = ver_rot;
kernels(:,:,end+1:end+size(hor_inv,3)) = flip(hor_inv,3);
ret = kernels;

%---------------VIZUALISATION---------------
if (visualize == true)
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

function kernel = generate_kernel(kernel_size)
kernel = zeros(kernel_size,kernel_size);
center = ceil(kernel_size/2);
kernel(center,:) = 2;
kernel(center - 1,:) = -1;
kernel(center + 1,:) = -1;
end

function se = kernel_to_strel(kernel)
se = kernel;
%se = (kernel == 2);
se(1,:) = [];
se(end,:) = [];
se(:,1) = [];
se(:,end) = [];
se = strel(se);
end

function ret = kernel_to_line(kernel)
ret = kernel;
ret(ret == -1) = 0;
ret(ret == 2) = 1;
ret(1,:) = [];
ret(end,:) = [];
ret(:,1) = [];
ret(:,end) = [];
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

function dist = euclidean_distance(point1, point2)
dist = sqrt((point1(1) - point2(1))^2 + (point1(2) - point2(2))^2);
end

function ret = diagonal_kernel(kernels, rotation)
    kernels_cnt = size(kernels,3);
    diagonal = kernels_cnt/2 + rotation;
    if (diagonal > kernels_cnt)
       diagonal = diagonal - kernels_cnt;
    end
    ret = kernels(:,:,diagonal);
end

function degrees = rotation_to_degrees(kernels, rotation)
kernels_cnt = size(kernels,3);
rot_res = 180/kernels_cnt; %How many degrees per rotation
degrees_arr = 0:rot_res:180-rot_res;
degrees = degrees_arr(rotation);
end

function kernel = degrees_to_kernel(kernels, degrees)
if (degrees > 180)
   degrees = degrees - 180; 
end
kernels_cnt = size(kernels,3);
rot_res = 180/kernels_cnt; %How many degrees per rotation
degrees_arr = 0:rot_res:180-rot_res;
[~,degrees_ind] = min(abs(degrees_arr - degrees));
kernel = kernels(:,:,degrees_ind);
end









