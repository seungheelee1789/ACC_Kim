% Pupil Size Analysis

% You need the following source codes (.m / .mexw64 / .cpp) to run this code.
% fitellipse.m plotellipse.m

close all;clear;clc;

v2 = VideoReader('source_video.avi');
SR = 29.896;

%%%%% draw rectangle for ROI
for i=1000 
    img(:,:,1) = rgb2gray(read(v2,i));    
end

figure(1)
imshow(img*3)
% eye
title('draw total ROI')
roi = drawrectangle
x_coor = roi.Vertices(:,1); y_coor = roi.Vertices(:,2);
% inside pupil
title('draw inside pupil')
roi2 = drawrectangle
x_coor2 = roi2.Vertices(:,1); y_coor2 = roi2.Vertices(:,2);
% outside pupil
title('draw outside pupil')
roi3 = drawrectangle
x_coor3 = roi3.Vertices(:,1); y_coor3 = roi3.Vertices(:,2);

for i=1:50
    disp(i)
    clear currimg img_ori img_ori2 img_ori3
    currimg = read(v2,i*30);
    img_ori = 255-currimg(round(y_coor(1)):round(y_coor(2)),round(x_coor(1)):round(x_coor(4)));  
    img_ori2 = 255-currimg(round(y_coor2(1)):round(y_coor2(2)),round(x_coor2(1)):round(x_coor2(4)));  
    img_ori3 = 255-currimg(round(y_coor3(1)):round(y_coor3(2)),round(x_coor3(1)):round(x_coor3(4))); 
    
    pupil_inout(i,1) = mean(mean(img_ori2));
    pupil_inout(i,2) = mean(mean(img_ori3));   
    
    pupil_ROI(i,:) = reshape(img_ori,1,size(img_ori,1)*size(img_ori,2));
end

% plot pupil-in intensity vs. pupil-out intensity
sz = 5;
figure()
subplot(2,2,1)
scatter(pupil_inout(:,1),pupil_inout(:,2),sz,'k')
xlabel('pupil-in')
ylabel('pupil-out')
limcri2 = ceil(max(pupil_inout,[],'all'));
limcri1 = floor(min(pupil_inout,[],'all'));
xlim([limcri1 limcri2])
ylim([limcri1 limcri2])
hold on
line([limcri1 limcri2], [limcri1 limcri2], 'color', 'r')
pupil_ROI = double(pupil_ROI);

thresholdindex(1) = mean(pupil_inout(:,1));                 % mean pupil_in intensity
thresholdindex(2) = mean(pupil_inout(:,2));                 % mean pupil_out intensity
threshold = (2*thresholdindex(1)+thresholdindex(2))/3;        % threshold boundary

pupil_ROI_re = pupil_ROI - threshold + 128;
pupil_ROI_re(find(pupil_ROI_re < 0)) = 0;
img_ori_re = img_ori - threshold + 128;
img_ori_re(find(img_ori_re < 0)) = 0;

threshold_re = thresholdindex - threshold + 128;

clear histodata
limcri1 = 90;
limcri2 = 150;
edges = [limcri1:1:limcri2];
edges2 = [limcri1:10:limcri2];
borderline1 = floor(threshold_re(2)) - limcri1;
borderline2 = ceil(threshold_re(1)) - limcri1;


for i=1:50
    clear aa
    histodata(i,:) = 100*histcounts(pupil_ROI_re(i,:),edges)/size(pupil_ROI_re,2);
    cri = min(histodata(i,borderline1:borderline2));
    boundary(i,1) = min(find(histodata(i,borderline1:borderline2) == cri)+borderline1-1);    
end

newcri = mean(boundary);

subplot(2,2,3)
imagesc((histodata))
hold on
line([threshold_re(1)-limcri1 threshold_re(1)-limcri1], [0 100],'color','c','LineWidth',2)
hold on
line([threshold_re(2)-limcri1 threshold_re(2)-limcri1], [0 100],'color','c','LineWidth',2)
hold on
line([newcri newcri], [0 100],'color','w','LineWidth',3,'LineStyle',':')
xticks([1:10:length(edges)])
xticklabels(edges2)
colormap(gca,'hot')  
colorbar
xlabel('pixel-intensity')
ylabel('frame')

subplot(1,2,2)
b = bar3((histodata));
for i=1:length(b)
    zdata = b(i).ZData;
    b(i).CData = zdata;
    b(i).FaceColor = 'interp';
    b(i).FaceAlpha = 0.5;
    b(i).EdgeColor = 'k';
    b(i).EdgeAlpha = 0.5;
    b(i).LineWidth = 0.0001;
end
xticks([1:10:length(edges)])
xticklabels(edges2)
ylim([0 100])
zlim([0.1 inf])
caxis([0 3])
colormap(gca,'jet')  
colorbar;
xlabel('pixel-intensity')
ylabel('frame')
zlabel('number of pixel')
view([5 25])    
point1 = [newcri, 100, 0];
point2 = [newcri, 100, 4];
point3 = [newcri, 0, 4];
point4 = [newcri, 0, 0];
points = [point1' point2' point3' point4'];
hold on
h = fill3(points(1,:),points(2,:),points(3,:),'r');
h.FaceAlpha = 0.6;
h.EdgeColor = 'none';
grid off


prelength = 0;          % start point (based on frame)
%targetlength = 30000;    % length of analysis (based on frame)
targetlength = v2.Duration * SR; 
pixelno = targetlength - prelength;

tic
smoothingfactor = 10;
clear pupilsize
fd = 0; 
frame_drop = []; 

for i1=prelength+1:targetlength
    clear currimg img_ori img
    currstate = [num2str((i1-prelength)/pixelno*100) '%'];
    disp(currstate)
    currimg = read(v2,i1);
    img_ori = currimg(round(y_coor(1)):round(y_coor(2)),round(x_coor(1)):round(x_coor(4)));  
    
    %%%%%% original image ####################   
    img_ori_pre = img_ori;
    
    %%%%%% binarization based on threshold ###  
    img = 255-img_ori-threshold;   
    img(img<0)=0;     img(img>0)=1;
    
    %%%%%% hole filling ######################
    HF1 = imfill(img,'holes');
    
    %%%%%% major regionprop and convex hull ##
    CH0 = bwconvhull(HF1,'object',4);
    clear a aa aprime centroid PixelList centroid_posi regionprops_data main_area areasize
    a = regionprops(CH0,'area','centroid','PixelList');
    regionprops_data{i1,1} = a;
    if size(a,1) > 0
        if size(a,1) > 1
            aprime = struct2cell(a);
            for jj=1:size(a,1)
                areasize(jj,1) = aprime{1,jj};
                centroid_posi(jj,1) = aprime{2,jj}(1,1);
                centroid_posi(jj,2) = aprime{2,jj}(1,2);
                PixelList{jj,1}(1,:) = aprime{3,jj}(:,1);
                PixelList{jj,1}(2,:) = aprime{3,jj}(:,2);
            end
            main_area = find(areasize==max(areasize));
            for jj=1:size(a,1)
                if jj ~= main_area
                    if centroid_posi(jj,1) < centroid_posi(main_area,1) & centroid_posi(jj,2) > centroid_posi(main_area,2)                        
                        % converting 0 to 1 at non-main area
                        CH0(PixelList{jj,1}(2,:),PixelList{jj,1}(1,:)) = 0;                        
                    end
                end
            end
            %aa = regionprops(CH0,'area','centroid','PixelList');         
            %aa = bwconvhull(CH0,4);
            aa = bwconvhull(CH0,'union',4);
        end  
    end     
    
    %%%%%% 1st gaussian smoothing ############
    CH0_double = double(CH0);
    GS1 = imgaussfilt(CH0_double,smoothingfactor);
    
    %%%%%% 1st binarization ##################
    BN1 = imbinarize(GS1);  

    %%%%%% 1st convex hull ###################
    CH1 = bwconvhull(BN1,'objects',4); 
    
    %%%%%% convex hull making ################
    clear a aa aprime centroid PixelList CHE1_2 CHE1_3
    a = regionprops(CH1,'area','centroid','PixelList');
    
    if size(a,1) > 1
        aprime = struct2cell(a);
        for jj=1:size(a,1)
            areasize(jj,1) = aprime{1,jj};
            PixelList{jj,1}(1,:) = aprime{3,jj}(:,1);
            PixelList{jj,1}(2,:) = aprime{3,jj}(:,2);
        end
        main_area = find(areasize==max(areasize));
        for jj=1:size(a,1)
            if jj ~= main_area
                % converting 0 to 1 at non-main area
                CH1(PixelList{jj,1}(2,:),PixelList{jj,1}(1,:)) = 0;                        
            end
        end
        clear a
        a = regionprops(CH1,'area','centroid','PixelList');
    end 
    
    if size(a,1) == 0
        x_posi(i1, 1) = x_posi(i1-1, 1);
        y_posi(i1, 1) = y_posi(i1-1, 1);
        CHE1_1{i1,1} = CHE1_1{i1-1,1};
        CHE1_2 = find(CHE1_1{i1-1,1} == 1);
        fd = fd + 1; 
        frame_drop(fd, 1) = i1; 
        
    else
        x_posi(i1, 1) = a.Centroid(1); % center of mass (x)
        y_posi(i1, 1) = a.Centroid(2); % center of mass (y)
        CHE1_1{i1,1} = edge(CH1);
        CHE1_2 = find(CHE1_1{i1,1} == 1);
    end
    
    pupilsize(i1,1) = length(find(CH1 == 1));
    for z=1:length(CHE1_2)
        CHE1_3(z,1) = floor(CHE1_2(z)/size(CH1,1));
        CHE1_3(z,2) = rem(CHE1_2(z),size(CH1,1));
    end
    
    x = CHE1_3';
    % Find the least squares geometric estimate 
    % (ellipse fitting method 1)
    [zg, ag, bg, alphag] = fitellipse(x);
%     Plot the results
%     hold on
%     plotellipse(zg, ag, bg, alphag, 'c')
    
    % pupil size by ellipse fit
    pupilsize(i1,2) = pi*ag*bg;    
end
toc

save pupilsize.mat pupilsize

if length(frame_drop) > 0 
    save frame_drop.mat frame_drop
end

figure()
subplot(2,1,1)
plot(pupilsize(:,1),'k')
hold on
title('pupil size (original)')
xlabel('frames (time)')
ylabel('pupil size (# of pixels)')
subplot(2,1,2)
title('pupil size (ellipse fit)')
plot(pupilsize(:,2),'r')
hold on
title('pupil size (original)')
xlabel('frames (time)')
ylabel('pupil size (# of pixels)')
disp("done!")
