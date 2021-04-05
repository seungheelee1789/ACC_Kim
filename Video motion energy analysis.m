% Video Motion Energy Analysis 

close all; clear;clc;

v2 = VideoReader('source_video.avi');
SR = 29.896;


%%%%% draw rectangle for ROI
for i=1 
    img(:,:,1) = 2*rgb2gray(read(v2,i));    
end

timecri = v2.NumFrames-1000;

pre = 0;

figure(1)
imshow(img)
hold on
title('draw nose movement ROI')
roi1 = drawrectangle;
x_coor{1,1} = roi1.Vertices(:,1); y_coor{1,1} = roi1.Vertices(:,2);

figure(2)
imshow(img)
hold on
title('draw licking ROI')
roi2 = drawrectangle;
x_coor{2,1} = roi2.Vertices(:,1); y_coor{2,1} = roi2.Vertices(:,2);

figure(3)
imshow(img)
hold on
title('draw whisking ROI')
roi3 = drawrectangle;
x_coor{3,1} = roi3.Vertices(:,1); y_coor{3,1} = roi3.Vertices(:,2);

% making matrix of ROI
tic
for i=1:timecri
    disp(horzcat(num2str(i*100/timecri),'%'))
    clear currimg
    currimg = read(v2,pre+i);
    for j=1:size(x_coor,1)
        ROI{j,1}(:,:,i) = currimg(round(y_coor{j,1}(1)):round(y_coor{j,1}(2)),round(x_coor{j,1}(1)):round(x_coor{j,1}(4)));   
    end
end
toc

% video motion energy calculation
for i=1:size(ROI,1)
    diffROI{i,1}(:,:,:) = abs(diff(ROI{i,1},1,3));
end

% data averaging
for i=1:size(ROI,1)
    tempmean(:,:) = mean(diffROI{i,1},2);
    mean_diffROI(i,:) = mean(tempmean,1);
    clear tempmean
end

figure()
imshow(img)
hold on
for k=1:size(x_coor,1)
    for j=1:4
        if j ~= 4
            hold on
            plot(x_coor{k,1}(j:j+1),y_coor{k,1}(j:j+1),'r','LineWidth',2)
        elseif j == 4
            hold on
            line([x_coor{k,1}(1) x_coor{k,1}(j)], [y_coor{k,1}(1) y_coor{k,1}(j)],'color','r','LineWidth',2)
        end
    end
end

datalabel = {'nose movement','licking','whisking'};

figure()
for i=1:3
    subplot(3,1,i)
    hold on
    plot(mean_diffROI(i,:))
    hold on
    title(datalabel(i))
    ylabel('VME')
    xlabel('frame')
end

disp("done!")