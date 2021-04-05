% spike raster plot and PSTH/PETH (Fig. 1s)

clear;clc;
% example dataset
load ACC_MDec_example.mat

% input arguments
sz = 3;
binsize = 50;
binsize2 = 1000/binsize;
edges = [1:binsize:6001];

% excluding lick contamination trials
deltatime(:,1) = event2{1,3} - event2{1,5};
filtered = find(deltatime > 500);

for k=1:5
    clear event
    if k <= 2
        % spontaneous lick onset(1), offset(2)
        PSTH{1,k} = zeros(length(event2{1,k}),6001);    % zero matrix
        event = event2{1,k};    % event time
    else
        % perceptual lick onset(3), offset(4), stim onset(5)
        PSTH{1,k} = zeros(length(filtered),6001);   % zero matrix
        event = event2{1,k}(filtered);  % event time (wo lick contamination trials)
    end    
    for i=1:length(event)
        % making binary matrix
        PSTH{1,k}(i,:) = spk(event(i)-2000:event(i)+4000);    
    end
    event3{1,k} = event;
end

% spontaneous lick (offset - onset)
cri1 = -1*(event3{1,2} - event3{1,1});

% perceptual lick (stim onset - lick onset)
cri2 = -1*(event3{1,3} - event3{1,5});

% perceptual lick (lick offset - lick onset)
cri3 = -1*(event3{1,4} - event3{1,3});

[kk1 index1] = sortrows(cri1');
[kk2 index2] = sortrows(cri2');
[kk3 index3] = sortrows(cri3');
kk1 = -1*kk1;
kk2 = -1*kk2;
kk3 = -1*kk3;

% spontaneous lick onset (trial alignment according to lick onset-offset)
for i=1:size(PSTH{1,1},1)
	PSTH{2,1}(i,:) = PSTH{1,1}(index1(i),:);
end
% spontaneous lick offset (trial alignment according to lick onset-offset)
for i=1:size(PSTH{1,2},1)
	PSTH{2,2}(i,:) = PSTH{1,2}(index1(i),:);
end
% perceptual lick onset (trial alignment according to lick onset-offset)
for i=1:size(PSTH{1,3},1)
	PSTH{2,3}(i,:) = PSTH{1,3}(index3(i),:);
end
% perceptual lick offset (trial alignment according to lick onset-offset)
for i=1:size(PSTH{1,4},1)
	PSTH{2,4}(i,:) = PSTH{1,4}(index3(i),:);
end
% perceptual stim onset (trial alignment according to stim onset-lick onset)
for i=1:size(PSTH{1,5},1)
	PSTH{2,5}(i,:) = PSTH{1,5}(index2(i),:);
end

% PSTH/PETH binning and smoothing
for i=1:5
    for j=1:size(PSTH{1,i},1)
        for k=1:length(edges)-1
            PSTH{3,i}(j,k) = sum(PSTH{1,i}(j,edges(k):edges(k+1)-1));
        end
    end
    % binned PSTH/PETH (Mean)
    PSTH2(i,:) = binsize2*mean(PSTH{3,i},1);
    % binned PSTH/PETH (SEM)
    PSTH3(i,:) = binsize2*std(PSTH{3,i},1)/sqrt(size(PSTH{3,i},1));
    % binned & smoothed PSTH/PETH (Mean)
    PSTH4(i,:) = smoothdata(PSTH2(i,:),'gaussian',400/binsize);
    % binned & smoothed PSTH/PETH (SEM)
    PSTH5(i,:) = smoothdata(PSTH3(i,:),'gaussian',400/binsize);
end

plottarget = [1 3]; % 1: spontaneous lick onset, 3: perceptual lick onset
browncode = [139 69 19]/255;
orangecode = [255 140 0]/255;

% data plot
figure()
timescale = [-2:binsize/1000:4-binsize/1000];
for i=1:2
    subplot(2,2,i+2)
    clear aa bb
    aa = PSTH4(plottarget(i),:);
    bb = PSTH5(plottarget(i),:);
    errorshade(timescale,aa,aa+bb,aa-bb,[.9 .9 .9])
    hold on
    plot(timescale,PSTH4(plottarget(i),:),'k')
    xlim([-2 4])    
    ylim([0 4])
    xlabel('time (s)')
    ylabel('spikes/s')
    if i == 1
        line([0 0], [0 4], 'Color',browncode,'LineWidth',2)
    elseif i == 2
        line([0 0], [0 4], 'Color','m','LineWidth',2)
    end
    
    subplot(2,2,i)
    for j=1:size(PSTH{2,plottarget(i)},1)
        clear tempscatter tempscatter2
        tempscatter = find(PSTH{2,plottarget(i)}(j,:) > 0)/1000;
        if length(tempscatter) > 0
            tempscatter2 = j*ones(1,length(tempscatter));
            hold on
            scatter(tempscatter,tempscatter2,sz,'k','filled')            
            hold on
            if i == 1
                scatter(2,j,sz,browncode)
                hold on
                scatter(2+kk1(j)/1000,j,sz,orangecode)
            elseif i == 2
                scatter(2,j,sz,'m')
                hold on
                scatter(2+kk3(j)/1000,j,sz,orangecode)
            end
        end
    end
    xlim([0 6])
    ylim([0 j])
    axis off
    if i == 1
        title('spontaneous lick onset')
    elseif i == 2
        title('perceptual lick onset')
    end
end

    

        