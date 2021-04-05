% network activity analysis (Fig. 3b)
% ACC Sinc neurons (n = 56)

clear;clc;
load datasource.mat

% data_source{:,1}: background trials
% data_source{:,2}: PL trials (stim onset)
% data_source{:,3}: Miss trials
% data_source{:,4}: PL trials (lick onset)
% data_source{:,5}: SL trials

binsize = 50;
edges = [1:binsize:4001];
pval = 0.05;
ylimcri1 = -0.2;
ylimcri2 = 0.4;
xlimcri = 1;
colorcode =[0 1 0; 0 0 0; 1 0 1; 139/255 69/255 19/255];
colorcode2 =[0.8 1 0.8; 0.8 0.8 0.8; 1 0.8 1; 195/255 155/255 119/255];
tltledata={'PL(stim)';'Miss(stim)';'PL(lick)';'SL(stim)'};

% pre 2s, post 2s PSTH making
for i=1:size(data_source,1)
    spk = spike_time{i,1};
    for k=1:5
        tempmat = zeros(length(data_source{i,k}),4001);
        for z=1:length(data_source{i,k})
            tempmat(z,spk(find(spk <= (data_source{i,k}(z) + 2000) & spk > (data_source{i,k}(z) - 2000))) - data_source{i,k}(z) + 2000) = 1;
        end
        heatmapdata{i,k} = tempmat;
        clear tempmat
    end
    clear spk
end

% data averaging
for i=1:size(heatmapdata,1)
    for k=1:5
        heatmapdata1{1,k}(i,:) = mean(heatmapdata{i,k},1);
    end
end

% data binning and smoothing
for i=1:size(heatmapdata,1)
    for k=1:5
        for z=1:length(edges)-1
            heatmapdata1{2,k}(i,z) = sum(heatmapdata1{1,k}(i,edges(z):edges(z+1)-1));
        end
        heatmapdata1{3,k}(i,:) = (1000/binsize)*smoothdata(heatmapdata1{2,k}(i,:),'gaussian',200/binsize);
    end
end

for i=1:size(heatmapdata,1)
    for k=1:5
        for z=1:length(edges)-1
            heatmapdata1{4,k}(i,z) = (heatmapdata1{3,k}(i,z) - mean(heatmapdata1{3,1}(i,:))) / (heatmapdata1{3,k}(i,z) + mean(heatmapdata1{3,1}(i,:)));
        end
    end
end

for k=1:5
    for i=1:size(heatmapdata,1)
        heatmapdata_re{2,k}(i,1) = mean(heatmapdata1{4,k}((i),1:(1000/binsize)*2));
        heatmapdata_re{4,k}(i,:) = heatmapdata1{4,k}((i),:);
    end
    [bb, index] = sortrows(heatmapdata_re{2,k}(:,1));
    heatmapdata_re{2,k}(:,2) = index;      
    heatmapdata_re{5,k} = heatmapdata_re{4,k};      
end

for k=1:5
    for i=1:size(heatmapdata,1)
        for z=1:length(edges)-1
            if heatmapdata_re{5,k}(i,z) == -1
                heatmapdata_re{5,k}(i,z) = 0;
            end
        end
    end
end

figure()
set(gcf,'Position',[150 150 800 250])
sgtitle('ACC Sinc neurons (n=56)')
curr = 0;
for k=1:4
    clear psth1 psth2
    curr = curr + 1;
    psth1(1,:) = mean(heatmapdata_re{5,k+1},1);
    psth2(1,:) = std(heatmapdata_re{5,k+1},1)/sqrt(size(heatmapdata,1));
    psth11(1,:) = mean(heatmapdata_re{5,1},1);
    psth12(1,:) = std(heatmapdata_re{5,1},1)/sqrt(size(heatmapdata,1));
    tline = [-2:binsize/1000:2-binsize/1000];
    
    subplot(4,4,[4+k,8+k,12+k])    
    hold on
    errorshade(tline, psth11, psth11+psth12, psth11-psth12, [.8 .8 .8])
    hold on
    plot(tline, psth11,'color', [.6 .6 .6], 'LineWidth', 2,'LineStyle',':')
    hold on
    errorshade(tline, psth1, psth1+psth2, psth1-psth2, colorcode2(k,:))
    hold on
    plot(tline, psth1,'color', colorcode(k,:), 'LineWidth', 2)
    hold on
    title(tltledata(k))
    hold on
    line([0 0], [ylimcri1 ylimcri2],'color',colorcode(k,:),'LineWidth', 2,'LineStyle','--')
    xlim([-xlimcri xlimcri])
    ylim([ylimcri1 ylimcri2])
    xlabel('time (s)')
    ylabel('normalized rate')
end

curr = 0;
for i=1:4
    curr = curr + 1;
    ROI = heatmapdata_re{5,i+1};
    ROC = heatmapdata_re{5,1};
    for k=1:size(ROI,2)
        a1 = ROI(:,k)';
        b1 = ROC(:,k)';
        p = signrank(a1,b1);
        aa(curr,k) = p;
    end
    clear ROI ROC ROC2
    subplot(4,4,i)
    imagesc(aa(curr,20:60))
    colormap(gray)
    caxis([0 pval])
    axis off
end
