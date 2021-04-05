% hit vs. miss comparison heatmap (Fig. 2a)
% sorting by visaully-evoked responsiveness (hit-based)

clear; clc;
load ACC_sensory_neurons.mat

% input arguments
binsize = 50;
binsize2 = 1000/binsize;
edges = [1:binsize:6001];

% matching the number of trials from hit and miss trials
for i=1:size(heatmapdata,1)    
    clear refcri
    refcri = min(size(heatmapdata{i,1},1),size(heatmapdata{i,2},1));    
    if refcri == size(heatmapdata{i,1},1)
        % # of hit < # of miss
        heatmapdata_re{i,1} = heatmapdata{i,1}(1:refcri,:);
        heatmapdata_re{i,2} = heatmapdata{i,2}(end-refcri+1:end,:);
    else
        % # of hit > # of miss
        heatmapdata_re{i,1} = heatmapdata{i,1}(end-refcri+1:end,:);
        heatmapdata_re{i,2} = heatmapdata{i,2}(1:refcri,:);
    end        
end
%
for i=1:size(heatmapdata,1)
    for k=1:2
        % all trials
        heatmapdata1{1,k}(i,:) = mean(heatmapdata{i,k},1);
        % trial number matched
        heatmapdata1_re{1,k}(i,:) = mean(heatmapdata_re{i,k},1);
    end
end

% PSTH binning & z-scoring
for i=1:size(heatmapdata,1)
    for k=1:2
        % PSTH binning
        for z=1:length(edges)-1
            final_re{1,k}(i,z) = binsize2*sum(heatmapdata1{1,k}(i,edges(z):edges(z+1)-1));
            final_re{2,k}(i,z) = binsize2*sum(heatmapdata1_re{1,k}(i,edges(z):edges(z+1)-1));
        end
        % z-scoring
        for j=1:2
            clear baseline1 baseline2
            baseline1 = mean(final_re{j,k}(i,1:binsize2*2));
            baseline2 = std(final_re{j,k}(i,1:binsize2*2));
            if baseline2 == 0
                baseline2 = 1;
            end
            for z=1:length(edges)-1
                final_z{j,k}(i,z) = (final_re{j,k}(i,z) - baseline1) / baseline2;
            end
            z_index(i,2*(j-1)+k) = mean(final_z{j,k}(i,binsize2*2+1:binsize2*2.5));
        end
    end
end

% sorting order of neurons according to visual responsiveness during hit trials
for i=1:size(heatmapdata,1)
    sortingindex(i,1) = -1*mean(final_z{2,1}(i,binsize2*2+1:binsize2*2.5));
end
[aa index] = sortrows(sortingindex);

for j=1:2
    for k=1:2
        for i=1:size(heatmapdata,1)
            final_sorted{j,k}(i,:) = final_z{j,k}(index(i),:);
        end
    end
end

timescale = [-2:1/binsize2:4];
timescale2 = [-2:1:4];
datalabelset = {'Hit (z)','Miss (z)'};
datalabelset2 = {'- all trials','- trial no. matched'};

figure()
set(gcf,'Position',[150 150 600 500])
sz = 3;
limcri = 10;
limcri2 = 3;
for i=1:2
    for j=1:2
        subplot(2,2,2*(2-j)+i)
        imagesc(final_sorted{j,i})        
        colormap(gca,'jet')
        caxis([-6 6])
        colorbar
        title(horzcat(datalabelset(i),datalabelset2(j)),'FontWeight','bold','FontSize',10)
        ylabel('ACC sensory neuron #')
        xlabel('Time (s)')
        xticks([1:binsize2:length(timescale)])
        xticklabels(timescale2)
        xlim([binsize2*1 binsize2*3+1])
    end    
end

% hit miss difference statisticis (sign-rank test)
stat(1,1) = signrank(z_index(:,1),z_index(:,2));
stat(1,2) = signrank(z_index(:,3),z_index(:,4));
disp(stat(1))
disp(stat(2))

