% hit vs. miss comparison z-score (Fig. 2b)

clear; clc;

figure()
set(gcf,'Position',[150 150 900 500])

for zzzz = 1:2
    clearvars -except zzzz
    if zzzz == 1
        load ACC_Sinc_neurons.mat
    elseif zzzz == 2
        load ACC_Sdec_neurons.mat
    end
    
    % input arguments
    binsize = 25;
    binsize2 = 1000/binsize;
    edges = [1:binsize:6001];
    timescale = [-2:1/binsize2:4];
    timescale2 = [-2:1:4];
    
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
    
    % hit miss difference statistics (time-series)
    for i=1:2
        for z=1:size(final_sorted{1,1},2)
            stat_final(i,z) = signrank(final_sorted{i,1}(:,z),final_sorted{i,2}(:,z));
        end
    end
    
    for i=2
        subplot(3,2,zzzz)
        imagesc(stat_final(i,:))
        caxis([0 0.05])
        colormap(gca,'gray')
        xticks([1:binsize2:length(timescale)])
        xticklabels(timescale2)
        xlim([binsize2*1 binsize2*3+1])
        axis off
        hold on
        title('p-val')
        
        subplot(3,2,[zzzz+2,zzzz+4])
        hold on
        if zzzz == 1
            title('ACC Sinc neurons (n = 56)')
        elseif zzzz == 2
            title('ACC Sdec neurons (n = 40)')
        end
        clear aa bb
        aa = mean(final_sorted{i,1},1);
        bb = std(final_sorted{i,1},1)/sqrt(size(heatmapdata,1));
        errorshade(timescale(1:end-1),aa,aa+bb,aa-bb,[.9 .9 .9])
        hold on
        plot(timescale(1:end-1),mean(final_sorted{i,1},1),'k')
        hold on
        clear aa bb
        aa = mean(final_sorted{i,2},1);
        bb = std(final_sorted{i,2},1)/sqrt(size(heatmapdata,1));
        if zzzz == 1
            errorshade(timescale(1:end-1),aa,aa+bb,aa-bb,[1 0.8 0.8])
            hold on
            plot(timescale(1:end-1),mean(final_sorted{i,2},1),'r')
            ylim([-.5 6])
        elseif zzzz == 2
            errorshade(timescale(1:end-1),aa,aa+bb,aa-bb,[0.8 0.8 1])
            hold on
            plot(timescale(1:end-1),mean(final_sorted{i,2},1),'b')
            ylim([-3 2])
        end
        xlim([-1 1])
        xlabel('Time (s)')
        ylabel('z-score')
    end
end