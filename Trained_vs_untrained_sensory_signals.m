% trained vs. untrained ACC sensory signals (Extended Data Fig. 3b,c)

clear; clc;

figure()
set(gcf,'Position',[150 150 700 500])
    
for zzzz=1:2
    clearvars -except zzzz
    if zzzz == 2
        load ACC_trained.mat
        final_re = trained_data;
    elseif zzzz == 1
        load ACC_untrained.mat
        final_re = untrained_data;
    end
    
    % input arguments
    binsize = 50;
    binsize2 = 1000/binsize;
    edges = [1:binsize:6001];
    timescale = [-2:1/binsize2:4];
    timescale2 = [-2:1:4];
    datalabelset = {'Untrained','Trained'};
    sz = 3;
    
    % making PSTH (1ms) for spontaneous lick and stimuli onset
    for i=1:size(final_re{1,1},1)
        for k=1:2
            clear baseline1 baseline2
            if k == 1
                baseline1 = mean(final_re{1,k}(i,1:binsize2*1.5));
                baseline2 = std(final_re{1,k}(i,1:binsize2*1.5));
            elseif k == 2
                baseline1 = mean(final_re{1,k}(i,1:binsize2*2));
                baseline2 = std(final_re{1,k}(i,1:binsize2*2));
            end
            if baseline2 == 0
                baseline2 = 1;
            end
            for z=1:length(edges)-1
                final_z{1,k}(i,z) = (final_re{1,k}(i,z) - baseline1) / baseline2;
            end
            final_smooth{1,k}(i,:) = smoothdata(final_z{1,k}(i,:),'gaussian',300/binsize);
            if k == 1
                z_index(i,k) = mean(final_smooth{1,k}(i,binsize2*2+1:binsize2*3));
                z_index2(i,k) = abs(mean(final_smooth{1,k}(i,binsize2*2+1:binsize2*3)));
            elseif k == 2
                z_index(i,k) = mean(final_smooth{1,k}(i,binsize2*2+1:binsize2*2.5));
                z_index2(i,k) = abs(mean(final_smooth{1,k}(i,binsize2*2+1:binsize2*2.5)));
            end
        end
    end
    
    for i=1:size(final_re{1,1},1)
        sortingindex(i,1) = -1*mean(final_smooth{1,1}(i,binsize2*2+1:binsize2*3));
        sortingindex2(i,1) = -1*mean(final_smooth{1,2}(i,binsize2*2+1:binsize2*2.5));
    end
    [aa index] = sortrows(sortingindex);
    [aa2 index2] = sortrows(sortingindex2);
    
    for i=1:size(final_re{1,1},1)
        final_smooth_sorted{1,1}(i,:) = final_smooth{1,1}(index(i),:);
        final_smooth_sorted{1,2}(i,:) = final_smooth{1,2}(index2(i),:);
    end
    
    clear cumdist
    edges2 = [0:0.02:5];
    edges3 = [0:0.02:4.98];
    clear currdata aa1 bb1
    currdata = abs(z_index(:,1));
    currdata(find(currdata >= 5)) = 5;
    aa1 = histcounts(currdata,edges2);
    bb1 = cumsum(aa1)/length(currdata);
    cumdist{1,1} = bb1;
    clear currdata aa1 bb1
    currdata = abs(z_index(:,2));
    currdata(find(currdata >= 5)) = 5;
    aa1 = histcounts(currdata,edges2);
    bb1 = cumsum(aa1)/length(currdata);
    cumdist{1,2} = bb1;    
    
    for i=2
        subplot(2,2,zzzz)
        imagesc(final_smooth_sorted{1,i})
        colormap(gca,'jet')
        if i == 1
            caxis([-5 5])
        elseif i == 2
            caxis([-5 5])
        end
        colorbar
        title(datalabelset(zzzz),'FontWeight','bold','FontSize',10)
        ylabel('ACC neurons')
        xlabel('Time (s)')
        xticks([1:binsize2:length(timescale)])
        xticklabels(timescale2)
        xlim([1 binsize2*4+1])
        subplot(2,2,3)
        if zzzz == 1
            hold on
            plot(edges3,cumdist{1,i},'LineWidth',1,'color',[.8 .8 .8])
        elseif zzzz == 2
            plot(edges3,cumdist{1,i},'k','LineWidth',1)
        end     
        xlim([0 5])
        xlabel('Visual signals (|z|)')
        ylabel('CDF')
    end    
end
subplot(2,2,3)
legend(datalabelset,'Location','southeast')
legend('boxoff')