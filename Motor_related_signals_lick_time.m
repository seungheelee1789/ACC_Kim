% Motor related signals and lick time (Fig. 2l,m)

clear;clc;
load ACC_Mdec_neurons_licktime.mat

% TemporalCoding{i,6}(:,1): trial # (sorted according to lick latency)
% TemporalCoding{i,6}(:,2): lick latency (ms)
% TemporalCoding{i,6}(:,3): stimuli onset (ms)

% input arguments
binsize = 50;       
segment = 8;
sw = 400;               % sliding window
load colorbrewer.mat

edges = [1:binsize:6001];

target = [1:1:size(spktime,1)];

% trial segmentation according to lick times
for i=1:length(target)
    segmentcri = floor(size(TemporalCoding{i,6},1) / segment);
    for k=1:segment
        % collecting trials (stimuli onset) from each segment
        TemporalCoding2{i,k}(:,1) = TemporalCoding{i,6}(segmentcri*(k-1)+1:segmentcri*k,3);
        % collecting trials (lick latency) from each segment
        TemporalCoding2_rev{i,k}(:,1) = TemporalCoding{i,6}(segmentcri*(k-1)+1:segmentcri*k,2);
        % collecting trials (lick onset) from each segment
        TemporalCoding21{i,k}(:,1) = TemporalCoding{i,6}(segmentcri*(k-1)+1:segmentcri*k,3)+TemporalCoding{i,6}(segmentcri*(k-1)+1:segmentcri*k,2);
    end
    
    % miss trials (stimuli onset)
    TemporalCoding2{i,segment+1}(:,1) = TemporalCoding{i,5}(:,2);
    for za = 1:size(TemporalCoding{i,5},1)
        TemporalCoding21{i,segment+1}(za,1) = sum(TemporalCoding{i,5}(za,1:2));
    end
    clear segmentcri
end

% mean lick time of each segment from individual neurons
for i=1:length(target)
    for k=1:segment        
        reaction_time(i,k) = mean(TemporalCoding2_rev{i,k}(:,1));
    end
end
reaction_time_final(1,:) = mean(reaction_time,1);
reaction_time_final(2,:) = std(reaction_time,1)/sqrt(length(target));

% make PSTH/PETH of each segment from individual neurons
for i=1:length(target)
    tempcri = max(max(spktime{i,1}),max(licktime{i,1})) + 800000;
    clear tempcell
    tempcell = zeros(1,tempcri);
    for z=1:length(spktime{i,1})
        tempcell(spktime{i,1}(z)) = 1;
    end
    
    for j=1:segment+1
        % perceptual lick trials
        if j ~= segment + 1
            for k=1:size(TemporalCoding2{i,j},1)
                % PSTH (stim onset)
                TemporalCoding3{i,j}(k,:) = tempcell(TemporalCoding2{i,j}(k)-2000:TemporalCoding2{i,j}(k)+4000);        %% 1ms bin-size
                % PETH (lick onset)
                TemporalCoding31{i,j}(k,:) = tempcell(TemporalCoding21{i,j}(k)-3000:TemporalCoding21{i,j}(k)+3000);     %% 1ms bin-size
            end
        
        % miss or spontaneous lick trials
        elseif j == segment + 1
            % PSTH (stim onset) - miss trials
            for k=1:size(TemporalCoding2{i,j},1)
                TemporalCoding3{i,j}(k,:) = tempcell(TemporalCoding2{i,j}(k)-2000:TemporalCoding2{i,j}(k)+4000);     %% 1ms bin-size
            end
            % PETH (lick onset) - spontaneous lick trials
            for k=1:size(sponlick{i,1},2)
                TemporalCoding31{i,segment + 1}(k,:) = tempcell(sponlick{i,1}(k)-3000:sponlick{i,1}(k)+3000);     %% 1ms bin-size
            end
        end
    end
end

% data binning and averaging
for i=1:length(target)
    % data binning
    for j=1:segment+1
        for k=1:size(TemporalCoding3{i,j},1)
            for z=1:length(edges)-1
                TemporalCoding4{i,j}(k,z) = sum(TemporalCoding3{i,j}(k,edges(z):edges(z+1)-1));               %% 50ms bin-size
            end
        end
        for k=1:size(TemporalCoding31{i,j},1)
            for z=1:length(edges)-1
                TemporalCoding41{i,j}(k,z) = sum(TemporalCoding31{i,j}(k,edges(z):edges(z+1)-1));               %% 50ms bin-size
            end
        end
    end
    % data averaging
    for j=1:segment+1
        TemporalCoding5{i,1}(j,:) = (1000/binsize)*mean(TemporalCoding4{i,j},1);           
        TemporalCoding51{i,1}(j,:) = (1000/binsize)*mean(TemporalCoding41{i,j},1);         
    end
end

% data smoothing
for i=1:length(target)
    for j=1:segment+1
        TemporalCoding5{i,2}(j,:) = smoothdata(TemporalCoding5{i,1}(j,:),'gaussian',sw/binsize);
        TemporalCoding51{i,2}(j,:) = smoothdata(TemporalCoding51{i,1}(j,:),'gaussian',sw/binsize);
    end
end

% data reorganization
for i=1:length(target)
    for z=1:size(TemporalCoding5{i,1},1)
        % average PSTH,PETH
        final_tempcoding{1,z}(i,:) = TemporalCoding5{i,1}(z,:);
        final_tempcoding2{1,z}(i,:) = TemporalCoding51{i,1}(z,:);
        % smoothed average PSTH,PETH
        final_tempcoding{2,z}(i,:) = TemporalCoding5{i,2}(z,:);
        final_tempcoding2{2,z}(i,:) = TemporalCoding51{i,2}(z,:);
    end
end

for i=1:length(target)
    for z=1:size(TemporalCoding5{i,1},1)-1
        reactiontime{i,1}(z,1) = mean(TemporalCoding21{i,z} - TemporalCoding2{i,z});
    end
end

for k=1:2
    for z=1:segment+1
        for i=1:length(target)
            for j=1:size(final_tempcoding{k,z},2)
                if final_tempcoding{k,z}(i,j)*0 ~= 0
                    final_tempcoding{k,z}(i,j) = 0;
                end
                if final_tempcoding2{k,z}(i,j)*0 ~= 0
                    final_tempcoding2{k,z}(i,j) = 0;
                end
            end
        end
    end
end

for i=1:segment+1
    meanfinal_tempcoding111(i,:) = mean(final_tempcoding{2,i},1);
    meanfinal_tempcoding112(i,:) = mean(final_tempcoding2{2,i},1);
    
    for k=1:length(target)
        stdtakehome1(k,i) = mean(final_tempcoding2{2,i}(k,(1000/binsize)*2+1:(1000/binsize)*2.5));
    end
    stdtakehome2(i,1) = std(stdtakehome1(:,i))/sqrt(length(target));
end

% plot average data
t1 = [-2:binsize/1000:4-binsize/1000];
t2 = [-3:binsize/1000:3-binsize/1000];
for i=1:segment+1
    colorcri(segment+2-i,:) = colorbrewer.seq.YlOrRd{1,segment+1}(i,:)/255;
end
%
reaction_time_final(3,:) = ((reaction_time_final(1,:)/binsize))+2000/binsize;
scatteraxis = ones(1,segment)*3;

for i=1:segment+1
    meanfinal_revised1(i,:) = smoothdata(meanfinal_tempcoding111(i,:),'gaussian',sw/binsize);
    meanfinal_revised2(i,:) = smoothdata(meanfinal_tempcoding112(i,:),'gaussian',sw/binsize);
end

figure()
sgtitle('ACC MDec neurons (n = 44)')
set(gcf,'Position',[150 150 700 600])
subplot(2,2,1)
title('stimuli onset')
hold on
for i=1:segment+1
    if i ~= segment+1
        plot(t1,meanfinal_revised1(i,:),'Color',[colorcri(i,:)],'LineWidth',2)
        hold on
        scatter(reaction_time_final(1,i)/1000,scatteraxis(i),15,colorcri(i,:),'filled')
    elseif i == segment+1
        plot(t1,meanfinal_revised1(i,:),'Color','k','LineWidth',2,'LineStyle','--')
    end
end
hold on
line([0 0],[1 3],'Color','g','LineWidth',2.5)
xlim([-1 3])
ylim([1 3])
xlabel('time (s)')
ylabel('spikes/s')

subplot(2,2,2)
title('lick onset')
hold on
for i=1:segment+1
    if i ~= segment+1
        plot(t2,meanfinal_revised2(i,:),'Color',[colorcri(i,:)],'LineWidth',2)                
    elseif i == segment+1
        plot(t2,meanfinal_revised2(i,:),'Color','k','LineWidth',2)
    end
    xlim([-2 2])
end
hold on
line([0 0],[1 3],'Color','m','LineWidth',2.5)
xlim([-2 2])
ylim([1 3])
xlabel('time (s)')
ylabel('spikes/s')

reaction_time_final_real = reaction_time_final';

for i=1:length(target)
    for k=1:size(TemporalCoding51{i,2},2)
        [r,p] = corrcoef(reactiontime{i,1}(:,1),TemporalCoding5{i,2}(1:end-1,k));
        corrdata_r_stim(i,k) = r(1,2);
        corrdata_p_stim(i,k) = p(1,2);
        clear r p
        [r,p] = corrcoef(reactiontime{i,1}(:,1),TemporalCoding51{i,2}(1:end-1,k));
        corrdata_r_lick(i,k) = r(1,2);
        corrdata_p_lick(i,k) = p(1,2);
        clear r p
    end
end
%
for i=1:length(target)
    aa(i,1) = 1 * sum((corrdata_p_stim(i,(1000/binsize)*2:(1000/binsize)*3)));
    bb(i,1) = 1 * sum((corrdata_p_lick(i,(1000/binsize)*2:(1000/binsize)*3)));
end
[cc, index1] = sortrows(aa);
align1 = index1;
[dd, index2] = sortrows(bb);
align2 = index2;

for i=1:length(target)
    corrdata_p_stim_aligned(i,:) = corrdata_p_stim(align1(i),:);
    corrdata_p_lick_aligned(i,:) = corrdata_p_lick(align2(i),:);
    corrdata_r_stim_aligned(i,:) = corrdata_r_stim(align1(i),:);
    corrdata_r_lick_aligned(i,:) = corrdata_r_lick(align2(i),:);
end

corrdata_r_stim_aligned2 = corrdata_r_stim_aligned;
for i=1:size(corrdata_p_stim_aligned,1)
    for j=1:size(corrdata_p_stim_aligned,2)
        if corrdata_p_stim_aligned(i,j) >= 0.05
            corrdata_r_stim_aligned2(i,j) = 0;
        end
    end
end

subplot(2,2,3)
imagesc(corrdata_r_stim_aligned2)
hold on
title('Pearson correlation')
caxis([-1 1]) 
xlim([(1000/binsize)*1 (1000/binsize)*5])
colormap(redblue)
ylabel('ACC MDec neurons')
xlabel('time (s)')
xticks([0:20:120])
xticklabels([-2:1:4])
colorbar


