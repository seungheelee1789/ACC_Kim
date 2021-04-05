% Visually evoked signals and lick time (Fig. 2j)
% example neuron (#9)

clear;clc;
load ACC_Sinc_neurons_licktime.mat

% TemporalCoding{i,6}(:,1): trial # (sorted according to lick latency)
% TemporalCoding{i,6}(:,2): lick latency (ms)
% TemporalCoding{i,6}(:,3): stimuli onset (ms)

% input arguments
exampleneuron = 9;      
binsize = 50;       
segment = 8;
sw = 200;               % sliding window

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
        reaction_time2(i,k) = std(TemporalCoding2_rev{i,k}(:,1))/sqrt(length(TemporalCoding2_rev{i,k}));
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
        TemporalCoding5{i,2}(j,1:(2000/binsize)) = smoothdata(TemporalCoding5{i,2}(j,1:(2000/binsize)),'gaussian',400/binsize);
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

% plot average data
t1 = [-2+binsize/1000:binsize/1000:4];

for i=1:segment+1
    colorcri(segment+2-i,:) = colorbrewer.seq.YlOrRd{1,segment+1}(i,:)/255;
end

for k=1:9
    for z=1:size(TemporalCoding4{9,k},1)
        asdf{k,1}(z,1) = 1000/binsize * mean(TemporalCoding4{9,k}(z,(1000/binsize)*2+1:(1000/binsize)*2.5));
    end
    
    meanStat(k,1) = mean(asdf{k,1});
    stdStat(k,1) = std(asdf{k,1})/sqrt(length(asdf{k,1}));
end

figure()
sgtitle('Example ACC Sinc neuron')
set(gcf,'Position',[150 150 350 600])
subplot(2,1,1)
title('stimuli onset')
hold on
for i=1:segment+1
    if i == 1
        plot(t1,final_tempcoding{2,10-i}(exampleneuron,:),'Color','k','LineWidth',2,'LineStyle','--')
    elseif i ~= 1
        plot(t1,final_tempcoding{2,10-i}(exampleneuron,:),'Color',[colorcri(10-i,:)],'LineWidth',2)
    end
end
hold on
line([0 0],[0 20],'Color','g','LineWidth',2.5)
xlim([-.5 .5])
ylim([0 20])
xlabel('time (s)')
ylabel('spikes/s')

aa = reaction_time(exampleneuron,:)/1000;
bb = reaction_time2(exampleneuron,:)/1000;
cc = meanStat;
dd = stdStat;
sz = 50;
[r1 p1] = corrcoef(aa',cc(1:segment));
disp(r1(1,2))
disp(p1(1,2))
[p,S] = polyfit(aa',cc(1:segment),1);
[y_fit,delta] = polyval(p,aa,S);

subplot(2,1,2)
title(strcat('Pearson correlation: r = ',num2str(r1(1,2)),', p = ',num2str(round(p1(1,2),4))))
hold on
for i=1:segment+1
    if i ~= segment+1
        line([aa(i)-bb(i) aa(i)+bb(i)],[cc(i) cc(i)],'Color',colorcri(i,:),'LineWidth',2)
        line([aa(i) aa(i)],[cc(i)-dd(i) cc(i)+dd(i)],'Color',colorcri(i,:),'LineWidth',2)
        hold on
        scatter(aa(i),cc(i),sz,colorcri(i,:),'filled')
    elseif i == segment+1
        line([1.6 1.6],[cc(i)-dd(i) cc(i)+dd(i)],'Color','k','LineWidth',2)
        hold on
        scatter(1.6,cc(i),sz,'k','filled')
    end
end
hold on
line([0 0],[1 3],'Color','m','LineWidth',2.5)
hold on
plot(aa,y_fit,'k','LineStyle','--')
xlim([.4 1.6])
ylim([4 16])
xlabel('lick latency (s)')
ylabel('spikes/s')
