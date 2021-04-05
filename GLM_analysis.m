%% Code to plot GLM coefficient

clear;clc;
load dataset_GLM.mat

for z=1:size(dataset,1)
    clearvars -except dataset z pval1 beta1 pval2 beta2
    disp(z)
    
    event{1,1} = dataset{z,1};  % stim onset time
    event{2,1} = dataset{z,3};  % lick bout onset time
    
    spk = dataset{z,2};
    cell_1ms = zeros(1,max(event{1,1}(end),event{2,1}(end))+60000);
    cell_1ms(spk) = 1; 
    
    binsize = 50;
    edges = [1:binsize:floor(length(cell_1ms)/binsize)*binsize];
    
    for i=1:length(edges)-1
        cell_binned(1,i) = sum(cell_1ms(edges(i):edges(i+1)-1));
    end
    
    bin_no = length(cell_binned);
    
    for i=1:size(event,1)
        event_binned{i,1} = floor(event{i,1}/binsize);
    end
    
    binsize2 = 1000/binsize;
    event1 = zeros(bin_no,binsize2*(1.5)+1);
    event2 = zeros(bin_no,binsize2*2+1);
    
    % make stim onset time occurence matrix
    for i=1:length(event_binned{1,1})
        for k=1:binsize2*1.5+1
            event1(event_binned{1,1}(i)-5*binsize2/10+k,k) = 1;
        end
    end
    
    % make lick bout onset time occurence matrix
    for i=1:length(event_binned{2,1})
        for k=1:binsize2*2+1
            event2(event_binned{2,1}(i)-10*binsize2/10+k,k) = 1;
        end
    end    
    
    % event information
    X_source = horzcat(event1,event2);
    % spiking activity (binned)
    Y_source = cell_binned';
    
    for zz=1:bin_no
        discard = sum(X_source(i,:));
    end
    
    X_source(find(discard == 0),:) = [];
    Y_source(find(discard == 0),:) = [];       
    
    % main GLM computation
    tic
    mdl = fitglm(X_source,Y_source,'linear','Distribution','poisson');
    toc
    %
    beta1(:,z) = mdl.Coefficients.Estimate(2:binsize2*15/10+2);    
    beta2(:,z) = mdl.Coefficients.Estimate(binsize2*15/10+3:end);
end

% checking outlier
beta1(find(beta1 < -50)) = 0;
beta2(find(beta2 < -50)) = 0;
pval = .05;

for i=1:size(beta1,1)
    pval_re1(1,i) = signrank(beta1(i,:),zeros(1,size(dataset,1)));
end

for i=1:size(beta2,1)
    pval_re2(1,i) = signrank(beta2(i,:),zeros(1,size(dataset,1)));
end

limcri = 1.2;
timescale1 = [-.5:binsize/1000:1];
timescale2 = [-1:binsize/1000:1];
browncolor = [153 76 0]/255;
browncolor2 = [255 153 51]/255;
%
figure()
set(gcf,'Position',[300 50 500 600])
subplot(3,2,1)
title('stimuli onset (individual neurons)')
for z=1:size(dataset,1)
    hold on
    plot(timescale1,beta1(:,z),'color',[.5 1 .5],'LineWidth',.1)
end
hold on
plot(timescale1,mean(beta1,2),'g','LineWidth',1.5)
xlabel('time(s)')
ylabel('GLM beta coef')
xlim([-.5 1])

subplot(3,2,5)
title('stimuli onset (mean/SEM)')
hold on
stdshade(beta1',.2,'g',timescale1);
hold on
plot(timescale1,mean(beta1,2),'g','LineWidth',1.5)
hold on
line([0 0], [-limcri limcri],'color','k')
ylim([-limcri limcri])
xlabel('time(s)')    
ylabel('GLM beta coef')
xlim([-.5 1])

subplot(12,2,15)
imagesc(pval_re1)
hold on
title('p-val')
colormap(gray)
caxis([0 pval])
xticks([1:binsize2/2:size(beta2,1)])
xticklabels([-1:0.5:1])
hold on
xlabel('time(s)')    
axis off

subplot(3,2,2)
title('lick onset (individual neurons)')
hold on
for z=1:size(dataset,1)
    plot(timescale2,beta2(:,z),'color',[1 .5 1],'LineWidth',.1)
end
hold on
plot(timescale2,mean(beta2,2),'m','LineWidth',1.5)
xlabel('time(s)')    
ylabel('GLM beta coef')
xlim([-1 1])

subplot(3,2,6)
title('lick onset (mean/SEM)')
hold on
stdshade(beta2',.2,'m',timescale2);
hold on
plot(timescale2,mean(beta2,2),'m','LineWidth',1.5)
hold on
line([0 0], [-limcri limcri],'color','k')
ylim([-limcri limcri])
xlabel('time(s)')    
ylabel('GLM beta coef')
xlim([-1 1])

subplot(12,2,16)
imagesc(pval_re2)
hold on
title('p-val')
colormap(gray)
caxis([0 pval])
xticks([1:binsize2/2:size(beta2,1)])
xticklabels([-1:0.5:1])
hold on
xlabel('time(s)')    
axis off