% spike waveform analysis (similar to Fig. 4g)
% example opto-tagging recording session
% data from mouse 8, session 13 (see Supplementary Table 4)

clear;clc
load Waveform_dataset.mat

sr = 30;            % sampling rate (kHz, Blackrock)
pre1 = 2;           % wide range of spike waveform (pre, ms)
post1 = 3;          % wide range of spike waveform (post, ms)
pre2 = 0.5;         % narrow range of spike waveform (pre, ms)
post2 = 1.5;        % narrow range of spike waveform (post, ms)

% plotting results
figure()
set(gcf,'Position',[350 250 500 400])
sz = 5;
timescale1 = [-pre1:1/sr:post1];
timescale2 = [-pre2:1/sr:post2];
FS=0;RS=0;UC=0;
for k=1:size(spike_width,1)
    hold on
    title('spike waveform','FontSize',12,'FontWeight','bold')
    hold on
    if spike_width(k,1) >= 0.45
        RS = RS + 1;
        plot(timescale2,waveform_mean_final(k,:),'k')    
        hold on
        scatter(spike_width(k,1),rand(1)/6 + 0.80,sz,'k')
    elseif spike_width(k,1) <= 0.35
        FS = FS + 1;
        plot(timescale2,waveform_mean_final(k,:),'r')    
        hold on
        scatter(spike_width(k,1),rand(1)/6 + 0.80,sz,'r')
    else
        UC = UC + 1;
        plot(timescale2,waveform_mean_final(k,:),'color',[.8 .8 .8])    
        hold on
        scatter(spike_width(k,1),rand(1)/6 + 0.80,sz,'MarkerEdgeColor',[.8 .8 .8])
    end        
    ylim([-1 1])
    xlim([-pre2 post2])
    xlabel('time (ms)','FontSize',10,'FontWeight','bold')
    ylabel('norm. voltage','FontSize',10,'FontWeight','bold')
end

hold on
plot ([.45 .45],[-1 1],'k')
hold on
plot ([.35 .35],[-1 1],'k')
%
str1 = horzcat('RS (n=',num2str(RS),', ',num2str(round(RS*100/size(spike_width,1),1)),'%)');
str2 = horzcat('FS (n=',num2str(FS),', ',num2str(round(FS*100/size(spike_width,1),1)),'%)');
str3 = horzcat('UC (n=',num2str(UC),', ',num2str(round(UC*100/size(spike_width,1),1)),'%)');
text(.5,-.4,str1,'FontSize',12,'Color','k')
text(.5,-.6,str2,'FontSize',12,'Color','r')
text(.5,-.8,str3,'FontSize',12,'Color',[.8 .8 .8])

    
