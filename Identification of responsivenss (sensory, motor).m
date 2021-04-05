% purpose of this code
% sorting trials (lick vs. no lick trials; impulsive vs. perceptual lick trials)
% responsiveness statistics (sensory and motor)

clear;clc;
% load lick time information
load Lick_Time_example_session.mat
% load spike time information
load Spike_Time_example_session.mat
% load event(stimuli) time information
load Trial_type_ms_example_session.mat

licktime = LickonTime_final;

% visual detection task
% 2000: go stimuli (flash) onset
% 2001: 75% luminance stimuli (flash) onset
% 2002: 50% luminance stimuli (flash) onset
% 2003: 25% luminance stimuli (flash) onset
% 123456: background trial (no stim) onset

% input arguments
IL_PL_cri = 2000;           % no pre-lick period prior to stimuli onset (ms)
stimuli_offset = 500;       % stimuli offset to rule out lick contamination (ms)
resp_wind_offset = 3000;    % response window offset (ms)
stimuli_cases = 4;          % the number of stimuli
pre_wd = 2000;              % pre-period window for PSTH (ms)
post_wd = 4000;             % post-period window for PSTH (ms)
binsize = 50;               % binsize for data analysis/visualization (ms)
p_val_cri = 0.01;           % p-value for responsiveness statistics

lick_trial_no = [3,5,8];

BG_trial = TrialType{1,5};
allstimuli = vertcat(TrialType{1,1:stimuli_cases});
alltrials = vertcat(TrialType{1,1:end});

Cell_rev = zeros(size(Cell,1),max(max(alltrials),max(LickonTime_final))+500000);
Cell_rev(:,1:size(Cell,2)) = Cell;
clear Cell
Cell = Cell_rev;

% PSTH_total_rev{1,:} all stimuli
% PSTH_total_rev{2,:} perceptual lick w/ lick contam. (stim onset)
% PSTH_total_rev{3,:} perceptual lick w/ lick contam. (lick onset)
% PSTH_total_rev{4,:} perceptual lick w/o lick contam. (stim onset)
% PSTH_total_rev{5,:} perceptual lick w/o lick contam. (lick onset)
% PSTH_total_rev{6,:} miss
% PSTH_total_rev{7,:} -2 ~ .5s no lick (stim onset)
% PSTH_total_rev{8,:} spontaneous lick (lick onset)
% PSTH_total_rev{9,:} no stim trial (catch)

% sorting out perceptual lick trials & miss trials (go stimuli)
for zz = 1:stimuli_cases
    curr=1; curr1=1; curr2=1; curr3=1;
    for i=1:length(TrialType{1,zz})
        % perceptual lick trials (go stimuli)
        % no licking prior to stimuli onset
        if length(find(TrialType{1,zz}(i)-IL_PL_cri < licktime & TrialType{1,zz}(i) > licktime)) == 0
            % licking between stimuli onset and response window offset
            if length(find(TrialType{1,zz}(i) < licktime & TrialType{1,zz}(i)+resp_wind_offset > licktime)) > 0
                % stimuli onset
                TrialType{2,zz}(curr,1) = TrialType{1,zz}(i);
                Lick_latency{1,zz}(curr,1) = licktime(min(find(TrialType{1,zz}(i) < licktime & TrialType{1,zz}(i)+resp_wind_offset > licktime))) - TrialType{1,zz}(i);
                Lick_latency{2,zz}(curr,1) = licktime(min(find(TrialType{1,zz}(i) < licktime & TrialType{1,zz}(i)+resp_wind_offset > licktime)));
                % lick onset
                TrialType{3,zz}(curr,1) = Lick_latency{2,zz}(curr,1);
                curr = curr + 1;
            end
        end
        
        % perceptual lick trials + no lick contamination
        % no licking prior to stimuli onset + no licking during stimuli period
        if length(find(TrialType{1,zz}(i)-IL_PL_cri < licktime & TrialType{1,zz}(i)+stimuli_offset > licktime)) == 0
            
            % trials for statistics (no licking between IL_PL_cri prior to stimuli onset and stimuli_offset)
            TrialType{7,zz}(curr3,1) = TrialType{1,zz}(i);
            curr3 = curr3 + 1;            
            
            % licking during stimuli offest and response window offset
            if length(find(TrialType{1,zz}(i)+stimuli_offset < licktime & TrialType{1,zz}(i)+resp_wind_offset > licktime)) > 0
                % stimuli onset
                TrialType{4,zz}(curr2,1) = TrialType{1,zz}(i);
                % lick onset
                TrialType{5,zz}(curr2,1) = licktime(min(find(TrialType{1,zz}(i) < licktime & TrialType{1,zz}(i)+resp_wind_offset > licktime)));
                curr2 = curr2 + 1;
            end
        end
        
        % miss trials
        % no licking between stimuli onset and response window offset
        if length(find(TrialType{1,zz}(i) < licktime & TrialType{1,zz}(i)+resp_wind_offset > licktime)) == 0
            TrialType{6,zz}(curr1,1) = TrialType{1,zz}(i);
            curr1 = curr1 + 1;
        end
    end
end

for i=1:6001
    allstimuli_mat(:,i) = allstimuli-2001+i;
end

curr = 0;
for i=1:length(lickbout_on)
    if length(find(lickbout_on(i) < allstimuli_mat & lickbout_off(i) > allstimuli_mat)) == 0  
        if lickbout_on(i) ~= lickbout_off(i)
            if lickbout_on(i) > pre_wd
                curr = curr + 1;
                spon_lick(curr,1) = lickbout_on(i);
                spon_lick2(curr,1) = lickbout_off(i);
            end
        end
    end
end

% making PSTH matrix
currentwork = 'making PSTH matrix';
for i=1:size(Cell,1)
    disp(horzcat(currentwork,' - Neuron ',num2str(i),' / ',num2str(size(Cell,1))))
    for j=1:stimuli_cases
        % different trial types of stimuli trials
        for z=1:size(TrialType,1)
            for k=1:size(TrialType{z,j},1)
                PSTH_1ms_total{i,1}{z,j}(k,1:pre_wd+post_wd+1) = Cell(i,TrialType{z,j}(k)-pre_wd:TrialType{z,j}(k)+post_wd);
                if z == 7
                    % sensory pre-activity
                    sensory_statistics{i,j}(k,1) = sum(PSTH_1ms_total{i,1}{z,j}(k,pre_wd-stimuli_offset+1:pre_wd));
                    sensory_statistics{i,j}(k,2) = sum(PSTH_1ms_total{i,1}{z,j}(k,pre_wd+1:pre_wd+stimuli_offset));
                end                    
            end
        end
        % spontaneous lick trials
        for k=1:length(spon_lick)
            PSTH_1ms_total{i,1}{size(TrialType,1)+1,j}(k,1:pre_wd+post_wd+1) = Cell(i,spon_lick(k)-pre_wd:spon_lick(k)+post_wd);
            % motor pre-activity
            motor_statistics{i,j}(k,1) = sum(PSTH_1ms_total{i,1}{size(TrialType,1)+1,j}(k,pre_wd-1500+1:pre_wd-500));
            % motor post-activity
            motor_statistics{i,j}(k,2) = sum(PSTH_1ms_total{i,1}{size(TrialType,1)+1,j}(k,pre_wd+1:pre_wd+1000));
        end  
        % background(no stim) trials
        for k=1:length(BG_trial)
            PSTH_1ms_total{i,1}{size(TrialType,1)+2,j}(k,1:pre_wd+post_wd+1) = Cell(i,BG_trial(k)-pre_wd:BG_trial(k)+post_wd);
        end           
    end
end

% scaling down (1ms to binsize written in input arguments)
edges = [1:binsize:pre_wd+post_wd+1];
binsize2 = 1000/binsize;

currentwork = 'scaling down (binning)';
for i=1:size(Cell,1)
    disp(horzcat(currentwork,' - Neuron ',num2str(i),' / ',num2str(size(Cell,1))))
    for j=1:stimuli_cases
        for z=1:size(TrialType,1)+2        % different trial types of stimuli trials + 2 (spontaneous lick + background(no stim) trials)
            if size(PSTH_1ms_total{i,1}{z,j},1) > 0
                for k=1:size(PSTH_1ms_total{i,1}{z,j},1)
                    for h=1:length(edges)-1
                        PSTH_total_rev{i,1}{z,j}(k,h) = sum(PSTH_1ms_total{i,1}{z,j}(k,edges(h):edges(h+1)-1));
                    end
                end
                PSTH_total_mean{z,j}(i,:) = binsize2 * mean(PSTH_total_rev{i,1}{z,j},1);
            else
                PSTH_total_rev{i,1}{z,j} = [];
                PSTH_total_mean{z,j}(i,:) = zeros(1,120);
            end
        end
    end
    
    for j=1:stimuli_cases
        for z=1:size(TrialType,1)+2        % different trial types of stimuli trials + 2 (spontaneous lick + background(no stim) trials)
            if size(PSTH_1ms_total{i,1}{z,j},1) > 0
                for h=1:length(edges)-1
                    if sum(z == lick_trial_no) > 0
                        % lick onset aligned trials
                        PSTH_total_zscore{z,j}(i,h) = (PSTH_total_mean{z,j}(i,h) - mean(PSTH_total_mean{z,j}(i,1:binsize2*1.5))) / std(PSTH_total_mean{z,j}(i,1:binsize2*1.5));
                    else
                        % stimuli onset aligned trials
                        PSTH_total_zscore{z,j}(i,h) = (PSTH_total_mean{z,j}(i,h) - mean(PSTH_total_mean{z,j}(i,1:binsize2*2))) / std(PSTH_total_mean{z,j}(i,1:binsize2*2));
                    end                    
                end
                PSTH_total_smooth{z,j}(i,:) = smoothdata(PSTH_total_zscore{z,j}(i,:),'gaussian',400/binsize);
            else
                PSTH_total_zscore{z,j}(i,:) = zeros(1,120);
                PSTH_total_smooth{z,j}(i,:) = zeros(1,120);
            end
        end
    end
end

% statistical analysis
currentwork = 'statistical analysis';
for i=1:size(Cell,1)
    disp(horzcat(currentwork,' - Neuron ',num2str(i),' / ',num2str(size(Cell,1))))
    for j=1:stimuli_cases
        % pre versus post spiking activities
        % sign-rank test (non-parametric paired dataset)
        sensory_statistics_p_sr(i,j) = signrank(sensory_statistics{i,j}(:,1),sensory_statistics{i,j}(:,2));
        sensory_delta(i,j) = mean(sensory_statistics{i,j}(:,2) - sensory_statistics{i,j}(:,1));
        motor_statistics_p_sr(i,j) = signrank(motor_statistics{i,j}(:,1),motor_statistics{i,j}(:,2));        
        motor_delta(i,j) = mean(motor_statistics{i,j}(:,2) - motor_statistics{i,j}(:,1));        
        
        % bootstrap analysis (re-sampling w/ replacement)
        clear pre_s_bs post_s_bs delta_s_bs
        pre_s_bs = bootstrp(5000,@mean,sensory_statistics{i,j}(:,1));      % resampling of pre activity
        post_s_bs = bootstrp(5000,@mean,sensory_statistics{i,j}(:,2));     % resampling of post activity
        delta_s_bs = post_s_bs-pre_s_bs;
        if sensory_delta(i,j) >= 0            
            sensory_statistics_p_bs(i,j) = length(find(delta_s_bs<0))/5000;  % in case of increasing statistics
        elseif sensory_delta(i,j) < 0
            sensory_statistics_p_bs(i,j) = length(find(delta_s_bs>0))/5000;  % in case of increasing statistics
        end
        clear pre_m_bs post_m_bs delta_m_bs
        pre_m_bs = bootstrp(5000,@mean,motor_statistics{i,j}(:,1));      % resampling of pre activity
        post_m_bs = bootstrp(5000,@mean,motor_statistics{i,j}(:,2));     % resampling of post activity
        delta_m_bs = post_m_bs-pre_m_bs;
        if motor_delta(i,j) >= 0            
            motor_statistics_p_bs(i,j) = length(find(delta_m_bs<0))/5000;  % in case of increasing statistics
        elseif motor_delta(i,j) < 0
            motor_statistics_p_bs(i,j) = length(find(delta_m_bs>0))/5000;  % in case of increasing statistics
        end
        
        response_type{1,j}(i,1) = sensory_statistics_p_sr(i,j);     % sign-rank test
        response_type{1,j}(i,2) = sensory_statistics_p_bs(i,j);     % bootstrap
        response_type{1,j}(i,3) = sensory_delta(i,j);
        response_type{1,j}(i,4) = motor_statistics_p_sr(i,j);       % sign-rank test
        response_type{1,j}(i,5) = motor_statistics_p_bs(i,j);       % boostrap
        response_type{1,j}(i,6) = motor_delta(i,j);
    end
end


