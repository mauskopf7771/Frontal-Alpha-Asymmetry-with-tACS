%% Frontal Alpha Asymmetry Depression and tACs %%
% Plugins required;
% Signal Processing Toolbox
% EEGLAB
% EEGLAB - Neuroelectrics toolbox
% EEGLAB - MARA

clear all;
close all;
clc;
eeglab
% Get subject list nd turn that into a loop
subject_list = dir('*.easy')
subject_list.name
n = length(subject_list)
%% Import all files, rename, add locations and save

for i = 1:n;
EEG = pop_easy([subject_list(i).name(1:end-5) '.easy'],1,0,[]);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 0,'setname',subject_list(i).name,'gui','off'); 
eeglab redraw;

%% Pre-Processing %%
 % Add Locs files %
 EEG=pop_chanedit(EEG, ...
    'lookup','/Users/[REDACTED]/Documents/Matlab_Toolboxes/eeglab2021.1/plugins/dipfit/standard_BESA/standard-10-5-cap385.elp',...
    'changefield',{7,'theta','-54'},'changefield',{7,'radius','0.51111'},'changefield',{8,'theta','54'},...
    'changefield',{8,'radius','0.51111'}); 
 EEG = eeg_checkset( EEG );
 eeglab redraw;

 % Rereference to Cz
 EEG = pop_reref( EEG, 2);
 EEG = eeg_checkset( EEG );
 eeglab redraw;
 % Filtering high pass 2Hz, low pass of 50Hz 
 EEG = pop_eegfiltnew(EEG, 'locutoff',2,'plotfreqz',1);
 EEG = pop_eegfiltnew(EEG, 'hicutoff',50,'hicutoff',60,'revfilt',1);
 EEG = eeg_checkset( EEG );
 eeglab redraw;

 %% Trial Rejection
 %EEG = pop_clean_rawdata(EEG, 'FlatlineCriterion',5,'ChannelCriterion',0.8,...
 %    'LineNoiseCriterion',4,'Highpass','off','BurstCriterion',20,'WindowCriterion',0.25,...
 %    'BurstRejection','off','Distance','Euclidian','WindowCriterionTolerances',[-Inf 7] );
 % [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
 % EEG = eeg_checkset( EEG );
 %eeglab redraw;
 
 %% Epoch
 EEG = pop_editeventvals(EEG,'insert',{1,[],[],[]},'changefield',{1,'type',1},'changefield',{1,'latency',0},...
     'changefield',{1,'latency_ms',0});
 EEG = pop_epoch( EEG, {  '1'  }, [0  300], 'epochinfo', 'yes');
 EEG = eeg_checkset( EEG );
 eeglab redraw;
 
 %% Trial Rejection
 % Apparently not needed and an adequate alternative is using component
 % rejection
 %% ICA
 % EEG = pop_runica(EEG, 'icatype', 'runica', 'extended',1,'interrupt','on');
 %% Automated form of component rejection
 % pop_processMARA ( ALLEEG,EEG,CURRENTSET )
 %[ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
 %EEG = eeg_checkset( EEG );
 %eeglab redraw;
 
 %% Save dataset
 EEG = pop_saveset( EEG, 'filename', [subject_list(i).name(1:end-5) '_Pre_Processed.set'],'filepath',pwd); 
end
 
 %% Component Rejection %%
 % for i = 1:n
 %    EEG = pop_loadset('filename', [subject_list{s}  '_ICA.set'],'filepath',[pwd]);
 %    [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
 %    EEG = eeg_checkset( EEG );
 %    eeglab redraw;
 % end
 
  %% Extract Frontal Alpha Asymmetry Value using Tesar's FAA Toolbox %%
  % Moditfied to remove the GUI and save values as columns to add to an excel
  % dataset.
  % ==================================
  % Computes frontal alpha asymmetry from specified data - channels,
  % frequency and latency.
  % This book/study is a result of the research funded by the project 
  % Nr. LO1611 with a financial support from the MEYS under the NPU I program.
  % Michael Tesar <tesarm@jcu.cz>
  % Ceske Budejovice, 2016
  %  ____  ____  ____  ____  ____  ____  ____  ____  ____  ____  ____  ____ 
  % ||n ||||e ||||u ||||r ||||o ||||p ||||a ||||c ||||a ||||b ||||r ||||a ||
  % ||__||||__||||__||||__||||__||||__||||__||||__||||__||||__||||__||||__||
  % |/__\||/__\||/__\||/__\||/__\||/__\||/__\||/__\||/__\||/__\||/__\||/__\|

  for i = 1:n;
      EEG = pop_loadset('filename', [subject_list(i).name(1:end-5) '_Pre_Processed.set'],'filepath',[pwd]);
      [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 0,'setname',subject_list(i).name,'gui','off');
      eeglab redraw;

      tmpData = mean(EEG.data, 3);
      L = tmpData(3,:);
      R = tmpData(4,:);
      % Pick frequency range
      FREQ_1 = 8;
      FREQ_2 = 13;
      % Compute power spectrum for left channel
      WIND = hamming(floor(length(L))/2);   % Changed to 2 was originally 1.5 sec time windows
      OVER = floor((length(L))/1.5/2);        % 50% overlap: previously /1.5/2
      SIGN = L';                              % Get signal
      [s, freqs, t, power] = spectrogram(SIGN, WIND, OVER, [], EEG.srate);
      indFreqs = find(freqs>FREQ_1 & freqs<FREQ_2);
      POW_L = power(indFreqs);
      % Compute power spectrum for right channel
      WIND = hamming(floor(length(R))/2);   % Get 1.5 sec time windows - changed to 2
      OVER = floor((length(R))/1.5/2);        % 50% overlap: previously /1.5/2
      SIGN = R';                              % Get signal
      [s, freqs, t, power] = spectrogram(SIGN, WIND, OVER, [], EEG.srate);
      indFreqs = find(freqs>FREQ_1 & freqs<FREQ_2);
      POW_R = power(indFreqs);
      % Compute whole FAA
      FAA(i) = mean(abs(log(POW_R)-log(POW_L)))
      disp(FAA(i))
  end

%% Save values into a dataframe
data = FAA'
%% load the excel file %% 

%% Concatenate with the behavioural data %%

%% End Message %%
for done = 1:20
    display('######################### DONE! #########################')
end
