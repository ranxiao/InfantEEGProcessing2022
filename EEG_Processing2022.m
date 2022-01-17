%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function: Preprocess individual-session EEG and perform spectral analysis
% Ran Xiao, Duke University
% Original code:  1/19/2016
% revision: 1/17/2022, tweak to handle example data collected in 2022 and
% integrate automatic IC rejection in EEGLab 2021.1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialize directories
addpath(genpath('./eeglab2021.1/')); 
DataDir = './Data/';
if ~exist(DataDir,'dir')
    mkdir(DataDir);
end
ResultDir = './Results/';
if ~exist(ResultDir,'dir')
    mkdir(ResultDir);
end
% handle new data file
FileName = 'Activity9 - 2048Hz_export_EEG_2cameras_OPAL.txt';
data = readtable([DataDir,FileName]);
% select 32 channel EEG data, 2048 srate, last second of data is nan, so
% removed
EEGdata = table2array(data(1:122800,2:33))'; 

% load data into EEGlab
eeglab;
EEG = pop_importdata('dataformat','array','nbchan',32,'data','EEGdata','srate',2048,'pnts',0,'xmin',0);
EEG = eeg_checkset( EEG );
% downsample to 250hz
EEG = pop_resample( EEG, 250);
EEG = eeg_checkset( EEG );
% load channel file
EEG = pop_chanedit(EEG, 'lookup','./BioSemi_32Ch.ced','load',{'./BioSemi_32Ch.ced' 'filetype' 'autodetect'});
EEG = eeg_checkset( EEG );
% re-reference to T7-T8
EEG = pop_reref( EEG, [7 24] ,'keepref','on');
EEG = eeg_checkset( EEG );
% bandpass to 0.2-30 Hz
EEG = pop_eegfiltnew(EEG, 'locutoff',0.2,'hicutoff',30,'plotfreqz',1);
EEG = eeg_checkset( EEG );
pop_eegplot( EEG, 1, 1, 1);
eeglab redraw;

%% Mark and reject segments with large fluctuations
% (GUI operation needed: 1. Mark segments; 2. Run code block; 3. Click reject button)
% Tip: Set window length as 20 sec to select bad segments
g = get(gcf, 'userdata');
BadSegments = int64(g.winrej(:,[1 2]));
save(strcat(ResultDir,FileName,'_BadSegments.mat'), 'BadSegments');

%% Bad channel rejection and interpolation
% 1. Channel rejection by calculating Kurtosis index and reject if exceeding 5
EEG = pop_rejchan(EEG, 'elec',[1:32] ,'threshold',5,'norm','on','measure','kurt');
EEG = eeg_checkset( EEG );
% 2. Interpolate bad channels by surrounding channels
EEG = pop_interp(EEG, ALLEEG(1).chanlocs, 'spherical');
EEG = eeg_checkset( EEG );
pop_eegplot( EEG, 1, 1, 1);

%% (OPTIONAL)Visual inspection as a secondary approach to identify addional ones
% Caution: Use Channel Number, instead of channel labels!
prompt = {'Enter indices for bad channels'};
dlg_title = 'Visual selection';
Inputs = inputdlg(prompt,dlg_title,[1 32]);
BadCh_Visual= str2num(Inputs{:});
EEG = pop_interp(EEG, BadCh_Visual, 'spherical');
EEG = eeg_checkset( EEG );
pop_eegplot( EEG, 1, 1, 1);

%% Spatial filtering through common average reference
EEG = pop_reref( EEG, []);
EEG = eeg_checkset( EEG );
pop_eegplot( EEG, 1, 1, 1);
eeglab redraw;

%% Independent component analysis
% 1. Run ICA with PCA to tackle lost rank
dataRank = rank(EEG.data);
EEG = eeg_checkset( EEG );
EEG = pop_runica(EEG, 'icatype', 'runica', 'extended',1,'interrupt','on','pca',dataRank);
EEG = eeg_checkset( EEG );
% pop_topoplot(EEG, 0, [1:dataRank] ,' resampled',[5 6] ,0,'electrodes','on');
% pop_selectcomps(EEG, [1:dataRank] );
% run automatic IC labeling
EEG = pop_iclabel(EEG, 'default');
EEG = eeg_checkset( EEG );

% based on IC labeling results, automatically identify Bad ICs
[p,I]=max(EEG.etc.ic_classification.ICLabel.classifications,[],2);
% bad IC criteria: non-brain pattern with probability over 70%
BadIC = find(p>=0.7 & I~=1);

% 2. Save ICA results and EEG data before rejection for later reference
EEG = eeg_checkset( EEG );
EEG = pop_saveset( EEG, 'filename',[FileName '_BeforeICARej.set'],'filepath',ResultDir);
EEG = eeg_checkset( EEG );
pop_selectcomps(EEG, [1:dataRank] );
saveas(gcf,[ResultDir FileName '_ICA' '.jpg']);

% remove bad ICs and save data
EEG = pop_subcomp( EEG, BadIC, 0);
EEG = eeg_checkset( EEG );
pop_eegplot( EEG, 1, 1, 1);
eeglab redraw;
% 4. Save EEG data after bad IC rejection
EEG = eeg_checkset( EEG );
EEG = pop_saveset( EEG, 'filename',[FileName '_AfterICARej.set'],'filepath',ResultDir);
EEG = eeg_checkset( EEG );

%% Perform Spectral analysis
% Loading preprocessed dataset
EEG = pop_loadset('filename',[FileName '_AfterICARej.set'],'filepath',ResultDir);
EEG = eeg_checkset( EEG );
Data = EEG.data;
% Calculate spectral powers using Pwelch
Sampling = 250;
Win = 2;
Nfft = Win*Sampling;
Overlap = 0.5*Nfft;
NoCh = size(Data,1);
Power = zeros(NoCh, Sampling/2*Win+1);
RelPower = zeros(NoCh, Sampling/2*Win+1);
for i = 1:1:NoCh
    [pxx,f] = pwelch(Data(i,:),hann(Nfft),Overlap,Nfft,Sampling);
    Power(i,:) = pxx;
    RelPower(i,:)=pxx/sum(pxx(1:find(f==30)));
end
save(strcat(ResultDir,FileName,'_Pw'), 'Power');
save(strcat(ResultDir,FileName,'_RelPw'),'RelPower');

% Visuallization of spectral profiles
load ChanLabel;
figure;
h=plot(f,Power');
xlim([2 10]);
% ylim([0 60]);
xlabel('Frequency');
ylabel('Power');
title([FileName ' Power']);
set(h, {'color'}, num2cell(jet(NoCh), 2));
legend(ChanLabel{:});
set(gcf, 'Position', get(0,'Screensize')); % Maximize figure.
saveas(gcf,[ResultDir FileName '_Power' '.fig']);
saveas(gcf,[ResultDir FileName '_Power' '.jpg']);

figure;
h=plot(f,RelPower');
xlim([2 10]);
ylim([0 0.1]);
xlabel('Frequency');
ylabel('Relative Power');
title([FileName ' Relative Power']);
set(h, {'color'}, num2cell(jet(NoCh), 2));
legend(ChanLabel{:});
set(gcf, 'Position', get(0,'Screensize')); % Maximize figure.
saveas(gcf,[ResultDir FileName '_RelPower' '.fig']);
saveas(gcf,[ResultDir FileName '_RelPower' '.jpg']);
