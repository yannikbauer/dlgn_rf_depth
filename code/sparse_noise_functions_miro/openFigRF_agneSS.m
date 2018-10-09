% TO CHECK IF WE USE SERIES IN SIZE TUNING
u4 = load('~/Dropbox/clean/locomotion_layer/animals/agneReAnalysis/currentData/uInfo_PV-Cre_SzTun.mat');

%%%% CHECK 
%%%% PV_0054 %%%% (4series)
%%%% CHECK s07/08!!!
% no RFs for ... but log files say nice responses
pv54ser = fetch(data.Series('series_num in (7,8,10,11)') & data.Mice('mouse_id like "%PVCre_0054%"'));
% s7 - check why exps are not grouped, Tf tuning if exists is hi 
openfig('/Volumes/busse_lab/users/sinem/SparseNoise_agneSS/sparseNoiseFigs_MUAE/z_co8/sparseNoiseMapon_PVCre_0054_07_(2).fig');
openfig('/Volumes/busse_lab/users/sinem/SparseNoise_agneSS/sparseNoiseFigs_MUAE/z_co8/sparseNoiseMapon_PVCre_0054_07_(13,15).fig');
openfig('/Volumes/busse_lab/users/sinem/SparseNoise_agneSS/sparseNoiseFigs_MUAE/z_co8/sparseNoiseMapon_PVCre_0054_07_(5).fig');
openfig('/Volumes/busse_lab/users/sinem/SparseNoise_agneSS/sparseNoiseFigs_MUAE/z_co8/sparseNoiseMapon_PVCre_0054_07_(7).fig');
openfig('/Volumes/busse_lab/users/sinem/SparseNoise_agneSS/sparseNoiseFigs_MUAE/z_co8/sparseNoiseMapon_PVCre_0054_07_(10).fig');
plot(data.TempFreqTuning(pv54ser(1), 'exp_num = 4'));
plot(data.TempFreqTuning(pv54ser(1), 'exp_num = 1'));
% 1 unit in the size tuning (leave out)
plot(agne.SizeTuningOptoSegment(stripKey(u4.uInfo, 'data.Units'), pv54ser(1)));

% s8 -  check why exps are not grouped, Tf tuning is weird but not clear %%%%%(PULV?)%%%%% 
% leave out
openfig('/Volumes/busse_lab/users/sinem/SparseNoise_agneSS/sparseNoiseFigs_MUAE/z_co10/sparseNoiseMap_PVCre_0054_08_(1).fig');
openfig('/Volumes/busse_lab/users/sinem/SparseNoise_agneSS/sparseNoiseFigs_MUAE/z_co10/sparseNoiseMap_PVCre_0054_08_(3).fig');
openfig('/Volumes/busse_lab/users/sinem/SparseNoise_agneSS/sparseNoiseFigs_MUAE/z_co8/sparseNoiseMapon_PVCre_0054_08_(8).fig');
openfig('/Volumes/busse_lab/users/sinem/SparseNoise_agneSS/sparseNoiseFigs_MUAE/z_co8/sparseNoiseMapon_PVCre_0054_08_(12).fig');
plot(data.TempFreqTuning(pv54ser(2), 'exp_num = 4'));
% plot(agne.SizeTuningOptoSegment(stripKey(u4.uInfo, 'data.Units'), pv54ser(2)));
%%% 
% s10 (ok)
openfig('/Volumes/busse_lab/users/sinem/SparseNoise_agneSS/sparseNoiseFigs_MUAE/z_co10/sparseNoiseMapon_PVCre_0054_10_(1,2,12).fig')
% s11 - (ok)
openfig('/Volumes/busse_lab/users/sinem/SparseNoise_agneSS/sparseNoiseFigs_MUAE/z_co8/sparseNoiseMap_PVCre_0054_11_(1).fig');
openfig('/Volumes/busse_lab/users/sinem/SparseNoise_agneSS/sparseNoiseFigs_MUAE/z_co8/sparseNoiseMap_PVCre_0054_11_(11).fig');

%%%% PV_0041 %%%%(1series)
%%%% CHECK
pv41ser = fetch(data.Series('series_num =15') & data.Mice('mouse_id like "%PVCre_0041%"'));
% s15 -  RFS a bit big? TF both hi and low?? %%%%(PULV?)(LGN?)%%%%% 
openfig('/Volumes/busse_lab/users/sinem/SparseNoise_agneSS/sparseNoiseFigs_MUAE/z_co8/sparseNoiseMap_PVCre_0041_15_(4).fig')
openfig('/Volumes/busse_lab/users/sinem/SparseNoise_agneSS/sparseNoiseFigs_MUAE/z_co8/sparseNoiseMap_PVCre_0041_15_(18).fig')
openfig('/Volumes/busse_lab/users/sinem/SparseNoise_agneSS/sparseNoiseFigs_MUAE/z_co8/sparseNoiseMap_PVCre_0041_15_(1,10,16).fig')
plot(data.TempFreqTuning(pv41ser, 'exp_num = 2'));
plot(agne.SizeTuningOptoSegment(stripKey(u4.uInfo, 'data.Units'), pv41ser))


%%%% Ntsr1-Cre_0071 %%%%(5series)
%%% CHECK
ntsr71ser= fetch(data.Series('series_num in (3,4,5,6,7)') & data.Mice('mouse_id like "%Ntsr1-Cre_0071%"'));
% S05 has DII STAINNING!!! (only one except for v1, according to db and logfile)
% s05 - Rfs kinda big? TF tun kinda low? 
openfig('/Volumes/busse_lab/users/sinem/SparseNoise_agneSS/sparseNoiseFigs_MUAE/z_co8/sparseNoiseMap_Ntsr1-Cre_0071_05_(1,7,13).fig')
plot(data.TempFreqTuning(ntsr71ser(3)) & data.Stimuli('exp_file_name like "%tftun%"'));
% s03 -  check why exps are not grouped, see clear progression in exp1 but no tfTuning
openfig('/Volumes/busse_lab/users/sinem/SparseNoise_agneSS/sparseNoiseFigs_MUAE/z_co8/sparseNoiseMap_Ntsr1-Cre_0071_03_(1).fig')
openfig('/Volumes/busse_lab/users/sinem/SparseNoise_agneSS/sparseNoiseFigs_MUAE/z_co8/sparseNoiseMap_Ntsr1-Cre_0071_03_(3).fig')
plot(agne.SizeTuningOptoSegment(stripKey(u4.uInfo, 'data.Units'), ntsr71ser(2))); % no size tuning exps included
%%% 
% s07(ok) but tftun is weird
openfig('/Volumes/busse_lab/users/sinem/SparseNoise_agneSS/sparseNoiseFigs_MUAE/z_co8/sparseNoiseMap_Ntsr1-Cre_0071_07_(1,11).fig')
plot(data.TempFreqTuning(ntsr71ser(5)) & data.Stimuli('exp_file_name like "%tftun%"'))
% s06 (ok)- seems like Rs are low/offscreen and tf tun is hi
openfig('/Volumes/busse_lab/users/sinem/SparseNoise_agneSS/sparseNoiseFigs_MUAE/z_co8/sparseNoiseMap_Ntsr1-Cre_0071_06_(4).fig')
openfig('/Volumes/busse_lab/users/sinem/SparseNoise_agneSS/sparseNoiseFigs_MUAE/z_co8/sparseNoiseMap_Ntsr1-Cre_0071_06_(6).fig')
plot(data.TempFreqTuning(ntsr71ser(4)) & data.Stimuli('exp_file_name like "%tftun%"'))
% s04 (ok) - quite nice
openfig('/Volumes/busse_lab/users/sinem/SparseNoise_agneSS/sparseNoiseFigs_MUAE/z_co8/sparseNoiseMap_Ntsr1-Cre_0071_04_(1).fig')
openfig('/Volumes/busse_lab/users/sinem/SparseNoise_agneSS/sparseNoiseFigs_MUAE/z_co8/sparseNoiseMap_Ntsr1-Cre_0071_04_(8).fig')
plot(data.TempFreqTuning(ntsr71ser(2)) & data.Stimuli('exp_file_name like "%tftun%"'))


%%%% Ntsr1-Cre_0039 %%%% (4series)
%%% CHECK
ntsr39ser= fetch(data.Series('series_num in (4,5,7,8)') & data.Mice('mouse_id like "%Ntsr1-Cre_0039%"'));
% s08  %%%%(PULV?)(LGN?)%%%%% 
% no real RFs, mixed tf tuning 
% size tuning: clear that RFs are not being hit, 
% PULV possible bc series notes says anterior than before but!
% monitor pos(40deg) makes me think that prolly not so anterior as to be in pulvinar bc 
% RFs should move more central and s07 (nice RFs) also have mon pos ~ 40-45
openfig('/Volumes/busse_lab/users/sinem/SparseNoise_agneSS/sparseNoiseFigs_MUAE/z_co8/sparseNoiseMap_Ntsr1-Cre_0039_08_(1,9).fig')
plot(data.TempFreqTuning(ntsr39ser(4)) & data.Stimuli('exp_file_name like "%tftun%"'))
%plot(agne.SizeTuningOptoSegment(ntsr39ser(4), 'exp_num = 8', 'rsq_c > .8'))
plot(agne.SizeTuningOptoSegment(stripKey(u4.uInfo, 'data.Units'), ntsr39ser(4)))

% S05 is STAINED!!!
% s05 (ok?) - RFs hard to make out but TF hi-ish and szTun clear not hitting RF
openfig('/Volumes/busse_lab/users/sinem/SparseNoise_agneSS/sparseNoiseFigs_MUAE/z_co8/sparseNoiseMapon_Ntsr1-Cre_0039_05_(1,3).fig')
openfig('/Volumes/busse_lab/users/sinem/SparseNoise_agneSS/sparseNoiseFigs_MUAE/z_co8/sparseNoiseMapon_Ntsr1-Cre_0039_05_(8).fig')
plot(data.TempFreqTuning(ntsr39ser(2)) & data.Stimuli('exp_file_name like "%tftun%"'))
%plot(agne.SizeTuningOptoSegment(ntsr39ser(2), 'rsq_c > .8'))
plot(agne.SizeTuningOptoSegment(stripKey(u4.uInfo, 'data.Units'), ntsr39ser(2)))
%
% s04 (ok) - can see nice small RFs, TF hi-ish
openfig('/Volumes/busse_lab/users/sinem/SparseNoise_agneSS/sparseNoiseFigs_MUAE/z_co20/sparseNoiseMap_Ntsr1-Cre_0039_04_(1).fig')
plot(data.TempFreqTuning(ntsr39ser(1)) & data.Stimuli('exp_file_name like "%tftun%"'))
% s07 (ok) - very nice progression, small RFs, tfTuning mixed but towards hi
openfig('/Volumes/busse_lab/users/sinem/SparseNoise_agneSS/sparseNoiseFigs_MUAE/z_co20/sparseNoiseMap_Ntsr1-Cre_0039_07_(1).fig')
openfig('/Volumes/busse_lab/users/sinem/SparseNoise_agneSS/sparseNoiseFigs_MUAE/z_co8/sparseNoiseMap_Ntsr1-Cre_0039_07_(3,12).fig')
plot(data.TempFreqTuning(ntsr39ser(3)) & data.Stimuli('exp_file_name like "%tftun%"'))




%%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%%
%%%% %%%% %%%% %%%% BELOW: SEEMS OK %%%% %%%% %%%% %%%% %%%% %%%% 
%%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%%

%%%% PV_0058 %%%% (2series)
% s13 (ok)
openfig('/Volumes/busse_lab/users/sinem/SparseNoise_agneSS/sparseNoiseFigs_MUAE/z_co10/sparseNoiseMap_PVCre_0058_13_(1,3).fig')
openfig('/Volumes/busse_lab/users/sinem/SparseNoise_agneSS/sparseNoiseFigs_MUAE/z_co10/sparseNoiseMap_PVCre_0058_13_(5).fig')
% s15 (ok)
openfig('/Volumes/busse_lab/users/sinem/SparseNoise_agneSS/sparseNoiseFigs_MUAE/z_co10/sparseNoiseMap_PVCre_0058_15_(1).fig')


%%%% PV_0046 %%%% (3series)
% s08 (ok)
openfig('/Volumes/busse_lab/users/sinem/SparseNoise_agneSS/sparseNoiseFigs_MUAE/z_co8/sparseNoiseMap_PVCre_0046_08_(1).fig')
openfig('/Volumes/busse_lab/users/sinem/SparseNoise_agneSS/sparseNoiseFigs_MUAE/z_co8/sparseNoiseMap_PVCre_0046_08_(3).fig')
openfig('/Volumes/busse_lab/users/sinem/SparseNoise_agneSS/sparseNoiseFigs_MUAE/z_co8/sparseNoiseMap_PVCre_0046_08_(10).fig')
% s06 (ok)
openfig('/Volumes/busse_lab/users/sinem/SparseNoise_agneSS/sparseNoiseFigs_MUAE/z_co8/sparseNoiseMap_PVCre_0046_06_(1).fig')
openfig('/Volumes/busse_lab/users/sinem/SparseNoise_agneSS/sparseNoiseFigs_MUAE/z_co8/sparseNoiseMap_PVCre_0046_06_(4).fig')
openfig('/Volumes/busse_lab/users/sinem/SparseNoise_agneSS/sparseNoiseFigs_MUAE/z_co8/sparseNoiseMap_PVCre_0046_06_(13).fig')
% s05 (ok)
openfig('/Volumes/busse_lab/users/sinem/SparseNoise_agneSS/sparseNoiseFigs_MUAE/z_co8/sparseNoiseMap_PVCre_0046_05_(2,11).fig')


%%
return
% total number of series (in this file)
nhere = 4+1+5+4+2+3; % should be 19 

%% check series recordings
% get all the series that are used (same as in makeSparseNoiseMats_AgneAnalysis)
u1 = load('/Users/sinem/Dropbox/clean/locomotion_layer/animals/agneReAnalysis/currentData/uInfo_%Ntsr1-Cre%_sponOpto.mat');
u2 = load('/Users/sinem/Dropbox/clean/locomotion_layer/animals/agneReAnalysis/currentData/uInfo_%Ntsr1-Cre%_ChR2_stillSponOptoLGN.mat');
u3 = load('/Users/sinem/Dropbox/clean/locomotion_layer/animals/agneReAnalysis/currentData/uInfo_%PV-Cre%_ChR2_stillSponOptoLGN.mat');
u4 = load('/Users/sinem/Dropbox/clean/locomotion_layer/animals/agneReAnalysis/currentData/uInfo_PV-Cre_SzTun.mat');
tmp = cellfun(@(u) stripKey(u.uInfo, 'data.Series'), {u1,u2,u3,u4},'uniformoutput',0);
% 19 series
allSeries = fetch(data.Series(cat(2,tmp{:}))); 
myMice = fetchn(data.Mice(allSeries), 'mouse_id')';
% 
stainedSeries = fetch(data.Series(allSeries, 'series_description like "%did%" or series_description like "%dii%" '));
% but there is also one blank ntsr_39 s04 (which also has a staining)
data.Series(allSeries, 'series_description like ""')
% first two series are "problematic" above
% Ntsr1_Cre 71 s05 (DiD) - only staining for animal (no v1 either)
% Ntsr1_Cre 39 s05 (DiI); 
%              s04 (DiD);  % v1 is also stained (DiD) 
