function descriptors = makeGrooveDescriptors1(wavFilename,audioDir,dataDir);

% function descriptors = makeGrooveDescriptors1(wavFilename,audioDir,dataDir);
%
% top level function to compute descriptors for groove experiments.
% 
% this are a recoded version of those used by Fabien in the 2011 paper with Guy Madison


% make the full path to the wave file - currently this code is mac/linux specific
longAudioFilename = [audioDir,'/',wavFilename];

% now compute input feature 
[feature,featureFs] = inputfeature(longAudioFilename);

% read beats_I, beats_II, tempo and ticks
% strip .wav from the wavefile
beatFilename_I = [dataDir,'/beats_I/',wavFilename(1:end-4),'.beats'];
beats_I = textread(beatFilename_I);
beats_I(beats_I <= 0) = []; % remove any zero valued or negative beats  

beatFilename_II = [dataDir,'/beats_II/',wavFilename(1:end-4),'.beats'];
beats_II = textread(beatFilename_II);
beats_II(beats_II <= 0) = []; % remove any zero valued or negative beats  

% can do this... but actually we can work out the tempo directly from the beats
% tempoFilename = [dataDir,'/tempo/',wavFilename(1:end-4),'.tempo'];
% tempo_I = textread(tempoFilename);
%tempo_I = 60./median(diff(beats_I));
%% and likewise for beats_II
%tempo_II = 60./median(diff(beats_II));

% work out tempi based on median IBI for each beat sequence
tempoPeriod_I = median(diff(beats_I));
tempoPeriod_II = median(diff(beats_II));


% I don't think that these are used.
%tickFilename = [dataDir,'/tick/',wavFilename(1:end-4),'.tick'];
%tick = textread(tickFilename);

% now compute self-similarity descriptors
% parameters
preferredTempo = 100;

ssdescs = ssimDescs(feature, featureFs, tempoPeriod_I, preferredTempo, beats_I);
% should write these to an output file .ssimd

if (~exist([dataDir,'/SelfSimilarityDescs'])),
	mkdir([dataDir,'/SelfSimilarityDescs']);
end;

fid = fopen([dataDir,'/SelfSimilarityDescs/',wavFilename(1:end-4),'.ssimd'],'w');
fprintf(fid,'%f %f %f %f %f',ssdescs); % expect 5 numbers
fclose(fid);

% now compute micro timing descriptors
[mtdescs] = microtDescs(feature, featureFs, beats_I, wavFilename);

if (~exist([dataDir,'/MicroTimingDescs'])),
	mkdir([dataDir,'/MicroTimingDescs']);
end;

fid = fopen([dataDir,'/MicroTimingDescs/',wavFilename(1:end-4),'.microtd'],'w');
fprintf(fid,'%f %f %f %f',mtdescs); % expect 4 numbers
fclose(fid);


% now compute periodicity descriptors

% if only have one level, or only want to analyse one level then
onlyOneLevel = 0;
if onlyOneLevel
	tempoPeriod_II = -1;
	beats_II = -1;
end

[pddescs] = periodfuncDescs(feature, featureFs, tempoPeriod_I, beats_I, tempoPeriod_II, beats_II);

if (~exist([dataDir,'/PeriodFuncDescs'])),
	mkdir([dataDir,'/PeriodFuncDescs']);
end;

fid = fopen([dataDir,'/PeriodFuncDescs/',wavFilename(1:end-4),'.pfds2'],'w');
% the file extension is .pfds2 because there's another function, but apparently I don't need to implement that.
fprintf(fid,'%f %f %f %f %f %f',pddescs); % expect 6 numbers
fclose(fid);

% concatenate all the descriptors together
descriptors = [ssdescs(:)' mtdescs(:)' pddescs(:)'];
