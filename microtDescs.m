function [descs] = microtDescs(feature, featureFs,  beats, wavFilename, metLevel);

% function [mtdescs] = microtDescs(feature, featureFs, beats, wavFilename, metLevel);
%
% compute micro-timing features
%
% original version required preferredTempo but don't actually need it.

%% Determine relevant metrical level for microtiming deviations
%parse name --> genre

if nargin<5 % if the metrical level isn't specified as an input then we can try to estimate it from the file name

	genre = wavFilename(1:strfind(wavFilename,'_12')-3);
	switch genre
		case 'greek'
			metLevel = 1/16;
		case 'india'
			metLevel = 1/8;
		case 'jazz'
			metLevel = 1/8;
		case 'samba'
			metLevel = 1/16;
		case 'wafrica'
			metLevel = 1/16;
		otherwise
			disp('unknown genre')
			metLevel = 1/16; % default to this??
	end
end

if metLevel>=1
    disp('metLevel not valid')
    pause
end

tolerance = 0.05; %in s %%%SHOULD REALLY DEPEND ON IBI
% a lot of stuff wasn't needed - including tempo period
% now do analysis based on the beats
beats = eventPeakAlign(beats,feature,featureFs,tolerance);
nBeats = length(beats);

% resample each beat to have 40 points..  - why 40 this time and not 100?
% put each beat into a column of a matrix and then take the mean
NPointsPerBeat = 40;
% pre-allocate matrix
pattern = zeros(NPointsPerBeat, nBeats-1);

for i=1:nBeats-1
    tmp_segment = feature(beats(i):beats(i+1)-1);
	pattern(:,i) = resample(tmp_segment,NPointsPerBeat,length(tmp_segment));
    %normalize amplitude -- this doesn't happen any more, should it?
    clear tmp_segment
end

% this is the mean beat pattern
aveBeatPattern = mean(pattern,2);
% find maximum value
maxBeatPattern = max(aveBeatPattern);

% identify consistent deviations from quantized positions

% determine how to split up the beat - based on metrical levels
factor = (1/4)/metLevel;
% note this will give an error if metLevel = 1/4 or 1/2 ...
% so probably need a catch here to do something else in that case.

% determine some window around which we can find peaks in the aveBeatPattern
halfWindowSize = (NPointsPerBeat/factor)/2-1;
for i=1:factor-1
	% identify the quantised position
    quantPos(i) = i*NPointsPerBeat/factor+1;
	% select a region around the quantised position
    tmp = aveBeatPattern(quantPos(i)-halfWindowSize:quantPos(i)+halfWindowSize);

	% find index of maximum value in the tmp range 
	[peakHeight,peakIndex] = max(tmp);
	% adjust for location of window, i.e. make it centre on 0
	peakIndex = peakIndex - halfWindowSize - 1;

	% now calculate deviations
	dev2(i) = abs(peakIndex);
	dev1(i) = dev2(i) * peakHeight/maxBeatPattern;	

	clear tmp;	
end

dev1 = dev1/halfWindowSize; %normalized between [0 1]
dev2 = dev2/halfWindowSize; %normalized between [0 1]

descs(1) = max(dev1); %% SYSTEMATIC MICROTIMING
descs(2) = mean(dev1);
descs(3) = max(dev2);
descs(4) = mean(dev2);

