function [descs] = periodfuncDescs(feature, featureFs, tempoPeriod_I, beats_I, tempoPeriod_II, beats_II);

% function [descs] = periodfuncDescs(feature, featureFs, tempoPeriod, beats_I, tempoPeriod_II, beats_II);
%
% compute some extra periodicity based descriptors - including event density and unsystematic microtimingDescs

% all of this up to acfPeaks is taken from the ssimDescs function... so probably want a sub-function to do this processing

maxLag = 5; % presumably 5 seconds
corrSize = round(maxLag * featureFs);
% check that feature is at least maxLag seconds long, if not, zero pad
feature = [feature(:)' zeros(1,corrSize - (length(feature))+1)];

% find the autocorrelation function and set half wave rectify.
fullAcf = max(xcorr(feature, corrSize, 'coeff'),0);  
% extract the latter half, so zero lag is sample '1'
acf = fullAcf(corrSize+1:end);

% now pick the peaks of the acf using MIREX peak picker (I have a different one)
peakHalfWindow = 0.05;

% need to check code for this - maybe replace with another function
acfPeaks = MIREX05_pickpeaks(acf, floor(featureFs * peakHalfWindow), -0.001);
acfPeaks = acfPeaks(2:end); %1st is for lag 0

% now for some sorting based on whether tempo_I is faster than tempo_II

% convert the tempo period into samples consistent with the input feature
tau_I = (featureFs * tempoPeriod_I);

% find the closest peak to tau
[tmpval,tmpindex] = min(abs(acfPeaks - tau_I));
tau_I = acfPeaks(tmpindex);

% if second tempo is unspecified, then set a flag
if (tempoPeriod_II==-1)
    onlyOneLevel = 1;
else       
    onlyOneLevel = 0; % BUG - must set this to zero.
	% put tempo period II into feature sampling rate
    tau_II = featureFs*tempoPeriod_II; 
	% adjust tau_II to the nearest acfPeaks
	[tmpval2,tmpindex2] = min(abs(acfPeaks - tau_II));
	tau_II = acfPeaks(tmpindex2);

end

% note... not adjusting beats to peaks of feature... - how come?

% put beats_I into feature samples
beats_I = round(beats_I*featureFs);
% and remove any bad beats
beats_I(beats_I <=0) = [];
beats_I(beats_I > length(feature)) = [];

% if only one level assign beats_fast = beats_I
% if we have more than one level of beat annotations
% similarly fix beats_II
% and then assign beats_fast and beats_slow
if onlyOneLevel

	beats_fast = beats_I;

	else
		% put beats_II into feature samples
		beats_II = round(beats_II*featureFs);
		% and remove any bad beats
		beats_II(beats_II <=0) = [];
		beats_II(beats_II > length(feature)) = []; % fabien doesn't have this line, but i think it's needed

	% now decide which sequence is beats_fast and which is beats_slow
	if tau_I > tau_II % i.e. tau_I has a greater period, hence slower then beats_I are beats_slow
		beats_fast = beats_II;		
		beats_slow = beats_I;
	else % otherwise 
		beats_fast = beats_I;		
		beats_slow = beats_II;
	end

end


% mean of energy variance on beat segments gives a measure of event density (or a noisy feature?)
% measure the energy from the start of one beat to the next 
% should we go up to beats_fast(i+1)-1 ? 
% or also should the event density measurement be centred on the beats?? 
% calculate for both slow and fast beats
for i=1:length(beats_fast)-1
    v_fast(i) = var(feature(beats_fast(i):min(length(feature),beats_fast(i+1))));
end
descs(1) = mean(v_fast); %% EVENT DENSITY

keyboard

if onlyOneLevel
    descs(2) = 0;
else
    for i=1:length(beats_slow)-1
        v_slow(i) = var(feature(beats_slow(i):min(length(feature),beats_slow(i+1))));
    end
    descs(2) = mean(v_slow);
end



%% deviations around beats at the faster level, normalized by #beats
halfWindowSize = 4; % in this context it's fixed - different from systematic microtiming 
dev = [];
for i=2:length(beats_fast)-1
    tmp = feature(beats_fast(i)-halfWindowSize:beats_fast(i)+halfWindowSize);
    indx = [1:2*halfWindowSize+1];
    dev = [dev (sum(indx.*tmp)/sum(tmp))-(halfWindowSize+1)]; %recentered temp centroid
    clear tmp indx
end
dev = dev/featureFs;
descs(3) = sum(dev)/length(beats_fast); 
descs(4) = sum(abs(dev))/length(beats_fast); %% UNSYSTEMATIC MICROTIMING DEVIATIONS
clear dev 



%% max correlation coefficient of the autocorrelation of deviations around beats at *half* the faster level, normalized by #beats
% interpolate beats by 2.
beats_veryFast = doubleBeats(beats_fast);
%beats_veryFast = round(beats_veryFast(1:end-1)); % these should all be in the sensible range, once we remove the last beat
beats_veryFast = round(beats_veryFast);
beats_veryFast(beats_veryFast > length(feature)) = []; % remove any outliers

% same halfWindowSize as before
dev = [];
% calculate the temporal centroid of the feature around the very fast beats
for i=2:length(beats_veryFast)-1
    halfWindowSize = round((beats_veryFast(i+1)-beats_veryFast(i))/2);
    tmp = feature(beats_veryFast(i)-halfWindowSize:beats_veryFast(i)+halfWindowSize);
    indx = [1:2*halfWindowSize+1];
    dev = [dev (sum(indx.*tmp)/sum(tmp))-(halfWindowSize+1)]; %recentered temp centroid
    clear tmp indx 
end

% dev is the sequence of temporal centroids for each very fast beat...
% we are taking the autocorrelation function over a very short lag range
% apparently the strongest periodicity lag in the acf (NOTE no acfs are unbiased - should they be?)
% is a measure of the period of repeated deviations - not sure what this is supposed to correlate to.

% what does this line do?? --
limCorr = ceil(5*featureFs/length(beats_veryFast));
tmp = xcorr(dev,limCorr,'coeff');
acf = max(tmp(limCorr+1:end),0); % will this take only the maximum positive 
clear tmp
%value of the highest coeff of dev autocorrelation
descs(5) = max(acf(2:end));
clear dev 

% Syncopation at the fastest level - go to an even faster level??



beats_fastest = round(doubleBeats(beats_veryFast)); % or should that be round(doubleBeats(doubleBeats(beats_fast))) ??


onBeat = zeros(1,length(feature));
% should we define for the offBeat as well?? 
% the line for onBeat is repeated in the old function
for i=4:4:length(beats_fastest)-4
    onBeat(beats_fastest(i):beats_fastest(i+2)-1) = feature(beats_fastest(i):beats_fastest(i+2)-1);
    offBeat(beats_fastest(i+2):beats_fastest(i+4)-1) = feature(beats_fastest(i+2):beats_fastest(i+4)-1);
end
onBeatNRG = sum(onBeat.^onBeat)/length(beats_fast);
offBeatNRG = sum(offBeat.^offBeat)/length(beats_fast);

% syncopation is the ratio of off-beat energy to on-beat energy
descs(6) = offBeatNRG/onBeatNRG;


function beatsOut = doubleBeats(beatsIn);

beatsOut(1:2:2*length(beatsIn)) = beatsIn; 
for i=1:length(beatsIn)-1,
	beatsOut(2*i) = beatsIn(i) + 0.5*(beatsIn(i+1)-beatsIn(i));
end

