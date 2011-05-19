function [descs] = ssimDescs(feature, featureFs, tempoPeriod, preferredTempo, beats);

% function [descs] = ssimDescs(feature, featureFs, tempoPeriod, preferredTempo, beats);
%
% compute self-similarity features 

% convert the tempo period into samples consistent with the input feature
tau = (featureFs * tempoPeriod);

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

% find the closest peak to tau
[tmpval,tmpindex] = min(abs(acfPeaks - tau));
tau = acfPeaks(tmpindex);

preferredTau = round(featureFs * 60./preferredTempo); % in samples - do we need to round?? 

tauCandidates = [tau/3 tau/2 tau 2*tau 3*tau]; % look for 1/3, 1/2, 1, 2x and 3x the preferred period

% find the closest tau canditate to the preferred tempo
[tmpval2,tmpindex2] = min(abs(tauCandidates - preferredTau));
closestTau = tauCandidates(tmpindex2);

% then find the nearest peak in the acf to the closets tau candidate
[tmpval3,tmpindex3] = min(abs(acfPeaks - closestTau));
closestTau = acfPeaks(tmpindex3);

weight = (preferredTau - abs(closestTau - preferredTau))/preferredTau;

%simply height of acf at T_pt
descs(1) = acf(closestTau); %% BEAT SALIENCE

%same as previous, but weighting with distance to prefered tempo
descs(2) = descs(1)*weight;


beats(beats <= 0) = []; % remove any beats before the signal
beats(beats > length(feature)/featureFs) = []; % remove any beats after the end of the feature

% find the median IBI
medianIBI = median(diff(beats));
medianIBI = medianIBI*featureFs;

fac = round(medianIBI/closestTau); 
if fac<1
    fac = 1/round(closestTau/medianIBI); % what does this do exactly - is it sensible to invert it??
end
% switch statement on the beats

beats = beats(:); % needs to be a colum vector in case the switch on fac ~=1
switch fac
    case 1
    case 2
        %add 1 beat between each beat
        for i=1:length(beats)-1
            b_supp(i) = beats(i)+(beats(i+1)-beats(i))/2;
        end
        beats = sort([beats; b_supp']);
    case 3
        %add 2 beats between each beat
        for i=1:length(beats)-1
            b_supp_1(i) = beats(i)+(beats(i+1)-beats(i))/3;
            b_supp_2(i) = beats(i)+2*(beats(i+1)-beats(i))/3;
        end
        beats = sort([beats; b_supp_1'; b_supp_2']);
    case 1/2
        % use one beat over 2
        beats = beats(1:2:end); % what about beats(2:2:end)?
    case 1/3
        % use one beat over 3
        beats = beats(1:3:end); % what about beats(2:3:end) and beats(3:3:end)? 
    otherwise
        disp('invalid fac')
        pause
end

% recenter each beat on closest onset
tolerance = 0.05; %in s %%%SHOULD REALLY DEPEND ON IBI
beats = eventPeakAlign(beats,feature,featureFs,tolerance);
nBeats = length(beats);

%resample each beat to have 100 points.. and then take the autocorrelation function of that.
NPointsPerBeat = 100;
for i=1:nBeats-1
    tmp_segment = feature(beats(i):beats(i+1)-1);
    seq(NPointsPerBeat*(i-1)+1:NPointsPerBeat*i) = ...
        resample(tmp_segment,NPointsPerBeat,length(tmp_segment));
    %normalize amplitude
    clear tmp_segment
end
% do you want to smooth the feature before computing the autocorrelation?
% won't there be a peak (potentially) every 100 samples because of possible discontinuities?

%compute self-similarity - not cutting in half this time... (shortcut) - 
acf_2 = max(xcorr(seq, NPointsPerBeat+1, 'coeff'),0); % do as before by HWR the acf

% finding the strength of the acf at a lag of 100 samples
descs(3) = acf_2(end-1);

% then reapply the same weight as before with descs(2);
descs(4) = descs(3)*weight;


% self-similarity of rhythmic patterns at the 
% periodicity function peak rate closest to 100BPM, regardless of 1/4-note rate
[tmpval5,tmpindex5] = min(abs(acfPeaks - preferredTau));
descs(5) = acf(acfPeaks(tmpindex5)); % there was a max statement here, but it's not needed


