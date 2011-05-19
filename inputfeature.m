function [feature,featureFs] = inputfeature(fname);

% function feature = inputfeature(fname)
% this version is authentic the code from fabien.

% take a wave file and then make the input feature
% as described in Madison et al 2011 paper

% calculate the energy envelope on 11.6ms (non-overlapping?) frames
% across 8 sub-bands which are roughly logarithmically spaced

% [0     - 100]
% [100   - 216]
% [216   - 467]
% [467   - 1009]
% [1009  - 2183]
% [2183  - 4719]
% [4719  - 10200]
% [10200 - 22050] 

% using 6th order butterworth filters


% first step, read audio file
% or do we want to make stereo features?? - is there a difference in groove perception from stereo to mono?
[x fs] = wavread(fname);
% turn the stereo file into a mono one.
x = mean(x,2);

% if signal is not at 44khz, then resample it 
if (fs ~=44100),
	x = resample(x,44100,fs);
end

lenx = length(x);

% now make filters these are the cut-off frequencies
numBands = 8;
lowf = 100;
maxf = fs/2;

% logarithmically spaced bands
fco = [0 lowf * (maxf / lowf) .^ ([0:numBands-1]/(numBands-1)) ] / (maxf);


% store them in cell structure, 

filtOrder = 6; % from the paper - should be an even number... 
% make 1 low pass, then 6 band pass, then 1 high pass
for i=1,
	[b{i},a{i}]=butter(filtOrder, fco(i+1),'low');
end
% apparently we need to make the filter order half for the band pass filters
% i guess this is due to stability... not guaranteed stable with order other than 6 -> untested!
for i=2:numBands-1, % make the band pass filters in a loop
	[b{i},a{i}] = butter(filtOrder/2,[fco(i) fco(i+1)]);
end
% now make the high pass filter
for i=numBands,
	[b{i},a{i}]=butter(filtOrder, fco(i),'high');
end

% now make the filtered versions of the audio
% this is probably going to eat a lot of memory (if long files)
% but want to listen first to hear that it's working
% definite issues with stability of the filters
% could use auditory filters instead?

% downsampling factor
dsFactor = 32;
% store sampling rate of decimated signal
dsFs = round(fs/dsFactor);
% taken from Fabien's code.
% arbitrary filter order, and cut-off
[bb,aa] = butter(4, 2*20/dsFs ,'Low');


% parameters for overlap add of feature
winlen = 32;
step = 16;

% pre-allocate output feature
numFrames = ceil(lenx/ dsFactor/step );
feature = zeros(numBands,numFrames);

times = ceil(length(x) / dsFactor);


for i=1:numBands

	% why is this the arithmetic mean? i guess the centre frequency is still linear even if bands are log-spaced?? 
	% work out group delay 
    gd1{i} = grpdelay(b{i}, a{i}, [fco(i) mean(fco(i:i+1)) fco(i+1)]);
    subBand = filter(b{i}, a{i}, x(1+round(gd1{i}(2)):end)); % you start after the delayed samples rather than shifting back the others bands to compensate... 

	% does it make sense to HWR and square the sub-band audio?
	% now HWR and square...
	subBand = subBand .* (subBand > 0) .^2; % as implemented by Fabien but this doesn't square
%	subBand = ( subBand .* (subBand > 0) ) .^2; % this version squares the result

	% downsample by dsFactor -- could use 'resample' function as well - would need to specify fs as well
	dsSubBand = decimate(subBand,dsFactor);
	% dsSubBand is not guaranteed to be non-negative - so HWR again - but not in Fabien's code, so maybe comment it out for now...
 	% dsSubBand = dsSubBand .* (dsSubBand > 0);

	% now smooth it using secondary filter
	% working out group delay again... do the same sub-band cut off frequencies apply??
	% since the filter coefficients are different, maybe the delays will be smaller? 
	% something doesn't sit right here based on the fact that the audio has been decimated
	gd2{i} = grpdelay(bb,aa,[fco(i) mean(fco(i:i+1)) fco(i+1)]); 
	% now make smoothed version of the subBand
	smoothDsSubBand = filter(bb,aa,dsSubBand(1+round(gd2{i}(2)):end));
	% zeropad the smoothDsSubBand in case it's not the right length
	smoothDsSubBand = [smoothDsSubBand; zeros(times - size(smoothDsSubBand,1), 1)];

	% and overlap add.
	pin = 0;
	pend = length(smoothDsSubBand)-winlen;
	ct = 0;
	while (pin<pend)
		ct = ct+1;
		% probably want an 'abs' in there - definitely want one.. but not for the moment.
		feature(i,ct) = sum((smoothDsSubBand(pin+1:pin+winlen)));
		pin = pin+step;
	end

end

% now do Weber's law compression on the feature
mu = 100;
feature = max(0,feature);
logDiffFeature = (log(1+mu*feature))/(log(1+mu));
feature = max(0,diff(logDiffFeature,1,2));

feature = sum(feature,1);
%feature = feature/sum(feature); % normalise to sum to unity

featureFs = fs / dsFactor / step; 
