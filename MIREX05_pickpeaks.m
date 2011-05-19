function result = MIREX05_pickpeaks(x, bandwidth, threshold)
% function result = pickpeaks(x, bandwidth, threshold, fig)
%  returns indices of local maxima in x within a given (optional)
%  half-bandwidth (default 2 points per side), and above threshold (if given).
%  A negative threshold value is used to indicate a required (positive)
%  difference between the peak value and the mean in the window around the peak.

if nargin < 1
	error('Usage: pickpeaks(x, bandwidth, threshold)');
end
if nargin < 2
	bandwidth = 2;		% default
end
check_threshold = nargin > 2;

% Find peaks
result = [];
peak_count = 0;
mid = 1;
while mid <= length(x)
	lo = max(1,mid - bandwidth);
	hi = min(length(x), mid + bandwidth);
	[max_val, max_pos] = max(x(lo:hi));
	max_pos = max_pos+lo-1;
	if (check_threshold & ...					% Reject low peaks
		(((threshold >= 0) & (max_val <= threshold)) | ...
			((threshold < 0) & (max_val - mean(x(lo:hi)) <= -threshold))))
		mid = hi + 1;							%  - move window on
	elseif max_pos == mid						% Peak found
		peak_count = peak_count + 1;
		result(peak_count) = max_pos;			%  - store it
		mid = hi + 1;							%  - move window on
	elseif max_pos > mid						% Possible peak
		mid = max_pos;							%  - centre window on it
	else										% No peak
		mid = mid + bandwidth;					%  - move window on
	end
end

	
