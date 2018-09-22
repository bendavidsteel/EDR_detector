function SCR = findpeaks(data, data_deriv, times, offset, start_WT, end_WT, thres_low, theta)


% 		This function finds the peaks of an EDA signal and returns basic properties.
% 		Also, peak_end is assumed to be no later than the start of the next peak. (Is this okay??)
% 
% 		********* INPUTS **********
% 		data:        DataFrame with EDA as one of the columns and indexed by a datetimeIndex
% 		offset:      the number of rising samples and falling samples after a peak needed to be counted as a peak
% 		start_WT:    maximum number of seconds before the apex of a peak that is the "start" of the peak
% 		end_WT:      maximum number of seconds after the apex of a peak that is the "rec.t/2" of the peak, 50% of amp
% 		thres:       the minimum uS change required to register as a peak, defaults as 0 (i.e. all peaks count)
% 		sampleRate:  number of samples per second, default=8
% 
% 		********* OUTPUTS **********
% 		peaks:               list of binary, 1 if apex of SCR
% 		peak_start:          list of binary, 1 if start of SCR
% 		peak_start_times:    list of strings, if this index is the apex of an SCR, it contains datetime of start of peak
% 		peak_end:            list of binary, 1 if rec.t/2 of SCR
% 		peak_end_times:      list of strings, if this index is the apex of an SCR, it contains datetime of rec.t/2
% 		amplitude:           list of floats,  value of EDA at apex - value of EDA at start
% 		max_deriv:           list of floats, max derivative within 1 second of apex of SCR

    sampleRate = 15;

	peaks = zeros(length(data_deriv),1);
	peak_sign = sign(data_deriv);
    peak_times = zeros(length(data),1);
    peak_eda = zeros(length(data),1);
    
	for i = round(offset)+1: length(data_deriv) - offset - 1
		if((peak_sign(i) == 1) && (peak_sign(i + 1) < 1))
			peaks(i) = 1;
            peak_times(i) = times(i);
            peak_eda(i) = data(i);
			for j = 1:round(offset)
                if ((i - j > 0) && (i + j > length(peak_sign)))
                    if (peak_sign(i - j) < 1) || (peak_sign(i + j) > -1)
                        %if peak_sign[i-j]==-1 or peak_sign[i+j]==1:
                        peaks(i) = 0;
                        break
                    end
                end
            end
        end
    end

	% Finding start of peaks
	peak_start = zeros(length(data_deriv),1);
	peak_start_times = zeros(length(data),1);
	max_deriv = zeros(length(data),1);
	rise_time = zeros(length(data),1);

	for i = 1: length(peaks)
		if (peaks(i) == 1)
			temp_start = max(1, i - sampleRate);
			max_deriv(i) = max(data_deriv(temp_start:i));
			start_deriv = theta * max_deriv(i);

			found = false;
			find_start = i;
			% has to peak within start_WT seconds
			while (found == false && find_start > (i - start_WT * sampleRate) && find_start > 0)
				if (data_deriv(find_start) < start_deriv)
					found = true;
					peak_start(find_start) = 1;
					peak_start_times(i) = times(find_start);
					rise_time(i) = times(i) - peak_start_times(i);
                end

				find_start = find_start - 1;
            end

			% If we didn't find a start
			if (found == false)
                if (i - (start_WT * sampleRate) > 0)
                    peak_start(i - (start_WT * sampleRate)) = 1;
                    peak_start_times(i) = times(i - start_WT * sampleRate);
                    rise_time(i) = start_WT;
                else
                    peak_start(1) = 1;
                    peak_start_times(i) = times(1);
                    rise_time(i) = start_WT;
                end
            end

			% Check if amplitude is too small
			if (thres_low > 0 && (data(i) - data(find(times == peak_start_times(i)))) < thres_low)
				peaks(i) = 0;
				peak_start(i) = 0;
				peak_start_times(i) = '';
				max_deriv(i) = 0;
				rise_time(i) = 0;
            end
        end
    end

	% Finding the end of the peak, amplitude of peak
	peak_end = zeros(length(data),1);
	peak_end_times = zeros(length(data),1);
	amplitude = zeros(length(data),1);
	decay_time = zeros(length(data),1);
	half_rise = zeros(length(data),1);
	SCR_width = zeros(length(data),1);

	for i = 1: length(peaks)
		if (peaks(i) == 1)
			peak_amp = data(i);
			start_amp = data(find(times == peak_start_times(i)));
			amplitude(i) = peak_amp - start_amp;

			half_amp = (amplitude(i) * 0.5) + start_amp;

			found = false;
			find_end = i;
			% has to decay within end_WT seconds
			while ((found == false) && (find_end < (i + (end_WT * sampleRate))) && (find_end < length(peaks)))
				if (data(find_end) < half_amp)
					found = true;
					peak_end(find_end) = 1;
					peak_end_times(i) = times(find_end);
					decay_time(i) = peak_end_times(i) - times(i);

					% Find width
					find_rise = i;
					found_rise = false;
					while (found_rise == false)
						if (data(find_rise) < half_amp)
							found_rise = true;
							half_rise(i) = times(find_rise);
							SCR_width(i) = peak_end_times(i) - times(find_rise);
                        end
						find_rise = find_rise - 1;
                    end

                elseif (peak_start(find_end) == 1)
					found = true;
					peak_end(find_end) = 1;
					peak_end_times(i) = times(find_end);
                end
				find_end = find_end + 1;
            end
        

			% If we didn't find an end
			if (found == false)
                if (i + (end_WT * sampleRate) < length(data))
                    [~,min_index] = min(data(i: i + (end_WT * sampleRate)));
                else
                    [~,min_index] = min(data(i: length(data)));
                end
                
				peak_end(i + min_index - 1) = 1;
				peak_end_times(i) = times(i + min_index - 1);
            end
        end
    end
            

	max_deriv = max_deriv * sampleRate; % now in change in amplitude over change in time form (uS/second)

    SCR = cell(1);
    SCR{1} = 1;
    
    j = 1;
    for i = 1:length(data)
        if (peaks(i) == 1)
            SCR{j} = [peak_times(i), peak_eda(i), rise_time(i), max_deriv(i), amplitude(i)];
        end
    end

end 