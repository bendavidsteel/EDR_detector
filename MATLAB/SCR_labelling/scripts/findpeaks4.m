function SCR = findpeaks4(data, times, alpha, turn_offset, upturn_thres, downturn_thres, deriv_thres, peak_dist, base_dist, amp_thres, offset)
    SCR = zeros(1,4);
    
    data_filtered = data;
    
    a = alpha;
    for i = 2:length(data_filtered)
        data_filtered(i) = data_filtered(i-1) + a*(data_filtered(i) - data_filtered(i-1));
    end
    
    data_deriv = zeros(length(data)-1, 1);
    for i = 1:length(data_deriv)
        data_deriv(i) = data_filtered(i+1) - data_filtered(i);
    end
    
    data_deriv2 = zeros(length(data)-2, 1);
    for i = 2:length(data_deriv2)
        data_deriv2(i) = data_filtered(i+1) - (2*data_filtered(i)) + data_filtered(i-1);
    end

    [max_deriv2, max_index] = max(data_deriv2);
    [min_deriv2, min_index] = min(data_deriv2);
    
    
    if ((min_index - turn_offset > max_index) && (max_deriv2 > upturn_thres) && (min_deriv2 < downturn_thres))
        [max_deriv, ~] = max(data_deriv);
        if (max_deriv > deriv_thres)
            
            found_peak = false;
            found_base = false;
            
            for i = 1:peak_dist
                if (min_index+i < length(data))
                    if ((data_deriv(min_index+i) < 0) && (data_deriv(min_index+i-1) > 0))
                        peak_eda = data(i);
                        time_eda = times(i);
                        found_peak = true;
                        break
                    end
                end
            end
            
            for i = 1:base_dist
                if (max_index-i > 0)
                    if ((data_deriv(max_index-i) < 0) && (data_deriv(max_index-i+1) > 0))
                        base_eda = data(i);
                        base_time = times(i);
                        found_base = true;
                        break
                    end
                end
            end
            
            if ((found_peak == true) && (found_base == true) && (peak_eda - base_eda > amp_thres) && (time_eda - base_time > offset))
                SCR(1) = time_eda;
                SCR(2) = peak_eda;
                SCR(3) = base_time;
                SCR(4) = base_eda;
            end
        end
    end
end