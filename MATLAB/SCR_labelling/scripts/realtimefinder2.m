function peaks = realtimefinder2(eda_csv, len_data, overlap, turn_offset, upturn_thres, downturn_thres, deriv_thres, peak_dist, base_dist, amp_thres, offset, alpha, method)
k = 1;

peaks = 0;
no_points = length(eda_csv);

for i = 1:round(len_data/overlap):(no_points-len_data)
    times = eda_csv(i:i+len_data-1, 1);
    data = eda_csv(i:i+len_data-1, 2);
    data_deriv = eda_csv(i:i+len_data-1, 5);
    data_deriv2 = eda_csv(i:i+len_data-1, 8);
    
    if (method == 1)
        SCR = findpeaks3(data, data_deriv, data_deriv2, times, turn_offset, upturn_thres, downturn_thres, deriv_thres, peak_dist, base_dist, amp_thres, offset);
    elseif (method == 2)
        SCR = findpeaks4(data, times, alpha, turn_offset, upturn_thres, downturn_thres, deriv_thres, peak_dist, base_dist, amp_thres, offset);
    end
    
    if(sum(SCR) ~= 0)
        if(isempty(find(peaks(:,1) == SCR(1))))
            peaks(k,1) = SCR(1);
            peaks(k,2) = SCR(2);
            peaks(k,3) = SCR(3);
            peaks(k,4) = SCR(4);
            
            %peaks(k) = [SCR(3),SCR(1)];
            %peaks(k+1) = [SCR(4),SCR(2)];

            k = k + 1;
        end
    end
end
end