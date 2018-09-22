function peaks = realtimefinder(eda_csv, len_data, overlap, offset, start_WT, end_WT, thres_low, theta, alpha)
k = 1;

peaks = zeros(1,5);
no_points = length(eda_csv);

for i = len_data:round(len_data/overlap):(no_points-len_data)
    times = eda_csv(i-len_data+1:i, 1);
    data = eda_csv(i-len_data+1:i, 2);
    %p = polyfit(eda_csv(i-len_data+1:i, 1), eda_csv(i-len_data+1:i, 2), poly_degree);
    %p = polyfit(eda_csv(i-len_data+1:i, 1), (1:len_data)', poly_degree);
    %q = polyder(p);
    %data_deriv = polyval(q, times);
    
    SCRs = findpeaks2(data, times, offset, start_WT, end_WT, thres_low, theta, alpha);
%     if (method == 1)
%         SCRs = findpeaks(data, data_deriv, times, offset, start_WT, end_WT, thres_low);
%     elseif (method == 2)
%         SCRs = findpeaks2(data, times, offset, start_WT, end_WT, thres_low, theta, alpha);
%     end
    
    if(SCRs{1} ~= 1)
        for j = 1:length(SCRs)
            if(isempty(find(peaks(:,1) == SCRs{j}(1))))
                peaks(k,1) = SCRs{j}(1);
                peaks(k,2) = SCRs{j}(2);

                k = k + 1;
            end
        end
    end
end
end