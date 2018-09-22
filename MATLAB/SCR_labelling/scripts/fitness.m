function score = fitness(edr_targets, peaks, len_data, overlap)
    score = 0;
    
    matches = zeros(length(peaks),1);
    target_matches = zeros(length(edr_targets),1);
    
    if (sum(peaks) == 0)
        score = -100;
    else
        [no_peaks,~] = size(peaks);
        for i = 1:no_peaks
            matches(i) = 0;
            for j = 1:length(edr_targets)
                diff = abs(edr_targets(j,1) - peaks(i,1));
                if (diff < 100)
                    score = score + 10;
                    matches(i) = j;
                    target_matches(j) = 1;
                elseif (diff < 200)
                    score = score + 5;
                    matches(i) = j;
                    target_matches(j) = 1;
                elseif (diff < 300)
                    score = score + 2;
                    matches(i) = j;
                    target_matches(j) = 1;
                elseif (diff < 500)
                    score = score + 1;
                    matches(i) = j;
                    target_matches(j) = 1;
                end
            end

            if (matches(i) == 0) 
                score = score - 20;
            end
        end

        hist = zeros(length(edr_targets),1);
        for i = 1:length(matches)
            if (matches(i) ~= 0)
                hist(matches(i)) = hist(matches(i)) + 1;
            end
        end
        for i = 1:length(hist)
            if (hist(i) > 1)
                score = score - 10;
            end
        end
        
        score = score + 5*sum(target_matches);
        
        score = score - 1*overlap/len_data;
    end
end