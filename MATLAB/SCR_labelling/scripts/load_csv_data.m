path = 'C:/Projects/Matlab/SCR_labelling/';

file = 'eda_data/recording_config1_eda.csv';

filename = fullfile(path, file);

fileID = fopen(filename);

eda_csv_cell = textscan(fileID, '%s %q %q %q %q %q %q %q %q %q', 'delimiter', ';');

fclose(fileID);



file = '2,4,4,0.04peaks.txt';

filename = fullfile(path, file);

fileID = fopen(filename);

peaks_txt = textscan(fileID, '%s');

fclose(fileID);

len_eda_csv = length(eda_csv_cell{1});

eda_csv = zeros(len_eda_csv, 10);

for i = 1:10
    for j = 1:len_eda_csv
        eda_csv_cell{i}(j) = strip(eda_csv_cell{i}(j), 'both', '"');
        eda_csv(j,i) = double(string(cell2mat(eda_csv_cell{i}(j))));
    end
end

% a = 0.25;
% 
% eda_csv1 = eda_csv;
% 
% for i = 2:2184
%     eda_csv1(i,2) = eda_csv1(i-1,2) + a*(eda_csv1(i,2) - eda_csv1(i-1,2));
% end
% 
% a = 0.1;
% 
% eda_csv2 = eda_csv;
% 
% for i = 2:2184
%     eda_csv2(i,2) = eda_csv2(i-1,2) + a*(eda_csv2(i,2) - eda_csv2(i-1,2));
% end

j = 1;
for i = 1:length(peaks_txt{1})
    if ~contains(peaks_txt{1}(i), ':')
        peaks(j,1) = double(string(cell2mat(peaks_txt{1}(i))));
        indices = find(eda_csv(:,1) == peaks(j,1));
        peaks(j,2) = eda_csv(indices, 2);
        j = j + 1;
    end
end

%peaks = realtimefinder(eda_csv, 116, 2, 24, 32, 20, 0.042940806147176, 3.867153104220367e-04, 0.544226844100348);
%peaks = realtimefinder2(eda_csv, 30, 3, 0, 0, 1000, 0, 1000, 1000, 0, 0, 0, 2);

plot(eda_csv(:,1), eda_csv(:,2))
hold on
scatter(peaks(:,1), peaks(:,2))
hold off