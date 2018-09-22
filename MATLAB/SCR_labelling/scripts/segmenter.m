path = 'C:/Projects/Matlab/SCR_labelling/';

file = 'eda_data/TP001376_c_eda.csv';

filename = fullfile(path, file);
fileID = fopen(filename);

eda_csv_cell = textscan(fileID, '%s %s %s %s %s %s %s %s %s %s', 'delimiter', ',');

fclose(fileID);

len = length(eda_csv_cell{1});

eda_csv = zeros(len, 10);

for i = 1:10
    for j = 1:len
        eda_csv_cell{i}(j) = strip(eda_csv_cell{i}(j), 'both', '"');
        eda_csv(j,i) = double(string(cell2mat(eda_csv_cell{i}(j))));
    end
end

a = 0.1;

for i = 2:len
    eda_csv(i,2) = eda_csv(i-1,2) + a*(eda_csv(i,2) - eda_csv(i-1,2));
end

eda_samples = zeros(ceil((len/10))+1,50);

j = 1;
for i = 1:10:length(eda_csv)-50
    m = eda_csv(i:i+50,2)';
    eda_samples(j,:) = eda_csv(i:i+49,2)';
    
    plot(eda_csv(i:i+50,1), eda_csv(i:i+50,2))
    
    currkey=0;
    while currkey~=1
        pause;
        currkey=get(gcf,'CurrentKey'); 
        if strcmp(currkey, 'y')
            eda_peak(j) = 1;
            currkey = 1;
        elseif strcmp(currkey, 'n')
            eda_peak(j) = 0;
            currkey = 1;
        else
            currkey = 0;
        end
    end
    
    j = j + 1;
end