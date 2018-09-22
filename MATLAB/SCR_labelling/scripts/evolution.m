path = 'C:/Projects/Matlab/SCR_labelling/';
file = 'eda_data/recording_config1_eda.csv';

filename = fullfile(path, file);
fileID = fopen(filename);

eda_csv_cell = textscan(fileID, '%s %q %q %q %q %q %q %q %q %q', 'delimiter', ';');

fclose(fileID);

file = 'peak_locations.txt';

filename = fullfile(path, file);
fileID = fopen(filename);

edr_targets_txt = textscan(fileID, '%s');

fclose(fileID);

no_points = length(eda_csv_cell{1});

eda_csv = zeros(no_points , 10);

for i = 1:10
    for j = 1:no_points
        eda_csv_cell{i}(j) = strip(eda_csv_cell{i}(j), 'both', '"');
        eda_csv(j,i) = double(string(cell2mat(eda_csv_cell{i}(j))));
    end
end

no_targets = length(edr_targets_txt{1});
edr_targets = zeros(no_targets,2);

for i = 1:no_targets
    edr_targets(i,1) = double(string(cell2mat(edr_targets_txt{1}(i))));
    indices = find(eda_csv(:,1) == edr_targets(i,1));
    edr_targets(i,2) = eda_csv(indices, 2);
end

%evolution parameters
mutation_rate = 0.2;
population = 20;
no_genes = 8;
no_parents = 2;
max_score = 0;
greatest_score = 0;
generation = 0;

chromosomes = zeros(population, no_genes); 
scores = zeros(population, 1);

%initial population
for i = 1:population
    %Setting range that arguments can be randomly generated in
    len_data = round((200-20)*rand(1) + 20);
    %len_data = 40;
    overlap = round((10-1)*rand(1) + 1);
    offset = round((40-0)*rand(1) + 0);
    start_WT = round((40-0)*rand(1) + 0);
    end_WT = round((40-0)*rand(1) + 0);
    thres_low = (1-0)*rand(1) + 0;
    theta = 0.1*rand(1);
    alpha = (1-0)*rand(1);
    
    %storing randomly generated arguments
    chromosomes(i, :) = [len_data, overlap, offset, start_WT, end_WT, thres_low, theta, alpha];
    %running algorithm with a set of arguments
    peaks = realtimefinder(eda_csv, len_data, overlap, offset, start_WT, end_WT, thres_low, theta, alpha);
    %comparing algorithm results with targets
    scores(i) = fitness(edr_targets, peaks, len_data, overlap);
end

parents = zeros(no_parents, no_genes);

while true
    %choosing parents based on which algorithm performed best
    [scores_ranks, index_ranks] = sort(scores);
    %saving best performing argument sets
    for i = 1:no_parents
        parents(i,:) = chromosomes(index_ranks(population - i + 1), :);
    end
    
    %creating new population of algorithms with arguments either taken from
    %one of the parents, or re randomly generated
    for i = 1:population
        gene_choices = round((no_parents - 1)*rand(no_genes,1) + 1);
      
        %crossover or mutation
        if (rand(1) < mutation_rate)
            len_data = round((200-20)*rand(1) + 20);
            %len_data = 40;
        else
            len_data = parents(gene_choices(1),1);
        end
        
        if (rand(1) < mutation_rate)
            overlap = round((10-1)*rand(1) + 1);
        else
            overlap = parents(gene_choices(2),2);
        end
        
        if (rand(1) < mutation_rate)
            offset = round((40-0)*rand(1) + 0);
        else
            offset = parents(gene_choices(3),3);
        end
        
        if (rand(1) < mutation_rate)
            start_WT = round((40-0)*rand(1) + 0);
        else
            start_WT = parents(gene_choices(4),4);
        end
        
        if (rand(1) < mutation_rate)
            end_WT = round((40-0)*rand(1) + 0);
        else
            end_WT = parents(gene_choices(5),5);
        end
        
        if (rand(1) < mutation_rate)
            thres_low = (1-0)*rand(1) + 0;
        else
            thres_low = parents(gene_choices(6),6);
        end
        
        if (rand(1) < mutation_rate)
            theta = 0.1*rand(1);
        else
            theta = parents(gene_choices(7),7);
        end
        
        if (rand(1) < mutation_rate)
            alpha = round((1-0)*rand(1) + 0);
        else
            alpha = parents(gene_choices(8),8);
        end
    
        chromosomes(i, :) = [len_data, overlap, offset, start_WT, end_WT, thres_low, theta, alpha];
        %running algorithm
        peaks = realtimefinder(eda_csv, len_data, overlap, offset, start_WT, end_WT, thres_low, theta, alpha);
        %scoring results
        scores(i) = fitness(edr_targets, peaks, len_data, overlap);
    end
    
    [max_score, index] = max(scores);
    %if this population features the best scorer of all time, save
    %argument set
    if(max_score > greatest_score)
        greatest_score = max_score;
        greatest_animal = chromosomes(index, :);
    end
    %display score and generation number
    max_score
    
    generation = generation + 1
end



