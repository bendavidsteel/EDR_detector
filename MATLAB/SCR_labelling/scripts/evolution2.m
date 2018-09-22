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
mutation_rate = 0.1;
population = 20;
no_genes = 12;
no_parents = 1;
max_score = 0;
greatest_score = -100;
generation = 0;

chromosomes = zeros(population, no_genes); 
scores = zeros(population, 1);

%initial population
for i = 1:population
    len_data = round((200-20)*rand(1) + 20);
    overlap = round((10-1)*rand(1) + 1);
    turn_offset = round(20*rand(1));
    upturn_thres = 0.001*rand(1);
    downturn_thres = 0.001*rand(1);
    deriv_thres = 0.001*rand(1);
    peak_dist = round(30*rand(1));
    base_dist = round(30*rand(1));
    amp_thres = 1*rand(1);
    offset = round(30*rand(1));
    alpha = (1-0)*rand(1) + 0;
    method = round((2-1)*rand(1) + 1);
    
    chromosomes(i, :) = [len_data, overlap, turn_offset, upturn_thres, downturn_thres, deriv_thres, peak_dist, base_dist, amp_thres, offset, alpha, method];
    
    peaks = realtimefinder2(eda_csv, len_data, overlap, turn_offset, upturn_thres, downturn_thres, deriv_thres, peak_dist, base_dist, amp_thres, offset, alpha, method);
    
    scores(i) = fitness(edr_targets, peaks, len_data, overlap);
end

parents = zeros(no_parents, no_genes);

while true
    %choosing parents
    [scores_ranks, index_ranks] = sort(scores);
    
    for i = 1:no_parents
        parents(i,:) = chromosomes(index_ranks(population - i + 1), :);
    end
    
    %breeding next generation
    for i = 1:population
        gene_choices = round((no_parents - 1)*rand(no_genes,1) + 1);
      
        %crossover or mutation
        if (rand(1) < mutation_rate)
            len_data = round((200-20)*rand(1) + 20);
        else
            len_data = parents(gene_choices(1),1);
        end
        
        if (rand(1) < mutation_rate)
            overlap = round((10-1)*rand(1) + 1);
        else
            overlap = parents(gene_choices(2),2);
        end
        
        if (rand(1) < mutation_rate)
            turn_offset = round(20*rand(1));
        else
            turn_offset = parents(gene_choices(3),3);
        end
        
        if (rand(1) < mutation_rate)
            upturn_thres = 0.001*rand(1);
        else
            upturn_thres = parents(gene_choices(4),4);
        end
        
        if (rand(1) < mutation_rate)
            downturn_thres = 0.001*rand(1);
        else
            downturn_thres = parents(gene_choices(5),5);
        end
        
        if (rand(1) < mutation_rate)
            deriv_thres = 0.001*rand(1);
        else
            deriv_thres = parents(gene_choices(6),6);
        end
        
        if (rand(1) < mutation_rate)
            peak_dist = round(30*rand(1));
        else
            peak_dist = parents(gene_choices(7),7);
        end
        
        if (rand(1) < mutation_rate)
            base_dist = round(30*rand(1));
        else
            base_dist = parents(gene_choices(8),8);
        end
        
        if (rand(1) < mutation_rate)
            amp_thres = 1*rand(1);
        else
            amp_thres = parents(gene_choices(9),9);
        end
        
        if (rand(1) < mutation_rate)
            offset = round(30*rand(1));
        else
            offset = parents(gene_choices(10),10);
        end
        
        if (rand(1) < mutation_rate)
            alpha = (1-0)*rand(1) + 0;
        else
            alpha = parents(gene_choices(11),11);
        end
        
        if (rand(1) < mutation_rate)
            method = round((2-1)*rand(1) + 1);
        else
            method = parents(gene_choices(12),12);
        end
         
    
        chromosomes(i, :) = [len_data, overlap, turn_offset, upturn_thres, downturn_thres, deriv_thres, peak_dist, base_dist, amp_thres, offset, alpha, method];
    
        peaks = realtimefinder2(eda_csv, len_data, overlap, turn_offset, upturn_thres, downturn_thres, deriv_thres, peak_dist, base_dist, amp_thres, offset, alpha, method);
    
        scores(i) = fitness(edr_targets, peaks, len_data, overlap);
    end
    
    [max_score, index] = max(scores);
    
    if(max_score > greatest_score)
        greatest_score = max_score;
        greatest_animal = chromosomes(index, :);
    end
    
    max_score
    
    generation = generation + 1
end