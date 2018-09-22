len_data = 40;
overlap = 2;
turn_offset = 0;
upturn_thres = 0;
downturn_thres = 0;
deriv_thres = 0;
peak_dist = 0;
base_dist = 0;
amp_thres = 0;
offset = 0;
alpha = 0;
method = 1;

path = 'C:/Projects/Matlab/SCR_labelling/';
file = 'eda_data/recording_config1_eda.csv';

filename = fullfile(path, file);
fileID = fopen(filename);

eda_csv_cell = textscan(fileID, '%s %q %q %q %q %q %q %q %q %q', 'delimiter', ';');

fclose(fileID);

no_points = length(eda_csv_cell{1});

eda_csv = zeros(no_points , 10);

for i = 1:10
    for j = 1:no_points
        eda_csv_cell{i}(j) = strip(eda_csv_cell{i}(j), 'both', '"');
        eda_csv(j,i) = double(string(cell2mat(eda_csv_cell{i}(j))));
    end
end


fig = figure;
ax = axes('Parent',fig,'position',[0.13 0.39  0.77 0.54]);
peaks = realtimefinder2(eda_csv, len_data, overlap, turn_offset, upturn_thres, downturn_thres, deriv_thres, peak_dist, base_dist, amp_thres, offset, alpha, method);

data_plot = plot(eda_csv(:,1), eda_csv(:,2));
hold on 
edr_plot = scatter(eda_csv(1,1), eda_csv(1,2));
if (sum(peaks) ~= 0)
    edr_plot = scatter([peaks(:,3), peaks(:,1)], [peaks(:,4), peaks(:,2)]);
end


amp_thres_slider = uicontrol('Parent',fig,'Style','slider','Position',[81,54,419,23],...
              'value',amp_thres, 'min',0, 'max',1,...
              'callback', @(src,eventdata)ampThresSliderCallBack(src, eventdata, edr_plot, eda_csv, len_data, overlap, turn_offset, upturn_thres, downturn_thres, deriv_thres, peak_dist, base_dist, amp_thres, offset, alpha, method));
bgcolor = fig.Color;
amp_thres_slider1 = uicontrol('Parent',fig,'Style','text','Position',[50,54,23,23],...
                'String','0','BackgroundColor',bgcolor);
amp_thres_slider2 = uicontrol('Parent',fig,'Style','text','Position',[500,54,23,23],...
                'String','1','BackgroundColor',bgcolor);
amp_thres_slider3 = uicontrol('Parent',fig,'Style','text','Position',[240,25,100,23],...
                'String','Amplitude Threshold','BackgroundColor',bgcolor);
            
            function ampThresSliderCallBack(source, event, edr_plot, eda_csv, len_data, overlap, turn_offset, upturn_thres, downturn_thres, deriv_thres, peak_dist, base_dist, amp_thres, offset, alpha, method)
                amp_thres = source.Value;
                peaks = realtimefinder2(eda_csv, len_data, overlap, turn_offset, upturn_thres, downturn_thres, deriv_thres, peak_dist, base_dist, amp_thres, offset, alpha, method);
                delete(edr_plot)
                if (sum(peaks) ~= 0)
                    edr_plot = scatter([peaks(:,3), peaks(:,1)], [peaks(:,4), peaks(:,2)]);
                end
            end
            
            
            
