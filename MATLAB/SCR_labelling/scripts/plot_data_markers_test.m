%load('plot_data_markers_test_workspace.mat') %w2_red_p
load('C:/Projects/Matlab/SCR_labelling/plot/plot_data_markers_workspace_w2_orange_c.mat') %w2_orange_c

%% Correct time values
disp('Correcting start time');
% Use formula: newvariablename = (variablename - start time)/1000
eda_time_corr = (eda_time-start_time)/1000;
hr_time_corr = (hr_time-start_time)/1000;
temp_time_corr = (temp_time-start_time)/1000;

%% Find time markers associated with each data type
eda_markers = [];
hr_markers = [];
temp_markers = [];
for i = 1:length(event_type)
    if(event_type(i) == 1)
        %eda
        eda_markers = [eda_markers,markers(i)];
    elseif(event_type(i) == 2)
       %hr
       hr_markers = [hr_markers,markers(i)];
    else
        %temp
        temp_markers = [temp_markers,markers(i)];
    end
end

%% Plot Caregiver Physiological Data
disp('Plotting data subplots');
title = strcat(name,' Physiological Data'); 
figure('Name',title)

eda_subplot = subplot(3,1,1);
plot(eda_time_corr,eda,color);grid on
hold 'on'
%{
if (~isempty(eda_markers)) %insert markers
    (vline(eda_markers,'k:')); 
end
%}
eda_sqi=(max(eda)+0.5)+eda_sqi; %insert SQI
plot(eda_time_corr,eda_sqi);
hold off
ylim([0 5]); %properties
ylabel('EDA');
xticks(0:100:2500);

hr_subplot = subplot(3,1,2);
plot(hr_time_corr,hr, color);grid on
hold on
%{
if (~isempty(hr_markers)) %insert markers
    vline(hr_markers,'k:');
end 
%}
hold off
ylabel('HR');
xticks(0:100:2500);

temp_subplot = subplot(3,1,3); 
%Find Peaks
[peak_val,peak_idx] = findpeaks(temp,'MinPeakProminence',0.5,'MaxPeakWidth',10); 
peak_time = temp_time_corr(peak_idx);
%Find Valleys
temp_inverted = -temp;
[valley_val,valley_idx] = findpeaks(temp_inverted,'MinPeakProminence',0.5,'MaxPeakWidth',10);
valley_time = temp_time_corr(valley_idx);
%Modify SQI
temp_sqi = (max(temp)+0.5)+temp_sqi;
%Create Subplot 
hold on
plot(temp_time_corr,temp,color,...
    peak_time,peak_val,'o',...
    valley_time,abs(valley_val),'*');
plot(temp_time_corr,temp_sqi);
grid on
%{
if (~isempty(temp_markers)) %insert markers
    vline(temp_markers,'k:'); 
end
%}
hold off
%Properties
xlabel('Time(s)');
ylabel('Temp(C)'); 
ylim([33 37]);
xticks(0:100:2500);

linkaxes([eda_subplot,hr_subplot,temp_subplot],'x');
xlim([120 inf]);

disp('Program complete!');