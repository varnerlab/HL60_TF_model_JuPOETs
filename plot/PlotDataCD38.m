% Matlab script to plot data -

% Load the data -
data_array = load('scaled_data_CD38.dat');

% Load the measurements -
measurement_array = load('../data/CD38-WTHL60.txt');

% What is the size of the data?
[number_of_rows,number_of_cols] = size(data_array);

% initialize -
mean_array = [];
std_array = [];

for row_index = 1:number_of_rows

  % calculate the (mean,std) for this data -
  mean_value = mean(data_array(row_index,:));
  std_value = std(data_array(row_index,:));

  mean_array = [mean_array mean_value measurement_array(row_index,2)];
  std_array = [std_array std_value measurement_array(row_index,3)];
end

figure
hold on
bar(1:(2*number_of_rows),mean_array);
errorbar(1:(2*number_of_rows),mean_array,std_array,'.')
