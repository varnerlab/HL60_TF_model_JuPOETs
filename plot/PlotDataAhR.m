% Matlab script to plot data -

% Load the data -
data_array = load('scaled_data_AhR.dat');

% What is the size of the data?
[number_of_rows,number_of_cols] = size(data_array);

% initialize -
mean_array = [];
std_array = [];

for row_index = 1:number_of_rows

  % calculate the (mean,std) for this data -
  mean_value = mean(data_array(row_index,:));
  std_value = std(data_array(row_index,:));

  mean_array = [mean_array mean_value];
  std_array = [std_array std_value];
end

figure
hold on
bar(1:number_of_rows,mean_array);
errorbar(1:number_of_rows,mean_array,std_array,'.')
