clear all; close all; clc;
addpath('log_semiring');

% scrappie levels
load('../DNA_dataset/scrappie_table.mat');
f = (f-mean(f))/std(f);
sigma = 0.5;
read_num = 11;

for block_num = 1:113
block_num

% load signal
y = csvread(['../DNA_dataset/chopped_reads_',num2str(read_num),'_fwd/signal',num2str(block_num),'.csv']);
% y = (y-mean(y))/std(y);
s = csvread(['../DNA_dataset/chopped_reads_',num2str(read_num),'_fwd/states',num2str(block_num),'.csv']);
s = s+1;
m = length(s)-1;

% decoding
log_prob_tot = fwd1(y, m, f, A, sigma, s(1));
log_prob = fwd2(y, m, f, sigma, s(2:end));

log_post = log_prob - log_prob_tot;
H = -(1/m)*log_post/log(2);
I = 2 - H

% save
fid = fopen(['../DNA_dataset/chopped_reads_',num2str(read_num),'_fwd/rates.txt'], 'a+');
fprintf(fid, '%d %0.2f %0.4f\n', block_num, sigma, I);
fclose(fid);

end