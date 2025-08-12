clear all; close all; clc;

load('scrappie_table.mat');
f = (f-mean(f))/std(f);

read_num = 10;

y = csvread(['../genomic_dataset/short_genomic_dataset_fwd2/signal',num2str(read_num),'.csv']);
y = (y-mean(y))/std(y);
s = csvread(['../genomic_dataset/short_genomic_dataset_fwd2/states',num2str(read_num),'.csv']);
m = length(s);
level_vec = f(s+1);
N = length(f);

length(y)/m


% alignment
[best_cost, dtw_matrix] = dtw_simple(level_vec,y);
best_cost/length(y)
[x_stretched_hat,ix] = dtw_best_path(level_vec,dtw_matrix);
plot(x_stretched_hat), hold on;plot(y)

% segmentation
segments = cell(1,m);
t_prev = 0;
lvl_samples = [];
durations = [];
level_mean_errs = [];
m_idx = 1;
for n = 1:length(y)
    lvl_samples = [lvl_samples, y(n)];
    if n==length(y) || ix(n) ~= ix(n+1)
        segments{1,m_idx} = lvl_samples;
        durations = [durations, length(lvl_samples)];
        level_mean_errs = [level_mean_errs, lvl_samples - level_vec(m_idx)];
        t_prev = n;
        m_idx = m_idx + 1;
        lvl_samples = [];
    end
end
jump_times = [0,cumsum(durations)];

%% split data
clr_len = 10;
m_block = 100;
N_reads = floor((m-2*clr_len)/m_block);
chopped_reads = cell(N_reads);

for bl_idx = 0:N_reads-1
    y1 = [];
    s1 = s(clr_len+bl_idx*m_block+1:clr_len+(bl_idx+1)*m_block);
    s1_s0 = s(clr_len+bl_idx*m_block);

    for i = clr_len+bl_idx*m_block+1:clr_len+(bl_idx+1)*m_block
        y1 = [y1,segments{1,i}];
    end
    
    csvwrite(['chopped_reads_',num2str(read_num),'_fwd2/signal',num2str(bl_idx+1),'.csv'], y1)
    csvwrite(['chopped_reads_',num2str(read_num),'_fwd2/states',num2str(bl_idx+1),'.csv'], [s1_s0,s1]);
end

csvwrite(['chopped_reads_',num2str(read_num),'_fwd2/dtw.txt'], best_cost/length(y));

%% plot signal

