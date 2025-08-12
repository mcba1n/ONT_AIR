import numpy as np
from scipy.io import savemat
from Bio import SeqIO
from ont_fast5_api.fast5_interface import get_fast5_file
import csv
from dtw import *
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import pysam
import h5py

## read in csv for Scrappie levels
f = {}
states_vec = []
f_vec = []
with open("scrappie_table.csv", 'r') as file:
    csvreader = csv.reader(file)
    header = next(csvreader)
    for row in csvreader:
        # make states an index for f_vec
        level = float(row[1])
        f[row[0]] = level
        f_vec.append(level)
        states_vec.append(row[0])

f_vec = (f_vec - np.mean(f_vec))/np.std(f_vec)

## input base sequences
tau = 5
samfile = pysam.AlignmentFile("calls2ref.bam", "rb")
n = 0
n_max = 40
lenfilt_reads = {}
lenfilt_ref_signals = {}
lenfilt_ref_states = {}
lenfilt_ref_isrev = {}

for x in samfile:
    read_id = x.query_name
    seq = x.query_sequence
    if seq is None or x.is_reverse==True:
        continue

    if x.reference_name=='chr1' and len(seq) >= 8000 and len(seq) <= 12000:

        # reference signal from Scrappie levels
        m = len(seq) - tau + 1
        level_seq = []
        state_seq = []
        for i in range(m):
            state = seq[i:i+tau]
            state_idx = states_vec.index(state)
            state_seq.append(state_idx)
            level_seq.append(f[state])

        # filter base distribution
        ct_vec = [0,0,0,0]
        for i in range(len(seq)):
            b = seq[i]
            if b=='A':
                ct_vec[0] += 1
            elif b=='T':
                ct_vec[1] += 1
            elif b=='C':
                ct_vec[2] += 1
            elif b=='G':
                ct_vec[3] += 1
        p = np.array(ct_vec)/len(seq)
        unif_dist_err = np.max(np.abs(p - 0.25))
        if unif_dist_err > 0.01:
            continue
    
        n += 1
        print(n)
        lenfilt_reads[read_id] = seq
        lenfilt_ref_signals[read_id] = np.array(level_seq)
        lenfilt_ref_states[read_id] = np.array(state_seq, dtype=int)

    if n == n_max:
        break

samfile.close()
print(len(lenfilt_reads), ' / ', n, 'reads left after filtering')
print('finished logging short reference sequences')

## output signals
file_path = 'fast5/all/'

# import signal from fast5 file (batch by batch in one chromosome)
N_batches = 37
signals = {}
k = 0
N_filt_reads = len(lenfilt_reads)
print('analysing', N_filt_reads, 'reads')

for batch_num in range(0,N_batches,1):
    print('reading batch', batch_num)
    fast5_file = 'batch' + str(batch_num) + '.fast5'
    with get_fast5_file(file_path + fast5_file, mode="r") as f5:
        for read in f5.get_reads():
            if read.read_id not in lenfilt_reads or len(read.get_raw_data())/len(lenfilt_ref_signals[read.read_id]) > 20:
               continue
            mean_qscore = read.handle['Analyses']['Basecall_1D_000']['Summary']['basecall_1d_template'].attrs['mean_qscore']
            print(read.read_id, mean_qscore)

           

            # read signal
            raw_sig = read.get_raw_data()
            ch = read.get_channel_info()
            ch_digitisation = ch['digitisation']
            ch_offset = ch['offset']
            ch_range = ch['range']
            ch_scale = ch_range / ch_digitisation
            sig = ch_scale * (raw_sig + ch_offset)
            sig = (sig - np.mean(sig))/np.std(sig)
            signals[read.read_id] = sig
            level_sig = lenfilt_ref_signals[read.read_id]

            # save files
            with open('short_genomic_dataset_fwd2/signal' + str(k+1) + '.csv', "ab") as f:
                numpy.savetxt(f, sig.reshape(1, sig.shape[0]), fmt='%1.3f', delimiter=",")
            with open('short_genomic_dataset_fwd2/states' + str(k+1) + '.csv', "ab") as f:
                state_seq = lenfilt_ref_states[read.read_id]
                numpy.savetxt(f, state_seq.reshape(1, state_seq.shape[0]), fmt='%d', delimiter=",")
            k += 1

print('finished generating short genomic dataset of', k, 'reads')

