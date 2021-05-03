#! /usr/bin/env python3

import argparse, sys, time

parser = argparse.ArgumentParser(
    description="This script split long reads into short reads based on the mean quality score of a sliding window, and then filter out short reads."	
)
parser.add_argument("fastq", metavar='<FASTQ>', type=str, help="the path to the reads fastq file/folder")
parser.add_argument("-w", "--window_size", metavar='<num>', type=int, help="the size of the sliding window (default: 250)", default=250)
parser.add_argument("-q", "--q_score", metavar='<num>', type=int, help="the average quality score of the sliding window (default: 18)", default=18)
parser.add_argument("-l", "--min_length", metavar='<num>', type=int, help="the minimum length of the sub-reads (default: 500)", default=500)

args = parser.parse_args(args=None if sys.argv[1:] else ['--help'])
fastq, window_size, q_score, min_length = args.fastq, args.window_size, args.q_score, args.min_length

#### create a dictionary that encode quality score
q_dic = dict(zip([chr(symbol) for symbol in range(33,127)], range(94)))

start_time = time.time()
#### import fastq file as a list of lists; [[name1,seq1,+,quality],[name2,seq2,+,quality]....]
with open(fastq, "r") as table:
     fastq_list = []
     fastq_seq = []
     for idx, line in enumerate(table):
         line = line.strip("\n")
         if idx%4 == 0 and idx != 0:
             fastq_list.append(fastq_seq)
             fastq_seq = []
         fastq_seq.append(line)
fastq_list.append(fastq_seq)    # save the last read

q_window = [25]*window_size
idx_range = [0,0] 

#### split
OUTPUT = open("after_split.fastq", "w")
for nano_read in fastq_list:
    range_start = 0    # to save the idx_range[1], which is the starting point of next split
    for idx, q_symbol in enumerate(nano_read[3]):
    
        #### keep the quality window updated 
        q_window.append(q_dic[q_symbol])
        del q_window[0]
        q_mean = sum(q_window)/window_size
        
        #### find the region of low quality
        ## 1) the region is (a, b)
        if q_mean < q_score and idx_range[0] == 0:
            idx_range[0]=idx
        if q_mean >= q_score and idx_range[0] > idx_range[1]:
            idx_range[1]=idx
        
        ## 2) no poor-quality region
        if idx == len(nano_read[3])-1:
            if idx_range == [0,0]:
                idx_range = [idx,idx+1]
                
        ## 3) the region is (a, end)
            else:
                idx_range[1] = idx+1
            
        ## export the region of high quality -- split at idx_range[0]
        if idx_range[0] < idx_range[1]:
            if idx_range[0] - range_start >= min_length:    # filter out short reads
                print(nano_read[0], "range: {}-{}".format(range_start, idx_range[0]), file=OUTPUT)
                print(nano_read[1][range_start:idx_range[0]+1], file=OUTPUT)
                print("+", file=OUTPUT)
                print(nano_read[3][range_start:idx_range[0]+1], file=OUTPUT)
            range_start = idx_range[1] 
            idx_range = [0,0]

OUTPUT.close()              
end_time = time.time()
print("running time: {} seconds".format(end_time-start_time))        


