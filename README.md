# NanoSplit
NanoSplit is intended to split long noisy reads into short reads based on the average quality score of a sliding window. The sub-reads with a length less than a certain threshold (500bp in default) will be filtered out. It takes **one** fastq file as input and produce one new fastq file in the current directory.

## Software requirement
* Python3

## Usage
```
usage: NanoSplit.py [-h] [-w <num>] [-q <num>] [-l <num>] <FASTQ>

This script split long reads into short reads based on the mean quality score of a sliding window, and then filter out
short reads.

positional arguments:
  <FASTQ>               the path to the reads fastq file/folder

optional arguments:
  -h, --help            show this help message and exit
  -w <num>, --window_size <num>
                        the size of the sliding window (default: 250)
  -q <num>, --q_score <num>
                        the average quality score of the sliding window (default: 18)
  -l <num>, --min_length <num>
                        the minimum length of the sub-reads (default: 500)
```

## Examples
```
NanoSplit all_barcode3_reads.fastq
NanoSplit.py -w 150 -q 20 -l 1000 all_barcode3_reads.fastq
```
