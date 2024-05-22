# kmerdemon

Our project develops a Python script similar to KmerGenie. It estimates optimal k-mer size and genome size from input reads of a haploid organism, accepting one or two fastQ files depending on read type. Users can specify k-mer length range, adjust k-mer sampling proportion as well as min and max values, and set an output file prefix. The tool outputs abundance histograms and estimates optimal k-mer size and genome size using a hash table for k-mer sampling and a dictionary for counting unique k-mers.

[Prerequisites](#prerequisites) | [Installation](#install) | [Basic Usage](#usage) 

<a name="prerequisites"></a>
## Prerequisites
`kmerdemon` requires the following python libraries to be installed:
- argparse
- random

Packages can be imported in Python
```
import argparse
import random
```

<a name="install"></a>
## Installation
To install `kmerdemon`, the following commands will be helpful

*add commands*

<a name="usage"></a>
## Basic Usage

To use `kmerdemon` the basic format is below
```
*FILL IN* kmerdemon [xyx]
```



- parse_fastq(): parses input file and outputs number of reads and writes reads to output file
- make_kmers(): given an input read and kmer size outputs a list of all kmers in the read
- abundance(): given a list of reads, selects a subset of reads for analysis and calculates raw abundance for each kmer in each read
- optimize(): find optimal kmer length given an input list of adundance data
- make_histogram(): given a dictionary of reads and their associates abundances, plot a histogram of unique kmers
- estimate_genome_size(): from analysis and input data, estimate genome size

