# kmerdemon

Our project develops a Python script similar to KmerGenie. It estimates optimal k-mer size and genome size from input reads of a haploid organism, accepting one or two fastQ files depending on read type. Users can specify k-mer length range, adjust k-mer sampling proportion as well as min and max values, and set an output file prefix. The tool outputs abundance histograms and estimates optimal k-mer size and genome size using a hash table for k-mer sampling and a dictionary for counting unique k-mers.

[Prerequisites](#prerequisites) | [Installation](#install) | [Basic Usage](#usage) 

<a name="prerequisites"></a>
## Prerequisites
`kmerdemon` requires the following python libraries to be installed:
- setuptools
- matplotlib

Packages can be imported in Python
```
import argparse
import random
import matplotlib.pyplot as plt
from collections import defaultdict
```

<a name="install"></a>
## Installation
To install `kmerdemon`, the following commands will be helpful

```bash
git clone https://github.com/sanaahrindha/kmerdemon.git
cd kmerdemon
pip install -r requirements.txt
pip install . # For user installation, add --user
```

<a name="usage"></a>
## Basic Usage

```bash
usage: kmerdemon [-h] [-o FILE] [-l int] [-k int] [-e float] input_files 

Tool to estimate optimal genome and k-mer size

positional arguments:
  input_files           Input FASTQ files

options:
  -h, --help            show this help message and exit
  -o FILE, --out FILE   Prefix for output files.Default: output
  -l int, --min_kmer_size int
                        Minimum k-mer length for analysis
  -k int, --max_kmer_size int
                        Maximum k-mer length for analysis
  -e float, --kmer_sampling_proportion float
                        Proportion of k-mers to sample
```



