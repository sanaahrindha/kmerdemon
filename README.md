# kmerdemon

Our project develops a Python script similar to KmerGenie. It estimates optimal k-mer size and genome size from input reads of a haploid organism, accepting one or two fastQ files depending on read type. Users can specify k-mer length range, adjust k-mer sampling proportion as well as min and max values, and set an output file prefix. The tool outputs abundance histograms and estimates optimal k-mer size and genome size using a hash table for k-mer sampling and a dictionary for counting unique k-mers.

[Prerequisites](#prerequisites) | [Installation](#install) | [Basic Usage](#usage) 

<a name="prerequisites"></a>
## Prerequisites
`kmerdemon` requires the following python libraries to be installed:
- argparse
- random

Packages can be imported in Python

<a name="install"></a>
## Installation
To install `kmerdemon`, the following commands will be helpful

*add commands*

<a name="usage"></a>
## Basic Usage

To use `kmerdemon` the basic format is below




