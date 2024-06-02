import argparse
import random
import matplotlib.pyplot as plt
from collections import defaultdict
import os

def parse_fastq(file_path, out_file): #parse a fastq file and return the number of reads and write reads to output file
    """
    Parse the input file to make it readable for the program

    Parameters
    ----------
    file_path : path
        path to the input fastq file

    out_file : path
        path to the output file used in by the program

    Returns
    -------
    num_reads : int
        number of reads counted in input file
    """
    with open(file_path, 'r') as file1, open(out_file, "w") as file2:
        num_reads = 0
        lines = file1.readlines()
        for line in lines:
            if line.startswith('@'):
                file2.write(line[1:])
                num_reads += 1
    return num_reads

def make_kmers(read, kmer_size):
    """
    Make a list of kmers from an input read

    Parameters
    ----------
    read : str
        single read from the input file

    kmer_size : int
        size of kmer

    Returns
    -------
    kmers : list
        list of strings of every kmer in the input string
    """
    kmers = [0]*(len(read)-kmer_size+1)
    for i in range(len(kmers)):
        kmers[i] = read[i:i+kmer_size]
    return kmers


def abundance(parsed_file, threshold, size): #given parsed file from parse_fastq() create subset of reads to analyze
    """
    Creates an abundance dictionary for all the kmers of specified size for a sample of the reads of the input file

    Parameters
    ----------
    parsed_file : path
        path to the parsed file

    threshold : float
        proportion of reads to be sampled

    size : int
        kmer size 

    Returns
    -------
    kmer_frequences : dictionary
        kmer strings as keys and counts of kmers as values
    """
    reads = []
    kmer_frequencies = {}
    with open(parsed_file, "r") as file1:
        lines = file1.readlines()
        for line in lines:
            rand = random.random()
            if rand < threshold:
                reads.append(line.strip()) 
    for read in reads:
        if (len(read)<size):
            raise ValueError("Maximum k-mer size cannot be larger than the length of the reads")
        kmers = make_kmers(read, size)
        for kmer in kmers:
            if kmer not in kmer_frequencies:
                kmer_frequencies[kmer] = 1
            else: 
                kmer_frequencies[kmer] += 1
    return kmer_frequencies 

def estimate_genome_size(num_unique_kmers, kmer_size, coverage):
    """
    Estimates genome size based on optimal number of unique kmers

    Parameters
    ----------
    num_unique_kmers : int
        maximum number of unique kmers

    kmer_size : int
        size of kmer

    coverage : int
        coverage of reads 

    Returns
    -------
    estimated_genome_size : int
        estimated genome size
    """ 
    estimated_genome_size = num_unique_kmers * kmer_size / coverage
    return estimated_genome_size


def create_histogram(kmer_counts, k):
    """
    Creates a histogram for an input dictionary

    Parameters
    ----------
    kmer_counts : dictionary
        path to the parsed file

    threshold : float
        proportion of reads to be sampled

    size : int
        kmer size 
    """
    histogram = defaultdict(int)
    for count in kmer_counts.values():
        histogram[count] += 1
    plot_histogram(histogram, k)
    
def plot_histogram(histogram, k):     
    """
    Using matplotlib, plots histogram dict with a k value

    Parameters
    ----------
    histogram : dictionary of dictionaries
        kmer length as keys and dictionaries as values
        for values: kmer strings as keys and counts of kmers as values

    k : int
        kmer size for which the histogram is being plotted

    Outputs
    -------
    Graph : plt
        creates the histogram with all required elements and saves it
        as a png to your directory
    """              
    abundances = sorted(histogram.keys())
    frequencies = [histogram[a] for a in abundances]
    plt.figure(figsize=(10,6))
    plt.bar(abundances, frequencies, width=1.0, edgecolor="black")
    plt.xlabel('K-mer Abundance')
    plt.ylabel('Frequency')
    plt.title(f'K-mer Abundance Histogram for k={k}')
    plt.yscale('log')
    plt.savefig(f'kmer_histogram_k{k}.png')
    plt.close()

def predict_best_k(histograms):
    """
    Gets the best kmer length given the abundance information for different kmer lengths

    Parameters
    ----------
    histograms : dictionary of dictionaries
        kmer length as keys and dictionaries as values
        for values: kmer strings as keys and counts of kmers as values

    Returns
    -------
    optimal_kmer_length : int
        kmer length that had the most unique kmers
    """    
    max_unique_kmers = 0
    for kmer_size, histogram in histograms.items():
        if len(histogram) > max_unique_kmers:
            max_unique_kmers = len(histogram)
            optimal_kmer_length = kmer_size
    return optimal_kmer_length


def main():
    parser = argparse.ArgumentParser(prog="kmerdemon", description="Tool to estimate optimal genome and k-mer size") #tool description 
    
    #Input 
    parser.add_argument("input_files", nargs="+", help="Input FASTQ files") #input file(s)
    
    #Output
    parser.add_argument("-o", "--out", default="output", help="Prefix for output files." "Default: output", metavar="FILE", type=str, required=False)

    #Options
    parser.add_argument("-l","--min_kmer_size", type=int, default=15, help="Minimum k-mer length for analysis", metavar="int", required=False) #optional, default set
    parser.add_argument("-k","--max_kmer_size", type=int, default=121, help="Maximum k-mer length for analysis", metavar="int", required=False) #optional, default set
    parser.add_argument("-e","--kmer_sampling_proportion", type=float, default=0.01, help="Proportion of k-mers to sample", metavar="float", required=False) #optional, default set
    
    #Parse args
    args = parser.parse_args() #process arguments

    files = args.input_files
    print(files)
    min_kmer_size = args.min_kmer_size
    print(min_kmer_size)
    max_kmer_size = args.max_kmer_size
    print(max_kmer_size)
    sampling_proportion = args.kmer_sampling_proportion
    output_prefix = args.out

    if min_kmer_size < 5:
        parser.error("Minimum k-mer size cannot be less than 5")

    if sampling_proportion <= 0 or sampling_proportion > 1:
        parser.error("Sampling proportion must be in the range (0,1]")

    if len(files) == 0 or len(files) > 2:
        parser.error("Requires 1 or 2 input files")

    num_reads = 0
    for file_path in files:
        print(file_path)
        if not os.path.exists(file_path):
            parser.error("File does not exist")
        file_name = os.path.basename(file_path)
        file_prefix, _ = os.path.splitext(file_name)
        output_file = f"{file_prefix}_parsed.txt"
        num_reads += parse_fastq(file_path, output_file)

    increment = int((max_kmer_size - min_kmer_size)/10)
    if increment == 0:
        increment = 1
    kmer_frequencies_by_size = {}
    # get abundance data for each kmer size in range from min to max, add to dictionary
    # keep track of optimal kmer length
    num_unique_kmers = 0
    optimal_kmer_length = -1
    for kmer_size in range(min_kmer_size, max_kmer_size, increment):
        kmer_frequencies_by_size[kmer_size] = abundance(output_file,sampling_proportion,kmer_size)
    optimal_kmer_length = predict_best_k(kmer_frequencies_by_size)
    num_unique_kmers = len(kmer_frequencies_by_size[optimal_kmer_length])

    if increment != 1:
        kmer_frequencies_by_size_better = {}
        for kmer_size in range(optimal_kmer_length-5, optimal_kmer_length+5, 2):
            kmer_frequencies_by_size_better[kmer_size] = abundance(output_file,sampling_proportion,kmer_size)
        optimal_kmer_length = predict_best_k(kmer_frequencies_by_size_better)
        num_unique_kmers = len(kmer_frequencies_by_size_better[optimal_kmer_length])

    print("Optimal k-mer length: ",optimal_kmer_length)
    print("Number of Unique k-mers: ",num_unique_kmers)
        
    output_file = f"{output_prefix}_kmer_{kmer_size}.txt"
    with open(output_file, 'w') as output_file:
        output_file.write(f"K-mer Size: {kmer_size}\n")
        output_file.write(f"Number of Unique K-mers: {num_unique_kmers}\n")
        output_file.write(f"Estimated Genome Size: {estimated_genome_size}\n")
            

if __name__ == "__main__":
    main()
    

#python kmer_estimator.py input1.fastq input2.fastq --min_kmer_size 15 --max_kmer_size 120 --output_prefix myoutput

