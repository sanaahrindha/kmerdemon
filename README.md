# kmerdemon

Our project entails the development of a tool similar to KmerGenie from lab2. This tool aims to estimate the optimal k-mer size and genome size based on a provided set of input reads. Implemented as a Python script, the tool will accept one or two fastQ files, depending on whether the reads are single end or paired end, respectively, from a sequence of a haploid organism.

Users will have the capability to specify both maximum and minimum k-mer lengths for analysis, adjust the proportion of sampled k-mers, and define the prefix for output files. Upon execution, the tool will generate a file containing information regarding the abundance histograms, showcasing the number of unique k-mers associated with various k-mer lengths, along with estimations for the optimal k-mer size and genome size.

Our tool will sample from the k-mers derived from the input reads. This sampling will be facilitated through a hash table, allowing the selection of a subset of k-mers for analysis. This approach not only enhances runtime efficiency but also optimizes the tool's performance. Subsequently, a dictionary will be utilized to tally the count of unique k-mers within the subset, aiding in the estimation of read coverage.

Ultimately, by deriving the count of unique reads from the subset across different potential k-mer lengths, our tool will enable the estimation of the optimal k-mer size, thereby offering valuable insights into genomic analysis.




