# Host CheckR

Host CheckR is a code written in R that performs host detection analysis on sequencing samples. It allows users to check the presence of host DNA in sequencing reads by comparing them to a reference database of host genomes. The code performs various steps such as database preparation, sample preparation, distance calculations, and result visualization.

## Dependencies

Host CheckR relies on the following R libraries:

- RColorBrewer
- kableExtra
- foreach
- doParallel
- ggplot2

Make sure you have these libraries installed before running the code.

## User-defined Variables

Before running the code, you need to set some user-defined variables based on your specific requirements. These variables are located at the beginning of the code:

- `make_database`: Determines whether to create the host database or not. Set it to "yes" if you want to create the database, or "no" if you already have the database.
- `sample_size`: Specifies the number of sequences to be checked in each sample.
- `prefix`: A prefix used for naming the output files.
- `sample_file`: The file containing information about the samples to be checked.

Make sure to modify these variables according to your needs.

## File Architecture

The code assumes a specific file architecture. It creates three main directories: "Host," "Sample," and "Result." If these directories don't exist, the code will create them automatically. Make sure you have the necessary permissions to create directories in the current working directory.

## Database Preparation

If `make_database` is set to "yes," the code will create the host database. The host database consists of reference genomes of different host organisms. The host organisms and their corresponding accessions are defined in the code. The code will generate a reference file for each host sample and perform the necessary steps to create the database, such as unzipping, combining fasta files, and making blast databases.

## Sample Preparation

The code reads the sample information from the specified `sample_file`. It determines the samples to be checked and prepares them for analysis. The code creates a directory for each sample, copies the fastq files, converts them to fasta format, and selects a subset of sequences based on the `sample_size` variable. Finally, it blasts the sample against each host database.

## Parallel Processing

To speed up the analysis, the code utilizes parallel processing. It detects the number of available CPU cores and creates a parallel cluster to run distance calculations in parallel. You can modify the `n.cores` variable to adjust the number of cores used for parallel processing.

## Result Visualization

After the analysis is completed, the code generates various visualizations to summarize the results. It creates a stacked barplot showing the proportion of reads assigned to different host organisms. It also generates a table summarizing the read proportions. Additionally, the code creates a PDF file containing the summary graph and a CSV file containing the detailed summary.

For the deer paper, the code compares the results to human VSP samples and generates an additional visualization specifically for this comparison.

## Output Files

The code generates several output files in the "Result" directory:

- `<prefix>_reference_sample.csv`: A CSV file containing the reference sample information.
- `<prefix>_summary-graph.pdf`: A PDF file containing the stacked barplot summarizing the read proportions.
- `<prefix>_summary-table.pdf`: A PDF file containing the table summarizing the read proportions.
- `<prefix>_summary.csv`: A CSV file containing the detailed summary of read proportions.
- `<prefix>_summary-graph_average-human.pdf`: A PDF file containing the comparison to human VSP samples.
- `<prefix>_host-read-origin.pdf`: A PDF file containing the grouped bar plot comparing human and deer samples.

Make sure to check these output files for the analysis results.

Please note that this readme file provides an overview of the code and its functionality. For detailed usage instructions and further customization, refer to the code comments and consult the relevant documentation for the R libraries used.
