##############
# MetaDamage #
##############

# Version 2.8
# 2022-01-02
# Rosie Everett, Becky Cribdon

The MetaDamage analysis estimates base substitutions of double-stranded libraries, focusing on C->T and G->A, by comparing them to reference sequences from GenBank. The MetaDamage also has the feature to assess single-stranded libraries. Each query sequence is aligned to its own reference.

The BLAST and reference-retrieval steps can take considerable time, so MetaDamage is best run on a server.

Three test FASTAs are included under Example/. test1.fasta has 6 reads, test2.fasta has 341, and test3.fasta has 0. They should produce the same output as in Example/Output/.

Pre-requisites
===========
The MetaDamage tool requires a conda environment to process. To install:

conda install -n perl -c bioconda perl-lwp-protocol-https

Quick start
===========
Query sequences must be present as one or more FASTA files. Collect all FASTAs in one directory alongside the MetaDamage scripts.

The meta-script will run the full analysis on all FASTA files in the current directory. It takes two arguments. The first, number of threads, is obligatory. The second, a path to a local BLAST database, is optional.

To run: 

1) enter new conda environment:

source activate perl

2) run MetaDamage:

perl MetaDamage_v2.8.pl -fastas [input_fasta] -threads [number of processing threads] 

Outputs:
[FASTA].aln.txt - Alignments for each query sequence and its reference.
[FASTA].blast.txt - BLAST output.
[FASTA].mismatches.txt - Proportions of base mismatches at the 5' and 3' ends across query sequences.
[FASTA].paired.txt - Contains each query sequence paired with its reference. Records are delimited by "@\n".
efetch_error_log.txt - If present, stores any warnings and errors produced by efetch during the reference-retrieval step.
MetaDamage_plots.pdf - Plots proportions of mismatches on the 5' and 3' ends across query sequences. Includes the 95% confidence intervals of C->T and G->A mismatches at the very end positions.
MetaDamage_CIs.txt - For those very end positions, lists the point estimates and 95% CI upper and lower bounds.


Flags
===========

USAGE: perl MetaDamage_v2.8.pl -flag [options]
FLAGS
-fastas : fasta input file, essential for all options
-threads : number of threads for blast algorithm to use
-blast_db : database to search. Optional, will default to searching NCBI if no database specified. This is SLOW!
-batches : batch size of sequences to run in analysis to make processing more efficient for high numbers of queries, default batch size 100
-blastfile : option to input a prepared blast output and skip blast search
-single_strand : TRUE to generate mismatch profiles suitable for interpretation with single stranded library preps
-out : option for output files root name
-help : help option to generate this output.


Additional notes 
================

1 Running BLAST
---------------
Outputs:
[input].blast.txt

A very simple BLAST. Output is in tabular format with an additional field for query length.

Flag -blastfile is an option to input a prepared blast output and skip blast search. The input must have followed the following format to create a tabular output with the additional field for query length:

blastn -db  [nucleotide database] -num_threads [x] -query [input FASTA] -out [output BLAST] -max_target_seqs 1 -max_hsps 1 -outfmt "6 std qlen"


2 Retrieving reference sequences using either blastdbcmd or efetch
------------------------------------------------------------------
Outputs:
[input].paired.txt

If you provide a local BLAST database, this step will use blastdbcmd to retrieve reference sequences from there. This is faster and more reliable.

If you do not provide a local BLAST database, this step will use the NCBI efetch tool (part of Entrez Direct) to retrieve reference sequences from GenBank. This is prone to warnings, but should continue running. Any warnings and errors are stored in efetch_error_log.txt.

You may need to install Entrez Direct: https://www.ncbi.nlm.nih.gov/books/NBK179288/
You can also manually download edirect.tar.gz from the FTP: https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/
Or install through conda: https://anaconda.org/bioconda/entrez-direct


3 Aligning query and reference sequences
----------------------------------------
Outputs:
[input].aln.txt

This stage uses the Needleman-Wunsch algorithm (https://github.com/sestaton/sesbio/blob/master/phylogenetics/needleman-wunsch.pl) to align every pair of sequences in the paired.txt file.

If you wish to merge samples, you can do so after alignment. Simply concatenate their .aln files.


4 Summarising mismatches
------------------------
Outputs:
[input].mismatches.txt

For each of the first and last 25 bases of every alignment, this script calculates the proportion of sequences where the reference has one base, but the query has a different base. It calculates proportions for all twelve possible mismatches.


5 Plotting mismatches
---------------------
Outputs:
MetaDamage_plots.pdf
MetaDdamage_CIs.txt

The R script plots all mismatches for every sample in one PDF. For single stranded libraries the R script is 'Plot_mismatches_ss.R'. For double stranded libraries the R script is 'Plot_mismatches_ds.R' The C->T and G->A mismatches are highlighted in red and blue respectively. Upticks in these at the 5' and 3' ends respectively suggest age-related damage. Error bars are 95% confidence intervals. If the point estimate was zero, the lower bound of the confidence interval is also set to zero.

The plots are similar to the mis-incorporation plots produced by mapDamage to aid interpretation.

The R dependencies required for this are ggplot2, gridExtra, grid.
