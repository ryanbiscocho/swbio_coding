hotspots.py

## Table of Contents
1. Introduction
2. Preparation
3. Command Line Usage
4. Downstream Use

## 1. Introduction
hotspots.py is a python script that parses a BED file that contains transposable element (TE) annotation information for a species' reference genome. These BED files are generated from our lab's custom TE annotation pipeline. The script identifies potential hotspots or coldspots of TE insertion by calculating how many TEs a given interval contains and how much this differs from the expected mean TE content of the interval (by calculating a fold-change value). Hotspots and coldspots are of interest as they may identify genes in which TEs have affected or have been selectively excluded from. As such this script can identify potentially important loci for TE insertion.

The script calculates the mean number of transposons in a given genomic interval (default: 10,000bp) and calculates the fold-changes for each interval, outputting a tab-delimited file. The tab-delimited file contains the name of the chromosome or scaffold analysed, the start position of the interval, the end position, the calculated fold-change and the species analysed.

Output file:
Scaffold	Start	End	Fold Change	Species

Another output file is also created that contains the same information but is filtered to only output fold changes that are lesser than -3x or greater than 3x, compared to the mean TE content expected.

## 2. Preparation
The script takes an input file in the BED format that contains TE annotation information. The script can be used with reference genomes consisting of chromosome-level assemblies or unplaced scaffolds. Chromosome sequences and scaffold sequences are treated equivalently in this script.

The input file should follow the following format:
Scaffold	Start	End	TE Type	Score	Strand

The file must first be sorted by scaffold number and then by start position. If not already sorted the following command can be used:

bedtools sort -i INPUT_FILE > OUTPUT_FILE

The script also requires further arguments such as the species to which the annotation file belongs to, the size or sequence length (in bp) of the genome assembly that the annotation file corresponds to and (optionally) the size of the interval that the user wants to look at for hotspots. If no interval is provided, the default interval used is 10,000bp.

## 3. Command Line Usage
The script can be used from the command line using the following syntax:
python3 hotspots.py INPUT_FILE GENUS SPECIES GENOME_SIZE INTERVAL

e.g. with the given files, the script can be called as follows:
python3 hotspots.py sortedSinonovaculaConstricta.bed Sinonovacula constricta 1220848272 100000

The script outputs two files named according to the species name provided. E.g. for the above command two output files, 'SinonovaculaConstrictaHotspots' and 'SinonovaculaConstrictaHotspots2' are created. The former is the raw output and the latter is the filtered output.

## 4. Downstream Use
The output can also be used to produce a graph to visualise hotspots or coldspots along a specific chromosome or scaffold. This involves filtering the output file for information of one chromosome/scaffold only and then plotting the graph in an R environment.

First use grep to filter only for the chromosome/scaffold you want to visualise:
grep CHROMOSOME_NAME FILE_NAME > OUTPUT_FILE
e.g. grep CM017555.1 SinonovaculaconstrictaHotspots2 > ScChrom1Hotspots

The graph can then be plotted in an R environment by importing the above output as a tab-delimited file and creating a bar chart using ggplot2 (or another data visualisation tool). Below is an example for chromosome CM017555.1 for Sinonovacula constricta.

library(ggplot2)
hotspots <- read.table(file.choose(), sep='\t', header= FALSE)
ggplot(data= hotspots, aes(x= V2, y= V4)) + geom_col(position= "dodge", aes(colour= V4)) + xlab("Chromosome Position") + ylab("Fold Change") + ggtitle("Sinonovacula constricta Chromosome CM017555.1") + scale_colour_gradient2(low="navyblue", mid="white", high= "red")

The output should be a bar chart with Chromosome Position (x axis) against Fold Change (y axis), where positive fold changes are coloured red and low ones blue.
