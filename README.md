# HiC-fast:

This repository hosts an automated data analysis pipeline, created using Snakemake, for HiC and HiChIP sequencing data. In summary, this pipeline trims and aligns HiC reads, then prepares the samples for downstream analysis by creating tag directories (used in Homer) and '.hic' files (used by Juicer tools). This pipeline uses HISAT2 to align reads, thus increasing the alignment speed compared to other pipelines. 

## Dependency:
- [Anaconda](https://conda.io/docs/user-guide/install/linux.html) 

## Installation:
Clone this repository and change into the cloned HiC-fast directory. 

To create an environment using the environment.yaml file, type the following:

`conda env create -f environment.yaml`

This will create a conda environment called hic_fast.

You will also need to install homer and download the juicer_tools jar file. It does not matter where you place Homer and the juicer_tools jar file (current version as of commit date is juicer_tools_1.19.02.jar), as you will specify the location of both your Homer install and the juicer_tools jar file in the 'config.yaml' file. 

To install homer: http://homer.ucsd.edu/homer/introduction/install.html 

The juicer_tools jar file can be [downloaded here](https://github.com/aidenlab/juicer/wiki/Download).

When analyzing your samples for TADS and loops, Homer suggests blacklisting problematic genomic regions. Instructions for extracting those regions are [provided here](http://homer.ucsd.edu/homer/interactions2/HiCTADsAndLoops.html). For convenience, the 'mm10_superdups.bed' file is provided in this repository.

## Usage note:

You must manually activate the conda environment prior to running the sh files. Type the following to activate the environment:

`conda activate hic_fast`

The reason for this requirement is a failure of the conda environment to successfully activate from within a shell script.

## Usage on an LSF cluser (Recommended usage!):

Due to the large library sizes generated by HiC (and other similar techniques, e.g. HiChIP), processing HiC samples will be most efficient on an LSF system, employing parallel processing.

Copy the config.yaml, run\_snakemake\_cluster.sh, cluster.json and Snakefile to your HiC project directory. This directory should also contain a directory called 'fastq' wherein all fastq files are placed. Make sure the project directory structure is as follows:
```
.
├── cluster.json
├── config.yaml
├── fastq
│   ├── negD1-WC-40_S2_L003_R1_001.fastq.gz
│   └── negD1-WC-40_S2_L003_R2_001.fastq.gz
├── run_snakemake_cluster.sh
└── Snakefile (required for all analyses)
```
Make the required changes to the config.yaml and cluster.json file.

Next, type `nohup sh run_snakemake_cluster.sh &` (to run in background).

## Usage on a local machine (will likely be slow!):

Copy the config.yaml, run\_snakemake.sh and Snakefile to your HiC project directory. This directory should also contain a directory called 'fastq' wherein all the fastq files are placed. Make sure the project directory structure is as follows:
```
.
├── config.yaml
├── fastq
│   ├── D1-WC_S2_L003_R1_001.fastq.gz
│   └── D1-WC_S2_L003_R2_001.fastq.gz
├── run_snakemake.sh
└── Snakefile
```
Make the required changes to the config.yaml file. Pay particular attention to the number of threads allowed for paralled processing ("max_threads"), as your local CPU may not be able to handle more than 4 threads. 

Finally, type `sh run_snakemake.sh` followed by the maximum number of CPU cores to be used by snakemake. For example, type `sh run_snakemake.sh 2` for 2 CPU cores. You can also type `nohup sh run_snakemake.sh 2 &` to run the pipeline in background.

## Alignment: Two-step or single-step alignment?

## Helpful tips for achieving better performance on LSF system
Increase "max_threads" in your config.yaml file
Adjust the walltimes and/or memory in the cluster.json file as needed to suit your samples and your particular LSF environment. Excessive walltimes and memory demands may increase the amount of time your jobs spend waiting in a queue, while insufficent walltimes and memory may cause your job to run out of memory or fail to run completely.


## Steps in Two-step alignment pipeline:

 ![ScreenShot](/dag/dag_twotier.png)

## Steps in Single-step alignment pipeline:

 ![ScreenShot](/dag/dag_singletier.png)

## Output directory structure:
```
.
├── cluster.json
├── config.yaml
├── fastq
│   ├── negD1-WC-40_S2_L003_R1_001.fastq.gz
│   └── negD1-WC-40_S2_L003_R2_001.fastq.gz
├── run_snakemake_cluster.sh
├── Snakefile
└── output
    ├── multiqc_report.html
    ├── trim_fastq
    	├── negD1-WC-40_S2_L003_R1.fq.gz_trimming_report.txt
    	└── negD1-WC-40_S2_L003_R2.fq.gz_trimming_report.txt
    ├── logs
    ├── fastqc
    	├── negD1-WC-40_S2_L003_R1_fastqc.html
    	├── negD1-WC-40_S2_L003_R1_fastqc.zip
    	├── negD1-WC-40_S2_L003_R2_fastqc.html
    	└── negD1-WC-40_S2_L003_R2_fastqc.zip
    ├── bam
    	├── negD1-WC-40_S2_L003_R1_dupsremoved.bam
    	└── negD1-WC-40_S2_L003_R2_dupsremoved.bam 
    └── homer
     	├── Merged.(tad/loop).2D.bed
    	└── negD1-WC-40_S2_L003 (Homer tag directory)
     		└── negD1-WC-40_S2_L003.hic

```


## File naming requirements

It is best practice to ensure that sequencing files are properly paired and are named with consistent fastq file endings (e.g. ".R1.fq.gz" and ".R2.fq.gz") prior to analyzing sequencing reads. However, the Snakemake pipeline is robust to mixtures of fastq file endings. It does this by detecting the most common forward read file ending (e.g. "\*.R1.fq.gz"), then renaming any files that do not conform to the most common fastq file ending.

Please note that input fastq file names that fail to conform to any of the expected fastq filename endings (for example, "\_R1\_001.fastq.gz", and ".R1.fq.gz" are examples of allowed fastq file endings) will be ignored by Snakemake (e.g. files named "Treatment.gz" or "Treatment.txt" or "Treatment" will be ignored). Only gzipped fastq files are allowed.

## Test dataset

A small example dataset, composed of one previously published mouse HiC dataset from ENCODE, is available [on Google Drive](https://drive.google.com/open?id=1ApjBYup9mOZMySgMIIDKSaTCJ3cbqeWm).

## Additional Snakemake options:

You can also customize the run\_snakemake.sh and run\_snakemake_cluster.sh scripts according to your own needs. You might wish to change the number of cores snakemake uses. Or you might want to do a dryrun. To explore additional options available in snakemake, type:

`source activate hic_fast`

followed by 

`snakemake --help`





