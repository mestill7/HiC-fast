import os
import sys
from os import listdir
from os.path import isfile, join
import re
from Bio import SeqIO
import gzip
from collections import Counter

configfile: "config.yaml"
myfastqpath = "fastq/"
sys.path.append(config["ngs_path"]) # needed to correctly find helper script
import ngs_helper.ngs_helper as ngs_helper #provides the helper functions

# HELPER FUNCTIONS
# Create function for creating rule sets

def choose_rule_all(config):
    """
    Selects the input needed for 'rule all' in snakemake pipeline

    Input Parameter: 
    config (dict): Dictionary derived from config.yaml and any additional 
    key:value pairs added during the file preperation steps. 

    Returns:
    list: List of required input files for 'rule all'. 
    """
    myout = []
    myout.append("output/multiqc_report.html")
    myout.append("output/logs/Compiled_alignment_metrics.txt")  
    myout.append(
        expand('output/bam/{sample}_R{readnum}_dupsremoved.bam', 
            sample=SAMPLES,readnum=[1,2]))
    ## Ensure that homer directory is made and hicexplorer is run
    myout.append(
        expand('output/homer/{sample}/{sample}.hic', 
            sample=SAMPLES))
    myout.append(
        expand('output/homer/{sample}_TAD/{sample}_50k_fdr_domains.bed', 
            sample=SAMPLES))
    myout.append(
        expand('output/homer/{sample}_10k_resolution.cool', 
            sample=SAMPLES))
    return(myout)
# End helper functions

# Retrieve the list of fastq files
onlyfiles, gzfiles = ngs_helper.getfilelist(myfastqpath)

# Raise exception if no fastq files present
if len(onlyfiles) == 0:
    raise NameError(
        "You do not seem to have any fastq files present to process. Exiting.")

# Raise exception if fastq files are a mixture of gzipped and non-gzipped files
if len(gzfiles) > 0 and len(gzfiles) != len(onlyfiles):
    myinput = "You have a mixture of gzipped files and non-gzipped files\n \
                Only {} of total {} files are gzipped!"
    raise NameError(print(myinput.format(len(gzfiles), len(onlyfiles))))

# Unify fastq file endings and return the final ending to be used.
if len(gzfiles) > 0:
    R1_file_ending, R2_file_ending, onlyfiles, gzfiles = \
            ngs_helper.fix_input_files(".gz", gzfiles, myfastqpath)
    suffix = ".gz"


sample_string = myfastqpath + "{sample}" + R1_file_ending + suffix
SAMPLES, = glob_wildcards(sample_string)

# Check the file pairing
# Raise exception for non-paired PE files
len_r1 = len([i for i in onlyfiles if i.endswith(R1_file_ending + suffix)])
if len_r1*2 != len(onlyfiles):
    myinput = "One or more samples do not have a read pair!\nIf using \
        paired-end samples, please ensure each sample has read 1 and \
        read 2 files\nAborting..."
    raise NameError(myinput)  # Raise exception to break workflow


# Retrieve read lengths for use as parameters
config["read_length"]= \
ngs_helper.check_readlength(suffix, gzfiles, R1_file_ending, myfastqpath)

# Calculate restriction site length and generate adaptor string
rs_len = max([len(i) for i in config["restriction_seq"].split(sep=",")])
rs_set = ["-a " + str(x) + "X " for x in config["restriction_seq"].split(sep=",")]
rs_set.extend(["-A " + str(x) + "X " for x in config["restriction_seq"].split(sep=",")])

# Generate input rule for Snakemake
rule all:
    input:
        choose_rule_all(config)

#Rules 
rule trim_fastq_galore:
    input:
        pair1 = (join(myfastqpath, "{sample}") + expand(
            "{ending}.gz", ending=R1_file_ending)[0]),
        pair2 = (join(myfastqpath, "{sample}") + expand(
            "{ending}.gz", ending=R2_file_ending)[0])
    output:
        trimmed_pair1 = temp(
            "output/trim_fastq/{sample}_R1_val_1.fq.gz"),
        trimmed_pair2 = temp(
            "output/trim_fastq/{sample}_R2_val_2.fq.gz"),
        original_r1 = temp("output/temp_dir/{sample}_R1.fq.gz"),
        original_r2 = temp("output/temp_dir/{sample}_R2.fq.gz")
    log:
        "output/logs/{sample}.trim_galore.log"
    params:
        read_length = config["read_length"],
        max_threads = config["max_threads"],
        rs_set = rs_set,
        read_length_fortrim = config["read_length_adaptormin"]
    run :
        shell("mkdir -p output/temp_dir")
        shell("cp {input.pair1} \
            output/temp_dir/{wildcards.sample}_R1.fq.gz")
        shell("cp {input.pair2} \
            output/temp_dir/{wildcards.sample}_R2.fq.gz")
        shell("trim_galore -j {params.max_threads} --paired -q 20 \
            -o output/trim_fastq \
        output/temp_dir/{wildcards.sample}_R1.fq.gz \
        output/temp_dir/{wildcards.sample}_R2.fq.gz")


#Use this condition for running the restriction site step.
rule trim_fastq_cutadapt:
    input:
        pair1 = "output/trim_fastq/{sample}_R1_val_1.fq.gz",
        pair2 = "output/trim_fastq/{sample}_R2_val_2.fq.gz"
    output:
        trimmed_pair1 = temp("output/trim_fastq/{sample}_R1_trimmed.fastq.gz"),
        trimmed_pair2 = temp("output/trim_fastq/{sample}_R2_trimmed.fastq.gz")
    log:
        "output/logs/{sample}.trim_cutadapt.log"
    params:
        read_length_fortrim = config["read_length_adaptormin"],
        trim_min = config["read_length_trim_min"],
        read_length = config["read_length"],
        max_threads = config["max_threads"],
        rs_set = rs_set
    run:
        if (config["read_length"]-config["read_length_trim_min"])<rs_len :
            # restriction sites entirely encompassed by the trim or using a two-stage system.
            shell("cutadapt -j {params.max_threads} {params.rs_set} --minimum-length {params.read_length_fortrim} \
                -o {output.trimmed_pair1} -p {output.trimmed_pair2} {input.pair1} {input.pair2}") 
        else: 
            shell("trimmomatic SE -threads {params.max_threads} {input.pair1} {output.trimmed_pair1} CROP:{params.trim_min}")
            shell("trimmomatic SE -threads {params.max_threads} {input.pair2} {output.trimmed_pair2} CROP:{params.trim_min}")

rule fastqc:
    input:
        pair1 = "output/trim_fastq/{sample}_R1_trimmed.fastq.gz",
        pair2 = "output/trim_fastq/{sample}_R2_trimmed.fastq.gz"  
    output:  
        fastqc_zipfile1 = "output/fastqc/{sample}_R1_trimmed_fastqc.zip",
        fastqc_zipfile2 = "output/fastqc/{sample}_R2_trimmed_fastqc.zip"
    params:
        max_threads = config["max_threads"]
    run:
        shell("fastqc -t {params.max_threads} -o output/fastqc -d output/temp_dir {input.pair1}")
        shell("fastqc -t {params.max_threads} -o output/fastqc -d output/temp_dir {input.pair2}")

rule alignment_R1:
    input:
        trimmed_reads = "output/trim_fastq/{sample}_R1_trimmed.fastq.gz"    
    output:
        bam = temp("output/bam/{sample}_R1.bam")
    log:
        "output/logs/{sample}.align_R1.log"
    params:
        max_threads = config["max_threads"],
        index = config["index"]
    run:
        shell("hisat2 -p {params.max_threads} -k 1 --no-spliced-alignment --max-seeds 10 --no-unal \
            -x {params.index} -U {input.trimmed_reads} -S output/bam/{wildcards.sample}_R1.sam")
        shell("samtools view -b output/bam/{wildcards.sample}_R1.sam > {output.bam}")
        shell("rm output/bam/{wildcards.sample}_R1.sam")

rule alignment_R2:
    input:
        trimmed_reads = "output/trim_fastq/{sample}_R2_trimmed.fastq.gz"    
    output:
        bam = temp("output/bam/{sample}_R2.bam")
    log:
        "output/logs/{sample}.align_R2.log"
    params:
        max_threads = config["max_threads"],
        index = config["index"]
    run:
        shell("hisat2 -p {params.max_threads} -k 1 --no-spliced-alignment --max-seeds 10 --no-unal \
            -x {params.index} -U {input.trimmed_reads} -S output/bam/{wildcards.sample}_R2.sam")
        shell("samtools view -b output/bam/{wildcards.sample}_R2.sam > {output.bam}")
        shell("rm output/bam/{wildcards.sample}_R2.sam")

rule rmdup:
    input: 
        pair1 = "output/bam/{sample}_R1.bam",
        pair2 = "output/bam/{sample}_R2.bam"
    output:
        r1 = "output/bam/{sample}_R1_dupsremoved.bam",
        r2 = "output/bam/{sample}_R2_dupsremoved.bam"
    log:
        "output/logs/{sample}.rmdup.log"
    run:
        shell("picard SortSam I={input.pair1} O=output/bam/{wildcards.sample}_R1_sorted.bam SORT_ORDER=coordinate")
        shell("picard SortSam I={input.pair2} O=output/bam/{wildcards.sample}_R2_sorted.bam SORT_ORDER=coordinate")
        shell("picard MarkDuplicates REMOVE_DUPLICATES=TRUE I=output/bam/{wildcards.sample}_R1_sorted.bam \
            O=output/bam/{wildcards.sample}_R1_dupsremoved.bam \
            M=output/logs/{wildcards.sample}_R1_dup_metrics.txt")
        shell("picard MarkDuplicates REMOVE_DUPLICATES=TRUE I=output/bam/{wildcards.sample}_R2_sorted.bam \
            O=output/bam/{wildcards.sample}_R2_dupsremoved.bam \
            M=output/logs/{wildcards.sample}_R2_dup_metrics.txt")
        shell("rm output/bam/{wildcards.sample}_R1_sorted.bam")
        shell("rm output/bam/{wildcards.sample}_R2_sorted.bam")
        shell("samtools index output/bam/{wildcards.sample}_R1_dupsremoved.bam")
        shell("samtools index output/bam/{wildcards.sample}_R2_dupsremoved.bam")
        # if config["keep_fastq"]=="FALSE":
        #     shell("rm output/trim_fastq/{wildcards.sample}_R1_trimmed.fastq.gz \
        #         output/trim_fastq/{wildcards.sample}_R2_trimmed.fastq.gz")


rule collect_stats:
    input:
        stage1_r1_fastq = "output/trim_fastq/{sample}_R1_trimmed.fastq.gz",
        stage1_r1_bam = "output/bam/{sample}_R1.bam",
        stage1_r1_dupsremoved_bam = "output/bam/{sample}_R1_dupsremoved.bam",
        stage1_r2_fastq = "output/trim_fastq/{sample}_R2_trimmed.fastq.gz",
        stage1_r2_bam = "output/bam/{sample}_R2.bam",
        stage1_r2_dupsremoved_bam = "output/bam/{sample}_R2_dupsremoved.bam",
        original_r1 = "output/temp_dir/{sample}_R1.fq.gz",
        original_r2 = "output/temp_dir/{sample}_R2.fq.gz"
    output:
        r1 = "output/logs/{sample}_R1_metrics.txt",
        r2 = "output/logs/{sample}_R2_metrics.txt"
    run:
        shell("unpigz -c {input.original_r1} | wc -l | awk '{{print $1/4}}' > output/logs/{wildcards.sample}_R1_metrics.txt")
        shell("unpigz -c {input.stage1_r1_fastq} | wc -l | awk '{{print $1/4}}' >> output/logs/{wildcards.sample}_R1_metrics.txt")
        shell("samtools view -c {input.stage1_r1_bam} >> output/logs/{wildcards.sample}_R1_metrics.txt")
        shell("samtools view -c {input.stage1_r1_dupsremoved_bam} >> output/logs/{wildcards.sample}_R1_metrics.txt")
        shell("unpigz -c {input.original_r2} | wc -l | awk '{{print $1/4}}' > output/logs/{wildcards.sample}_R2_metrics.txt")
        shell("unpigz -c {input.stage1_r2_fastq} | wc -l | awk '{{print $1/4}}' >> output/logs/{wildcards.sample}_R2_metrics.txt")
        shell("samtools view -c {input.stage1_r2_bam} >> output/logs/{wildcards.sample}_R2_metrics.txt")
        shell("samtools view -c {input.stage1_r2_dupsremoved_bam} >> output/logs/{wildcards.sample}_R2_metrics.txt")

rule compile_stats:
    input:
        expand("output/logs/{sample}_R{readnum}_metrics.txt", sample=SAMPLES,readnum=[1,2])
    output:
        "output/logs/Compiled_alignment_metrics.txt"
    run:
        import pandas as pd
        import re
        dict_of_counts = {}
        for file in input:
            sample=re.sub(r'_metrics.txt', '', re.sub(r'.*/', '', file)) #equivalent of 'basename' command
            sample
            with open(file, "r") as infile:
                mytext=infile.read()
                mytext=list(filter(None, mytext.split("\n"))) #remove empty strings
                dict_of_counts[sample] = mytext

        dataframe = pd.DataFrame(dict_of_counts)
        dataframe.index = ["starting_reads","trimmed_reads","aligned_reads","dedup_reads"]
        dataframe.to_csv(output[0], sep = '\t')



# Introduce homer 
rule create_tagdir:
    input:
        r1 = "output/bam/{sample}_R1.bam",
        r2 = "output/bam/{sample}_R2.bam"
    output:
        required = "output/homer/{sample}/tagLengthDistribution.txt",
        timecheck = temp("output/homer/{sample}/Test_tagdir_completion.txt"),
        r1 = temp("output/bam/{sample}_R1.chr.bam"),
        r2 = temp("output/bam/{sample}_R2.chr.bam")
    log:
        "output/logs/{sample}.tagdir.log"
    params:
        homer_bin = config["homer_bin"],
        homer_genome = config["homer_genome"],
        main_chr = config["main_chr"]
    run:
        shell("picard SortSam I={input.r1} O=output/bam/{wildcards.sample}_R1_sorted.bam SORT_ORDER=coordinate")
        shell("picard SortSam I={input.r2} O=output/bam/{wildcards.sample}_R2_sorted.bam SORT_ORDER=coordinate")
        shell("samtools index output/bam/{wildcards.sample}_R1_sorted.bam")
        shell("samtools index output/bam/{wildcards.sample}_R2_sorted.bam")
        shell("samtools view -b output/bam/{wildcards.sample}_R1_sorted.bam {params.main_chr} > {output.r1}")
        shell("samtools view -b output/bam/{wildcards.sample}_R2_sorted.bam {params.main_chr} > {output.r2}")
        restriction = config["restriction_seq"].split(sep=",")[0]
        shell("samtools index {output.r1}")
        shell("samtools index {output.r2}")
        shell("{params.homer_bin}/makeTagDirectory output/homer/{wildcards.sample} {output.r1},{output.r2} -tbp 1 -genome {params.homer_genome}")
        shell("{params.homer_bin}/makeTagDirectory output/homer/{wildcards.sample} -update -genome {params.homer_genome}\
            -removePEbg -restrictionSite {restriction} -both -removeSelfLigation -removeSpikes 10000 5 2> output/logs/{wildcards.sample}_tag_log.txt")
        shell("touch {output.timecheck}")

rule create_hic:
    input:
        required = "output/homer/{sample}/tagLengthDistribution.txt",
        timecheck = "output/homer/{sample}/Test_tagdir_completion.txt"
    output:
        "output/homer/{sample}/{sample}.hic"
    log:
        "output/logs/{sample}.createhic.log"
    params:
        homer_bin = config["homer_bin"],
        juicer_jar = config["juicer_jar"],
        homer_genome = config["homer_genome"]
    run:
        shell("PATH={params.homer_bin}:$PATH; \
            tagDir2hicFile.pl output/homer/{wildcards.sample}/ -juicer auto -genome {params.homer_genome} -p 1 \
        -juicerExe \"java -jar {params.juicer_jar}\"")

## Rule for pairing reads
rule pair_reads:
    input: 
        r1 = "output/bam/{sample}_R1.bam",
        r2 = "output/bam/{sample}_R2.bam"
    output:
        rname="output/bam/{sample}_common_readnames",
        r1="output/bam/{sample}_R1_paired.bam",
        r2="output/bam/{sample}_R2_paired.bam"
    log:
        "output/logs/{sample}.pair_reads.log"
    run:
        shell("picard SortSam I=output/bam/{wildcards.sample}_R1.bam O=output/bam/{wildcards.sample}_R1_sorted.bam SORT_ORDER=queryname")
        shell("picard SortSam I=output/bam/{wildcards.sample}_R2.bam O=output/bam/{wildcards.sample}_R2_sorted.bam SORT_ORDER=queryname")
        shell("samtools view output/bam/{wildcards.sample}_R1_sorted.bam  | awk '{{print $1}}' > output/bam/{wildcards.sample}_R1_readnames")
        shell("samtools view output/bam/{wildcards.sample}_R2_sorted.bam  | awk '{{print $1}}' > output/bam/{wildcards.sample}_R2_readnames")
        shell("comm --nocheck-order -12 output/bam/{wildcards.sample}_R1_readnames output/bam/{wildcards.sample}_R2_readnames > output/bam/{wildcards.sample}_common_readnames")
        shell("samtools view -b -N output/bam/{wildcards.sample}_common_readnames output/bam/{wildcards.sample}_R1_sorted.bam > output/bam/{wildcards.sample}_R1_paired.bam")
        shell("samtools view -b -N output/bam/{wildcards.sample}_common_readnames output/bam/{wildcards.sample}_R2_sorted.bam > output/bam/{wildcards.sample}_R2_paired.bam")
        shell("samtools index output/bam/{wildcards.sample}_R1_paired.bam")
        shell("samtools index output/bam/{wildcards.sample}_R2_paired.bam")

rule generate_rest_sites:
    params:
        homer_bin = config["genome_fasta"],
        rest = config["restriction_seq"]
    output: 
        "output/homer/rest_site_positions.bed"
    run:
        shell("hicFindRestSite --fasta {params.homer_bin} --searchPattern {params.rest} -o {output}")

rule build_matrix: 
    input: 
        r1 = "output/bam/{sample}_R1_paired.bam",
        r2 = "output/bam/{sample}_R2_paired.bam",
        rest = "output/homer/rest_site_positions.bed"
    output:
        h5name="output/homer/{sample}_5k_resolution.h5",
        coolname="output/homer/{sample}_10k_resolution.cool"
    params:
        chrom_sizes = config["chrom_sizes"]
    log:
        "output/logs/{sample}.build_matrix.log"
    run:
        shell("hicBuildMatrix -s {input.r1} {input.r2} -o output/homer/{wildcards.sample}_10k_resolution.cool \
          --binSize 10000 --QCfolder output/homer/{wildcards.sample}_QC --chromosomeSizes {params.chrom_sizes} \
          -rs {input.rest} -seq GATC --danglingSequence GATC")
        shell("hicBuildMatrix -s {input.r1} {input.r2} -o output/homer/{wildcards.sample}_5k_resolution.h5 \
          --binSize 5000 --QCfolder output/homer/{wildcards.sample}_QC --chromosomeSizes {params.chrom_sizes} \
          -rs {input.rest} -seq GATC --danglingSequence GATC")

rule merge_matrix: 
    input: 
        "output/homer/{sample}_5k_resolution.h5"
    output: 
        k50="output/homer/{sample}_50k_resolution.h5",
        k25="output/homer/{sample}_25k_resolution.h5",
        k10="output/homer/{sample}_10k_resolution.h5"
    log:
        "output/logs/{sample}.merge_matrix.log"
    run:
        shell("hicMergeMatrixBins -m {input}  -o output/homer/{wildcards.sample}_50k_resolution.h5 -nb 10")
        shell("hicMergeMatrixBins -m {input}  -o output/homer/{wildcards.sample}_25k_resolution.h5 -nb 5")
        shell("hicMergeMatrixBins -m {input}  -o output/homer/{wildcards.sample}_10k_resolution.h5 -nb 2")

rule normalize_matrix: 
    input: 
        k50="output/homer/{sample}_50k_resolution.h5",
        k25="output/homer/{sample}_25k_resolution.h5",
        k10="output/homer/{sample}_10k_resolution.h5",
        k5="output/homer/{sample}_5k_resolution.h5"
    output: 
        k50="output/homer/{sample}_50k_resolution_norm.h5",
        k25="output/homer/{sample}_25k_resolution_norm.h5",
        k10="output/homer/{sample}_10k_resolution_norm.h5",
        k5="output/homer/{sample}_5k_resolution_norm.h5"  
    log:
        "output/logs/{sample}.norm_matrix.log"
    run:
        shell("hicNormalize -m {input.k5} {input.k10} {input.k25} {input.k50} \
         --normalize norm_range -o {output.k5} {output.k10} {output.k25} {output.k50}")

rule correct_matrix: 
    input: 
        k50="output/homer/{sample}_50k_resolution_norm.h5",
        k25="output/homer/{sample}_25k_resolution_norm.h5",
        k10="output/homer/{sample}_10k_resolution_norm.h5",
        k5="output/homer/{sample}_5k_resolution_norm.h5"
    output:  
        k50="output/homer/{sample}_50k_resolution_corrected.h5",
        k25="output/homer/{sample}_25k_resolution_corrected.h5",
        k10="output/homer/{sample}_10k_resolution_corrected.h5",
        k5="output/homer/{sample}_5k_resolution_corrected.h5"
    run:
        shell("hicCorrectMatrix correct -m {input.k5} --filterThreshold -1.5 5 -o {output.k5}")
        shell("hicCorrectMatrix correct -m {input.k10} --filterThreshold -1.5 5 -o {output.k10}")
        shell("hicCorrectMatrix correct -m {input.k25} --filterThreshold -1.5 5 -o {output.k25}")
        shell("hicCorrectMatrix correct -m {input.k50} --filterThreshold -1.5 5 -o {output.k50}")

rule calculate_TADs: 
    input: 
        k50="output/homer/{sample}_50k_resolution_corrected.h5",
        k25="output/homer/{sample}_25k_resolution_corrected.h5",
        k10="output/homer/{sample}_10k_resolution_corrected.h5"
    output: 
        k50="output/homer/{sample}_TAD/{sample}_50k_fdr_domains.bed",
        k25="output/homer/{sample}_TAD/{sample}_25k_fdr_domains.bed",
        k10="output/homer/{sample}_TAD/{sample}_10k_fdr_domains.bed"
    log:
        "output/logs/{sample}.calc_TADs.log"
    run:
        shell("mkdir -p output/homer/{wildcards.sample}_TAD")
        shell("hicFindTADs -m {input.k10} --outPrefix output/homer/{wildcards.sample}_TAD/{wildcards.sample}_10k_fdr \
            --numberOfProcessors 16 --correctForMultipleTesting fdr")
        shell("hicFindTADs -m {input.k25} --outPrefix output/homer/{wildcards.sample}_TAD/{wildcards.sample}_25k_fdr \
            --numberOfProcessors 16 --correctForMultipleTesting fdr")      
        shell("hicFindTADs -m {input.k50} --outPrefix output/homer/{wildcards.sample}_TAD/{wildcards.sample}_50k_fdr \
            --numberOfProcessors 16 --correctForMultipleTesting fdr")

rule calculate_loops: 
    input: 
        k10="output/homer/{sample}_10k_resolution.cool"
    output: 
        k10="output/homer/{sample}_loop/{sample}_10k_loops.bedgraph"
    log:
        "output/logs/{sample}.calc_TADs.log"
    run:
        shell("mkdir output/homer/{wildcards.sample}_loop")
        shell("hicDetectLoops -m {input.k10} -o {output.k10} \
            --maxLoopDistance 2000000 --windowSize 10 --peakWidth 6 \
            --pValuePreselection 0.05 --pValue 0.05")

# Create multiqc report
rule run_multiqc:
    input:
        tad = expand('output/homer/{sample}_TAD/{sample}_50k_fdr_domains.bed', 
            sample=SAMPLES),
        cool = expand('output/homer/{sample}_10k_resolution.cool', 
            sample=SAMPLES),
        k10 = expand('output/homer/{sample}_loop/{sample}_10k_loops.bedgraph', 
            sample=SAMPLES)
    output:
        multiqc_report = "output/multiqc_report.html"
    params:
        multiqc_config = config["multiqc_yaml"]
    shell:
        "multiqc . -f --config {params.multiqc_config}"


