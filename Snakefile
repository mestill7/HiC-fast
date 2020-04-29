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
    if config["two_tier"] == "FALSE":     
        myout.append(
            expand('output/bam/{sample}_R{readnum}_dupsremoved.bam', 
                sample=SAMPLES,readnum=[1,2]))
    else :
        myout.append(
            expand('output/bam/{sample}_R{readnum}_merged.bam',sample=SAMPLES,readnum=[1,2]))
    if config["to_homer"] == "TRUE" :
        myout.append(
            expand('output/homer/{sample}/{sample}.tad.2D.bed', 
                sample=SAMPLES))
        myout.append(
            expand('output/homer/{sample}/{sample}.hic', 
                sample=SAMPLES))
        if config["to_merge"] == "TRUE": #create merged loops/tads and use all samples in scoring
            myout.append('output/homer/Merged.tad.scores.txt')
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

#Rules for two-tier and single-tier
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
        # if ({params.read_length} - {params.trim_min} < rs_len) or config["two_tier"] == "TRUE"  : 
        if config["two_tier"] == "TRUE" or (config["read_length"]-config["read_length_trim_min"])<rs_len :
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

# Create multiqc report
rule run_multiqc:
    input:
        qc = expand("output/fastqc/{sample}_R{readnum}_trimmed_fastqc.zip", sample=SAMPLES,readnum=[1,2])
    output:
        multiqc_report = "output/multiqc_report.html"
    params:
        multiqc_config = config["multiqc_yaml"]
    shell:
        "multiqc . -f --config {params.multiqc_config}"

if config["two_tier"] == "FALSE" :
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


#Add additional rules for Homer, two-tier approach
##Two-tier approach
if config["two_tier"] == "TRUE" :
    rule generate_unmaplist:
        input:
            pair1 = "output/trim_fastq/{sample}_R1_trimmed.fastq.gz",
            pair2 = "output/trim_fastq/{sample}_R2_trimmed.fastq.gz",
            bam1 = "output/bam/{sample}_R1.bam",
            bam2 = "output/bam/{sample}_R2.bam"
        output:  
            r1 = temp("output/bam/{sample}_R1_unal_reads.txt"),
            r2 = temp("output/bam/{sample}_R2_unal_reads.txt")
        params:
            max_threads = config["max_threads"]
        run:
            #Get aligned read names from bam file
            shell("samtools view output/bam/{wildcards.sample}_R1.bam | cut -f1 > output/bam/{wildcards.sample}_R1_aligned_reads.txt")
            shell("samtools view output/bam/{wildcards.sample}_R2.bam | cut -f1 > output/bam/{wildcards.sample}_R2_aligned_reads.txt")
            #Get all read names from the trimmed fastq files
            shell("zcat {input.pair1} | sed -n '1~4p' | sed 's/ .*//' | sed 's/@//g' > output/bam/{wildcards.sample}_R1_all_reads.txt")
            shell("zcat {input.pair2} | sed -n '1~4p' | sed 's/ .*//' | sed 's/@//g' > output/bam/{wildcards.sample}_R2_all_reads.txt")

            #Sort read names
            shell("sort -T output/temp_dir --parallel={params.max_threads} \
            --output=output/bam/{wildcards.sample}_R1_aligned_reads_sorted.txt  \
            output/bam/{wildcards.sample}_R1_aligned_reads.txt")
            shell("sort -T output/temp_dir --parallel={params.max_threads} \
            --output=output/bam/{wildcards.sample}_R1_all_reads_sorted.txt  \
            output/bam/{wildcards.sample}_R1_all_reads.txt")
            shell("sort -T output/temp_dir --parallel={params.max_threads} \
            --output=output/bam/{wildcards.sample}_R2_aligned_reads_sorted.txt \
             output/bam/{wildcards.sample}_R2_aligned_reads.txt")
            shell("sort -T output/temp_dir --parallel={params.max_threads} \
            --output=output/bam/{wildcards.sample}_R2_all_reads_sorted.txt  \
            output/bam/{wildcards.sample}_R2_all_reads.txt")
            # Compare lists
            shell("comm -23 output/bam/{wildcards.sample}_R1_all_reads_sorted.txt \
            output/bam/{wildcards.sample}_R1_aligned_reads_sorted.txt > \
            output/bam/{wildcards.sample}_R1_unal_reads.txt")
            shell("comm -23 output/bam/{wildcards.sample}_R2_all_reads_sorted.txt \
            output/bam/{wildcards.sample}_R2_aligned_reads_sorted.txt > \
            output/bam/{wildcards.sample}_R2_unal_reads.txt")
            shell("rm output/bam/{wildcards.sample}_R1_aligned_reads*.txt \
                output/bam/{wildcards.sample}_R1_all_reads*.txt")

    rule split_unmap:
        input:
            r1 = "output/bam/{sample}_R1_unal_reads.txt",
            r2 = "output/bam/{sample}_R2_unal_reads.txt"
        output:  
            r1 = temp("output/temp_dir/{sample}_R1_splitlist_part1"),
            r2 = temp("output/temp_dir/{sample}_R2_splitlist_part1")
        run:  
            shell("split -l 50000 output/bam/{wildcards.sample}_R1_unal_reads.txt output/temp_dir/{wildcards.sample}_R1_unalsplit.")
            shell("split -l 50000 output/bam/{wildcards.sample}_R2_unal_reads.txt output/temp_dir/{wildcards.sample}_R2_unalsplit.")
            shell("ls output/temp_dir/{wildcards.sample}_R1_unalsplit.* > output/temp_dir/{wildcards.sample}_R1_splitlist_part1")
            shell("ls output/temp_dir/{wildcards.sample}_R2_unalsplit.* > output/temp_dir/{wildcards.sample}_R2_splitlist_part1")

    rule rename_unmap:
        input:
            r1 = "output/temp_dir/{sample}_R1_splitlist_part1",
            r2 = "output/temp_dir/{sample}_R2_splitlist_part1"  
        output:  
            r1 = temp("output/temp_dir/{sample}_R1_splitlist"),
            r2 = temp("output/temp_dir/{sample}_R2_splitlist") 
        run:
            import os
            for file in input:
                with open(file, "r") as infile:
                    for lines in infile:
                        lines=lines.strip()
                        os.rename(lines,lines+".txt")
            shell("ls output/temp_dir/{wildcards.sample}_R1_unalsplit.*.txt > output/temp_dir/{wildcards.sample}_R1_splitlist")
            shell("ls output/temp_dir/{wildcards.sample}_R2_unalsplit.*.txt > output/temp_dir/{wildcards.sample}_R2_splitlist")

    rule filter_unmap:
        input:
            r1 = "output/temp_dir/{sample}_R1_splitlist",
            r2 = "output/temp_dir/{sample}_R2_splitlist"
        output:  
            r1 = temp("output/temp_dir/{sample}_R1_unmapped.fastq.gz"),
            r2 = temp("output/temp_dir/{sample}_R2_unmapped.fastq.gz")
        run:
            import re 
            with open(input.r1, "r") as infile:
                for lines in infile:
                    lines=lines.strip()
                    mybase=re.sub(r'.*/', '', lines) #equivalent of 'basename' command
                    splitsuffix=re.sub(r'.txt', '', re.sub(r'.*_unalsplit.', '', mybase))
                    shell("filterbyname.sh in=output/trim_fastq/{wildcards.sample}_R1_trimmed.fastq.gz \
                    out=output/temp_dir/{wildcards.sample}_R1_{splitsuffix}.fastq.gz \
                    names=output/temp_dir/{wildcards.sample}_R1_unalsplit.{splitsuffix}.txt include=true ow=t")

            with open(input.r2, "r") as infile:
                for lines in infile:
                    lines=lines.strip()
                    mybase=re.sub(r'.*/', '', lines) #equivalent of 'basename' command
                    splitsuffix=re.sub(r'.txt', '', re.sub(r'.*_unalsplit.', '', mybase))
                    shell("filterbyname.sh in=output/trim_fastq/{wildcards.sample}_R2_trimmed.fastq.gz \
                    out=output/temp_dir/{wildcards.sample}_R2_{splitsuffix}.fastq.gz \
                    names=output/temp_dir/{wildcards.sample}_R2_unalsplit.{splitsuffix}.txt include=true ow=t")

            shell("cat output/temp_dir/{wildcards.sample}_R1_*.fastq.gz > output/temp_dir/{wildcards.sample}_R1_unmapped.fastq.gz")
            shell("cat output/temp_dir/{wildcards.sample}_R2_*.fastq.gz > output/temp_dir/{wildcards.sample}_R2_unmapped.fastq.gz")

            with open(input.r1, "r") as infile:
                for lines in infile:
                    lines=lines.strip()
                    mybase=re.sub(r'.*/', '', lines) #equivalent of 'basename' command
                    splitsuffix=re.sub(r'.txt', '', re.sub(r'.*_unalsplit.', '', mybase))
                    shell("rm output/temp_dir/{wildcards.sample}_R1_{splitsuffix}.fastq.gz")
                    shell("rm output/temp_dir/{wildcards.sample}_R1_unalsplit.{splitsuffix}.txt")
            with open(input.r2, "r") as infile:
                for lines in infile:
                    lines=lines.strip()
                    mybase=re.sub(r'.*/', '', lines) #equivalent of 'basename' command
                    splitsuffix=re.sub(r'.txt', '', re.sub(r'.*_unalsplit.', '', mybase))
                    shell("rm output/temp_dir/{wildcards.sample}_R2_{splitsuffix}.fastq.gz")
                    shell("rm output/temp_dir/{wildcards.sample}_R2_unalsplit.{splitsuffix}.txt")                    

    rule trim_stage2:
        input:
            pair1 = "output/temp_dir/{sample}_R1_unmapped.fastq.gz",
            pair2 = "output/temp_dir/{sample}_R2_unmapped.fastq.gz"
        output:  
            trimmed_pair1 = temp("output/trim_fastq/{sample}_R1_unmapped_trimmed.fastq.gz"),
            trimmed_pair2 = temp("output/trim_fastq/{sample}_R2_unmapped_trimmed.fastq.gz")
        params:
            max_threads = config["max_threads"],
            trim_min = config["read_length_trim_min"]
        run:
            shell("trimmomatic SE -threads {params.max_threads}  \
                {input.pair1} {output.trimmed_pair1} CROP:{params.trim_min}")
            shell("trimmomatic SE -threads {params.max_threads}  \
                {input.pair2} {output.trimmed_pair2} CROP:{params.trim_min}")        

    rule alignment_R1_stage2:
        input:
            trimmed_reads = "output/trim_fastq/{sample}_R1_unmapped_trimmed.fastq.gz"    
        output:
            bam = temp("output/bam/{sample}_R1_stage2.bam")
        log:
            "output/logs/{sample}.align_stage2_R1.log"
        params:
            max_threads = config["max_threads"],
            index = config["index"]        
        run:
            shell("hisat2 -p {params.max_threads} -k 1 --no-spliced-alignment --max-seeds 10 --no-unal -x {params.index} \
            -U {input.trimmed_reads} -S output/bam/{wildcards.sample}_R1.sam")
            shell("samtools view -b output/bam/{wildcards.sample}_R1.sam > {output.bam}")
            shell("rm output/bam/{wildcards.sample}_R1.sam")

    rule alignment_R2_stage2:
        input:
            trimmed_reads = "output/trim_fastq/{sample}_R2_unmapped_trimmed.fastq.gz"    
        output:
            bam = temp("output/bam/{sample}_R2_stage2.bam")
        log:
            "output/logs/{sample}.align_stage2_R2.log"
        params:
            max_threads = config["max_threads"],
            index = config["index"]        
        run:
            shell("hisat2 -p {params.max_threads} -k 1 --no-spliced-alignment --max-seeds 10 --no-unal -x {params.index} \
            -U {input.trimmed_reads} -S output/bam/{wildcards.sample}_R2.sam")
            shell("samtools view -b output/bam/{wildcards.sample}_R2.sam > {output.bam}")
            shell("rm output/bam/{wildcards.sample}_R2.sam")

    rule rmdup_stage2:
        input: 
            pair1 = "output/bam/{sample}_R1_stage2.bam",
            pair2 = "output/bam/{sample}_R2_stage2.bam"
        output:
            r1 = "output/bam/{sample}_R1_stage2_dupsremoved.bam",
            r2 = "output/bam/{sample}_R2_stage2_dupsremoved.bam"
        log:
            "output/logs/{sample}.rmdup.log"
        run:
            shell("picard SortSam I={input.pair1} O=output/bam/{wildcards.sample}_R1_sorted.bam SORT_ORDER=coordinate")
            shell("picard SortSam I={input.pair2} O=output/bam/{wildcards.sample}_R2_sorted.bam SORT_ORDER=coordinate")
            shell("picard MarkDuplicates REMOVE_DUPLICATES=TRUE I=output/bam/{wildcards.sample}_R1_sorted.bam \
                O=output/bam/{wildcards.sample}_R1_stage2_dupsremoved.bam \
                M=output/logs/{wildcards.sample}_R1_stage2_dup_metrics.txt")
            shell("picard MarkDuplicates REMOVE_DUPLICATES=TRUE I=output/bam/{wildcards.sample}_R2_sorted.bam \
                O=output/bam/{wildcards.sample}_R2_stage2_dupsremoved.bam \
                M=output/logs/{wildcards.sample}_R2_stage2_dup_metrics.txt")
            shell("rm output/bam/{wildcards.sample}_R1_sorted.bam")
            shell("rm output/bam/{wildcards.sample}_R2_sorted.bam")

    rule collect_stats:
        input:
            stage1_r1_fastq = "output/trim_fastq/{sample}_R1_trimmed.fastq.gz",
            stage1_r1_bam = "output/bam/{sample}_R1.bam",
            stage1_r1_dupsremoved_bam = "output/bam/{sample}_R1_dupsremoved.bam",
            stage2_r1_fastq = "output/trim_fastq/{sample}_R1_unmapped_trimmed.fastq.gz",
            stage2_r1_bam = "output/bam/{sample}_R1_stage2.bam",
            stage2_r1_dupsremoved_bam = "output/bam/{sample}_R1_stage2_dupsremoved.bam",
            stage1_r2_fastq = "output/trim_fastq/{sample}_R2_trimmed.fastq.gz",
            stage1_r2_bam = "output/bam/{sample}_R2.bam",
            stage1_r2_dupsremoved_bam = "output/bam/{sample}_R2_dupsremoved.bam",
            stage2_r2_fastq = "output/trim_fastq/{sample}_R2_unmapped_trimmed.fastq.gz",
            stage2_r2_bam = "output/bam/{sample}_R2_stage2.bam",
            stage2_r2_dupsremoved_bam = "output/bam/{sample}_R2_stage2_dupsremoved.bam",
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
            shell("unpigz -c {input.stage2_r1_fastq} | wc -l | awk '{{print $1/4}}' >> output/logs/{wildcards.sample}_R1_metrics.txt")
            shell("samtools view -c {input.stage2_r1_bam} >> output/logs/{wildcards.sample}_R1_metrics.txt")
            shell("samtools view -c {input.stage2_r1_dupsremoved_bam} >> output/logs/{wildcards.sample}_R1_metrics.txt")

            shell("unpigz -c {input.original_r2} | wc -l | awk '{{print $1/4}}' > output/logs/{wildcards.sample}_R2_metrics.txt")
            shell("unpigz -c {input.stage1_r2_fastq} | wc -l | awk '{{print $1/4}}' >> output/logs/{wildcards.sample}_R2_metrics.txt")
            shell("samtools view -c {input.stage1_r2_bam} >> output/logs/{wildcards.sample}_R2_metrics.txt")
            shell("samtools view -c {input.stage1_r2_dupsremoved_bam} >> output/logs/{wildcards.sample}_R2_metrics.txt")
            shell("unpigz -c {input.stage2_r2_fastq} | wc -l | awk '{{print $1/4}}' >> output/logs/{wildcards.sample}_R2_metrics.txt")
            shell("samtools view -c {input.stage2_r2_bam} >> output/logs/{wildcards.sample}_R2_metrics.txt")
            shell("samtools view -c {input.stage2_r2_dupsremoved_bam} >> output/logs/{wildcards.sample}_R2_metrics.txt")

    rule merge_bam:
        input:
            stage1_r1 = "output/bam/{sample}_R1_dupsremoved.bam",
            stage2_r1 = "output/bam/{sample}_R1_stage2_dupsremoved.bam",
            stage1_r2 = "output/bam/{sample}_R2_dupsremoved.bam",
            stage2_r2 = "output/bam/{sample}_R2_stage2_dupsremoved.bam"
        output:
            r1 = "output/bam/{sample}_R1_merged.bam",
            r2 = "output/bam/{sample}_R2_merged.bam"
        params:
            max_threads = config["max_threads"],
            index = config["index"]
        run:
            shell("samtools merge -@ {params.max_threads} {output.r1} {input.stage1_r1} {input.stage2_r1}"),
            shell("samtools index {output.r1}"),
            shell("samtools merge -@ {params.max_threads} {output.r2} {input.stage1_r2} {input.stage2_r2}"),
            shell("samtools index {output.r2}")

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
                with open(file, "r") as infile:
                    mytext=infile.read()
                    mytext=list(filter(None, mytext.split("\n"))) #remove empty strings
                    dict_of_counts[sample] = mytext

            dataframe = pd.DataFrame(dict_of_counts)
            dataframe.index = ["starting_reads","trimmed_reads","aligned_reads","dedup_reads","unmapped_stage1","aligned_stage2","dedup_stage2"]
            dataframe.to_csv(output[0], sep = '\t')

# Introduce homer 
if config["two_tier"] == "TRUE" :
    rule create_tagdir:
        input:
            r1 = "output/bam/{sample}_R1_merged.bam",
            r2 = "output/bam/{sample}_R2_merged.bam"
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
            shell("samtools view -b {input.r1} {params.main_chr} > {output.r1}")
            shell("samtools view -b {input.r2} {params.main_chr} > {output.r2}")
            restriction = config["restriction_seq"].split(sep=",")[0]
            shell("{params.homer_bin}/makeTagDirectory output/homer/{wildcards.sample} {output.r1},{output.r2} -tbp 1 -genome {params.homer_genome}")
            shell("{params.homer_bin}/makeTagDirectory output/homer/{wildcards.sample} -update -genome {params.homer_genome}\
                -removePEbg -restrictionSite {restriction} -both -removeSelfLigation -removeSpikes 10000 5 2> output/logs/{wildcards.sample}_tag_log.txt")
            shell("touch {output.timecheck}")

else :
    rule create_tagdir:
        input:
            r1 = "output/bam/{sample}_R1_dupsremoved.bam",
            r2 = "output/bam/{sample}_R2_dupsremoved.bam"
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
            shell("samtools view -b {input.r1} {params.main_chr} > {output.r1}")
            shell("samtools view -b {input.r2} {params.main_chr} > {output.r2}")
            restriction = config["restriction_seq"].split(sep=",")[0]
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

rule find_tadloops:
    input:
        required = "output/homer/{sample}/tagLengthDistribution.txt",
        timecheck = "output/homer/{sample}/Test_tagdir_completion.txt"
    output:
        tad = "output/homer/{sample}/{sample}.tad.2D.bed"
    log:
        "output/logs/{sample}.find_tadloops.log"
    params:
        homer_bin = config["homer_bin"],
        homer_genome = config["homer_genome"],
        homer_mask = config["homer_mask"]
    run:
        shell("PATH={params.homer_bin}:$PATH; {params.homer_bin}/findTADsAndLoops.pl find  output/homer/{wildcards.sample} \
            -cpu 1 -res 3000 -window 15000 -genome {params.homer_genome} -p {params.homer_mask}")

rule merge_tadloops:
    input:
        tad = expand("output/homer/{sample}/{sample}.tad.2D.bed", sample=SAMPLES),
        loop = expand("output/homer/{sample}/{sample}.loop.2D.bed", sample=SAMPLES)
    output:
        "output/homer/Merged.tad.scores.txt"
    params:
        homer_bin = config["homer_bin"],
        homer_genome = config["homer_genome"]
    log:
        "output/logs/merge_tadloops.log"
    run:
        #Assemble the 
        myfiles=' '.join(input.tad)
        shell("PATH={params.homer_bin}:$PATH; merge2Dbed.pl {myfiles} -tad > output/homer/Merged.tad.2D.bed")
        myfiles=' '.join(input.loop)
        shell("PATH={params.homer_bin}:$PATH; merge2Dbed.pl {myfiles} -loop > output/homer/Merged.loop.2D.bed")
        myfiles=re.sub(r'/[^/]*.loop.2D.bed', '', myfiles)
        shell("PATH={params.homer_bin}:$PATH; findTADsAndLoops.pl score -tad output/homer/Merged.tad.2D.bed \
            -loop output/homer/Merged.loop.2D.bed -o \"output/homer/Merged\" -d {myfiles} \
            -cpu 1 -res 3000 -window 15000 -genome {params.homer_genome}")
