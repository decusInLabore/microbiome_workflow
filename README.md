# Microbiome Workflow


* Creeate a project directory
* cd into that project directory
* Create required directories
```
mkdir scripts workdir basedata FASTQ_files

cd scripts

git clone git@github.com:decusInLabore/microbiome_workflow.git

cd bulkRNAseq_workflow/analyses/Main_Analysis

```

# Alignment using the nf-core ampliseq pipeline
The nf-core ampliseq pipeline is described [here](https://nf-co.re/ampliseq).

To get started, fastq files need to be listed with their file paths and biological sample names in the samplesheet.tsv file. This needs to be a tab-deliminated file and it needs to be saved with the .tsv extension. 

## First round of running ampliseq to determine quality cut-offs
In a first instance the ampliseq pipelines is run on an paired-end sequencing experiment like this:
```
nextflow run nf-core/ampliseq \
    -r 2.3.2 \
    -profile crick \
    --input '../scripts/samplesheet.tsv' \
    --FW_primer GTGYCAGCMGCCGCGGTAA \
    --RV_primer GGACTACNVGGGTWTCTAAT \
    --outdir "../workdir/results" \
    --ignore_failed_trimming \
    --untilQ2import \
    -resume
 ```
 This run mainly is to determine good values for the --trunclenf and --trunclenr parameters in the next step. Leaving away the --untilQ2import flag will run the pipeline up until the data visualization module. 
 
## Second round of running the ampliseq pipeline to run up until the data visualization step
The second run of the pipeline uses the above parameters:
```
nextflow run nf-core/ampliseq \
    -r 2.3.2 \
    -profile crick \
    --input '../scripts/samplesheet.tsv' \
    --FW_primer GTGYCAGCMGCCGCGGTAA \
    --RV_primer GGACTACNVGGGTWTCTAAT \
    --outdir "../workdir/results" \
    --ignore_failed_trimming \
    --trunclenf 200 \
    --trunclenr 180 \
    -resume
    
```

## The third stage of running the pipeline will include the data visualization module
In order to activate the data visualization module, a metadata.tsv file needs to be provided. It is important that this file contains only sample names that have been successfully processed by the ampliseq pipeline. 

```
nextflow run nf-core/ampliseq \
    -r 2.3.2 \
    -profile crick \
    --input '../scripts/samplesheet.tsv' \
    --FW_primer GTGYCAGCMGCCGCGGTAA \
    --RV_primer GGACTACNVGGGTWTCTAAT \
    --outdir "../workdir/results" \
    --ignore_failed_trimming \
    --metadata '../scripts/metadata.tsv' \
    -resume
```



# Downstream Data Viualization and Analysis

