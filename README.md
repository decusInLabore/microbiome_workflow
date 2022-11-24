# bulkRNAseq_workflow


* Creeate a project directory
* cd into that project directory
* Create required directories
```
mkdir scripts workdir basedata FASTQ_files

cd scripts

git clone git@github.com:decusInLabore/bulkRNAseq_workflow.git

cd bulkRNAseq_workflow/analyses/Main_Analysis

```
Start R

```
module purge;source /camp/stp/babs/working/software/modulepath_new_software_tree_2018-08-13;module load pandoc/2.2.3.2-foss-2016b;ml R/4.0.3-foss-2020a;R
```

Adjust paths in file Part_0_prepare_alignment_and_input_files_for_Crick_RNA_Seq.Rmd


# Quickstart
* Minimal input file requirement: design file, (RSEM-) count file, res(dd) DEseq2 outputs for LRT and differential gene expression comparisons (either generated in this package or from outside sources), a PCA file (optional) and a file containing the most variable genes (optional). 
* Create in the projectdirectory a directory with the name ```mkdir wordir```. Clone into that this github repo ``` git clone git@github.com:decusInLabore/bulkRNAseq_workflow.git ```. Move into that directory ```cd bulkRNAseq_workflow```
* Review the <a href="https://github.com/decusInLabore/bulkRNAseq_workflow/blob/main/Part_0_prepare_input_files.Rmd" target="_blank">Part_0_prepare_input_files.Rmd</a> script
* Edit the <a href="https://github.com/decusInLabore/bulkRNAseq_workflow/blob/main/PartA_Automatic_Setup.Rmd" target="_blank">PartA_Automatic_Setup.Rmd</a> setup script. Change the parameters of the biologicSeqTools2::assembleBiologicProject() function according to your needs. This is the only place in which project-specific parameters need to be specified. 
* Run the workflow by doing 
```
sh runWorkflow.sh
```
* Additional html files (please remove all dots (.) from the file name) can be uploaded to the web-server by doing
```
cp /path/to/your/additional/html_file.html /camp/stp/babs/www/boeings/bioLOGIC_external/data/[project_id]/html/
```

* Powerpoint slide shows can be added to the project as follows: In powerpoint, export the slideshow as JPEG file and save these in a slide_folder. Then copy the slide folder to the webserver like so:
```
cp -r /path/to/your/slide_folder /camp/stp/babs/www/boeings/bioLOGIC_external/data/[project_id]/slides
```


This workflow can be run in three modes: 
* Option 1: Starting from raw fastq files, performing the alignment, differential gene expression analysis and interactive visualization
* Option 2: Starting from a (RSEM-) read count matrix performing differential gene expression analysis and interactive visualization 
* Option 3: Starting from DEseq2 output files to create an interactive data visulization

Preparation of the input files is lined out in more detail in the <a href="https://github.com/decusInLabore/bulkRNAseq_workflow/blob/main/Part_0_prepare_input_files.Rmd" target="_blank">Part_0_prepare_input_files.Rmd</a> script. 

# Required input files
## Required for all options

## Only required for option 1
### Sample specification sheet
The easiest is create a basedesign file with the following columns:
* |sampleID| sequencing facility sample id, optional
* |sample.id|  suggestion for a format: [dataseries] _ [sample.group] _ [replicate]
* |sample.group|
* |dataseries|
* |comp_1|comp_2|...|comp_N|
* |LRT_Treatment|LRT_...|
* |f_experimental_factor_1|f_experimental_factor_2|...
* |dataseries_colors| may be specified in a dataseries_color column with a unique
 hex code (#FF0000) for each dataseries
* sample.group_colors may be specified in a sample.group_color column with one entry per 
* sample.group

Once done, save this file in 

This file can be saved in projectFolder/data/base.design.txt

### Meta data sheet

# Option 2: Starting from preprossesed outside data

#### Input Format External DESeq2 Analysis
One possible input format for the DESeq2 result files can be reviewed in the *example_DESeq2_inputs* folder. Create one folder for differential gene expression type res(dds) outputs and one folder for LRT-type res(dds) outputs. You can either save the res(dds) outputs right away (with row.names = TRUE) or add a gene_id column carrying the gene alignment gene identifiers. The file name of the dds output file will be used for the data visualization. 

#### Load DESeq2 input files from Example example_DESeq2_inputs (see also file vis_project_partA.r)

*(1) Design file*

All sample.ids given in the sample.id column of the design file have to be present as column names in the TPM/normalizedCounts file. In this example that is the dfTPM file. 
```
dfDesign <- read.delim(
  "example_DESeq2_inputs/design.txt", 
  sep = "\t",
  stringsAsFactors = F
)

head(dfDesign)
```
*(2) Contrast Table (e.g. log-fold changes)*
Deposit one DESeq2 result file per comparison and name the file A_condition_vs_B_condition.txt. This will be the name of the comparison. 
```
DEseq2resultDir <- "example_DESeq2_inputs/DESeq2"
allfiles <- paste0(DEseq2resultDir, "/", list.files(DEseq2resultDir))
allfiles <- allfiles[grep(".txt", allfiles)]
contrastNames <- gsub(paste0(DEseq2resultDir, "/"), "", allfiles)
contrastNames <- gsub(".txt", "", contrastNames)
primaryAlignmentGeneID <- Obio@parameterList$primaryAlignmentGeneID
primaryAlignmentGeneID <- "ENSMUSG"

View example for a single comparison file
res <- read.delim(allfiles[1], header = T, sep="\t")
head(res)

for (i in 1:length(allfiles)){
  colName <- contrastNames[i]
  res <- read.delim(allfiles[i], header = T, sep="\t")
  names(res) = paste(names(res), colName, sep="_")
  res[[primaryAlignmentGeneID]] = rownames(res)
  
  
  names(res) = gsub("log2FoldChange", "logFC", names(res))
  names(res) = gsub(
    "logFC",
    paste("contrast_", i, "_logFC", sep=""),
    names(res)
  )
  
  names(res) = gsub(
    "padj",
    paste("contrast_", i, "_padj", sep=""),
    names(res)
  )
  
  names(res) = gsub(
    "stat",
    paste("contrast_", i, "_stat", sep=""),
    names(res)
  )
  
  res$baseMean <- log2(res$baseMean)
  names(res) = gsub(
    "baseMean",
    paste("contrast_", i, "_lg2BaseMean", sep=""),
    names(res)
  )
  
  #Remove all rows without a padj
  padj.col = grep("padj", names(res))[1]
  res[,padj.col][is.na(res[,padj.col])] = ""
  res = res[res[,padj.col] != "", ]
  res[,padj.col] <- as.numeric(res[,padj.col])
  
  ## Add log10p column ##
  padj  <- names(res)[grep("_padj_", names(res))]
  lg10p <- gsub("padj", "lg10p", padj)
  
  for (z in 1:length(padj)){
    preprocess <- as.numeric(res[,padj[z]])
    minNum <- min(preprocess[preprocess != 0])
    preprocess[preprocess == 0] <- minNum
    
    # if (length(grep("padj_LRT", padj[i])) > 0){
    #     preprocess <- as.numeric(res[,padj[z]])
    #     minNum <- min(preprocess[preprocess != 0])
    #     preprocess[preprocess == 0] <- minNum
    # } else {
    #     preprocess <- as.numeric(res[,padj[z]])
    # }
    
    temp <- -1*log10(preprocess)
    #temp[temp >= 50] = 50
    res[,lg10p[z]] <- temp
  }
  
  col.vector = c(
    primaryAlignmentGeneID,
    names(res)[grep("contrast", names(res))]
  )
  
  res = res[,col.vector]
  
  ## Make all numeric columns numeric ##
  res[,grep("contrast_", names(res))] <- apply(res[,grep("contrast_", names(res))], 2, as.numeric)
  
  if (i == 1){
    dfContrastTable <- res
  } else {
    dfContrastTable <- merge(
      dfContrastTable,
      res,
      by.x = primaryAlignmentGeneID,
      by.y = primaryAlignmentGeneID,
      all = TRUE
    )
    dfContrastTable[is.na(dfContrastTable)] <- 0
  }
}


head(dfContrastTable)
```

(3) Sample-level normalized count or TPM values
```
dfTPM <- read.delim(
  "example_DESeq2_inputs/dfTPM.txt", 
  sep = "\t",
  stringsAsFactors = F
)

head(dfTPM)
```

(4) PCA Table
```
dfPCA <- read.delim(
    "example_DESeq2_inputs/dfPCA.txt", 
    header = T,
    sep ="\t",
    stringsAsFactors = F
)

head(dfPCA)

## Now the most variable genes ##
FNvar <- "example_DESeq2_inputs/most.variable.features.txt"
dfVar <- read.delim(
    FNvar,
    header = T,
    sep = "\t",
    stringsAsFactors = F
)

head(dfVar)
```


## Custom Creation of Input Files
#### Create Contrast Input
Create a folder that contains all DEseq2 outputs in the NFcore_bulkRNAseq_vis folder, e.g. NFcore_bulkRNAseq_vis/DEseq2, within the project. In this folder, save all the DESeq2 outputs as follows (with dds being the DESeq2 object):

```
library(DESeq2)
res <- results(dds, contrast = contrast.vector)
res = data.frame(res)
write.table(res, "../DESeq2/A_Conditon_vs_B_Condition.txt", sep="\t")
```

Start running script 
vis_project_partA.r

#### Create PCA Input

Use either 
```
rld <- rlog(dds)
```

for few samples or 

```
rld <- vst(dds)
```

for many samples.

Then select the number of most variable genes to consider - 500 by default

```
Ntop4pca <- 500
```
Then get the PCA coordinates

```
rv <- rowVars(assay(rld))
select <- order(rv, decreasing = TRUE)[seq_len(Ntop4pca)]
pcaSelectionVec <- row.names(assay(rld)[select, ])
pca = prcomp(t(assay(rld)[pcaSelectionVec, ]))

PCApercentVar <- pca$sdev^2/sum(pca$sdev^2)

## Add percent variation plot ##
PercentVariation <- round(100*PCApercentVar,1)
PCdimension <- paste0("PC", 1:length(PercentVariation))
df <- data.frame(
    PercentVariation,
    PCdimension
)

df <- df[df$PercentVariation > 0,]


selVec <- c("sample.id", "sample.group", "sample.group_color")
selVec <- selVec[selVec %in% names(dfDesign)]
df.design.pca <- unique(dfDesign[,selVec])
df.pca = data.frame(pca$x)
df.pca[["sample.id"]] <- row.names(df.pca)

df.pca <- merge(
    df.design.pca,
    df.pca,
    by.x = "sample.id",
    by.y = "sample.id"
)

df.pca <- df.pca[order(df.pca$sample.id),]
names(df.pca) <- gsub("[.]", "_", names(df.pca))

dfPCA <- df.pca
```


## Create R-object for the project
Open the bulkRNAseq_workflow/PartA_Automatic_Setup.Rmd script in a R-editor such as RStudio.

Edit the entries according to your project. 

To run the script after editing, use the following R-version, if possible:
```
module purge;source /camp/stp/babs/working/software/modulepath_new_software_tree_2018-08-13;module load pandoc/2.2.3.2-foss-2016b;ml R/4.0.3-foss-2020a
```

Once done, save and run the script.

Section 1 of this script will create the *bulkRNAseq_workflow/design/biologic.settings.file.csv* file
This file contains all settings for the default data analysis and visualization. 

Section 2 of this script will create a design file for a Crick experiment. 


