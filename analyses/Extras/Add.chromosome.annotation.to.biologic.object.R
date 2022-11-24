---
title: "Add Chromosome Annotation"
author: "Stefan Boeing stefan.boeing@crick.ac.uk"
date: 'Compiled: `r format(Sys.Date(), "%B %d, %Y")`'
abstract: "This module provides a correlation analysis for a bulk RNA-Sequencing Experiment."
output: 
    html_document:
        code_folding: hide
        df_print: tibble
        toc: true
        toc_depth: 5
        toc_float: true
        css: src/style/style.css

always_allow_html: yes

---

```{css setup_css, echo=FALSE}


.table{
  width:auto;
  font-size: 10px;
}

```

```{r setup, include=FALSE}
knitr::opts_chunk$set(
    tidy = TRUE,
    tidy.opts = list(width.cutoff = 120),
    message = FALSE,
    warning = FALSE
)

# module purge;source /camp/stp/babs/working/software/modulepath_new_software_tree_2018-08-13;module load pandoc/2.2.3.2-foss-2016b;ml R/4.0.3-foss-2020a

if (!require("remotes")){
  install.packages("remotes")
}

remotes::install_github("rstudio/renv")

if (!file.exists("renv.lock")){
    renv::init()
} else {
    renv::restore(prompt = FALSE)
}

urlString <- "biologic.crick.ac.uk"
```



```{r set_directories, eval=T}
## Setup plot collection object
library(knitr)
library(ggplot2)
library(ggpubr)
library(DT)
library(biologicSeqTools2)

addCorCatsToLabDb <- FALSE
figureCount <- 1
chnkVec <- as.vector(NULL, mode = "character")

VersionPdfExt <- paste0(".V", gsub("-", "", Sys.Date()), ".pdf")

if (dir.exists("/Volumes/babs/working/boeings/")){
    hpc.mount <- "/Volumes/babs/working/boeings/"
} else if (dir.exists("Y:/working/boeings/")){
    hpc.mount <- "Y:/working/boeings/"
} else if (dir.exists("/camp/stp/babs/working/boeings/")){
    hpc.mount <- "/camp/stp/babs/working/boeings/"
} else {
    hpc.mount <- ""
}





## Heatmap setup ##
## Select Heatmap samples ##
FN <- paste0(hpc.mount, "Projects/reference_data/pwd_folder/babs.txt")
dbTable <- read.delim(
  FN,
  header = F,
  sep = "\t",
  stringsAsFactors = F
)

db.pwd <- as.vector(dbTable[1,1])
#setwd(Obio@parameterList$localWorkDir)

# The active biologic data object is expected to be found in ../../../../data/biologic_active_object/
source("load.biologic.robj.R")

# ObioFN <- paste0("../",list.files("..")[grep(".bioLOGIC.Robj", list.files(".."))])
# 
# load(ObioFN)


Obio <- setMountingPoint(Obio)
Obio <- setAnalysisPaths(Obio)
Obio <- setCrickGenomeAndGeneNameTable(Obio)
Obio <- createAnalysisFolders(
    Obio
)
Obio <- setDataBaseParameters(Obio)


Obio@parameterList[["reportFigDir"]] <- paste0(Obio@parameterList$localWorkDir,Obio@parameterList$project_id, "/report_figures/")

## Create outputfolders ##
if (!dir.exists(paste0(Obio@parameterList$localWorkDir,Obio@parameterList$project_id))){
    dir.create(paste0(Obio@parameterList$localWorkDir,Obio@parameterList$project_id))
}

if (!dir.exists(Obio@parameterList$reportFigDir)){
    dir.create(Obio@parameterList$reportFigDir)
}

figureCount <- 1

```
