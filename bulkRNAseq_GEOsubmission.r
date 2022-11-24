###############################################################################
## GEO Submission Metadata     ################################################
###############################################################################

###############################################################################
# Schaefer lab                                                                #
# Experiment: Worm and Chemical Infection                                     #
# Primary data anlalysis                                                      #
# Ensembl reference                                                           #
###############################################################################


rm(list = ls())

###############################################################################
## Initializing                                                              ##

#library(bioLOGIC)

if (dir.exists("/Volumes/babs/working/boeings/")){
    hpc.mount <- "/Volumes/babs/working/boeings/"
} else if (dir.exists("Y:/working/boeings/")){
    hpc.mount <- "Y:/working/boeings/"
} else if (dir.exists("/camp/stp/babs/working/boeings/")){
    hpc.mount <- "/camp/stp/babs/working/boeings/"
} else {
    hpc.mount <- ""
}

# source(
#     paste0(
#         hpc.mount,
#         "Stefan/protocol_files/github/boeings/packages/cDev.package.SBwebtools.V7.r"
#     )
# )

source("assets/SBwebtools.pckg.r")

###############################################################################
## Set libpaths                                                              ##

if (length(.libPaths()) > 2){
    .libPaths(.libPaths()[2:3])
}

## Done setting Stefan's libpaths                                            ##
###############################################################################

###############################################################################
## Create S4 object for this analysis                                        ##

## Check parameters before proceeding ##



###############################################################################
## Set database password                                                     ##
# if (exists("db.pwd")){
#     print("Database password is set.")
# } else {
#     ## Set database password ##
#     library(rstudioapi)
#     db.pwd <- rstudioapi::askForPassword("Please enter database password")
# }

FN <- paste0(hpc.mount, "Projects/reference_data/documentation/BC.parameters.txt")
dbTable <- read.delim(
    FN,
    sep = "\t",
    stringsAsFactors = F
)

db.pwd <- as.vector(dbTable[1,1])
##                                                                           ##
###############################################################################

###############################################################################
## Add Crick parameters to report                                            ##

## Determine sequencer for report
## Crick Sequencing
## more /camp/stp/babs/working/data/crick_runs_instruments.csv | grep "190205_K00102_0309_AH37TNBBXY"
## CRUK Sequencing
## more /camp/stp/babs/working/data/lif_runs_instruments.csv | grep "190205_K00102_0309_AH37TNBBXY"
## Gathering flow cell details
## https://github.com/10XGenomics/supernova/blob/master/tenkit/lib/python/tenkit/illumina_instrument.py#L12-L45

## Determine read length of unknown fastq files
## zcat SRR2926047_pass_1.fastq.gz | awk '{if(NR%4==2) print length($1)}' | head

## Done                                                                      ##
###############################################################################



###############################################################################
## Set Project Parameters                                                    ##

## Load Obio object ##
ObioFN <- paste0("../",list.files("..")[grep(".bioLOGIC.Robj", list.files(".."))])

load(ObioFN)


# checkFile = paste0(
#     Obio@parameterList$project_id,
#     ".bioLOGIC.Robj"
# )
# 
# if (ObioFN != checkFile){
#     exit()
# }

Obio <- setMountingPoint(Obio)
Obio <- setAnalysisPaths(Obio)
Obio <- setCrickGenomeAndGeneNameTable(Obio)
Obio <- createAnalysisFolders(
    Obio,
    baseDir="/camp/stp/babs/working/boeings/Projects/",
    localBaseDir = paste0(hpc.mount, "Projects/")
)
Obio <- setDataBaseParameters(Obio)

Obio@parameterList$DEseq2Dir <- paste0(
    Obio@parameterList$localWorkDir,
    "DESeq2/"
)

## Set Project Parameters                                                    ##
###############################################################################

###############################################################################
## Add annotation file                                                       ##

Obio <- addGeneAnnotation(Obio)

## Done adding annotation file                                               ##
###############################################################################



# df.design <- read.delim(
#     design.file, 
#     header = TRUE, 
#     sep = "\t",
#     stringsAsFactors = FALSE
# )

df.design <- Obio@dfDesign



df.design$FASTQ <- df.design$sample.id
Sample.name <- df.design$FASTQ
title <-  df.design$FASTQ
source.name <-  df.design$FASTQ
organism <- rep(Obio@parameterList$species, length(Sample.name))
characteristics.tissue <- rep("", length(Sample.name))
characteristics.tag <- rep("", length(Sample.name))
molecule <- rep("mRNA", length(Sample.name))
description <- rep("", length(Sample.name))
description2 <- rep("", length(Sample.name))
processed.data.file <- rep(paste0(Obio@parameterList$RSEMcountDataFile), length(Sample.name))
original.file.location <- df.design$original.NGS
raw.file <- df.design$NGS
#raw.file <- gsub("195", "195A", raw.file)
raw.file <- gsub(Obio@parameterList$fastqDir, "", as.vector(raw.file))

raw.file2 <- df.design$original.NGS


dfGEO <- data.frame(
    Sample.name,
    title,
    source.name,
    organism,
    characteristics.tissue,
    characteristics.tag,
    molecule,
    description,
    description2,
    processed.data.file,
    raw.file
)

## Copy count file into folder
GEOdir <- gsub("workdir", paste0(Obio@parameterList$project_id, "_GEO_submission"), Obio@parameterList$localWorkDir)

if (!exists(GEOdir)){
    dir.create(GEOdir)    
}


cmd <- paste0( "cp ", Obio@parameterList$localWorkDir, "RSEM/",Obio@parameterList$RSEMcountDataFile, " ", GEOdir )
system(cmd)

## Copy fastq files into folder       
cmd <- paste0("cp ", df.design$NGS, " ", GEOdir)
sapply(cmd, function(x) system(x))

## Create md5sum
# md5sum * > checkfile.txt
md5FN <- paste0(
    GEOdir, 
    "checkfile.txt"
)


cmd <- paste0("md5sum ", GEOdir, "*.* > ", md5FN)
system(cmd)

# cp ../RSEM/sll357.count.data.txt .
# md5sum sll357.count.data.txt

dfMD5 <- read.delim(
    md5FN,
    header=F,
    sep = " ",
    stringsAsFactors = F
)
dfMD5$V2 <- NULL


names(dfMD5) <- c("md5sum", "FN")
dfMD5$FN <- gsub(GEOdir, "", as.vector(dfMD5$FN))

dfGEO <- merge(
    dfGEO,
    dfMD5,
    by.x = "raw.file",
    by.y = "FN",
    all = T
)

write.table(
    dfGEO,
    paste0(GEOdir, "GEO.outputfile.txt"),
    row.names = FALSE,
    sep = "\t"
)

###############################################################################
## Upload to GEO

ml purge
ml ncftp/3.2.6-foss-2016b
ncftpput -F -R -z -u geoftp -p "rebUzyi1" ftp-private.ncbi.nlm.nih.gov ./uploads/sjboulton_WFS5ymYo ./DRIPSeq_GEO

ml purge
ml ncftp/3.2.6-foss-2016b
ncftpput -F -R -z -u geoftp -p "rebUzyi1" ftp-private.ncbi.nlm.nih.gov ./uploads/schaefer_lab_KPsRsW2K /camp/stp/babs/working/boeings/Projects/schaefera/tobias.ackels/320_ASL_TA_bulk_RNAseq_mouse_olfactory_bulb_projection_RN19134/asl320_GEO_submission


