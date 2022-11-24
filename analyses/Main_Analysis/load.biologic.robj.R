## This is to load the biologic Robj.
# The biologic R-object is expected in data/biologic_active_object

###############################################################################
## Recommended R-version                                                     ##

# module purge;source /camp/stp/babs/working/software/modulepath_new_software_tree_2018-08-13;module load pandoc/2.2.3.2-foss-2016b;ml R/4.0.3-foss-2020a

# Check processes
# lsof /dev/pts/*

## Done                                                                      ##
###############################################################################


###############################################################################
## Set the environment                                                       ##

# if (!require("remotes")){
#   install.packages("remotes")
# }
#
# remotes::install_github("rstudio/renv")

if (!file.exists("renv.lock")){
  renv::init()
} else {
  renv::restore(
    prompt=FALSE
  )
}



#renv::install("bioc::DESeq2")
#renv::install("bioc::clusterProfiler")
#renv::install("decusInLabore/biologicSeqTools2")
#renv::install("jokergoo/ComplexHeatmap")

## Done                                                                      ##
###############################################################################

###############################################################################
## Load biologic object                                                      ##

# Scripts will run for now in projectDir/scripts/bulkrnaseq_workflow
check <- list.files("../../../../data/biologic_active_object/")
check <- check[grep("bioLOGIC.Robj$", check)]

if (length(check) == 0){
    stop(paste0("Check if a biologic object has been initiated and is stored in [projectDir]/data/biologic_active_object/. "))
} else if (length(check) > 1){
    stop(paste0("More than one .biologic.Robj file is present in [projectDir]/data/biologic_active_object/. Please move the outdated biologic object into another folder"))
} else {
    library(biologicSeqTools2)
    base::load(paste0("../../../../data/biologic_active_object/", check))
}

## Done loading biologic object                                              ##
###############################################################################
