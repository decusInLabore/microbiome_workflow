###############################################################################
## Save Obio Object                                                          ##

biologic_active_object_dir <- paste0(
    Obio@parameterList$folder,
    "data/biologic_active_object/"
)

if (!dir.exists(biologic_active_object_dir )){
    dir.create(biologic_active_object_dir, recursive = T)
}

save(Obio,
     file = paste0(
       biologic_active_object_dir,
       Obio@parameterList$project_id,
       ".bioLOGIC.Robj"
     )
)

print(paste0("R bioLOGIC object saved in", biologic_active_object_dir,"."))

##                                                                           ##
###############################################################################


