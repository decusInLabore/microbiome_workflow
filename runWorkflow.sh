wait_on_lsf() { ## wait on jobs{
    sleep 300
    n=`squeue --name=$project  | wc -l`
    while [ $n -ne 1 ]
    do
    n=`squeue --name=$project  | wc -l`
    ((z=$n-1))
    #number of running
    echo "$project jobs running: $z"
    #number of pending
    sleep 300
    done
}

project_id=bio

project=rA_$project_id
sbatch --time=06:00:00 --wrap "module purge;source /camp/stp/babs/working/software/modulepath_new_software_tree_2018-08-13;module load pandoc/2.2.3.2-foss-2016b;ml R/4.0.3-foss-2020a;Rscript runRmd.r PartA_Automatic_Setup.Rmd" --job-name=$project --mem=7G -o ../$project.slurm >> ../commands.txt

wait_on_lsf

project=rB_$project_id
sbatch --time=06:00:00 --wrap "module purge;source /camp/stp/babs/working/software/modulepath_new_software_tree_2018-08-13;module load pandoc/2.2.3.2-foss-2016b;ml R/4.0.3-foss-2020a;Rscript runRmd.r PartB_Analysis.Rmd" --job-name=$project --mem=200G -o ../$project.slurm >> ../commands.txt

wait_on_lsf

project=rC_$project_id
sbatch --time=06:00:00 --wrap "module purge;source /camp/stp/babs/working/software/modulepath_new_software_tree_2018-08-13;module load pandoc/2.2.3.2-foss-2016b;ml R/4.0.3-foss-2020a;Rscript runRmd.r PartC_Database_Upload.Rmd" --job-name=$project --mem=200G -o ../$project.slurm >> ../commands.txt
