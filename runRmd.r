args = commandArgs(trailingOnly = T)

rmarkdown::render(args[1], output_dir ="../html_local")
