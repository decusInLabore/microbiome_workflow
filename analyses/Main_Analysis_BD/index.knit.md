---
output:
  bookdown::html_document2:
    includes:
      in_header: header.html
  bookdown::gitbook:
    includes:
      in_header: header.html

---


# Bulk RNA-Seq Analysis

  This is a _sample_ book written in **Markdown**. You can use anything that Pandoc's Markdown supports, e.g., a math equation $a^2 + b^2 = c^2$.

The **bookdown** package can be installed from CRAN or Github:


```r
install.packages("bookdown")
# or the development version
# devtools::install_github("rstudio/bookdown")
```

Remember each Rmd file contains one and only one chapter, and a chapter is defined by the first-level heading `#`.

To compile this example to PDF, you need XeLaTeX. You are recommended to install TinyTeX (which includes XeLaTeX): <https://yihui.name/tinytex/>.




```r
## Try to retrieve project data from db ##
db.pwd2 <- "zU3ufd9L"
db.user2 <- "reader"
host2 <- "clvd1-db-u-p-17.thecrick.org"
projectParams <- Obio@documentationParams

tryCatch({
    dbDB = DBI::dbConnect(
        drv = RMySQL::MySQL(),
        user = db.user2,
        password = db.pwd2,
        host = host2,
        dbname = "clarity_shadow"
    )
    dfProposal <-  DBI::dbGetQuery(
        dbDB,
        paste0("SELECT * FROM clarify_asf_proposals WHERE project_name ='",Obio@projectDetailList$lims.id,"'")
    )
    dbDisconnect(dbDB)
}, error = function(x) {
    message("Project Database could not be reached or has no entry in Obio@parameterList$lims.id for this analysis.")
})
```

```
## Project Database could not be reached or has no entry in Obio@parameterList$lims.id for this analysis.
```

```r
###############################################################################
## Helper
firstup <- function(x) {
    substr(x, 1, 1) <- toupper(substr(x, 1, 1))
    x
}
##
###############################################################################


if (exists("dfProposal")){
    if (!is.na(dfProposal[1,"ProjectAlias"]) & dfProposal[1,"ProjectAlias"] != ""){
        projectParams[["title"]] = paste0(dfProposal[1,"ProjectAlias"], " - ", dfProposal[1,"project_name"])
    }

    if (!is.na(dfProposal[1,"project_user"]) & dfProposal[1,"project_user"] != ""){
        labString <- firstup(dfProposal[1,"user_lab"])
        labString <- substr(labString, 1, (nchar(labString) - 1))

        projectParams[["subtitle"]] = paste0(labString, " Lab - ", dfProposal[1,"project_user"])
        projectParams[["subtitle"]] <- gsub("^ Lab - ", "", projectParams[["subtitle"]])

    }

    if (!is.na(dfProposal[1,"proposal_text"]) & dfProposal[1,"proposal_text"] != ""){
        projectParams[["abstract"]] = dfProposal[1,"proposal_text"]

    }
}

## Escape all special characters
projectParams <- lapply(
  projectParams, function(x)
  #gsub("([.|()\\^{}+$*?]|\\[|\\])", "\\\\\1", x)
  gsub("([.|()/\\^{}+$*?]|\\[|\\])", " ", x)
)

projectParams <- lapply(
  projectParams, function(x)
  #gsub("([.|()\\^{}+$*?]|\\[|\\])", "\\\\\1", x)
  gsub("\\\n", " ", x)
)
```



---
title: "Transcriptomic analysis of Sh2d2a knockout regulatory T cells following a time-course of anti-CD3 stimulation"
subtitle:  "David M Briscoe Lab; Literature Bulk-rna Seq GSE134515"
author:
    - Bioinformatics: Stefan Boeing^[The Francis Crick Institute, stefan.boeing@crick.ac.uk]
date: 'Compiled: July 26, 2022'

abstract: |
    "CD4 CD25high regulatory T cells were FACS-sorted from Sh2d2a knockout and wild type mice and activated with 1 µg ml anti-CD3  clone: 145-2C11  for 2, 8 and 24 hours  Experiments were performed in duplicate conditions "

description: "This is a minimal example of using the bookdown package to write a book. The output format for this example is bookdown::gitbook."
---
