load.packages = function(noupdate=TRUE)
{
  if("BiocGenerics" %in% rownames(installed.packages()) == FALSE) {
    biocLite <- NULL; rm(biocLite)
    source("http://bioconductor.org/biocLite.R")
    biocLite(suppressUpdates=noupdate , ask=FALSE)
  }
  if("biomaRt" %in% rownames(installed.packages()) == FALSE) {
    biocLite <- NULL; rm(biocLite)
    source("http://bioconductor.org/biocLite.R")
    biocLite("biomaRt", suppressUpdates=noupdate , ask=FALSE)
  }
  if("sm" %in% rownames(installed.packages()) == FALSE) {install.packages("sm")}
  #library(sm,warn.conflicts = FALSE, quietly = TRUE,verbose =FALSE)
  if("VennDiagram" %in% rownames(installed.packages()) == FALSE) {install.packages("VennDiagram")}
  #library(VennDiagram,warn.conflicts = FALSE, quietly = TRUE,verbose =FALSE)
  if("knitr" %in% rownames(installed.packages()) == FALSE) {install.packages("knitr")}
  #requireNamespace("knitr")
  #library(knitr,warn.conflicts = FALSE, quietly = TRUE,verbose =FALSE)
  if("rmarkdown" %in% rownames(installed.packages()) == FALSE) {install.packages("rmarkdown")}
  #library(rmarkdown,warn.conflicts = FALSE, quietly = TRUE,verbose =FALSE)
  
  if("DESeq2" %in% rownames(installed.packages()) == FALSE) {
    biocLite <- NULL; rm(biocLite)
    source("http://bioconductor.org/biocLite.R")
    biocLite("DESeq2",suppressUpdates=noupdate , ask=FALSE)
  }
  #library(DESeq2,warn.conflicts = FALSE, quietly = TRUE,verbose =FALSE)
  
  #if("MESS" %in% rownames(installed.packages()) == FALSE) {install.packages("MESS")}
  #requireNamespace("MESS")
  #library(MESS,warn.conflicts = FALSE, quietly = TRUE,verbose =FALSE)
  if("scatterplot3d" %in% rownames(installed.packages()) == FALSE) {
    biocLite <- NULL; rm(biocLite)
    source("http://bioconductor.org/biocLite.R")
    biocLite("scatterplot3d",suppressUpdates=noupdate , ask=FALSE)
  }
  #library(scatterplot3d,warn.conflicts = FALSE, quietly = TRUE,verbose =FALSE)
  
  if("stringi" %in% rownames(installed.packages()) == FALSE) {install.packages("stringi")}
  #requireNamespace("stringi")
  #library(stringi,warn.conflicts = FALSE, quietly = TRUE,verbose =FALSE)
  
  if("xlsx" %in% rownames(installed.packages()) == FALSE) {install.packages("xlsx")}
  #library(xlsx,warn.conflicts = FALSE, quietly = TRUE,verbose =FALSE)
  
  if("caRpools" %in% rownames(installed.packages()) == FALSE) {install.packages("caRpools")}
  #library(caRpools,warn.conflicts = FALSE, quietly = TRUE,verbose =FALSE)
  
  return(TRUE)
}
