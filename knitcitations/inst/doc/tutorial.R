## ------------------------------------------------------------------------
library("knitcitations")
cleanbib()

## ------------------------------------------------------------------------
options("citation_format" = "pandoc")

## ----echo=FALSE, include=FALSE, results="hide"---------------------------
## Dummy block to stop R check failing when it tries to tangle and source the vignette
citep(citation("knitr"))
citet("10.1098/rspb.2013.1372")

## ---- message=FALSE------------------------------------------------------
write.bibtex(file="references.bib")

