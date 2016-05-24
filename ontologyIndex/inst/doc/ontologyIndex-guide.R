## ------------------------------------------------------------------------
suppressPackageStartupMessages(library(ontologyIndex))
data(hpo)

## ----eval=FALSE----------------------------------------------------------
#  ontology <- get_ontology(file, qualifier)

## ---- echo=FALSE---------------------------------------------------------
data.frame(property=names(hpo), class=sapply(hpo, class), stringsAsFactors=FALSE, row.names=NULL)

## ------------------------------------------------------------------------
hpo$name["HP:0001873"]
hpo$ancestors[["HP:0001873"]]
hpo$name[hpo$ancestors[["HP:0001873"]]]

## ------------------------------------------------------------------------
minimal_set(hpo, c("HP:0001871", "HP:0001873", "HP:0011877"))

## ------------------------------------------------------------------------
get_ancestors(hpo, c("HP:0001873", "HP:0011877"))

