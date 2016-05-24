## ----compile-settings, include=FALSE-------------------------------------
library("methods")
library("knitr")
opts_chunk$set(tidy = FALSE, warning = FALSE, message = FALSE, 
               cache = FALSE, comment = NA, verbose = TRUE)
basename <- gsub(".Rmd", "", knitr:::knit_concord$get('infile')) 


## ------------------------------------------------------------------------
library('RNeXML')
data(bird.orders)

## ------------------------------------------------------------------------
birds <- add_trees(bird.orders)
birds <- add_basic_meta(
  title = "Phylogeny of the Orders of Birds From Sibley and Ahlquist",

  description = "This data set describes the phylogenetic relationships of the
     orders of birds as reported by Sibley and Ahlquist (1990). Sibley
     and Ahlquist inferred this phylogeny from an extensive number of
     DNA/DNA hybridization experiments. The ``tapestry'' reported by
     these two authors (more than 1000 species out of the ca. 9000
     extant bird species) generated a lot of debates.

     The present tree is based on the relationships among orders. The
     branch lengths were calculated from the values of Delta T50H as
     found in Sibley and Ahlquist (1990, fig. 353).",

  citation = "Sibley, C. G. and Ahlquist, J. E. (1990) Phylogeny and
     classification of birds: a study in molecular evolution. New
     Haven: Yale University Press.",

  creator = "Sibley, C. G. and Ahlquist, J. E.",
	nexml=birds)

## ------------------------------------------------------------------------
birds <- add_basic_meta(citation = citation("ape"), nexml = birds)

## ----message=FALSE, results='hide'---------------------------------------
birds <- taxize_nexml(birds, "NCBI")

## ------------------------------------------------------------------------
prefixes <- get_namespaces(birds)
prefixes["dc"]

## ------------------------------------------------------------------------
modified <- meta(property = "prism:modificationDate", content = "2013-10-04")

## ------------------------------------------------------------------------
birds <- add_meta(modified, birds) 

## ------------------------------------------------------------------------
history <- meta(property = "skos:historyNote",
  content = "Mapped from the bird.orders data in the ape package using RNeXML")

## ------------------------------------------------------------------------
birds <- add_meta(history, 
                birds, 
                level = "trees",
                namespaces = c(skos = "http://www.w3.org/2004/02/skos/core#"))

## ------------------------------------------------------------------------
birds

## ------------------------------------------------------------------------
meta <- get_metadata(birds) 

## ------------------------------------------------------------------------
names(meta)[1:4]

## ------------------------------------------------------------------------
meta[["dc:title"]]

## ------------------------------------------------------------------------
prefixes <- get_namespaces(birds)
prefixes["dc"]

## ----get_citation--------------------------------------------------------
get_citation(birds)

## ----get_taxa------------------------------------------------------------
get_taxa(birds)

## ------------------------------------------------------------------------
otu_meta <- get_metadata(birds, level="otu")
otu_meta[1:4]

