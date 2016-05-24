## ----supplement-compile-settings, include=FALSE--------------------------
library("methods")
library("knitr")
opts_chunk$set(tidy = FALSE, warning = FALSE, message = FALSE, 
               cache = FALSE, comment = NA, verbose = TRUE)
basename <- 'S4' 


## ----include=FALSE-------------------------------------------------------
library("RNeXML")

## ------------------------------------------------------------------------
nex <- new("nexml")

## ------------------------------------------------------------------------
slotNames(nex)

## ------------------------------------------------------------------------
nex@version
nex@namespaces

## ------------------------------------------------------------------------
class(nex@characters)

## ------------------------------------------------------------------------
length(nex@characters)

## ------------------------------------------------------------------------
nex@characters <- new("ListOfcharacters", list(new("characters")))

## ------------------------------------------------------------------------
length(nex@characters)

## ------------------------------------------------------------------------
class(nex@characters[[1]])

## ------------------------------------------------------------------------
f <- system.file("examples", "trees.xml", package="RNeXML")
nex <- nexml_read(f)

## ------------------------------------------------------------------------
nex@otus[[1]]@otu[[1]]@label

## ------------------------------------------------------------------------
nex@otus[[1]]@otu[[1]]@meta <- 
  c(meta("skos:note", 
          "This species was incorrectly identified"),
         nex@otus[[1]]@otu[[1]]@meta)

