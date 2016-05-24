## ----supplement-compile-settings, include=FALSE--------------------------
library("methods")
library("knitr")
opts_chunk$set(tidy = FALSE, warning = FALSE, message = FALSE, 
               cache = FALSE, comment = NA, verbose = TRUE, eval=require("rrdf"))
basename <- 'sparql' 



## ----include=FALSE-------------------------------------------------------
#  library("RNeXML")

## ------------------------------------------------------------------------
#  library("rrdf")
#  library("XML")
#  library("phytools")
#  library("RNeXML")

## ------------------------------------------------------------------------
#  nexml <- nexml_read(system.file("examples/primates.xml", package="RNeXML"))

## ------------------------------------------------------------------------
#  rdf <- get_rdf(system.file("examples/primates.xml", package="RNeXML"))
#  tmp <- tempfile()  # so we must write the XML out first
#  saveXML(rdf, tmp)
#  graph <- load.rdf(tmp)

## ------------------------------------------------------------------------
#  root <- sparql.rdf(graph,
#  "SELECT ?uri WHERE {
#      ?id <http://rs.tdwg.org/ontology/voc/TaxonConcept#rank> <http://rs.tdwg.org/ontology/voc/TaxonRank#Order> .
#      ?id <http://rs.tdwg.org/ontology/voc/TaxonConcept#toTaxon> ?uri
#  }")

## ------------------------------------------------------------------------
#  get_name <- function(id) {
#    max <- length(nexml@otus[[1]]@otu)
#    for(i in 1:max) {
#      if ( nexml@otus[[1]]@otu[[i]]@id == id ) {
#        label <- nexml@otus[[1]]@otu[[i]]@label
#        label <- gsub(" ","_",label)
#        return(label)
#      }
#    }
#  }

## ------------------------------------------------------------------------
#  recurse <- function(node){
#  
#      # fetch the taxonomic rank and id string
#      rank_query <- paste0(
#          "SELECT ?rank ?id WHERE {
#              ?id <http://rs.tdwg.org/ontology/voc/TaxonConcept#toTaxon> <",node,"> .
#              ?id <http://rs.tdwg.org/ontology/voc/TaxonConcept#rank> ?rank
#            }")
#      result <- sparql.rdf(graph, rank_query)
#  
#      # get the local ID, strip URI part
#      id <- result[2]
#      id <- gsub("^.+#", "", id, perl = TRUE)
#  
#      # if rank is terminal, return the name
#      if (result[1] == "http://rs.tdwg.org/ontology/voc/TaxonRank#Species") {
#          return(get_name(id))
#      }
#  
#      # recurse deeper
#      else {
#          child_query <- paste0(
#              "SELECT ?uri WHERE {
#                  ?id <http://www.w3.org/2000/01/rdf-schema#subClassOf> <",node,"> .
#                  ?id <http://rs.tdwg.org/ontology/voc/TaxonConcept#toTaxon> ?uri
#              }")
#          children <- sparql.rdf(graph, child_query)
#  
#          return(paste("(",
#                       paste(sapply(children, recurse),
#                             sep = ",", collapse = "," ),
#                       ")",
#                       get_name(id), # label interior nodes
#                       sep = "", collapse = ""))
#      }
#  }
#  

## ------------------------------------------------------------------------
#  newick <- paste(recurse(root), ";", sep = "", collapse = "")
#  tree <- read.newick(text = newick)
#  collapsed <- collapse.singles(tree)
#  plot(collapsed,
#       type='cladogram',
#       show.tip.label=FALSE,
#       show.node.label=TRUE,
#       cex=0.75,
#       edge.color='grey60',
#       label.offset=-9)

