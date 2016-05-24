#' extract all phylogenetic trees in ape format
#' 
#' extract all phylogenetic trees in ape format
#' @param nexml a representation of the nexml object from  which the data is to be retrieved
#' @return returns a list of lists of multiphylo trees, even if all trees are in the same `trees` node (and hence the outer list will be of length 1) or if there is only a single tree (and hence the inner list will also be of length 1.  This guarentees a consistent return type regardless of the number of trees present in the nexml file, and also preserves any heirarchy/grouping of trees.  
#' @export
#' @import plyr 
#' @examples
#' comp_analysis <- system.file("examples", "comp_analysis.xml", package="RNeXML")
#' nex <- nexml_read(comp_analysis)
#' get_trees_list(nex)
#' @seealso \code{\link{get_trees}} \code{\link{get_flat_trees}} \code{\link{get_item}}
get_trees_list <- function(nexml) as(nexml, "multiPhyloList")



#' extract a phylogenetic tree from the nexml
#' 
#' extract a phylogenetic tree from the nexml
#' @param nexml a representation of the nexml object from  which the data is to be retrieved
#' @return an ape::phylo tree, if only one tree is represented.  Otherwise returns a list of lists of multiphylo trees.  To consistently recieve the list of lists format (preserving the heriarchical nature of the nexml), use \code{\link{get_trees_list}} instead.  
#' @export
#' @seealso \code{\link{get_trees}} \code{\link{get_flat_trees}} \code{\link{get_item}}
#' @examples
#' comp_analysis <- system.file("examples", "comp_analysis.xml", package="RNeXML")
#' nex <- nexml_read(comp_analysis)
#' get_trees(nex)
get_trees <- function(nexml) as(nexml, "phylo")


#' get_flat_trees  
#' 
#' extract a single multiPhylo object containing all trees in the nexml
#' @details Note that this method collapses any heirachical structure that may have been present as multiple `trees` nodes in the original nexml (though such a feature is rarely used).  To preserve that structure, use \code{\link{get_trees}} instead.  
#' @return a multiPhylo object (list of ape::phylo objects).  See details.  
#' @param nexml a representation of the nexml object from  which the data is to be retrieved
#' @export
#' @seealso \code{\link{get_trees}} \code{\link{get_trees}} \code{\link{get_item}} 
#' @examples
#' comp_analysis <- system.file("examples", "comp_analysis.xml", package="RNeXML")
#' nex <- nexml_read(comp_analysis)
#' get_flat_trees(nex)
get_flat_trees <- function(nexml) flatten_multiphylo(get_trees_list(nexml))



####### Coercion methods ######### 

setAs("nexml", "multiPhyloList", function(from){
   map <- get_otu_maps(from) 
   unname(lapply(from@trees, 
           function(X){
             out <- unname(lapply(X@tree,  toPhylo, map[[X@otus]]))
             class(out) <- "multiPhylo"
             out
           }))
})


# Always collapses all trees nodes into a multiphylo
setAs("nexml", "multiPhylo", function(from){
   map <- get_otu_maps(from) 
   out <- unname(lapply(from@trees, 
           function(X){
             out <- unname(lapply(X@tree,  toPhylo, map[[X@otus]]))
             class(out) <- "multiPhylo"
             out
           }))
  flatten_multiphylo(out)
})


#' Flatten a multiphylo object
#' 
#' @details NeXML has the concept of multiple <trees> nodes, each with multiple child <tree> nodes.
#' This maps naturally to a list of multiphylo  objects.  Sometimes
#' this heirarchy conveys important structural information, so it is not discarded by default. 
#' Occassionally it is useful to flatten the structure though, hence this function.  Note that this
#' discards the original structure, and the nexml file must be parsed again to recover it.  
#' @param object a list of multiphylo objects 
#' @export
flatten_multiphylo <- function(object){
  out <- unlist(object, FALSE, FALSE)
  class(out) <- "multiPhylo"
  out
}


setAs("nexml", "phylo", function(from){ 
    if(length(from@trees[[1]]@tree) == 1){
      maps <- get_otu_maps(from)
      otus_id <- from@trees[[1]]@otus
      out <- toPhylo(from@trees[[1]]@tree[[1]], maps[[otus_id]])
    } else { 
      warning("Multiple trees found, Returning multiPhylo object")
      out <- as(from, "multiPhylo") 
    }
    if(length(out) == 1)
     out <- flatten_multiphylo(out)
   out 
  })



########### Main internal function for converting nexml to phylo ########



#' nexml to phylo 
#' 
#' nexml to phylo coercion 
#' @param tree an nexml tree element 
#' @param otus a character string of taxonomic labels, named by the otu ids.  
#' e.g. (from get_otu_maps for the otus set matching the relevant trees node. 
#' @return phylo object.  If a "reconstructions" annotation is found on the 
#' edges, return simmap maps slot as well.  
toPhylo <- function(tree, otus){
  otu <- NULL # Avoid CRAN NOTE as per http://stackoverflow.com/questions/8096313/no-visible-binding-for-global-variable-note-in-r-cmd-check
  ## Extract the nodes list
  nodes <- sapply(unname(tree@node), 
                  function(x) 
                    c(node = unname(x@id), 
                      otu = missing_as_na(x@otu)))
  
  # If any edges have lengths, use this routine
  if(any(sapply(tree@edge, function(x) length(x@length) > 0)))
    edges <- sapply(unname(tree@edge), 
                    function(x) 
                      c(source = unname(x@source), 
                        target = unname(x@target), 
                        length = if(identical(x@length, numeric(0)))
                                  NA 
                                 else 
                                   unname(x@length), 
                        id = unname(x@id)))
  else # no edge lengths, use this routine 
    edges <- sapply(unname(tree@edge), 
                    function(x) 
                      c(source = unname(x@source), 
                        target = unname(x@target), 
                        id = unname(x@id)))

  nodes <- data.frame(t(nodes), stringsAsFactors=FALSE)
  names(nodes) <- c("node", "otu")

## Identifies tip.label based on being named with OTUs while others are NULL
## FIXME Should instead decide that these are tips based on the edge labels?
  nodes <- cbind(plyr::arrange(nodes, otu), id = 1:dim(nodes)[1])  # Also warns because arrange isn't quoting the column name.  

## NB: these ids are the ape:id numbers by which nodes are identified in ape::phylo
## Arbitrary ids are not supported - ape expecs the numbers 1:n, starting with tips. 

##  nodes$node lists tip taxa first (see arrange fn above), since
##  APE expects nodes numbered 1:n_tips to be to correspond to tips.
  source_nodes <- match(edges["source",], nodes$node)
  target_nodes <- match(edges["target",], nodes$node)

  ##  Define elements of a phylo class object ##
  #--------------------------------------------#

  ## define edge matrix
  edge <- unname(cbind(source_nodes, target_nodes))
  if("length" %in% rownames(edges))
    edge.length <- as.numeric(edges["length",])       
  else
    edge.length <- NULL

  ## define tip labels
  tip_otus <- as.character(na.omit(nodes$otu))   
  tip.label <- otus[tip_otus]

  # Count internal nodes (assumes bifurcating tree. Does ape always assume this?) 
  # FIXME use a method that does not assume bifurcating tree... 
  Nnode <- length(tip.label) - 1 

  # assemble the phylo object, assign class and return.  
  phy = list(edge=edge, 
             tip.label = unname(tip.label), 
             Nnode = Nnode)
  if(!is.null(edge.length))
    phy$edge.length = edge.length # optional fields
  class(phy) = "phylo"



  ## Check for simmap  


  phy
}


## Helper function
missing_as_na <- function(x){
  if(length(x) == 0)
    NA
  else
    unname(x)
}






