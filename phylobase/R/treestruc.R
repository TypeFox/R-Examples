
##' Test trees for polytomies, inline nodes (singletons), or reticulation
##' 
##' Methods to test whether trees have (structural) polytomies, inline
##' nodes (i.e., nodes with a single descendant), or reticulation
##' (i.e., nodes with more than one ancestor). \code{hasPoly} only
##' check for structural polytomies (1 node has more than 2
##' descendants) and not polytomies that result from having edges with
##' a length of 0.
##' 
##' @aliases hasSingle
##' @param object an object inheriting from class \code{phylo4}
##' @return Logical value
##' @note Some algorithms are unhappy with structural polytomies (i.e., >2
##' descendants from a node), with single-descendant nodes, or with
##' reticulation; these functions check those properties.  We haven't bothered
##' to check for zero branch lengths: the consensus is that it doesn't come up
##' much, and that it's simple enough to test \code{any(edgeLength(x) == 0)} in
##' these cases.  (Single-descendant nodes are used e.g. in OUCH, or in other
##' cases to represent events occurring along a branch.)
##' @author Ben Bolker
##' @rdname treeStructure-methods
##' @export
##' @keywords misc
##' @examples
##' 
##' tree.owls.bis <- ape::read.tree(text="((Strix_aluco:4.2,Asio_otus:4.2):3.1,Athene_noctua:7.3);")
##' owls4 <- as(tree.owls.bis, "phylo4")
##' hasPoly(owls4)
##' hasSingle(owls4)
##'
setGeneric("hasSingle", function(object) {
    standardGeneric("hasSingle")
})

##' @rdname treeStructure-methods
##' @aliases hasSingle,phylo4-method
setMethod("hasSingle", signature(object="phylo4"),
  function(object) {
      if (nEdges(object) == 0) {
          return(FALSE)
      }
      ## This is about 3 times slower than using the C++
      ## function tabulateTips
      ##  degree <- tabulate(edges(object, drop.root=TRUE)[, 1])
      degree <- tabulateTips(object@edge[, 1])
      any(degree == 1)
})

##' @rdname treeStructure-methods
##' @aliases hasRetic
##' @export
setGeneric("hasRetic", function(object) {
    standardGeneric("hasRetic")
})

##' @rdname treeStructure-methods
##' @aliases hasRetic,phylo4-method
setMethod("hasRetic", signature(object="phylo4"), function(object) {
    if (nEdges(object)==0) {
        return(FALSE)
    }
    ## this is about the same (slightly faster on 10,000 tips)
    ##  than using the C++ function
    ancest <- tabulate(edges(object)[, 2])
    any(ancest > 1)    
})

##' @rdname treeStructure-methods
##' @aliases hasPoly
##' @export
setGeneric("hasPoly", function(object) {
    standardGeneric("hasPoly")
})

##' @rdname treeStructure-methods
##' @aliases hasPoly,phylo4-method
setMethod("hasPoly", signature(object="phylo4"), function(object) {
  if (nEdges(object)==0) {
      return(FALSE)
  }
  ## This is about 3 times slower than using the C++
  ## function tabulateTips
  ## degree <- tabulate(edges(object, drop.root=TRUE)[, 1])
  degree <- tabulateTips(object@edge[, 1])
  any(degree > 2)
})
