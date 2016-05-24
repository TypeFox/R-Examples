## classes for holding multiple tree objects

##' multiPhylo4 and extended classes
##' 
##' Classes for lists of phylogenetic trees.  These classes and methods are
##' planned for a future version of \code{phylobase}.
##' 
##' 
##' @name multiPhylo-class
##' @aliases multiPhylo-class multiPhylo4-class multiPhylo4d-class tbind
##' @docType class
##' @keywords classes
## @export
setClass("multiPhylo4", representation(phylolist = "list", 
    tree.names = "character"), prototype = list(phylolist = list(), 
    tree.names = character(0)))

setClass("multiPhylo4d", representation(tip.data = "data.frame"), 
    contains = "multiPhylo4")

setMethod("initialize", "multiPhylo4", function(.Object, ...) {
    message("multiPhylo and multiphylo4d not yet implemented", 
            "Try using a list of phylo4(d) objects and lapply().")
})

##' multiPhylo4 and extended classes
##' 
##' Classes for lists of phylogenetic trees.  These classes and methods are
##' planned for a future version of \code{phylobase}.
##' 
##' 
##' @name multiPhylo-class
##' @aliases multiPhylo-class multiPhylo4-class multiPhylo4d-class tbind
##' @docType class
##' @keywords classes
setAs("multiPhylo", "multiPhylo4", function(from, to) {
    trNm <- names(from)
    if(is.null(trNm)) trNm <- character(0)
    newobj <- new("multiPhylo4", phylolist = lapply(from, function(x)
                                 as(x, "phylo4")),
                  tree.names = trNm)
    newobj
})


setAs("multiPhylo4", "multiPhylo", function(from, to) {
    y <- lapply(from@phylolist, function(x) as(x, "phylo"))
    names(y) <- from@tree.names
    if (hasTipData(from))
        warning("discarded tip data")
    class(y) <- "multiPhylo"
    y
})
