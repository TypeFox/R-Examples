
##' Methods to test, access (and modify) the root of a phylo4 object.
##'
##' @rdname root-methods
##' @aliases isRooted
##' @docType methods
##' @param x a \code{phylo4} or \code{phylo4d} object.
##' @param value a character string or a numeric giving the new root.
##' @return \describe{
##'   \item{isRooted}{logical whether the tree is rooted}
##'   \item{rootNode}{the node corresponding to the root}
##' }
##' @include phylo4-class.R phylo4-methods.R phylo4-accessors.R
##' @export
##' @author Ben Bolker, Francois Michonneau
##' @examples
##' data(geospiza)
##' isRooted(geospiza)
##' rootNode(geospiza)
setGeneric("isRooted", function(x) {
    standardGeneric("isRooted")
})

##' @rdname root-methods
##' @aliases isRooted,phylo4-method
setMethod("isRooted", signature(x="phylo4"),
 function(x) {
    ## hack to avoid failure on an empty object
    if(nTips(x) == 0) return(FALSE)
    any(edges(x)[, 1] == 0)
})

##' @rdname root-methods
##' @aliases rootNode
##' @export
setGeneric("rootNode", function(x) {
    standardGeneric("rootNode")
})

##' @rdname root-methods
##' @aliases rootNode,phylo4-method
setMethod("rootNode", signature(x="phylo4"),
 function(x) {
    if (!isRooted(x))
        return(NA)
    rootnd <- unname(edges(x)[which(edges(x)[, 1] == 0), 2])
    getNode(x, rootnd)
})

##' @rdname root-methods
##' @aliases rootNode<-
##' @export
setGeneric("rootNode<-", function(x, value) {
    standardGeneric("rootNode<-")
})

##' @name rootNode<-
##' @rdname root-methods
##' @aliases rootNode<-,phylo4-method
setReplaceMethod("rootNode", signature(x="phylo4"),
 function(x, value) {
    stop("Root node replacement not implemented yet")
})

