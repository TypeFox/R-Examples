
##' Adding data to a phylo4 or a phylo4d object
##'
##' \code{addData} adds data to a \code{phylo4} (converting it in a
##' \code{phylo4d} object) or to a \code{phylo4d} object
##'
##' Rules for matching data to tree nodes are identical to those used
##' by the \code{\link{phylo4d-methods}} constructor.
##'
##' If any column names in the original data are the same as columns
##' in the new data, ".old" is appended to the former column names and
##' ".new" is appended to the new column names.
##'
##' The option \code{pos} is ignored (silently) if \code{x} is a
##' \code{phylo4} object. It is provided for compatibility reasons.
##'
##' @param x a phylo4 or a phylo4d object
##' @param tip.data a data frame (or object to be coerced to one)
##' containing only tip data
##' @param node.data a data frame (or object to be coerced to one)
##' containing only node data
##' @param all.data a data frame (or object to be coerced to one)
##' containing both tip and node data
##' @param merge.data if both \code{tip.data} and \code{node.data} are
##' provided, it determines whether columns with common names will be
##' merged together (default TRUE). If FALSE, columns with common
##' names will be preserved separately, with ".tip" and ".node"
##' appended to the names. This argument has no effect if
##' \code{tip.data} and \code{node.data} have no column names in
##' common.
##' @param pos should the new data provided be bound \code{before} or
##' \code{after} the pre-existing data?
##' @param \dots additional arguments to control how matching between
##' data and tree (see Details section of
##' \code{\link{phylo4d-methods}} for more details).
##' @return \code{addData} returns a \code{phylo4d} object.
##' @author Francois Michonneau
##' @seealso \code{\link{tdata}} for extracting or updating data and
##' \code{\link{phylo4d-methods}} constructor.
##' @keywords methods
##' @rdname addData-methods
##' @include phylo4d-class.R
##' @export
##' @examples
##'   data(geospiza)
##'   nDt <- data.frame(a=rnorm(nNodes(geospiza)), b=1:nNodes(geospiza),
##'                     row.names=nodeId(geospiza, "internal"))
##'   t1 <- addData(geospiza, node.data=nDt)
setGeneric("addData", function(x, ...) {
    standardGeneric("addData")
})

##' @rdname addData-methods
##' @aliases addData-methods addData,phylo4-method
setMethod("addData", signature(x="phylo4d"),
  function(x, tip.data=NULL, node.data=NULL, all.data=NULL,
           merge.data=TRUE, pos=c("after", "before"), ...) {

    pos <- match.arg(pos)

    ## apply formatData to ensure data have node number rownames and
    ## correct dimensions
    tip.data <- formatData(phy=x, dt=tip.data, type="tip", ...)
    node.data <- formatData(phy=x, dt=node.data, type="internal", ...)
    all.data <- formatData(phy=x, dt=all.data, type="all", ...)
    ## combine data as needed
    new.data <- .phylo4Data(x=x, tip.data=tip.data, node.data=node.data,
        all.data=all.data, merge.data=merge.data)

    if (all(dim(new.data) == 0)) {
        return(x)
    }
    if (all(dim(x@data) == 0)) {
        x@data <- new.data
        return(x)
    }

    if (identical(pos, "after")) {
        new.data <- merge(x@data, new.data, by=0, all=TRUE,
            sort=FALSE, suffixes=c(".old", ".new"))
    } else {
        new.data <- merge(new.data, x@data, by=0, all=TRUE,
            sort=FALSE, suffixes=c(".new", ".old"))
    }
    row.names(new.data) <- new.data[["Row.names"]]
    x@data <- new.data[, -match("Row.names", names(new.data)), drop = FALSE]

    x
})

##' @rdname addData-methods
##' @aliases addData,phylo4d-method
setMethod("addData", signature(x="phylo4"),
  function(x, tip.data=NULL, node.data=NULL, all.data=NULL,
           merge.data=TRUE, pos=c("after", "before"), ...) {
    phylo4d(x, tip.data=tip.data, node.data=node.data, all.data=all.data,
            merge.data=merge.data, ...)
})
