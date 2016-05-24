
##' Tests for presence of data associated with trees stored as phylo4d objects
##' 
##' Methods that test for the presence of data associated with trees stored as
##' \code{phylo4d} objects.
##' 
##' \code{nData} tests for the presence of data associated with the object.
##'
##' \code{hasTipData} and \code{hasNodeData} tests for the presence of
##' data associated with the tips and the internal nodes
##' respectively. The outcome of the test is based on row names of the
##' data frame stored in the \code{data} slot. If no rows have names
##' from the set \code{nodeId(x, "tip")}, then \code{hasTipData}
##' returns FALSE.  Likewise, if no rows have names from the set
##' \code{nodeId(x, "internal")}, then \code{hasNodeData} returns
##' FALSE.
##' 
##' @param x a \code{phylo4d} object
##' @return \describe{
##' 
##'  \item{\code{nData}}{returns the number of datasets (i.e.,
##' columns) associated with the object.}
##' 
##'  \item{\code{hasTipData}, \code{hasNodeData}}{return \code{TRUE}
##' or \code{FALSE} depending whether data associated with the
##' tree are associated with either tips or internal nodes respectively.}}
##' @section Methods: \describe{ \item{hasNodeData}{\code{signature(object =
##' "phylo4d")}: whether tree has internal node data}
##' \item{hasTipData}{\code{signature(object = "phylo4d")}: whether tree has
##' data associated with its tips} }
##' @author Ben Bolker, Thibault Jombart, Francois Michonneau
##' @seealso \code{\link{phylo4d-methods}} constructor and
##' \code{\linkS4class{phylo4d}} class.
##' @rdname phylo4d-accessors
##' @aliases hasTipData
##' @keywords methods
##' @docType methods
##' @include phylo4d-class.R phylo4d-methods.R
##' @export
##' @examples
##'   data(geospiza)
##'   nData(geospiza)       ## 5
##'   hasTipData(geospiza)  ## TRUE
##'   hasNodeData(geospiza) ## FALSE
##'
setGeneric("hasTipData", function(x) {
    standardGeneric("hasTipData")
})

##' @rdname phylo4d-accessors
##' @aliases hasTipData-method,phylo4d-method
setMethod("hasTipData", signature(x="phylo4d"),
 function(x) {
    ncol(tdata(x, type="tip", empty.columns=FALSE)) > 0
})

##' @rdname phylo4d-accessors
##' @aliases hasNodeData-methods
##' @export
setGeneric("hasNodeData", function(x) {
    standardGeneric("hasNodeData")
})

##' @rdname phylo4d-accessors
##' @aliases hasNodeData,phylo4d-method
setMethod("hasNodeData", signature(x="phylo4d"),
 function(x) {
    ncol(tdata(x, type="internal", empty.columns=FALSE)) > 0
})

##' @rdname phylo4d-accessors
##' @aliases nData
##' @export
setGeneric("nData", function(x) {
     standardGeneric("nData")
})

##' @rdname phylo4d-accessors
##' @aliases nData,phylo4d-method
setMethod("nData", signature(x="phylo4d"), function(x) {
    ncol(x@data)
})
