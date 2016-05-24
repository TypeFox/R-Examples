##' Generic method for extracting children or parent elements of object (e.g. a graph)
##'
##' @title Extract children or parent elements of object
##' @export
##' @aliases children parents
##' @param object Object
##' @param \dots Additional arguments
##' @author Klaus K. Holst
"children" <- function(object,...) UseMethod("children")
##' @export
"parents" <- function(object,...) UseMethod("parents")

##' @export
parents.lvmfit <- function(object,...) parents(Model(object),...)

##' @export
children.lvmfit <- function(object,...) children(Model(object),...)

##' @export
parents.lvm <- function(object,var,...) {
  A <- index(object)$A
  if (missing(var)) {
    return(rownames(A))
  }
  if (inherits(var,"formula"))
    var <- all.vars(var)
  res <- lapply(var, function(v) rownames(A)[A[,v]!=0])
  unique(unlist(res))
}

##' @export
children.lvm <- function(object,var,...) {
  A <- index(object)$A
  if (missing(var)) {
    return(rownames(A))
  }
  if (inherits(var,"formula"))
    var <- all.vars(var)
  res <- lapply(var, function(v) rownames(A)[A[v,]!=0])
  unique(unlist(res))

}
