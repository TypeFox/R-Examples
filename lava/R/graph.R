##' Extract graph
##'
##' Extract or replace graph object
##'
##'
##' @aliases Graph Graph<-
##' @usage
##'
##' Graph(x, ...)
##'
##' Graph(x, ...) <- value
##'
##' @param x Model object
##' @param value New \code{graphNEL} object
##' @param \dots Additional arguments to be passed to the low level functions
##' @author Klaus K. Holst
##' @seealso \code{\link{Model}}
##' @keywords graphs models
##' @export
##' @examples
##'
##' m <- lvm(y~x)
##' Graph(m)
##'
##' @export
`Graph` <-
function(x,...) UseMethod("Graph")

##' @export
`Graph.lvm` <-
function(x,add=FALSE,...) {
  if ((is.null(x$graph) || length(x$graph)==0) & add) {
    m <- Model(x)
    return(plot(m,noplot=TRUE))
  }
  else return(x$graph)
}

##' @export
`Graph.lvmfit` <- function(x,...) Graph.lvm(x,...)

##' @export
"Graph<-" <- function(x,...,value) UseMethod("Graph<-")

##' @export
"Graph<-.lvmfit" <- function(x,...,value) { x$graph <- value; return(x) }

##' @export
"Graph<-.lvm" <- function(x,...,value) { x$graph <- value; return(x) }
