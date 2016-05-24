##' Extract model
##'
##' Extract or replace model object
##'
##'
##' @aliases Model Model<-
##' @usage
##'
##' Model(x, ...)
##'
##' Model(x, ...) <- value
##'
##' @param x Fitted model
##' @param value New model object (e.g. \code{lvm} or \code{multigroup})
##' @param \dots Additional arguments to be passed to the low level functions
##' @return Returns a model object (e.g. \code{lvm} or \code{multigroup})
##' @author Klaus K. Holst
##' @seealso \code{\link{Graph}}
##' @keywords models
##' @examples
##'
##' m <- lvm(y~x)
##' e <- estimate(m, sim(m,100))
##' Model(e)
##'
##' @export
`Model` <- function(x,...) UseMethod("Model")


##' @export
`Model.default` <- function(x,...) x

##' @export
`Model.lvm` <- function(x,...) x

##' @export
`Model.lvmfit` <- function(x,...) x$model

##' @export
`Model.multigroup` <- function(x,...) x$lvm

##' @export
`Model.multigroupfit` <- function(x,...) x$model

##' @export
"Model<-" <- function(x,...,value) UseMethod("Model<-")

##' @export
"Model<-.lvm" <- function(x,...,value) { x <- value; return(x) }
##' @export
"Model<-.lvmfit" <- function(x,...,value) { x$model <- value; return(x) }
##' @export
"Model<-.multigroup" <- function(x,...,value) { x$lvm <- value; return(x) }
##' @export
"Model<-.multigroupfit" <- function(x,...,value) { x$model <- value; return(x) }
