##' Defines variables as binary
##' 
##' Defines variables of a latent variable model (\code{lvm}-object from the
##' \code{lava} package) as binary.
##' 
##' This function sets the status of variables to 'binary'. For use with the
##' \code{estimate} method a probit-link will be assumed. Used with the
##' \code{sim} method normal distributed data will be simulated followed by
##' thresholding at 0. To simulate data where the dichotomous variable has a
##' direct effect on the outcome the \code{distribution} method can be used, e.g.
##' \code{distribution(m,~x) <- binomial.lvm("probit")}.
##' 
##' @aliases binary binary<- binary.lvm binary<-.lvm
##' @param x \code{lvm}-object
##' @param var Formula or vector of variable names
##' @param value Formula or vector of variable names
##' @param \dots Additional arguments parsed to lower-level functions
##' @usage
##' \method{binary}{lvm}(x,var=NULL, ...)
##' \method{binary}{lvm}(x, ...) <- value
##' @return \code{lvm}-object (or vector of variable names if called without any
##' arguments)
##' @author Klaus K. Holst
##' @keywords models regression
##' @export
`binary` <- function(x,...) UseMethod("binary")

##' @export
"binary<-" <- function(x,...,value) UseMethod("binary<-")

##' @S3method binary<- lvm
"binary<-.lvm" <- function(x,...,value) {
  if (class(value)[1]=="formula") {
    return(binary(x,all.vars(value),...))
  }
  binary(x, value, ...)
}

##' @S3method binary lvm
`binary.lvm` <-
function(x,var=NULL, ...) {
  if (is.null(var)) {
    ## binidx <- tryCatch(unlist(nodeData(Graph(x), attr="binary")),error=function(e) NULL)
    binidx <- unlist(x$attributes$binary)
    if (length(binidx)>0)
      return(names(binidx)[binidx])
    else
      return(NULL)
  }
  ##  if (is.null(nodeDataDefaults(Graph(x))$binary)) {
  ##    nodeDataDefaults(Graph(x),"binary") <- FALSE
  ##  } 
  
  ##  x <- addattr(x,attr="shape",var=var,val="box")
  x$attributes$binary[var] <- TRUE
  x$attributes$normal[var] <- FALSE
  ## nodeData(Graph(x), var, attr="binary") <- TRUE
  ## nodeData(Graph(x), var, attr="normal") <- FALSE
  covfix(x,var,NULL) <- 1
  ##  distribution(x, var) <- probit.lvm
  return(x)
}

