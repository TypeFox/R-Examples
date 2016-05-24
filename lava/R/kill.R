##' Generic method for removing elements of object
##'
##' @title Remove variables from (model) object.
##' @aliases kill kill<- rmvar rmvar<-
##' @param x Model object
##' @param value Vector of variables or formula specifying which nodes to
##' remove
##' @param \dots additional arguments to lower level functions
##' @usage
##' kill(x, ...) <- value
##' @seealso \code{cancel}
##' @author Klaus K. Holst
##' @keywords models regression
##' @export
##' @examples
##'
##' m <- lvm()
##' addvar(m) <- ~y1+y2+x
##' covariance(m) <- y1~y2
##' regression(m) <- c(y1,y2) ~ x
##' ### Cancel the covariance between the residuals of y1 and y2
##' cancel(m) <- y1~y2
##' ### Remove y2 from the model
##' rmvar(m) <- ~y2
##'
"kill" <- function(x, ...) UseMethod("kill")

##' @export
"rmvar" <- function(x, ...) UseMethod("rmvar")

##' @export
"kill<-" <- function(x, ..., value) UseMethod("kill<-")
##' @export
"rmvar<-" <- function(x, ..., value) UseMethod("rmvar<-")


##' @export
"kill<-.lvm" <- function(x, ..., value) {
  kill(x,value)
}

##' @export
"rmvar<-.lvm" <- get("kill<-.lvm")

##' @export
"kill.lvm" <- function(x, value, ...) {
  if (inherits(value,"formula")) value <- all.vars(value)
  idx <- which(names(x$exfix)%in%value)
  if (length(idx)>0) {
      x$attributes$parameter[idx] <- x$expar[idx] <- x$exfix[idx] <- NULL
      if (length(x$exfix)==0) {
          x$exfix <- x$expar <- x$attributes$parameter <- NULL
      }
      index(x) <- reindex(x)
  }
  idx <- which(vars(x)%in%value)
  if (length(idx)==0)
    return(x)
  vv <- vars(x)[idx]
  keep <- setdiff(seq_along(vars(x)),idx)
  x$M <- x$M[keep,keep,drop=FALSE]
  x$par <- x$par[keep,keep,drop=FALSE]
  x$fix <- x$fix[keep,keep,drop=FALSE]
  x$covpar <- x$covpar[keep,keep,drop=FALSE]
  x$covfix <- x$covfix[keep,keep,drop=FALSE]
  x$cov <- x$cov[keep,keep,drop=FALSE]
  x$mean <- (x$mean)[-idx]
  x$exogenous <- setdiff(exogenous(x),vv)
  x$latent[vv] <- NULL
  index(x) <- reindex(x)
  return(x)
}

##' @export
"rmvar.lvm" <- get("kill.lvm")
