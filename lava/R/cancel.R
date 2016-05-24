##' Generic cancel method
##'
##' @title Generic cancel method
##' @param x Object
##' @param \dots Additioal arguments
##' @author Klaus K. Holst
##' @aliases cancel<-
##' @export
"cancel" <- function(x,...) UseMethod("cancel")

##' @export
"cancel<-" <- function(x,...,value) UseMethod("cancel<-")

##' @export
"cancel<-.lvm" <- function(x, ..., value) {
  cancel(x,value,...)
}


##' @export
cancel.lvm <- function(x,value,...) {
  if (inherits(value,"formula")) {
      ##      yx <- all.vars(value)
    lhs <- getoutcome(value)
    if (length(lhs)==0) yy <- NULL else yy <- decomp.specials(lhs)
    xf <- attributes(terms(value))$term.labels
    if(identical(all.vars(value),xf))
      return(cancel(x,xf))
    res <- lapply(xf,decomp.specials)
    xx <- unlist(lapply(res, function(z) z[1]))
    for (i in yy) {
      for (j in xx)
        cancel(x) <- c(i,j)
      }
    index(x) <- reindex(x)
  return(x)
  }

  for (v1 in value)
    for (v2 in value)
      if (v1!=v2)
        {
          if (all(c(v1,v2)%in%vars(x))) {
            x$M[v1,v2] <- 0
            x$par[v1,v2] <- x$fix[v1,v2] <-
              x$covpar[v1,v2] <- x$covfix[v1,v2] <- NA
            x$cov[v1,v2] <- 0
          }
        }
  x$parpos <- NULL
  index(x) <- reindex(x)
  return(x)
}
