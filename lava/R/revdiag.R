##' Create/extract 'reverse'-diagonal matrix or off-diagonal elements
##' @title Create/extract 'reverse'-diagonal matrix or off-diagonal elements
##' @aliases revdiag revdiag<- offdiag offdiag<-
##' @usage
##' revdiag(x,...)
##' offdiag(x,type=0,...)
##'
##' revdiag(x,...) <- value
##' offdiag(x,type=0,...) <- value
##' @param x vector
##' @param value For the assignment function the values to put in the diagonal
##' @param type 0: upper and lower triangular, 1: upper triangular, 2: lower triangular, 3: upper triangular + diagonal, 4: lower triangular + diagonal
##' @param \dots additional arguments to lower level functions
##' @author Klaus K. Holst
##' @export
revdiag <- function(x,...) {
    if (NCOL(x)==1) {
      res <- matrix(0,length(x),length(x))
      revdiag(res) <- x
      return(res)
    }
    n <- max(ncol(x),nrow(x))
    x[cbind(rev(seq(n)),seq(n))]
  }

##' @export
"revdiag<-" <- function(x,...,value) {
  n <- max(ncol(x),nrow(x))
  x[cbind(rev(seq(n)),seq(n))] <- value
  x
}


##' @export
offdiag <- function(x,type=0,...) {
    ##if (NCOL(x)==1) return(NULL)
    if (type%in%c(1,3)) {
        ii <- which(upper.tri(x,diag=(type==3)))
    } else if (type%in%c(2,4)) {
        ii <- which(lower.tri(x,diag=(type==4)))
    } else {
        ii <- c(which(lower.tri(x,diag=FALSE)),which(upper.tri(x,diag=FALSE)))
    }
    res <- x[ii]
    class(res) <- c("offdiag",class(res))
    attributes(res) <-
        c(attributes(res),list(type=type,dimension=dim(x),index=ii))
    return(res)
  }

##' @export
"offdiag<-" <- function(x,type=0,...,value) {
    if (type%in%c(1,3)) {
        ii <- which(upper.tri(x,diag=(type==3)))
    } else if (type%in%c(2,4)) {
        ii <- which(lower.tri(x,diag=(type==4)))
    } else {
        ii <- c(which(lower.tri(x,diag=FALSE)),which(upper.tri(x,diag=FALSE)))
    }
    x[ii] <- value
    return(x)
}

##' @export
print.offdiag <- function(x,...) {
    type <- attr(x,"type")
    nn <- attr(x,"dimension")
    M <- matrix(NA,nn[1],nn[2])
    M[attr(x,"index")] <- x
    print(M,na.print="")
}
