# some collapsing functions for ff
# (c) 2007 Jens Oehlschägel
# Licence: GPL2
# Provided 'as is', use at your own risk
# Created: 2007-10-09
# Last changed: 2007-10-25

# source("d:/mwp/eanalysis/ff/R/CFUN.R")

# xx TODO: extend this for weighted means, weighted median etc., see http://tolstoy.newcastle.edu.au/R/help/02a/1073.html and http://tolstoy.newcastle.edu.au/R/help/02a/1060.html or google "Re: [R] Weighted median"

#! \name{CFUN}
#! \alias{CFUN}
#! \alias{ccbind}
#! \alias{crbind}
#! \alias{cfun}
#! \alias{cquantile}
#! \alias{csummary}
#! \alias{cmedian}
#! \alias{clength}
#! \alias{csum}
#! \alias{cmean}
#! \title{ Collapsing functions for batch processing }
#! \description{
#!   These are used in aggregating the chunks resulting from batch processing. They are usually called via \command{\link{do.call}}
#! }
#! \usage{
#! ccbind(\dots)
#! crbind(\dots)
#! cfun(\dots, FUN, FUNARGS = list())
#! cquantile(\dots, probs = seq(0, 1, 0.25), na.rm = FALSE, names = TRUE, type = 7)
#! csummary(\dots, na.rm = "ignored")
#! cmedian(\dots, na.rm = FALSE)
#! clength(\dots, na.rm = FALSE)
#! csum(\dots, na.rm = FALSE)
#! cmean(\dots, na.rm = FALSE)
#! }
#! \arguments{
#!   \item{\dots}{ \code{\dots} }
#!   \item{FUN}{ a aggregating function }
#!   \item{FUNARGS}{ further arguments to the aggregating function }
#!   \item{na.rm}{ TRUE to remove NAs }
#!   \item{probs}{ see \code{\link{quantile}} }
#!   \item{names}{ see \code{\link{quantile}} }
#!   \item{type}{ see \code{\link{quantile}} }
#! }
#! \details{
#!  \tabular{lll}{
#!   \strong{CFUN}        \tab \strong{FUN}              \tab \strong{comment} \cr
#!   \command{ccbind}     \tab \command{\link{cbind}}    \tab like \command{cbind} but respecting names \cr
#!   \command{crbind}     \tab \command{\link{rbind}}    \tab like \command{rbind} but respecting names \cr
#!   \command{cfun}       \tab                           \tab \command{crbind} the input chunks and then apply 'FUN' to each column \cr
#!   \command{cquantile}  \tab \command{\link{quantile}} \tab \command{crbind} the input chunks and then apply 'quantile' to each column \cr
#!   \command{csummary}   \tab \command{\link{summary}}  \tab \command{crbind} the input chunks and then apply 'summary' to each column \cr
#!   \command{cmedian}    \tab \command{\link{median}}   \tab \command{crbind} the input chunks and then apply 'median' to each column \cr
#!   \command{clength}    \tab \command{\link{length}}   \tab \command{crbind} the input chunks and then determine the number of values in each column\cr
#!   \command{csum}       \tab \command{\link{sum}}      \tab \command{crbind} the input chunks and then determine the sum values in each column\cr
#!   \command{cmean}      \tab \command{\link{mean}}     \tab \command{crbind} the input chunks and then determine the (unweighted) mean in each column\cr
#!   }
#!   In order to use CFUNs on the result of \code{\link{lapply}} or \code{\link{ffapply}} use \code{\link{do.call}}.
#! }
#! \note{
#!  Currently - for command line convenience - we map the elements of a single list argument to \dots, but this may change in the future.
#! }
#! \section{ff options}{
#!   xx TODO: extend this for weighted means, weighted median etc., \cr
#!   see http://tolstoy.newcastle.edu.au/R/help/02a/1073.html and http://tolstoy.newcastle.edu.au/R/help/02a/1060.html or google "Re: [R] Weighted median"
#! }
#! \value{
#!   depends on the CFUN used
#! }
#! \author{ Jens Oehlschlägel }
#! \seealso{ \code{\link{ffapply}}, \code{\link{do.call}}, \code{\link{na.count}} }
#! \examples{
#!    X <- lapply(split(rnorm(1000), 1:10), summary)
#!    do.call("crbind", X)
#!    do.call("csummary", X)
#!    do.call("cmean", X)
#!    do.call("cfun", c(X, list(FUN=mean, FUNARGS=list(na.rm=TRUE))))
#!    rm(X)
#! }
#! \keyword{ manip }
#! \keyword{ list }

ccbind <- function(...){
  l <- list(...); if (length(l)==1 && is.list(l[[1]])) l <- l[[1]]
  nl <- length(l)
  if (nl){
    checkit <- sapply(l, function(x)if (is.matrix(x)) c(n=nrow(x), nonam=is.null(rownames(x))) else c(n=length(x), nonam=is.null(names(x))) )
    need.names <- length(unique(checkit["n",])) != 1
    if (need.names){
      if (any(checkit["nonam",]))
        stop("either need equal lengths or names")
      nam <- unique(do.call("c", lapply(l, function(x)if(is.matrix(x)) rownames(x) else names(x))))
      mat <- do.call("cbind", lapply(l, function(x)if (is.matrix(x)) {
        m <- matrix(NA, length(nam), ncol(x), dimnames=list(nam, colnames(x)))
        m[rownames(x),] <- x
        m
      } else {
        v <- rep(NA, length(nam))
        names(v) <- nam
        v[names(x)] <- x
        v
      }))
      mat
    }else{
      do.call("cbind", l)
    }
  }else{
    return(NULL)
  }
}

crbind <- function(...){
  l <- list(...); if (length(l)==1 && is.list(l[[1]])) l <- l[[1]]
  nl <- length(l)
  if (nl){
    checkit <- sapply(l, function(x)if (is.matrix(x)) c(n=ncol(x), nonam=is.null(colnames(x))) else c(n=length(x), nonam=is.null(names(x))) )
    need.names <- length(unique(checkit["n",])) != 1
    if (need.names){
      if (any(checkit["nonam",]))
        stop("either need equal lengths or names")
      nam <- unique(do.call("c", lapply(l, function(x)if(is.matrix(x)) colnames(x) else names(x))))
      mat <- do.call("rbind", lapply(l, function(x)if (is.matrix(x)) {
        m <- matrix(NA, nrow(x), length(nam), dimnames=list(rownames(x), nam))
        m[,colnames(x)] <- x
        m
      } else {
        v <- rep(NA, length(nam))
        names(v) <- nam
        v[names(x)] <- x
        v
      }))
      mat
    }else{
      do.call("rbind", l)
    }
  }else{
    return(NULL)
  }
}



cfun <- function(..., FUN, FUNARGS=list()){
  l <- list(...); if (length(l)==1 && is.list(l[[1]])) l <- l[[1]]
  nl <- length(l)
  if (nl){
    l <- do.call("crbind", l)
    ret <- do.call("apply", c(list(l, 2, FUN), FUNARGS))
    if (is.list(ret)){
      do.call("ccbind", ret)
    }else{
      ret
    }
  }else{
    return(NULL)
  }
}


cmedian <- function(..., na.rm=FALSE){
  cfun(..., FUN=median, FUNARGS=list(na.rm=na.rm))
}

csummary <- function(..., na.rm="ignored"){
  cfun(..., FUN=summary)
}

cquantile <- function(..., probs = seq(0, 1, 0.25), na.rm = FALSE, names = TRUE, type = 7){
  cfun(..., FUN=quantile, FUNARGS=list(probs = probs, na.rm = na.rm, names = names, type = type))
}

crange <- function(..., na.rm="ignored"){
  cfun(..., FUN=range)
}



# faster implementations of clength, csum, cmean
clength <- function(..., na.rm=FALSE){
  l <- list(...); if (length(l)==1 && is.list(l[[1]])) l <- l[[1]]
  nl <- length(l)
  if (nl){
    checkit <- sapply(l, function(x)c(n=length(x), nonam=is.null(names(x))))
    need.names <- length(unique(checkit["n",])) != 1
    if (need.names){
      if (any(checkit["nonam",]))
        stop("either need equal lengths or names")
      nam <- unique(do.call("c", lapply(l, names)))
      if (na.rm){
        n <- as.integer(!is.na(l[[1]][nam]))
      }else{
        return(rep(nl, length.out=length(l[[1]][nam])))
      }
      if (nl>1){
        for (i in 2:nl){
            n <- n + !is.na(l[[i]][nam])
        }
      }
      return(n)
    }else{
      if (na.rm){
        n <- as.integer(!is.na(l[[1]]))
      }else{
        return(rep(nl, length.out=length(l[[1]])))
      }
      if (nl>1){
        for (i in 2:nl){
            n <- n + !is.na(l[[i]])
        }
      }
      return(n)
    }
  }else{
    return(NULL)
  }
}

csum <- function(..., na.rm=FALSE){
  l <- list(...); if (length(l)==1 && is.list(l[[1]])) l <- l[[1]]
  nl <- length(l)
  if (nl){
    checkit <- sapply(l, function(x)c(n=length(x), nonam=is.null(names(x))))
    need.names <- length(unique(checkit["n",])) != 1
    if (need.names){
      if (any(checkit["nonam",]))
        stop("either need equal lengths or names")
      nam <- unique(do.call("c", lapply(l, names)))
      if (na.rm)
        s <- ifelse(is.na(l[[1]][nam]), 0L, l[[1]][nam])
      else{
        s <- l[[1]][nam]
      }
      if (nl>1){
        for (i in 2:nl){
          if (na.rm)
            s <- ifelse(is.na(l[[i]][nam]), s, s + l[[1]][nam])
          else{
            s <- s + l[[i]][nam]
          }
        }
      }
      s
    }else{
      if (na.rm)
        s <- ifelse(is.na(l[[1]]), 0L, l[[1]])
      else{
        s <- l[[1]]
      }
      if (nl>1){
        for (i in 2:nl){
          if (na.rm)
            s <- ifelse(is.na(l[[i]]), s, s + l[[1]])
          else{
            s <- s + l[[i]]
          }
        }
      }
      s
    }
  }else{
    return(NULL)
  }
}

cmean <- function(..., na.rm=FALSE){
  l <- list(...); if (length(l)==1 && is.list(l[[1]])) l <- l[[1]]
  nl <- length(l)
  if (nl){
    checkit <- sapply(l, function(x)c(n=length(x), nonam=is.null(names(x))))
    need.names <- length(unique(checkit["n",])) != 1
    if (need.names){
      if (any(checkit["nonam",]))
        stop("either need equal lengths or names")
      nam <- unique(do.call("c", lapply(l, names)))
      if (na.rm){
        na <- is.na(l[[1]][nam])
        s <- ifelse(na, 0L, l[[1]][nam])
        n <- as.integer(!na)
      }else{
        s <- l[[1]][nam]
      }
      if (nl>1){
        for (i in 2:nl){
          if (na.rm){
            na <- is.na(l[[i]][nam])
            s <- ifelse(na, s, s + l[[1]][nam])
            n <- n + !na
          }else{
            s <- s + l[[i]][nam]
          }
        }
      }
      if (na.rm)
        s/n
      else
        s/nl
    }else{
      if (na.rm){
        na <- is.na(l[[1]])
        s <- ifelse(na, 0L, l[[1]])
        n <- as.integer(!na)
      }else{
        s <- l[[1]]
      }
      if (nl>1){
        for (i in 2:nl){
          if (na.rm){
            na <- is.na(l[[i]])
            s <- ifelse(na, s, s + l[[1]])
            n <- n + !na
          }else{
            s <- s + l[[i]]
          }
        }
      }
      if (na.rm)
        s/n
      else
        s/nl
    }
  }else{
    return(NULL)
  }
}
