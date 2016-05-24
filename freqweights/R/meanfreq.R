## Emilio Torres Manzanera
## University of Oviedo
## Time-stamp: <2014-10-09 10:31 emilio on emilio-despacho>
## ============================================================



## ============================================================
##
## ============================================================

##' Descriptive statistics of a frequency table.
##'
##' Computes the descriptive statistics of a frequency table.
##'
##' These functions compute various weighted versions of standard estimators.
##'
##' \code{meanfreq}, \code{sdfreq}, \code{quantilefreq}, \code{covfreq},
##' \code{corfreq} estimate the mean, standard desviation, quantiles,
##' covariances and correlation matrix, respectively. In this last two cases,
##' resulst are equals to the \code{pairwise.complete.obs} option of \code{cov}
##' and \code{cor} of the desaggregated data, respectively.
##'
##' Missing values or cases with non-positive frequency weights are
##' automatically removed.
##'
##' If \code{freq} is not null, the data set must contain a column with that name. These
##' variable are removed from the data set in order to calculate the
##' descriptive statistics.
##'
##' The dot versions are intented to be used when programing. The \code{tfq} may be a \code{tablefreq} object or a matrix or a data frame with the last column being the frequency weights.
##'
##' The algorithm of \code{quantilefreq} are based on
##' \code{\link[Hmisc]{wtd.quantile}}.
##'
##' The intern functions are for programming purpose. It does not check the data.
##' @name statsfreq
##' @rdname statsfreq
##' @aliases meanfreq quantilefreq covfreq corfreq sdfreq
##' @param data any object that can be processed by \code{link{tablefreq}}.
##' @param freq a single name of the variable specifying frequency weights. 
##' @param probs A vector of quantiles to compute. Default is 0 (min), .25, .5,
##' .75, 1 (max).
##' @param tfq a \code{tablefreq} object, or a matrix, data frame with the last column being the frequency wweights
##' @return \code{meanfreq} and \code{sdfreq} return vectors.
##' \code{quantilefreq} returns a vector or matrix.
##' \code{covfreq} and \code{corfreq} the estimated covariance matrix and
##' correlation matrix, respectively. \code{scalefreq} return a data frame or matrix
##' @seealso
##' \code{\link{tablefreq}}, \code{\link[Hmisc]{wtd.quantile}}
##' @references Andrews, Chris,
##' \url{https://stat.ethz.ch/pipermail/r-help/2014-March/368350.html}
##' @note The author would like to thank Prof. Frank E. Harrell Jr. who allowed the reutilisation of part of his code.
##' @keywords univar
##' @examples
##' if(require(hflights)) {
##'   meanfreq(hflights[,c("ArrDelay","DepDelay")])
##'   sdfreq(hflights[,c("ArrDelay","DepDelay")])
##'   corfreq(hflights[,c("ArrDelay","DepDelay")])
##' }
##'
##' tfq <- tablefreq(iris$Sepal.Length)
##' tfq
##' 
##' meanfreq(iris$Sepal.Length)
##' meanfreq(tfq,freq="freq")
##' .meanfreq(tfq)
##' 
##' dat <- iris[,1:4]
##' quantilefreq(dat)
##' corfreq(dat)
##' 
##' tfq <- tablefreq(dat)
##' .meanfreq(tfq)
##' .quantilefreq(tfq)
##' .corfreq(tfq)
##'
##' ## dplyr integration
##' library(dplyr)
##' tfq  %>% 
##'   summarise( mean = .meanfreq(cbind(Sepal.Length,freq)),
##'             sd = .sdfreq(cbind(Sepal.Length,freq)))
##'
##' tfq <- tablefreq(iris)
##' tfq %>% group_by(Species) %>% 
##'   summarise( mean = .meanfreq(cbind(Sepal.Length,freq)),
##'             sd = .sdfreq(cbind(Sepal.Length,freq)))
##' @name statsfreq
##' @rdname statsfreq
NULL

##' @rdname statsfreq
##' @export
meanfreq <- function(data, freq=NULL){
  tfq <- tablefreq(data,freq=freq)
  .meanfreq(tfq)
}

##' @rdname statsfreq
##' @export
.meanfreq <- function(tfq){
  if(any(class(tfq) %in% c('matrix', 'double', 'numeric')))
     tfq <- tbl_df(as.data.frame(tfq))
  nms <- tbl_vars(tfq)
  freq <- nms[length(nms)]
  nms <- nms[-length(nms)]
  weigthedsum <- function(inms){
    argsfilter <- paste( "!is.na(", inms,"), ",freq," > 0",", !is.na(",freq,")")
    argssum    <- paste( inms," = sum( ", freq ,"* ", inms, ") / sum( ",freq," )")
    meantbl <- tfq %>% evaldp(filter, argsfilter) %>% evaldp(summarise, argssum)
    meantbl[,1]
  }
  result <- plyr::laply(nms,weigthedsum)
  if(length(result)> 1) {
    names(result) <- nms
    return(result)} else{
      return(as.numeric(result))
  }
}


##' @rdname statsfreq
##' @export
quantilefreq <-  function(data, probs = c(0, 0.25, 0.5, 0.75, 1), freq = NULL){
   tfq <- tablefreq(data, freq=freq)
  .quantilefreq(tfq, probs)
}


##' @rdname statsfreq
##' @export
.quantilefreq <- function(tfq, probs = c(0, 0.25, 0.5, 0.75, 1)) {
   if (any(probs < 0 | probs > 1))
    stop("quantilefreq: Probabilities must be between 0 and 1 inclusive")
  nams <- paste(format(round(probs * 100,
                             if (length(probs) > 1)
                             2 - log10(diff(range(probs)))
                             else 2)),
                "%", sep = "")
  result <- laply(seq_len(ncol(tfq)-1),function(i){
    ok <- unlist(!is.na(tfq[,i]) & tfq[,ncol(tfq)]>0)
    x <- unlist(tfq[ok,i])
    w <- unlist(tfq[ok,ncol(tfq)])
    n <- sum(w)
    order <- 1 + (n - 1) * probs
    low <- pmax(floor(order), 1)
    high <- pmin(low + 1, n)
    order <- order%%1
    allq <- approx(cumsum(w), x,
                   xout = c(low, high), method = "constant",
                   f = 1, rule = 2)$y
    k <- length(probs)
    quantiles <- (1 - order) * allq[1:k] + order * allq[-(1:k)]
    names(quantiles) <- nams
    return(quantiles)
  })
  rownames(result) <- colnames(tfq)[-ncol(tfq)]
  result
}

## ============================================================
##
## ============================================================


covfreqcompletecases <- function(x){
  x <- as.matrix(x)
  if (!all(is.finite(x)))
    stop("'x' must contain finite values only")
  w <- x[,ncol(x)]
  x <- x[,-ncol(x),drop=FALSE]
  n <- sum(w)
  center <-  colSums(w * x) / n
  xcw <- sqrt(w) * sweep(x, 2, center, check.margin = TRUE)
  r <- crossprod(xcw)
  if( n == 1) {
    r <- 0
  } else {
    r <- r/(n-1)
  }
  rownames(r) <- colnames(x)
  colnames(r) <- colnames(x)
  return(r)
}





##' @rdname statsfreq
##' @export
covfreq <- function(data, freq = NULL) {
  tfq <- tablefreq(data, freq=freq)
  ##print(paste("covfreq a ",dim(x)))
  .covfreq(tfq)
}


##' @rdname statsfreq
##' @export
.covfreq <- function(tfq){
  if( all(complete.cases(tfq)) ) {
    return(covfreqcompletecases(tfq))
  } else  {
    ##print(dim(x))
    ncy <- ncx <- ncol(tfq) - 1
    r <- matrix(0, nrow = ncx, ncol = ncy)
    for (i in seq_len(ncx)) {
      for (j in seq_len(i)) {
        cols <-  c(i,j,ncol(tfq))
        ok <- complete.cases(tfq[,cols])
        covb <- covfreqcompletecases(tfq[ok,cols, drop=FALSE])
        r[i,j] <- r[j,i] <- covb[1,2]
      }
    }
    rownames(r) <- colnames(tfq)[-ncol(tfq)]
    colnames(r) <- colnames(tfq)[-ncol(tfq)]
    return(r)
  }
}


##' @rdname statsfreq
##' @export
sdfreq <- function(data, freq=NULL){
  tfq <- tablefreq(data, freq=freq)
  .sdfreq(tfq)
}


##' @rdname statsfreq
##' @export
.sdfreq <- function(tfq){
  if(any(class(tfq) %in% c('matrix', 'double', 'numeric')))
    tfq <- tbl_df(as.data.frame(tfq))
  nms <- tbl_vars(tfq)
  freq <- nms[length(nms)]
  nms <- nms[-length(nms)]
  weigthedsd <- function(inms){
    argsfilter <- paste( "!is.na(", inms,"), ",freq," > 0",", !is.na(",freq,")")
    argssum    <- paste0( inms," = sum( ", freq ," * ", inms, ") / sum( ",freq," )")
    tblinms <- tfq %>% evaldp(filter, argsfilter)
    meantbl <-  tblinms %>% evaldp(summarise, argssum)
    meanvar <- meantbl[,1]
    argssd <- paste0( inms," = (sum( ", freq ," * (", inms," - ", meanvar, ")^2) / (sum( ",freq," ) -1 ))^(0.5)")
    result <- tblinms %>% evaldp(summarise, argssd)
    result[,1]
  }
  result <- plyr::laply(nms,weigthedsd)
  if(length(result)> 1) {
    names(result) <- nms
    return(result)} else{
      return(as.numeric(result))
  }
  ## result <- laply(seq_len(ncol(tfq)-1),function(i){
  ##   ok <- !is.na(tfq[,i]) & tfq[, ncol(tfq)] > 0
  ##   x <- tfq[ok, i]
  ##   w <- tfq[ok, ncol(tfq)]
  ##   xbar <- sum( x * w) / sum(w)
  ##   if( sum(w) != 1 ) {
  ##     return(sum( w *(x-xbar)^2)/(sum(w)-1))
  ##   } else {
  ##     return(0)
  ##   }
  ## })
  ## if(length(result)> 1) {
  ##    return(result)} else{
  ##     return(as.numeric(result))
  ## }
}


##' @rdname statsfreq
##' @export
scalefreq <- function(data,freq=NULL){
  tfq <- tablefreq(data, freq=freq)
  .scalefreq(tfq)
}


##' @rdname statsfreq
##' @export
.scalefreq <- function(tfq){
  center <- .meanfreq(tfq)
  sd <- .sdfreq(tfq)
  x <- sweep(x=tfq[, -ncol(tfq),drop=FALSE], MARGIN=2L, STATS=center,FUN="-", check.margin=FALSE)
  x <- sweep(x=x, MARGIN=2L, STATS=sd,FUN="/", check.margin=FALSE)
  newtbl <- cbind(x, tfq[,ncol(tfq)])
  colnames(newtbl) <- colnames(tfq)
  attr(newtbl,"freq") <- attr(tfq, "freq")
  attr(newtbl,"colweights") <- attr(tfq,"colweights")
  newtbl
}


##' @rdname statsfreq
##' @export
corfreq <- function(data, freq = NULL){
  tfq <- tablefreq(data, freq=freq)
  .corfreq(tfq)
}


##' @rdname statsfreq
##' @export
.corfreq <- function(tfq){
  ncy <- ncx <- ncol(tfq) -1
  r <- matrix(0, nrow = ncx, ncol = ncy)
  for (i in seq_len(ncx)) {
    for (j in seq_len(i)) {
      cols <-  c(i,j,ncol(tfq))
      ok <- complete.cases(tfq[,cols])
      covb <- covfreqcompletecases(tfq[ok,cols, drop=FALSE])
      r[i,j] <- r[j,i] <- covb[1,2] / sqrt(covb[1,1] * covb[2,2])
    }
  }
  rownames(r) <- colnames(tfq)[-ncol(tfq)]
  colnames(r) <- colnames(tfq)[-ncol(tfq)]
  return(r)
}
