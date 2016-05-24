#' Fast squared summing in different bins
#' 
#' \code{binned_sum} implements fast squared summing for given bins by calling c-code,
#' which can be used to calculate variance and standard deviation
#' Please note that incorrect use of this function may crash your R-session.
#' the values of \code{bins} must be in between 1:\code{nbins} and \code{bin} may not 
#' contain \code{NA}
#' @useDynLib ffbase
#' @param x \code{numeric} vector with the data to be summed squared
#' @param mean \code{numeric} vector with an optional mean to be subtracted from the data to be summed and squared
#' @param bin \code{integer} vector with the bin number for each observation
#' @param nbins \code{integer} maximum bin number 
#' @param ... will be passed on to the implementation. 
#' @return \code{numeric} matrix where each row is a bin
#' @export
binned_sumsq <- function (x, mean=rep(0, nbins), bin, nbins=max(bin), ...){
  UseMethod("binned_sumsq")
}

#' @return \code{numeric} matrix where each row is a bin
#' @rdname binned_sumsq
#' @method binned_sumsq default
#' @export
#' @export binned_sumsq.default
binned_sumsq.default <- function (x, mean=rep(0, nbins), bin, nbins=max(bin), ...){
   stopifnot(length(x)==length(bin))
   if (is.factor(bin)){
     bins <- levels(bin)
     nbins <- length(bins)
   } else {
     bins <- seq_len(nbins)
   }
  stopifnot(nbins==length(mean))
  res <- matrix(0, nrow=nbins, ncol=3, dimnames=list(bin=bins, c("count", "sumsq", "<NA>")))
   .Call("binned_sumsq", as.numeric(x), as.numeric(mean), as.integer(bin), as.integer(nbins), res, PACKAGE = "ffbase")
   res
}

#' @return \code{numeric} matrix where each row is a bin
#' @rdname binned_sumsq
#' @method binned_sumsq ff
#' @export
#' @export binned_sumsq.ff
binned_sumsq.ff <- function(x, mean=rep(0, nbins), bin, nbins=max(bin), ...){
  INDEX <- list(...)$INDEX
  if (!is.null(INDEX)){
    bins <- seq_len(nbins)
    res <- matrix(0, nrow=nbins, ncol=3, dimnames=list(bin=bins, c("count", "sumsq", "<NA>")))
    for (i in chunk(INDEX, ...)){
      Log$chunk(i)
      bin <- seq.int(i[1], i[2]) / ((length(INDEX)+1)/nbins) + 1
      .Call("binned_sumsq", as.numeric(x[INDEX[i]]), as.numeric(mean), as.integer(bin), as.integer(nbins), res, PACKAGE = "ffbase")
    }
    return(res)
  }
  
  if (is.factor.ff(bin)){
    bins <- levels(bin)
    nbins <- length(bins)
  } else {
    bins <- seq_len(nbins)
  }
  res <- matrix(0, nrow=nbins, ncol=3, dimnames=list(bin=bins, c("count", "sumsq","<NA>")))
  for (i in chunk(x, ...)){
    Log$chunk(i)
    .Call("binned_sumsq", as.numeric(x[i]), as.numeric(mean), as.integer(bin[i]), as.integer(nbins), res, PACKAGE = "ffbase")
  }
  res
}


##### quick testing code ######
# x <- as.numeric(1:100000)
# bin <- as.integer(runif(length(x), 1, 101))
# x[1] <- NA
# 
# binned_sumsq(1:10, rep(1, 10), 1:10, nbins=10)
# binned_sumsq(c(1000,NA), 1:2, 1:2, nbins=2L)
