#' Fast summing in different bins
#' 
#' \code{binned_sum} implements fast summing for given bins by calling c-code.
#' @useDynLib ffbase
#' @param x \code{numeric} vector with the data to be summed
#' @param bin \code{integer} vector with the bin number for each data point
#' @param nbins \code{integer} maximum bin number 
#' @param ... used by binned_sum.ff
#' @return \code{numeric} matrix where each row is a bin
#' @export
binned_sum <- function (x, bin, nbins=max(bin), ...){
  UseMethod("binned_sum")
}

#' @rdname binned_sum
#' @method binned_sum default
#' @export
#' @export binned_sum.default
binned_sum.default <- function (x, bin, nbins=max(bin), ...){
   stopifnot(length(x)==length(bin))
   if (is.factor(bin)){
     bins <- levels(bin)
     nbins <- length(bins)
   } else {
     bins <- seq_len(nbins)
   }
   res <- matrix(0, ncol=3, nrow=nbins, dimnames=list(bin=bins, c("count", "sum", "NA")))
   .Call("binned_sum", as.numeric(x), as.integer(bin), as.integer(nbins), res, PACKAGE = "ffbase")
   res
}

#' @rdname binned_sum
#' @export
#' @export binned_sum.ff
binned_sum.ff <- function(x, bin, nbins=max(bin), ...){
  INDEX <- list(...)$INDEX
  if (!is.null(INDEX)){
    bins <- seq_len(nbins)
    res <- matrix(0, nrow=nbins, ncol=3, dimnames=list(bin=bins, c("count", "sum", "NA")))
    for (i in chunk(INDEX, ...)){
      Log$chunk(i)
      bin <- seq.int(i[1], i[2]) / ((length(INDEX)+1)/nbins) + 1
      .Call("binned_sum", as.numeric(x[INDEX[i]]), as.integer(bin), as.integer(nbins), res, PACKAGE = "ffbase")
    }
    return(res)
  }
  
  stopifnot(length(x)==length(bin))
  if (is.factor.ff(bin)){
    bins <- levels(bin)
    nbins <- length(bins)
  } else {
    bins <- seq_len(nbins)
  }
  res <- matrix(0, nrow=nbins, ncol=3, dimnames=list(bin=bins, c("count", "sum", "NA")))
  for (i in chunk(x, ...)){
    Log$chunk(i)
    .Call("binned_sum", as.numeric(x[i]), as.integer(bin[i]), as.integer(nbins), res, PACKAGE = "ffbase")
  }
  res
}

##### quick testing code ######
# x <- as.numeric(1:100000)
# bin <- as.integer(runif(length(x), 1, 101))
# x[1] <- NA
# 
# binned_sum(1:10, 1:10, nbins=10)
# binned_sum(c(1000,NA), factor(c("M","V")), nbins=2L)
# system.time({
#   replicate(50, {tapply(x, bin, function(i){c(sum=sum(i, na.rm=TRUE), na=sum(is.na(i)))})})
# })
#    
# system.time({
#   replicate(50, {binned_sum(x, bin, nbins=100L)})
# })
# 
# x <- ff(1:1e5)
# b <- ff(as.factor(rep(c("M","V"), 1e5/2)))
# o <- ff(20:1)
# binned_sum.ff(x, b, nbins=5, INDEX=o, by=1e2)
