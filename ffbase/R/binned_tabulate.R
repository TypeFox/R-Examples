#' Fast tabulating in different bins
#' 
#' \code{binned_sum} implements fast tabulating for given bins by calling c-code. 
#' It also returns the number of NA's per bin.
#' Please note that incorrect use of this function may crash your R-session.
#' the values of \code{bins} must be between \code{1} and \code{nbins} and may not contain \code{NA}.
#' The values of \code{x} must be between \code{1} and \code{nlevels}.
#' @useDynLib ffbase
#' @param x \code{factor} or \code{integer} vector with the data to be tabulated
#' @param bin \code{integer} vector with the bin number for each data point
#' @param nbins \code{integer} maximum bin number 
#' @param nlevels \code{integer} number of levels used in x
#' @param ... used by binned_tabulate.ff
#' @return \code{numeric} matrix where each row is a bin and each column a level
#' @export
binned_tabulate <- function (x, bin, nbins=max(bin), nlevels=nlevels(x), ...){
  UseMethod("binned_tabulate")
}

#' @rdname binned_tabulate
#' @method binned_tabulate default
#' @export
#' @export binned_tabulate.default
binned_tabulate.default <- function (x, bin, nbins=max(bin), nlevels=nlevels(x), ...){
   stopifnot(length(x)==length(bin))
   if (is.factor(bin)){
     bins <- levels(bin)
     nbins <- length(bins)
   } else {
     bins <- seq_len(nbins)
   }
   lev <- if (nlevels(x)) c("na", levels(x))
          else c("na", 1:nlevels)
   res <- matrix(0L, nrow=nbins, ncol=length(lev), dimnames=list(bin=bins, level=lev))
   .Call("binned_tabulate", as.integer(x), as.integer(bin), as.integer(nbins), as.integer(nlevels), res, PACKAGE = "ffbase")
   res
}

#' @rdname binned_tabulate
#' @method binned_tabulate ff
#' @export
#' @export binned_tabulate.ff
binned_tabulate.ff <- function(x, bin, nbins=max(bin), nlevels=nlevels(x), ...){
  lev <- if (nlevels(x)) c("na", levels(x))
         else c("na", 1:nlevels)

  INDEX <- list(...)$INDEX
  if (!is.null(INDEX)){
    bins <- seq_len(nbins)
    res <- matrix(0L, nrow=nbins, ncol=length(lev), dimnames=list(bin=bins, level=lev))
    for (i in chunk(INDEX, ...)){
      Log$chunk(i)
      bin <- seq.int(i[1], i[2]) / ((length(INDEX)+1)/nbins) + 1
      .Call("binned_tabulate", as.integer(x[INDEX[i]]), as.integer(bin), as.integer(nbins), as.integer(nlevels), res, PACKAGE = "ffbase")
    }
    return(res)
  }
  
  if (is.factor.ff(bin)){
    bins <- levels(bin)
    nbins <- length(bins)
  } else {
    bins <- seq_len(nbins)
  }
  res <- matrix(0L, nrow=nbins, ncol=length(lev), dimnames=list(bin=bins, level=lev))
  for (i in chunk(x, ...)){
    Log$chunk(i)
    .Call("binned_tabulate", as.integer(x[i]), as.integer(bin[i]), as.integer(nbins), as.integer(nlevels), res, PACKAGE = "ffbase")
  }
  res
}

####### quick test ###################
# size <- 1e5
# x <- sample(c(1:4,NA), size=size, replace=TRUE)
# bin <- sample(1:100, size=size, replace=TRUE)
# nbins <- max(bin, na.rm=TRUE)
# nlevels <- max(x, na.rm=TRUE)
# 
# binned_tabulate(x, bin, nbins, nlevels)
# 
# 
# system.time(
#     replicate( 50
#              , binned_tabulate(x, bin, nbins, nlevels)
#              )
#            )
# 
# size <- 1e5
# x <- ff(sample(c(1:4,NA), size=size, replace=TRUE), vmode="integer")
# bin <- ff(sample(1:100, size=size, replace=TRUE))
# nbins <- max(bin, na.rm=TRUE)
# nlevels <- max(x, na.rm=TRUE)
# o <- ff(as.integer(200:1))
# binned_tabulate.ff(x, bin, nbins, nlevels, INDEX=o, by=1e2)
