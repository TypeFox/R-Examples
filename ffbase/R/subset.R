#' Subsetting a ff vector or ffdfdata frame
#'
#' @export
#' @export subset.ff
#' @aliases subset.ff subset.ffdf
#' @method subset ff
#' @param x \code{ff} vector or \code{ffdf} data.frame to be subset
#' @param subset an expression, \code{ri}, \code{bit} or logical \code{ff} vector that can be used to index x
#' @param ... not used
#' @return a new ff vector containing the subset, data is physically copied
subset.ff <- function(x, subset, ...){
  if (missing(subset)){
    return(ff::clone(x))
  }
  ss <- as.expression(substitute(subset))
  try(ss <- subset, silent=TRUE)
  .ffdf <- ffdf(x=x)
  idx <- ffwhich(.ffdf, ss)
	x[idx]
}

#' @export
#' @export subset.ffdf
subset.ffdf <- function(x, subset, select, drop = FALSE, ...){
  # remove rownames otherwise we have errors...
  rownames(x) <- NULL
  if (missing(subset)){
    idx = ffseq_len(nrow(x))
  } else {  
    ss <- as.expression(substitute(subset))
    #try(ss <- subset, silent=TRUE)
    idx <- ffwhich.ffdf(x, ss, envir=parent.frame())
  }
  if (missing(select)){
    select <- names(x)
  } else {
    nl <- as.list(seq_along(x))
    names(nl) <- names(x)
    select <- eval(substitute(select), nl, parent.frame())
  }
  x[idx, select, drop=drop]
}


# quick testing
# x <- as.ffdf(iris)
# log <- x$Species == "setosa"
# 
# subset(x, Species=="setosa")
#ffwhich(x, log)
