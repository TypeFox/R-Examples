plot.linear.bayesx <- function(x, diagnostics = FALSE, ...) 
{
  if(is.null(x)) {
    warning("there is nothing to plot!")
    return(invisible(NULL))
  }
  args <- list(...)
  if(is.null(attr(x, "specs")$is.factor))
    attr(x, "specs")$is.factor <- FALSE
  if(!attr(x, "specs")$is.factor) {
    for(j in 2L:ncol(x))
      x[, j] <- x[, j] - mean(x[, j], na.rm = TRUE)
    args$x <- x
    args$diagnostics <- diagnostics
    do.call("plot.sm.bayesx", args)
  } else {
    if(!diagnostics) {
      args$x <- x
      do.call("plotblock", args)
    } else {
      if(!is.null(attr(x[[1]], "sample")) && diagnostics < 2L) {
        lx <- length(x)
        nx  <- length(attr(x[[1L]], "sample"))
        args$x <- matrix(0, nx, lx)
        for(j in 1L:lx)
          args$x[,j] <- attr(x[[j]], "sample")
        args$selected <- attr(x, "specs")$term
        args$var <- FALSE
        do.call("plotsamples", args)
      } else warning("term with no variance samples to plot!")
    }
  }

  return(invisible(NULL))
}

