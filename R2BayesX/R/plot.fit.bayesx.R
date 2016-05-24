plot.fit.bayesx <-
function(x, ...)
{
  if(length(x) >= 1L) {
    if(any(grepl("bayesx", class(x[[1L]]))))
      x <- list(x)
    nx <- length(x)
    pval <- list()
    for(i in 1L:nx)
      pval[[i]] <- list(effects = x[[i]])
    class(pval) <- "bayesx"
    plot(pval, ...)
  } else stop("there is nothing to plot!")

  return(invisible(NULL))
}

