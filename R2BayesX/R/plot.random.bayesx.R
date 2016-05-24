plot.random.bayesx <-
function(x, diagnostics = FALSE, ...)
{	
  if(is.null(x)) {
    warning("there is nothing to plot!")
    return(invisible(NULL))
  }
  args <- list(...)
  args$x <- x
  args$diagnostics <- diagnostics
  do.call("plot.mrf.bayesx", args)

  return(invisible(NULL))
}

