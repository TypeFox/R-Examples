plot.varying.bayesx <-
function(x, ...)
{
  args <- list(...)
  if(is.null(args$ask))
    ask <- FALSE
  else
    ask <- args$ask
  nx <- length(x)
  if(!ask)
    setmfrow(nx)
  else
    par(ask = TRUE)
  for(i in 1L:nx) {
    args$x <- x[[i]]
    do.call("plot", args)
  }

  return(invisible(NULL))
}

