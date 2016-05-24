plot.bnd <-
function(x, ...)
{
  if(is.null(x))
    return(invisible(NULL))
  args <- list(...)
  args$map <- x
  do.call("plotmap", args)

  return(invisible(NULL))
}

