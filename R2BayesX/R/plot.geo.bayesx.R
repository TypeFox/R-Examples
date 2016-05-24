plot.geo.bayesx <-
function(x, diagnostics = FALSE, ...)
{
  if(!is.null(x)) {
    args <- list(...)
    if(diagnostics == FALSE) {
      xattr <- attributes(x)
      if(is.null(args$map)) {
        plottype <- "plot3d"
        args$x <- x[, 2L:ncol(x)]
      } else {
        plottype <- "plot.mrf.bayesx"
        args$x <- x[,c(1L, 4L:ncol(x))]
        class(args$x) <- "mrf.bayesx"
      }
      an <- names(xattr)
      for(att in an)
        if(!att %in% c("dim", "dimnames", "names")) {
            attr(args$x, att) <- xattr[[att]]
        }
      do.call(plottype, args)
    } else coeffun(x, args, diagnostics)
  } else warning("there is nothing to plot!")

  return(invisible(NULL))
}

