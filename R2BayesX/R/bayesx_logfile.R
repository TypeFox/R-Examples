bayesx_logfile <- function(x, model = 1L)
{
  x <- get.model(x, model)
  bayesx.log <- NULL
  if(inherits(x, "bayesx")) {
    bayesx.log <- x[[model]]$bayesx.run$log
    if(length(bayesx.log) < 2L || class(bayesx.log) == "integer") {
      if(!is.null(x[[model]]$logfile))
        bayesx.log <- x[[model]]$logfile
    }
    if(!is.null(bayesx.log)) {
      if(!is.character(bayesx.log))
        print(bayesx.log)
      else
        writeLines(bayesx.log)
    }
  }
  if(is.null(bayesx.log))
    warning("logfile is not available!")

  return(invisible(bayesx.log))
}

