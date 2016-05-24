bayesx_prgfile <- function(x, model = 1L)
{
  x <- get.model(x, model)
  bayesx.prg <- NULL
  if(inherits(x, "bayesx")) {
    bayesx.prg <- x[[model]]$bayesx.prg$prg
      if(!is.null(bayesx.prg))
        cat(bayesx.prg)
  }
  if(is.null(bayesx.prg))
    warning("program file is not available!")

  return(invisible(bayesx.prg))
}

