bayesx_runtime <- runtime <- function(x, model = 1L)
{
  x <- get.model(x, model)
  runtime <- NULL
  if(inherits(x, "bayesx")) {
    runtime <- x[[model]]$bayesx.run$runtime
      if(!is.null(runtime))
        print(runtime)
  }
  if(is.null(runtime))
    warning("runtime information is not available!")

  return(invisible(runtime))
}

