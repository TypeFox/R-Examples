c.bayesx <- function(...)
{
  b <- list(...)
  x <- NULL
  for(i in 1L:length(b))
    x <- c(x, b[i])
  Call <- match.call()
  names(x) <- as.character(Call[-1L])
  class(x) <- "bayesx"

  return(x)
}

