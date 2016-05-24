s4dim <- function(x)
{
  dim <- 1L
  if(!is.null(dim(x))) {
    if(ncol(x) == 10L)
      dim <- 2L
    if(ncol(x) == 11L)
      dim <- 2L
    if(ncol(x) == 16L)
      dim <- 2L
  }
  return(dim)
}

