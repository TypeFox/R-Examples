setmfrow <- function(np, ask = NULL) 
{
  if(np > 12L) {
    if(is.null(ask)) ask <- TRUE
    par(ask = ask)
    par(mfrow = c(4L, 3L))
  } else par(mfrow = n2mfrow(np))

  return(invisible(NULL))
}

