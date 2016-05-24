unique.id <- function(x)
{
  rval <- .Call("unique_id",
    as.numeric(x),
    as.numeric(unique(x)))

  return(rval)
}

