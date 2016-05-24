cpos <-
function(p, np) 
{
  rval <- .Call("cpos",
    as.numeric(p),
    as.integer(np),
    as.numeric(c(0, 0)))

  return(rval)
}

