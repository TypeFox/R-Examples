Boutrp <-
function (x, cmat) 
{
  ngroup <- ncol(x$t)
  f <- x$strata
  ni <- unlist(lapply(split(f, f = f), length))
  gnames <- names(x$t0)
  names(ni) <- gnames
  if (any(ni < 5)) {
    warning("For sample sizes les than 5 this function hardly makes sense!")
  }
  chains <- x$t
  out <- CCdrp(x = chains, cmat = cmat)
  return(out)
}

