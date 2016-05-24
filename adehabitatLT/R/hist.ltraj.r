"hist.ltraj" <-
function(x, which="dx/sqrt(dt)", ...)
{
  if (!inherits(x,"ltraj"))
    stop("x should be of class ltraj")
  opar <- par(mfrow=n2mfrow(length(x)))
  on.exit(par(opar))
  toto <- lapply(x, function(i) {
      if (!is.null(attr(i, "infolocs")))
          i <- cbind(i, attr(i, "infolocs"))
      ex<- parse(text=which)
      coin <- eval(ex, envir=i)
      tutu <- hist(coin, main=attr(i,"burst"), xlab=ex,...)
      return(tutu)
  })
  invisible(toto)
}

