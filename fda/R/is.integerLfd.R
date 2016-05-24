is.integerLfd <- function(Lfdobj)
{
  #  check whether Lfd object is a simple differential operator

  #  Last modified 9 February 2007

  nderiv  <- Lfdobj$nderiv
  bintwrd <- TRUE
  if (nderiv > 0) {
    bwtlist <- Lfdobj$bwtlist
      if (!is.null(bwtlist)) {
    	  nderiv <- Lfdobj$nderiv
    	  for (j in 1:nderiv) {
          bfdj <- bwtlist[[j]]
          if (any(bfdj$coefs != 0.0)) bintwrd <- FALSE
    	  }
	}
    }
  bintwrd
}
