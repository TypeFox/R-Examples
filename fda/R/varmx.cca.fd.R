varmx.cca.fd <- function(ccafd, nx=201)
{
#  VARMX_CCA  Apply varimax rotation to CCA weight functions
#             and scores
#  Arguments:
#  CCAFD ... An object of the CCA.FD class produced by a call to
#            function CCA.FD.
#  CCAWTFD2 ... A functional parameter object for the canonical weight
#                functions for the second set of functions.
#  CCAVAR1  ... Canonical variate scores for first  set of functions.
#  CCAVAR2  ... Canonical variate scores for second set of functions.
#  Return:  Arguments after rotation.

#  last modified 6 December 2012

  ccawtfd1 <- ccafd[[1]]
  ccawtfd2 <- ccafd[[2]]
  ccacorr  <- ccafd[[3]]
  ccavar1  <- ccafd[[4]]
  ccavar2  <- ccafd[[5]]

  wtcoef1 <- ccawtfd1$coefs
  wtcoef2 <- ccawtfd2$coefs

  basisobj <- ccawtfd1$basis
  rangex   <- basisobj$rangeval
  x        <- seq(rangex[1], rangex[2], len=nx)
  ccawtmat1 <- eval.fd(x, ccawtfd1)
  ccawtmat2 <- eval.fd(x, ccawtfd2)
  #  If fdmat is a 3-D array, stack into a matrix
  ccawtmat  <- rbind(ccawtmat1, ccawtmat2)
  #  compute rotation matrix for varimax rotation of ccawtmat
  rotmat <- varmx(ccawtmat)
  #  rotate coefficients and scores
  wtcoef1 <- wtcoef1 %*% rotmat
  wtcoef2 <- wtcoef2 %*% rotmat
  #  rotate ccawt objects
  ccawtfd1$coefs <- wtcoef1
  ccawtfd2$coefs <- wtcoef2
  #  rotate cca scores
  ccavar1 <- ccavar1 %*% rotmat
  ccavar2 <- ccavar2 %*% rotmat
  ccavard <- dim(ccavar1)
  canvarvalues      <- array(0, c(ccavard[1], ccavard[2], 2))
  canvarvalues[,,1] <- ccavar1
  canvarvalues[,,2] <- ccavar2

  ccavarmxlist <- list(ccawtfd1, ccawtfd2, ccacorr, canvarvalues, rotmat)
  names(ccavarmxlist) <- c("ccawtfd1", "ccawtfd2", "ccacorr",
                           "ccavar1",  "ccavar2")
  
  class(ccavarmxlist) <- "cca.fd"
  
  return(ccavarmxlist)

}
