lines.fdSmooth <- function(x, Lfdobj=int2Lfd(0), nx=201, ...){
  lines(x$fd, Lfdobj=Lfdobj, nx=nx, ...)
}

lines.fd <- function(x, Lfdobj=int2Lfd(0), nx=201, ...)
{
  #  Plot a functional data object FD using lines in a pre-existing plot.
  #  If there are multiple variables, each curve will appear in the same plot.
  #  The remaining optional arguments are the same as those available
  #     in the regular "lines" function.
  
# Last modified 2008.07.05 by Spencer Graves
  #  previously modified 2007.05.03 and 1 October 2005
  fdobj <- x
  
  if (!inherits(fdobj,  "fd"))  stop(
		"First argument is not a functional data object.")
  if (!inherits(Lfdobj, "Lfd")) stop(
      "Second argument is not a linear differential operator.")

  coef   <- fdobj$coefs
  coefd  <- dim(coef)
  ndim   <- length(coefd)
  nbasis <- coefd[1]
  nrep   <- coefd[2]
  if (ndim > 3) nvar <- coefd[3] else nvar <- 1
  crvnames <- fdobj$fdnames[[2]]
  varnames <- fdobj$fdnames[[3]]

  basisobj <- fdobj$basis
#
  xlim <- par('usr')[1:2]
  if(par('xlog')) xlim <- 10^xlim
# 
  rngx <- basisobj$rangeval 
  xmin <- max(rngx[1], xlim[1])
  xmax <- min(rngx[2], xlim[2])
  x.        <- seq(xmin, xmax, length=nx)
  fdmat    <- eval.fd(x.,fdobj,Lfdobj)

  if (length(dim(coef)) < 2) {
    lines (x.,fdmat,...)
  }
  if (length(dim(coef)) ==2 ) {
    matlines (x.,fdmat,...)
  }
  if (length(dim(coef)) == 3) {
    for (ivar in 1:nvar) {
      matlines (x.,fdmat[,,ivar],type="l",lty=1,
                main=varnames[ivar],...)
    }
  }
  invisible()
}
