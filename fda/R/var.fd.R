#var <- function (x, ...) UseMethod("var")

#var.default <- function (x, y = NULL, na.rm = FALSE, use, ...) {
## First preprocess 'use'
#  if(missing(use)){
#    use <- {
#      if (na.rm) "complete.obs"
#      else "all.obs"
#    }
#  }
## Now use the official non-generic.
#  stats::var (x, y = NULL, na.rm = FALSE, use)
#}

# Seemed to work but did not pass R CMD check:
#* checking for missing documentation entries ... WARNING
#Undocumented code objects:
#  sd sd.default var var.default
#All user-level objects in a package should have documentation

# I'm not sure how to fix, so revert

#var.fd <- function(x, y=x, ...)
var.fd <- function(fdobj1, fdobj2=fdobj1)
{
#  fdobj1 <- x
#  fdobj2 <- y
  #  compute the variance and covariance functions for functional observations

  #  Last modified 26 October 2005

  if (!(inherits(fdobj1, "fd"))) stop(
		"Argument FDX not a functional data object.")
  if (!(inherits(fdobj2, "fd"))) stop(
		"Argument FDY not a functional data object.")

  coefx   <- fdobj1$coefs
  coefy   <- fdobj2$coefs
  coefdobj1  <- dim(coefx)
  coefdobj2  <- dim(coefy)
  basisx  <- fdobj1$basis
  basisy  <- fdobj2$basis
  nbasisx <- basisx$nbasis
  nbasisy <- basisy$nbasis

  if (coefdobj1[2] != coefdobj2[2]) stop(
    	"Number of replications are not equal.")
  if (length(coefdobj1) == 2) {
    	if(length(coefdobj2) == 2) {
      		coefvar   <- var(t(coefx),t(coefy))
      		coefnames <- list(dimnames(coefx)[[1]], dimnames(coefy)[[1]])
      		varbifd   <- bifd(coefvar, basisx, basisy, coefnames)
    	} else stop("Both arguments must be univariate.")
  } else {
    	nvar    <- coefdobj1[3]
    	npair   <- nvar*(nvar+1)/2
    	coefvar <- array(0,c(nbasisx,nbasisx,1,npair))
       varnames <- fdobj1$fdnames[[3]]
    	m <- 0
       bivarnames <- vector("list",npair)
    	for (i in 1:nvar) for (j in 1:i) {
      		m <- m + 1
      		coefvar[,,1,m] <- var(t(coefx[,,i]),t(coefx[,,j]))
          bivarnames[m] <- paste(varnames[i],"vs",varnames[j])
    	}
       bifdnames = fdobj1$fdnames
       bifdnames[[3]] <- bivarnames
    	varbifd <- bifd(coefvar, basisx, basisx, bifdnames)
  }
  return(varbifd)
}
