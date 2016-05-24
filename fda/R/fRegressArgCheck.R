fRegressArgCheck <- function(yfdPar, xfdlist, betalist, wt=NULL) 
{
#  FREGRESS_ARGCHECK checks the first four arguments for the functions
#  for function regression, including FREGRESS.

#  --------------------  Check classes of arguments  --------------------

#  check YFDPAR and compute sample size N

if (inherits(yfdPar, "fd")) yfdPar <- fdPar(yfdPar)

if (!(inherits(yfdPar, "fdPar") || is.numeric(yfdPar))) stop(
  "First argument is not of class 'fd', 'fdPar' or 'numeric'.")

if (inherits(yfdPar, "fdPar")) {
    yfd   <- yfdPar$fd
    ycoef <- yfd$coefs
    N     <- dim(ycoef)[2]
}

if (is.numeric(yfdPar)) {
    N <- length(yfdPar)
}

#  check that xfdlist is a list object

#  check XFDLIST

if (inherits(xfdlist, "fd") || inherits(xfdlist, "numeric")) 
    xfdlist <- list(xfdlist)

if (!inherits(xfdlist, "list")) stop(
	"Argument XFDLIST is not a list object.")

#  get number of independent variables p

p <- length(xfdlist)

#  check BETALIST

if (inherits(betalist, "fd")) betalist <- list(betalist)

if (!inherits(betalist, "list")) stop(
	"Argument BETALIST is not a list object.")

if (length(betalist) != p)  {
	cat(paste("\nNumber of regression coefficients does not match\n",
		       "number of independent variables."))
	stop("")
}

#  check that the regression is functional, and extract the range

if (inherits(yfdPar, "fdPar")) {
    rangeval <- yfdPar$fd$basis$rangeval
} else {
    allscalar <- TRUE
    for (j in 1:p) {
        if (inherits(xfdlist[[j]], "fd")) {
            rangeval <- xfdlist[[j]]$basis$rangeval            
            allscalar <- FALSE
            break
        }
    }
    if (allscalar) stop(
        paste("The dependent variable and all the independent",   
              "variables are scalar."))
}


#  --------------------  check contents of arguments  -------------------

#  XFDLIST:
 
#  If the object is a vector of length N,
#  it is converted to a functional data object with a
#  constant basis

onebasis <- create.constant.basis(rangeval)
onesfd   <- fd(1,onebasis)

xerror <- FALSE
for (j in 1:p) {
    xfdj <- xfdlist[[j]]
    if (inherits(xfdj, "fd")) {
        xcoef <- xfdj$coefs
        if (length(dim(xcoef)) > 2) stop(
            paste("Covariate",j,"is not univariate."))
        #  check size of coefficient array
        Nj <- dim(xcoef)[2]
        if (Nj != N) {
            print(
               paste("Incorrect number of replications in XFDLIST",
                     "for covariate",j))
            xerror = TRUE
        }
    } 
    if (inherits(xfdj, "numeric")) {
        if (!is.matrix(xfdj)) xfdj = as.matrix(xfdj)
	  Zdimj <- dim(xfdj)
        if (Zdimj[1] != N && Zdimj != 1) {
            print(paste("Vector in XFDLIST[[",j,"]] has wrong length."))
				    xerror = TRUE 
	  } 
        if (Zdimj[2] != 1) {
            print(paste("Matrix in XFDLIST[[",j,"]] has more than one column."))
				    xerror = TRUE 
	  } 
        xfdlist[[j]] <- fd(matrix(xfdj,1,N), onebasis)
    } 
    if (!(inherits(xfdlist[[j]], "fd") || 
          inherits(xfdlist[[j]], "numeric"))) {
      print(paste("XFDLIST[[",j,"]] is neither an FD object nor numeric."))
      xerror = TRUE
    }
}
    
#  BETALIST:

berror <- FALSE
for (j in 1:p) {
	betafdParj <- betalist[[j]]
	if (inherits(betafdParj, "fd") || inherits(betafdParj, "basisfd")) {
		betafdParj    <- fdPar(betafdParj)
		betalist[[j]] <- betafdParj
	}
	if (!inherits(betafdParj, "fdPar")) {
		print(paste("BETALIST[[",j,"]] is not a FDPAR object."))
		berror <- TRUE
	}
}

if (xerror || berror) stop(
    "An error has been found in either XFDLIST or BETALIST.")

#  check weights

if (is.null(wt)) wt = rep(1,N)
if (length(wt) != N) stop("Number of weights not equal to N.")
if (any(wt < 0))     stop("Negative weights found.")

return(list(yfdPar=yfdPar, xfdlist=xfdlist, betalist=betalist, wt=wt))

}