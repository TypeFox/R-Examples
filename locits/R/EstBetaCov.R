EstBetaCov <-
function (x, nz, filter.number = 1, family = "DaubExPhase", smooth.dev = var, 
    AutoReflect = TRUE, WPsmooth.type = "RM", 
    binwidth = 0, mkcoefOBJ, ThePsiJ, Cverbose = 0, verbose=0, OPLENGTH = 10^5, 
    ABB.tol=0.1, ABB.plot.it=FALSE, ABB.verbose=0,
    ABB.maxits=10, do.init=TRUE, ...) 
{
    if (missing(nz))
	stop("You need to specify nz")
    if (missing(x))
	stop("You need to specify x")
    TT <- length(x)

    if (nz < 1)
	stop("nz has to be >= 1")
    else if (nz > TT)
	stop(paste("nz has to be <= ", TT))

    if (do.init==TRUE)	{
	    ans <- .C("initThmStore", error=as.integer(0), PACKAGE="locits")

	    if (ans$error != 0)
		stop(paste("initThmStore: returned ", ans$error))
	    }

    if (binwidth==0)	{
	binwidth <- AutoBestBW(x=x, filter.number=filter.number,
		family=family, smooth.dev=smooth.dev,
		AutoReflect=AutoReflect, tol=ABB.tol, plot.it=ABB.plot.it,
		verbose=ABB.verbose, maxits=ABB.maxits) 
	if (verbose > 0)
	    cat("Choosing auto bandwidth of ", binwidth, "\n")
	}

#
# LACF computes the estimates of S and beta that we need using running mean
# smoothing with bandwidth given by binwidth (which might be auto-chosen)
#

	EWS <- ewspec3(x=x, filter.number=filter.number, family=family,
		smooth.dev=smooth.dev, AutoReflect=AutoReflect,
		WPsmooth.type=WPsmooth.type, binwidth=binwidth, ...)
		

    J <- nlevelsWT(EWS$S)

#
# Construct smoothed EWS estimator
#
    Smat <- matrix(EWS$S$D, nrow = J, ncol = length(x), 
        byrow = TRUE)
    Svec <- as.vector(Smat)
#
# Construct smoothed wavelet periodogram
#
    SmoothWmat <- matrix(EWS$SmoothWavPer$D, nrow = J, ncol = length(x), 
        byrow = TRUE)
    SmoothWvec <- as.vector(SmoothWmat)
#
# Generate wavelet quantities if they are not already supplied
#
    if (missing(mkcoefOBJ)) {
        mkcoefOBJ <- mkcoef(-J, filter.number = filter.number, 
            family = family)
    }
    if (missing(ThePsiJ)) {
        ThePsiJ <- PsiJ(-J, filter.number = filter.number, family = family, 
            OPLENGTH = OPLENGTH)
    }
    PsiJvec <- as.numeric(unlist(ThePsiJ))
    lPsiJ <- length(PsiJvec)
    lvPsiJ <- as.numeric(unlist(lapply(ThePsiJ, length)))
    linPsiJ <- c(0, cumsum(lvPsiJ)[1:(length(lvPsiJ) - 1)])
    lvPsiJ <- (lvPsiJ - 1)/2

#
# Set up space for covariance matrix
#
    Sigma <- matrix(0, nrow = J, ncol = J)
#
# Compute covariance matrix using the CstarIcov function
#
    for (ell in (1:J)) for (j in ell:J) {

	if (verbose>0)
		cat("ell: ", ell, " j: ", j, "\n")
        psil <- mkcoefOBJ[[ell]]
        lpsil <- length(psil)
        psij <- mkcoefOBJ[[j]]
        lpsij <- length(psij)
        ans <- .C("CstarIcov", ell = as.integer(ell), j = as.integer(j), 
            nz = as.integer(nz), s = as.integer(binwidth), TT = as.integer(TT), 
            IIvec = as.double(SmoothWvec), Svec = as.double(Svec), 
            J = as.integer(J), PsiJ = as.double(PsiJvec), lPsiJ = as.integer(lPsiJ), 
            linPsiJ = as.integer(linPsiJ), lvPsiJ = as.integer(lvPsiJ), 
            psil = as.double(psil), lpsil = as.integer(lpsil), 
            psij = as.double(psij), lpsij = as.integer(lpsij), 
            verbose = as.integer(Cverbose), ans = as.double(0), 
            error = as.integer(0), PACKAGE="locits")
        if (ans$error != 0) 
            return
        Sigma[ell, j] <- ans$ans
	Sigma[j, ell] <- ans$ans
    }

betahat <- SmoothWmat[, nz]

ll <- list(betahat = betahat, Sigma=Sigma)

return(ll)
}
