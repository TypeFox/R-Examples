Rvarlacf <-
function (x, nz, filter.number = 1, family = "DaubExPhase", smooth.dev = var, 
    AutoReflect = TRUE, lag.max = NULL, WPsmooth.type = "RM", 
    binwidth = 0, mkcoefOBJ, ThePsiJ, Cverbose = 0, verbose=0, OPLENGTH = 10^5, 
    var.lag.max = 3, ABB.tol=0.1, ABB.plot.it=FALSE, ABB.verbose=0,
    ABB.maxits=10, ...) 
{
    if (missing(nz))
	stop("You need to specify nz")
    if (missing(x))
	stop("You need to specify x")
    TT <- length(x)
    ans <- .C("initThmStore", error=as.integer(0), PACKAGE="locits")

    if (ans$error != 0)
	stop(paste("initThmStore: returned ", ans$error))

    if (binwidth==0)	{
	binwidth <- AutoBestBW(x=x, filter.number=filter.number,
		family=family, smooth.dev=smooth.dev,
		AutoReflect=AutoReflect, tol=ABB.tol, plot.it=ABB.plot.it,
		verbose=ABB.verbose, maxits=ABB.maxits) 
	if (verbose > 0)
	    cat("Choosing auto bandwidth of ", binwidth, "\n")
	}
    the.lacf <- lacf(x = x, filter.number = filter.number, family = family, 
        smooth.dev = smooth.dev, AutoReflect = AutoReflect, lag.max = lag.max, 
        WPsmooth.type = WPsmooth.type, binwidth = binwidth, ...)
    J <- the.lacf$J
    Smat <- matrix(the.lacf$S$D, nrow = J, ncol = length(x), 
        byrow = TRUE)
    Svec <- as.vector(Smat)
    SmoothWmat <- matrix(the.lacf$SmoothWP$D, nrow = J, ncol = length(x), 
        byrow = TRUE)
    SmoothWvec <- as.vector(SmoothWmat)
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
    Sigma <- matrix(0, nrow = J, ncol = J)
    for (ell in (1:(J - 1))) for (j in ell:J) {

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
    }
    A <- ipndacw(-J, filter.number = filter.number, family = family)
    sA <- solve(A)
    PMat <- PsiJmat(-J, filter.number = filter.number, family = family, 
        OPLENGTH = OPLENGTH)
    pd <- ncol(PMat)
    mid <- (pd - 1)/2
    cvar <- rep(0, var.lag.max + 1)
    for (tau in 0:var.lag.max) {
        kappa <- sA %*% PMat[, mid + tau]
        secsum <- 0
        for (ell in (1:(J - 1))) for (j in (ell + 1):J) {
            secsum <- secsum + kappa[ell] * kappa[j] * Sigma[ell, 
                j]
        }
        bigsum <- 0
        for (ell in 1:J) bigsum <- bigsum + (kappa[ell]^2) * 
            Sigma[ell, ell]
        bigsum <- bigsum + 2 * secsum
        cvar[tau + 1] <- bigsum
    }
    cvar <- pmax(cvar, 0)
    l <- list(lag = 0:var.lag.max, cvar = cvar, the.lacf = the.lacf, nz=nz)
    class(l) <- "lacfCI"
    l
}
