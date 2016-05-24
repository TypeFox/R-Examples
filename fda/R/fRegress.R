fRegress <- function(y, xfdlist, betalist, wt=NULL,
                     y2cMap=NULL, SigmaE=NULL, returnMatrix=FALSE)
{

#  FREGRESS  Fits a functional linear model using multiple
#  functional independent variables with the dependency being
#  pointwise or concurrent.
#  The case of a scalar independent variable is included by treating
#  it as a functional independent variable with a constant basis
#  and a unit coefficient.
#
#  Arguments:
#  Y        ... an object for the dependent variable,
#               which may be:
#                   a functional data object,
#                   a functional parameter (fdPar) object, or
#                   a vector
#  XFDLIST  ... a list object of length p with each list
#               containing an object for an independent variable.
#               the object may be:
#                   a functional data object or
#                   a vector
#               if XFDLIST is a functional data object or a vector,
#               it is converted to a list of length 1.
#  BETALIST ... a list object of length p with each list
#               containing a functional parameter object for
#               the corresponding regression function.  If any of
#               these objects is a functional data object, it is
#               converted to the default functional parameter object.
#               if BETALIST is a functional parameter object
#               it is converted to a list of length 1.
#  WT       ... a vector of nonnegative weights for observations
#  Y2CMAP   ... the matrix mapping from the vector of observed values
#               to the coefficients for the dependent variable.
#               This is output by function SMOOTH_BASIS.  If this is
#               supplied, confidence limits are computed, otherwise not.
#  SIGMAE   ... Estimate of the covariances among the residuals.  This
#               can only be estimated after a preliminary analysis
#               with FREGRESS.
#  RETURNMATRIX ... If False, a matrix in sparse storage model can be returned
#               from a call to function BsplineS.  See this function for
#               enabling this option.
#
#  Returns FREGRESSLIST  ... A list containing seven members with names:
#    yfdPar      ... first  argument of FREGRESS
#    xfdlist     ... second argument of FREGRESS
#    betalist    ... third  argument of FREGRESS
#    betaestlist ... estimated regression functions
#    yhatfdobj   ... functional data object containing fitted functions
#    Cmat        ... coefficient matrix for the linear system defining
#                    the regression coefficient basis vector
#    Dmat        ... right side vector for the linear system defining
#                    the regression coefficient basis vector
#    Cmatinv     ... inverse of the coefficient matrix, needed for
#                    function FREGRESS.STDERR that computes standard errors
#    wt          ... weights for observations
#    df          ... degrees of freedom for fit

#  Last modified 28 December 2012 by Jim Ramsay

arglist <- fRegressArgCheck(y, xfdlist, betalist, wt)

yfdPar   <- arglist$yfdPar
xfdlist  <- arglist$xfdlist
betalist <- arglist$betalist
wt       <- arglist$wt

p <- length(xfdlist)
N <- dim(xfdlist[[1]]$coef)[1]

wtconst <- var(wt) == 0

#  --------------------------------------------------------------------------
#  branch depending on whether the dependent variable is functional or scalar
#  --------------------------------------------------------------------------

if (inherits(yfdPar, "fdPar") || inherits(yfdPar, "fd")) {
    if (inherits(yfdPar, "fd")) yfdPar = fdPar(yfdPar)

    #  ----------------------------------------------------------------
    #           YFDPAR is a functional parameter object
    #  ----------------------------------------------------------------

    #  extract dependent variable information

    yfdobj    <- yfdPar$fd
    ylambda   <- yfdPar$lambda
    yLfdobj   <- yfdPar$Lfd
    ycoef     <- yfdobj$coefs
    ycoefdim  <- dim(ycoef)
    N         <- ycoefdim[2]
    ybasisobj <- yfdobj$basis
    rangeval  <- ybasisobj$rangeval
    ynbasis   <- ybasisobj$nbasis
    onesbasis <- create.constant.basis(rangeval)
    onesfd    <- fd(1,onesbasis)

    if (length(ycoefdim) > 2) stop("YFDOBJ from YFDPAR is not univariate.")

    #  --------  set up the linear equations for the solution  -----------

    #  compute the total number of coefficients to be estimated

    ncoef <- 0
    for (j in 1:p) {
        betafdParj <- betalist[[j]]
        if (betafdParj$estimate) {
        	ncoefj     <- betafdParj$fd$basis$nbasis
        	ncoef      <- ncoef + ncoefj
        }
    }

    Cmat <- matrix(0,ncoef,ncoef)
    Dmat <- rep(0,ncoef)

    #  loop through rows of CMAT

    mj2 <- 0
    for (j in 1:p) {
        betafdParj <- betalist[[j]]
        if (betafdParj$estimate) {
            betafdj    <- betafdParj$fd
            betabasisj <- betafdj$basis
            ncoefj     <- betabasisj$nbasis
            #  row indices of CMAT and DMAT to fill
            mj1    <- mj2 + 1
            mj2    <- mj2 + ncoefj
            indexj <- mj1:mj2
            #  compute right side of equation DMAT
            xfdj <- xfdlist[[j]]
            if (wtconst) {
                xyfdj <- xfdj*yfdobj
            } else {
                xyfdj <- (xfdj*wt)*yfdobj
            }
            wtfdj <- sum(xyfdj)
            Dmatj <- inprod(betabasisj,onesfd,0,0,rangeval,wtfdj)
            Dmat[indexj] <- Dmatj
            #  loop through columns of CMAT
            mk2 <- 0
            for (k in 1:j) {
                betafdPark <- betalist[[k]]
                if (betafdPark$estimate) {
                    betafdk    <- betafdPark$fd
                    betabasisk <- betafdk$basis
                    ncoefk     <- betabasisk$nbasis
                    #  column indices of CMAT to fill
                    mk1 <- mk2 + 1
                    mk2 <- mk2 + ncoefk
                    indexk <- mk1:mk2
                    #  set up two weight functions
                    xfdk <- xfdlist[[k]]
                    if (wtconst) {
                        xxfdjk <- xfdj*xfdk
                    } else {
                        xxfdjk <- (xfdj*wt)*xfdk
                    }
                    wtfdjk <- sum(xxfdjk)
                    Cmatjk <- inprod(betabasisj, betabasisk, 0, 0,
                                     rangeval, wtfdjk)
                    Cmat[indexj,indexk] <- Cmatjk
                    Cmat[indexk,indexj] <- t(Cmatjk)
                }
            }
            #  attach penalty term to diagonal block
            lambdaj <- betafdParj$lambda
            if (lambdaj > 0) {
                Rmatj <- betafdParj$penmat
                if (is.null(Rmatj)) {
                    Lfdj  <- betafdParj$Lfd
                    Rmatj <- eval.penalty(betabasisj, Lfdj)
                }
                Cmat[indexj,indexj] <- Cmat[indexj,indexj] +
                                       lambdaj*Rmatj
            }
        }
    }

    Cmat <- (Cmat+t(Cmat))/2

    #  check Cmat for singularity

    eigchk(Cmat)


    #  solve for coefficients defining BETA

    Lmat    <- chol(Cmat)
    Lmatinv <- solve(Lmat)
    Cmatinv <- Lmatinv %*% t(Lmatinv)

    betacoef <- Cmatinv %*% Dmat

    #  set up fdPar objects for reg. fns. in BETAESTLIST

    betaestlist <- betalist
    mj2 <- 0
    for (j in 1:p) {
        betafdParj <- betalist[[j]]
        if (betafdParj$estimate) {
            betafdj <- betafdParj$fd
            ncoefj  <- betafdj$basis$nbasis
            mj1     <- mj2 + 1
            mj2     <- mj2 + ncoefj
            indexj  <- mj1:mj2
            coefj   <- betacoef[indexj]
            betafdj$coefs <- as.matrix(coefj)
            betafdParj$fd <- betafdj
        }
        betaestlist[[j]] <- betafdParj
    }

    #  set up fd objects for predicted values in YHATFDOBJ

    nfine     <- max(501,10*ynbasis+1)
    tfine     <- seq(rangeval[1], rangeval[2], len=nfine)
    yhatmat <- matrix(0,nfine,N)
    for (j in 1:p) {
	  xfdj       <- xfdlist[[j]]
        xmatj      <- eval.fd(tfine, xfdj, 0, returnMatrix)
        betafdParj <- betaestlist[[j]]
        betafdj    <- betafdParj$fd
        betavecj   <- eval.fd(tfine, betafdj, 0, returnMatrix)
        yhatmat    <- yhatmat + xmatj*as.vector(betavecj)
    }
    yhatfdobj <- smooth.basis(tfine, yhatmat, ybasisobj)$fd

    df <- NA

    #  -----------------------------------------------------------------------
    #        Compute pointwise standard errors of regression coefficients
    #               if both y2cMap and SigmaE are supplied.
    #  -----------------------------------------------------------------------

    if (!(is.null(y2cMap) || is.null(SigmaE))) {

        #  check dimensions of y2cMap and SigmaE

        y2cdim = dim(y2cMap)
        if (y2cdim[1] != ynbasis || y2cdim[2] != dim(SigmaE)[1]) {
             stop("Dimensions of Y2CMAP not correct.")
        }

        ybasismat = eval.basis(tfine, ybasisobj, 0, returnMatrix)

        deltat    = tfine[2] - tfine[1]

        #  compute BASISPRODMAT

        basisprodmat = matrix(0,ncoef,ynbasis*N)

        mj2 = 0
        for (j in 1:p) {
            betafdParj = betalist[[j]]
            betabasisj = betafdParj$fd$basis
            ncoefj     = betabasisj$nbasis
            bbasismatj = eval.basis(tfine, betabasisj, 0, returnMatrix)
            xfdj       = xfdlist[[j]]
            tempj      = eval.fd(tfine, xfdj, 0, returnMatrix)
            #  row indices of BASISPRODMAT to fill
            mj1    = mj2 + 1
            mj2    = mj2 + ncoefj
            indexj = mj1:mj2
            #  inner products of beta basis and response basis
            #    weighted by covariate basis functions
            mk2 = 0
            for (k in 1:ynbasis) {
                #  row indices of BASISPRODMAT to fill
                mk1    = mk2 + 1
                mk2    = mk2 + N
                indexk = mk1:mk2
                tempk  = bbasismatj*ybasismat[,k]
                basisprodmat[indexj,indexk] =
                                 deltat*crossprod(tempk,tempj)
            }
        }

        #  compute variances of regression coefficient function values

        c2bMap    = solve(Cmat,basisprodmat)
        VarCoef   = y2cMap %*% SigmaE %*% t(y2cMap)
        CVariance = kronecker(VarCoef,diag(rep(1,N)))
        bvar      = c2bMap %*% CVariance %*% t(c2bMap)
        betastderrlist = vector("list", p)
        mj2 = 0
        for (j in 1:p) {
            betafdParj = betalist[[j]]
            betabasisj = betafdParj$fd$basis
            ncoefj     = betabasisj$nbasis
            mj1 	     = mj2 + 1
            mj2 	     = mj2 + ncoefj
            indexj     = mj1:mj2
            bbasismat  = eval.basis(tfine, betabasisj, 0, returnMatrix)
            bvarj      = bvar[indexj,indexj]
            bstderrj   = sqrt(diag(bbasismat %*% bvarj %*% t(bbasismat)))
            bstderrfdj = smooth.basis(tfine, bstderrj, betabasisj)$fd
            betastderrlist[[j]] = bstderrfdj
        }
    } else {
        betastderrlist = NULL
        bvar           = NULL
        c2bMap         = NULL
    }


    #  -------------------------------------------------------------------
    #                       Set up output list object
    #  -------------------------------------------------------------------

    fRegressList <-
       	list(yfdPar         = yfdPar,
                 xfdlist        = xfdlist,
                 betalist       = betalist,
                 betaestlist    = betaestlist,
                 yhatfdobj      = yhatfdobj,
                 Cmat           = Cmat,
                 Dmat           = Dmat,
                 Cmatinv        = Cmatinv,
                 wt             = wt,
                 df             = df,
                 y2cMap         = y2cMap,
                 SigmaE         = SigmaE,
                 betastderrlist = betastderrlist,
                 bvar           = bvar,
                 c2bMap         = c2bMap)

 }

 #  -------------------------------------------------------------------

 if (inherits(yfdPar,"numeric")) {

    #  ----------------------------------------------------------------
    #                   YFDPAR is scalar or multivariate
    #  ----------------------------------------------------------------

    ymat <- as.matrix(yfdPar)
    N    <- dim(ymat)[1]

    Zmat  <- NULL
    Rmat  <- NULL
    pjvec <- rep(0,p)
    ncoef <- 0
    for (j in 1:p) {
        xfdj       <- xfdlist[[j]]
        xcoef      <- xfdj$coefs
        xbasis     <- xfdj$basis
        betafdParj <- betalist[[j]]
        bbasis     <- betafdParj$fd$basis
        bnbasis    <- bbasis$nbasis
        pjvec[j]   <- bnbasis
        Jpsithetaj <- inprod(xbasis,bbasis)
        Zmat       <- cbind(Zmat,crossprod(xcoef,Jpsithetaj))
        if (betafdParj$estimate) {
            lambdaj    <- betafdParj$lambda
            if (lambdaj > 0) {
                Lfdj  <- betafdParj$Lfd
                Rmatj <- lambdaj*eval.penalty(bbasis,Lfdj)
            } else {
                Rmatj <- matrix(0,bnbasis,bnbasis)
            }
            if (ncoef > 0) {
                zeromat <- matrix(0,ncoef,bnbasis)
                Rmat    <- rbind(cbind(Rmat,       zeromat),
                                 cbind(t(zeromat), Rmatj))
            } else {
                Rmat  <- Rmatj
                ncoef <- ncoef + bnbasis
            }
        }
    }

    #  -----------------------------------------------------------
    #          set up the linear equations for the solution
    #  -----------------------------------------------------------

    #  solve for coefficients defining BETA

    if (any(wt != 1)) {
        rtwt   <- sqrt(wt)
        Zmatwt <- Zmat*rtwt
        ymatwt <- ymat*rtwt
        Cmat   <- t(Zmatwt) %*% Zmatwt + Rmat
        Dmat   <- t(Zmatwt) %*% ymatwt
    } else {
        Cmat <- t(Zmat) %*% Zmat + Rmat
        Dmat <- t(Zmat) %*% ymat
    }

    eigchk(Cmat)

    Cmatinv  <- solve(Cmat)

    betacoef <- Cmatinv %*% Dmat


    #  compute and print degrees of freedom measure

    df <- sum(diag(Zmat %*% Cmatinv %*% t(Zmat)))

    #  set up fdPar object for BETAESTFDPAR

    betaestlist <- betalist
    mj2 <- 0
    for (j in 1:p) {
        betafdParj <- betalist[[j]]
        betafdj    <- betafdParj$fd
        ncoefj     <- betafdj$basis$nbasis
        mj1        <- mj2 + 1
        mj2        <- mj2 + ncoefj
        indexj     <- mj1:mj2
        betacoefj        <- betacoef[indexj]
        betaestfdj       <- betafdj
        betaestfdj$coefs <- as.matrix(betacoefj)
        betaestfdParj    <- betafdParj
        betaestfdParj$fd <- betaestfdj
        betaestlist[[j]] <- betaestfdParj
    }

    #  set up fd object for predicted values

    yhatmat <- matrix(0,N,1)
    for (j in 1:p) {
        xfdj <- xfdlist[[j]]
        if (inherits(xfdj, "fd")) {
            xbasis     <- xfdj$basis
            xnbasis    <- xbasis$nbasis
            xrng       <- xbasis$rangeval
            nfine      <- max(501,10*xnbasis+1)
            tfine      <- seq(xrng[1], xrng[2], len=nfine)
            deltat     <- tfine[2]-tfine[1]
            xmat       <- eval.fd(tfine, xfdj, 0, returnMatrix)
            betafdParj <- betaestlist[[j]]
            betafdj    <- betafdParj$fd
            betamat    <- eval.fd(tfine, betafdj, 0, returnMatrix)
            fitj       <- deltat*(crossprod(xmat,betamat) -
						                 0.5*(outer(xmat[1,    ],betamat[1,    ]) +
				                          outer(xmat[nfine,],betamat[nfine,])))
            yhatmat    <- yhatmat + fitj
        } else{
	          betaestfdParj <- betaestlist[[j]]
            betavecj      <- betaestfdParj$fd$coefs
            yhatmat       <- yhatmat + xfdj %*% t(betavecj)
        }
    }
    yhatfdobj <- yhatmat

    #  -----------------------------------------------------------------------
    #        Compute pointwise standard errors of regression coefficients
    #               if both y2cMap and SigmaE are supplied.
    #  -----------------------------------------------------------------------

    if (!(is.null(y2cMap) || is.null(SigmaE))) {

        #  check dimensions of y2cMap and SigmaE

        y2cdim <- dim(y2cMap)
        if (y2cdim[2] != dim(SigmaE)[1])  stop(
                         "Dimensions of Y2CMAP not correct.")


        #  compute linear mapping c2bMap takinging coefficients for
        #  response into coefficients for regression functions

        c2bMap <- Cmatinv %*% t(Zmat)
        y2bmap <- c2bMap
        bvar   <- y2bmap %*% as.matrix(SigmaE) %*% t(y2bmap)
        betastderrlist <- vector("list",p)
        mj2 <- 0
        for (j in 1:p) {
	      betafdParj <- betalist[[j]]
            betabasisj <- betafdParj$fd$basis
	      ncoefj     <- betabasisj$nbasis
            mj1        <- mj2 + 1
            mj2        <- mj2 + ncoefj
            indexj     <- mj1:mj2
            bvarj      <- bvar[indexj,indexj]
            betarng    <- betabasisj$rangeval
            nfine      <- max(c(501,10*ncoefj+1))
            tfine      <- seq(betarng[1], betarng[2], len=nfine)
            bbasismat  <- eval.basis(tfine, betabasisj, 0, returnMatrix)
            bstderrj   <- sqrt(diag(bbasismat %*% bvarj %*% t(bbasismat)))
            bstderrfdj <- smooth.basis(tfine, bstderrj, betabasisj)$fd
            betastderrlist[[j]] <- bstderrfdj
        }
    } else {
        betastderrlist = NULL
        bvar           = NULL
        c2bMap         = NULL
    }

    #  -----------------------------------------------------------------------
    #                  Set up output list object
    #  -----------------------------------------------------------------------

    fRegressList <-
       	list(yfdPar         = yfdPar,
                 xfdlist        = xfdlist,
                 betalist       = betalist,
                 betaestlist    = betaestlist,
                 yhatfdobj      = yhatfdobj,
                 Cmat           = Cmat,
                 Dmat           = Dmat,
                 Cmatinv        = Cmatinv,
                 wt             = wt,
                 df             = df,
                 y2cMap         = y2cMap,
                 SigmaE         = SigmaE,
                 betastderrlist = betastderrlist,
                 bvar           = bvar,
                 c2bMap         = c2bMap)
 }

 return(fRegressList)

}

#  ------------------------------------------------------------------------

eigchk <- function(Cmat) {

    #  check Cmat for singularity

    eigval <- eigen(Cmat)$values
    ncoef  <- length(eigval)
    if (eigval[ncoef] < 0) {
	     neig <- min(length(eigval),10)
        cat("\nSmallest eigenvalues:\n")
        print(eigval[(ncoef-neig+1):ncoef])
        cat("\nLargest  eigenvalues:\n")
        print(eigval[1:neig])
        stop("Negative eigenvalue of coefficient matrix.")
    }
    if (eigval[ncoef] == 0) stop("Zero eigenvalue of coefficient matrix.")
    logcondition <- log10(eigval[1]) - log10(eigval[ncoef])
    if (logcondition > 12) {
        warning("Near singularity in coefficient matrix.")
        cat(paste("\nLog10 Eigenvalues range from\n",
                  log10(eigval[ncoef])," to ",log10(eigval[1]),"\n"))
    }
}

