fRegress.stderr <- function(y, y2cMap, SigmaE, returnMatrix=FALSE, ...) {

#  FREGRESS.STDERR  computes standard error estimates for regression
#       coefficient functions estimated by function FREGRESS.
#
#  Arguments:
#
#  FREGRESSLIST ... a list object produced by function FREGRESS.  This
#                   is indicated by Y in the arguments since R syntax
#                   requires all of tghe fRegress family of functions to
#                   use this notation.
#  Y2CMAP       ... the matrix mapping from the vector of observed values
#                   to the coefficients for the dependent variable.
#                   This is output by function SMOOTH_BASIS.  If this is
#                   supplied, confidence limits are computed, otherwise not.
#  SIGMAE       ... Estimate of the covariances among the residuals.  This
#                   can only be estimated after a preliminary analysis
#                   with FREGRESS.
#
#  Returns:
#
#  BETASTDERRLIST ... a list object, each list containing a fdPar object
#                     for the standard error of a regression function.
#  BVAR           ... the symmetric matrix of sampling variances and
#                     covariances for the matrix of regression coefficients
#                     for the regression functions.  These are stored
#                     column-wise in defining BVARIANCE.
#  C2BMAP         ... the matrix mapping from response variable coefficients
#                     to coefficients for regression coefficients
#  RETURNMATRIX ... If False, a matrix in sparse storage model can be returned
#               from a call to function BsplineS.  See this function for
#               enabling this option.

#  Last modified 8 May 2012 by Jim Ramsay

#  generic argument y is actually fRegressList

fRegressList <- y

#  get number of independent variables

  xfdlist  <- fRegressList$xfdlist
  yfdPar   <- fRegressList$yfdPar
  betalist <- fRegressList$betalist
  Cmatinv  <- fRegressList$Cmatinv

  p <- length(xfdlist)

#  compute number of coefficients

  ncoef <- 0
  for (j in 1:p) {
    betaParfdj <- betalist[[j]]
    ncoefj     <- betaParfdj$fd$basis$nbasis
    ncoef      <- ncoef + ncoefj
  }

  if (inherits(yfdPar, "fdPar") || inherits(yfdPar, "fd")) {

    #  ----------------------------------------------------------------
    #           YFDPAR is functional for a functional parameter
    #  ----------------------------------------------------------------

    if (inherits(yfdPar, "fd")) yfdPar <- fdPar(yfdPar)

    #  get number of replications and basis information for YFDPAR

    yfd       <- yfdPar$fd
    N         <- dim(yfd$coefs)[2]
    ybasisobj <- yfdPar$fd$basis
    rangeval  <- ybasisobj$rangeval
    ynbasis   <- ybasisobj$nbasis
    nfine     <- max(501,10*ynbasis+1)
    tfine     <- seq(rangeval[1], rangeval[2], len=nfine)
    deltat    <- tfine[2] - tfine[1]

    ybasismat <- eval.basis(tfine, ybasisobj, 0, returnMatrix)

    #  compute BASISPRODMAT

    basisprodmat <- matrix(0,ncoef,ynbasis*N)

    mj2 <- 0
    for (j in 1:p) {
      betafdParj <- betalist[[j]]
      betabasisj <- betafdParj$fd$basis
      ncoefj     <- betabasisj$nbasis
      bbasismatj <- eval.basis(tfine, betabasisj, 0, returnMatrix)
      xfdj       <- xfdlist[[j]]
      tempj      <- eval.fd(tfine, xfdj, 0, returnMatrix)
        #  row indices of BASISPRODMAT to fill
      mj1    <- mj2 + 1
      mj2    <- mj2 + ncoefj
      indexj <- mj1:mj2
        #  inner products of beta basis and response basis
        #    weighted by covariate basis functions
      mk2 <- 0
      for (k in 1:ynbasis) {
            #  row indices of BASISPRODMAT to fill
        mk1    <- mk2 + 1
        mk2    <- mk2 + N
        indexk <- mk1:mk2
        tempk  <- bbasismatj*ybasismat[,k]
        basisprodmat[indexj,indexk] <- deltat*crossprod(tempk,tempj)
        }
    }

    #  check dimensions of Y2CMAP

    y2cdim <- dim(y2cMap)
    if (y2cdim[1] != ynbasis ||
        y2cdim[2] != dim(SigmaE)[1])
      stop("Dimensions of Y2CMAP not correct.")

    #  compute variances of regression coefficient function values

    c2bMap    <- Cmatinv %*% basisprodmat
    VarCoef   <- y2cMap %*% SigmaE %*% t(y2cMap)
    CVariance <- kronecker(VarCoef,diag(rep(1,N)))
    bvar      <- c2bMap %*% CVariance %*% t(c2bMap)
    betastderrlist <- vector("list",p)
    mj2 <- 0
    for (j in 1:p) {
      betafdParj <- betalist[[j]]
      betabasisj <- betafdParj$fd$basis
      ncoefj     <- betabasisj$nbasis
      mj1 	     <- mj2 + 1
      mj2 	     <- mj2 + ncoefj
      indexj 	   <- mj1:mj2
      bbasismat  <- eval.basis(tfine, betabasisj, 0, returnMatrix)
      bvarj      <- bvar[indexj,indexj]
      bstderrj   <- sqrt(diag(bbasismat %*% bvarj %*% t(bbasismat)))
      bstderrfdj <- smooth.basis(tfine, bstderrj, betabasisj)$fd
      betastderrlist[[j]] <- bstderrfdj
    }

  } else {

    #  ----------------------------------------------------------------
    #                   YFDPAR is scalar or multivariate
    #  ----------------------------------------------------------------

    ymat <- as.matrix(yfdPar)
    N    <- dim(ymat)[1]

    Zmat  <- NULL
    for (j in 1:p) {
      xfdj <- xfdlist[[j]]
      if (inherits(xfdj, "fd")) {
        xcoef      <- xfdj$coefs
        xbasis     <- xfdj$basis
        betafdParj <- betalist[[j]]
        bbasis     <- betafdParj$fd$basis
        Jpsithetaj <- inprod(xbasis,bbasis)
        Zmat       <- cbind(Zmat,t(xcoef) %*% Jpsithetaj)
      }
      else if (inherits(xfdj, "numeric")) {
        Zmatj <- xfdj
        Zmat  <- cbind(Zmat,Zmatj)
      }
    }

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
      ncoefj <- betabasisj$nbasis
      mj1    <- mj2 + 1
      mj2    <- mj2 + ncoefj
      indexj <- mj1:mj2
      bvarj  <- bvar[indexj,indexj]
      xfdj   <- xfdlist[[j]]
      if (inherits(xfdj,"fd")) {
        betarng    <- betabasisj$rangeval
        nfine      <- max(c(501,10*ncoefj+1))
        tfine      <- seq(betarng[1], betarng[2], len=nfine)
        bbasismat  <- eval.basis(tfine, betabasisj, 0, returnMatrix)
        bstderrj   <- sqrt(diag(bbasismat %*% bvarj %*% t(bbasismat)))
        bstderrfdj <- smooth.basis(tfine, bstderrj, betabasisj)$fd
      } else {
        bsterrj    <- sqrt(diag(bvarj))
        onebasis   <- create.constant.basis(betabasisj$rangeval)
        bstderrfdj <- fd(t(bstderrj), onebasis)
      }
      betastderrlist[[j]] <- bstderrfdj
    }
  }

  stderrList <-
    list(betastderrlist = betastderrlist,
         bvar           = bvar,
         c2bMap         = c2bMap)

  return(stderrList)

}

