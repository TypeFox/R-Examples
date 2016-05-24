fRegress.CV <- function(y, xfdlist, betalist, wt=NULL, CVobs=1:N,
                        returnMatrix=FALSE, ...)
{

# FREGRESS.CV computes cross-validated error sum of squares
# for scalar or functional responses. NOTE: ordinary and
# generalized cross validation scores are now returned by fRegress
# when scalar responses are used.

# last modified 28 July 2012 by Jim Ramsay

#  check the arguments

argList  <- fRegressArgCheck(y, xfdlist, betalist, wt)
yfdPar   <- argList$yfdPar
xfdlist  <- argList$xfdlist
betalist <- argList$betalist
wt       <- argList$wt

# extract dimensions of the data and analysis

p <- length(xfdlist)
N <- dim(xfdlist[[1]]$coef)[2]
M <- length(CVobs)

#  branch to either scalar or functional dependent variable

if (inherits(yfdPar, "numeric"))  {

    #  scalar dependent variable case

    yvec   <- yfdPar
    SSE.CV <- 0
    errfd  <- c()
    for (m in 1:M) {
      i        <- CVobs[m]
      #  eliminate case i from the weights
      wti <- wt[-i]
      xfdlisti <- vector("list",p)
      for (j in 1:p) {
        xfdj          <- xfdlist[[j]]
        if (inherits(xfdj, "numeric")) {
          betafdParj <- betalist[[j]]
          betafdj    <- betafdParj$fd
          basisj     <- betafdj$basis
          betarangej <- basisj$rangeval
          conbasisj  <- create.constant.basis(betarangej)
          xfdj       <- fd(matrix(xfdj,1,N), conbasisj)
        }
        basisj <- xfdj$basis
        coefj  <- xfdj$coefs
        if (dim(coefj)[1] == 1) coefj <- matrix(coefj[-i],1,N-1)
        else                    coefj <- as.matrix(coefj[,-i])
        xfdlisti[[j]] <- fd(coefj,basisj)
      }
      yveci         <- yvec[-i]
      fRegressListi <- fRegress(yveci, xfdlisti, betalist, wti)
      betaestlisti  <- fRegressListi$betaestlist
      yhati <- 0
      for (j in 1:p) {
        betafdParj <- betaestlisti[[j]]
        betafdj    <- betafdParj$fd
        xfdj       <- xfdlist[[j]]
        bbasisj    <- betafdj$basis
        rangej     <- bbasisj$rangeval
        nfine      <- max(501, bbasisj$nbasis*10+1)
        tfine      <- seq(rangej[1], rangej[2], len=nfine)
        delta      <- tfine[2]-tfine[1]
        betavec    <- eval.fd(tfine, betafdj, 0, returnMatrix)
        xveci      <- eval.fd(tfine, xfdj[i], 0, returnMatrix)
        yhati      <- yhati + delta*(sum(xveci*betavec) -
                                    0.5*( xveci[1]    *betavec[1] +
                                          xveci[nfine]*betavec[nfine] ))
      }
      errfd[i] = yvec[i] - yhati;
      SSE.CV <- SSE.CV + errfd[i]^2
    }
 } else {

    #  functional dependent variable case

    yfd      <- yfdPar$fd
    SSE.CV   <- 0
    errcoefs <- c()
    for(m in 1:length(CVobs)){
      # index of case to eliminate
      i <-  CVobs[m]
      # eliminate case i from the weights
      wti <- wt[-i]
      # eliminate case i from covariates
      txfdlist <- xfdlist
      for(k in 1:p){
        txfdlist[[k]] <- xfdlist[[k]][-i]
      }
      # eliminate case i from dependent variable
      yfdi <- yfd[-i]
      # carry out the functional regression analysis
      tres <- fRegress(yfdi,txfdlist,betalist,wti)
      #  extract the regression coefficient functions
      betaestlisti <- tres$betaestlist
      #  compute the fit to the data for case i
      yhatfdi <- 0
      for(k in 1:p){
        betafdPark = betaestlisti[[k]]
        betafdk    = betafdPark$fd
        xfdk       = xfdlist[[k]]
        xfdik      = xfdk[i]
        tempfd     = xfdik*betafdk
        yhatfdi <- yhatfdi + tempfd
      }
      #  compute the residual function
      errfdi   <- yfd[i] - yhatfdi
      #  increment the error sum of squares by the integral of the
      #  square of the residual function
      SSE.CV   <- SSE.CV + inprod(errfdi,errfdi)
      #  add the coefficients for the residual function
      errcoefs <- cbind(errcoefs,errfdi$coefs)
    }
    #  set up the functional data object for the residual fns
    errfd <- fd(errcoefs,errfdi$basis)
    names(errfd$fdnames)[[3]] <- "Xval Errors"
}
return(list(SSE.CV=SSE.CV,errfd.cv=errfd))
}


