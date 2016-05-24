landmarkreg <- function(fdobj, ximarks, x0marks=xmeanmarks,
                        WfdPar=NULL, monwrd=FALSE, ylambda=1e-10,
                        returnMatrix=FALSE)
{
#  Arguments:
#  FDOBJ   ... functional data object for curves to be registered
#  XIMARKS ... N by NL array of times of interior landmarks for
#                 each observed curve
#  XOMARKS ... vector of length NL of times of interior landmarks for
#                 target curve
#  WFDPAR  ... a functional parameter object defining a warping function.  
#                 If NULL, registration is done using linear interpolation
#                 of lamdmark times in XIMARKS plotted against corresponding 
#                 target times in X0MARKS.
#  MONWRD  ... If TRUE, warping functions are estimated by monotone smoothing,
#                 otherwise by regular smoothing.  The latter is faster, but
#                 not guaranteed to produce a strictly monotone warping
#                 function.  If MONWRD is 0 and an error message results
#                 indicating nonmonotonicity, rerun with MONWRD = 1.
#                 Default:  TRUE
#  YLAMBDA ... smoothing parameter to be used in computing the registered
#                 functions.  For high dimensional bases, local wiggles may be
#                 found in the registered functions or its derivatives that are
#                 not seen in the unregistered functions.  In this event, this
#                 parameter should be increased.
#  Returns:
#  FDREG   ... a functional data object for the registered curves
#  WARPFD  ... a functional data object for the warping functions
#  WFD     ... a functional data object for the W functions defining the
#              warping functions
#  RETURNMATRIX ... If False, a matrix in sparse storage model can be returned
#               from a call to function BsplineS.  See this function for
#               enabling this option.

 #  Last modified 14 November 2012 by Jim Ramsay

  #  check FDOBJ

  if (!(inherits(fdobj,  "fd"))) stop(
		"Argument fdobj  not a functional data object.")

  #  extract information from curve functional data object and its basis

  coef   <- fdobj$coefs
  coefd  <- dim(coef)
  ndim   <- length(coefd)
  ncurve <- coefd[2]
  if (ndim > 2) {
      nvar <- coefd[3]
  } else {
      nvar <- 1
  }

  basisobj  <- fdobj$basis
  nbasis    <- basisobj$nbasis
  rangeval  <- basisobj$rangeval
  fdParobj  <- fdPar(basisobj, 2, ylambda)

  #  check landmarks

  if (is.vector(ximarks) | is.data.frame(ximarks) ) ximarks = as.matrix(ximarks)
  ximarksd <- dim(ximarks)
  if (ximarksd[1] != ncurve) stop(
     "Number of rows of XIMARKS is incorrect.")
    if (any(ximarks <= rangeval[1]) || any(ximarks >= rangeval[2])) stop(
     "Some landmark values are not within the range.")
  nlandm <- dim(ximarks)[2]
  xmeanmarks <- apply(ximarks,2,mean)

  if (length(x0marks) != nlandm) stop(
     "Number of target landmarks not equal to number of curve landmarks.")

  #  set up default WfdPar

  if (is.null(WfdPar)) {
    basisobj  <- fdobj$basis
    rangex    <- basisobj$rangeval
    wnbasis   <- length(x0marks) + 2
    wbasis    <- create.bspline.basis(rangex, wnbasis)
    WfdParobj <- fdPar(wbasis)
  }

  #  check WFDPAR

  WfdPar <- fdParcheck(WfdPar)

  #  set up WFD0 and WBASIS

  Wfd0   <- WfdPar$fd
  wLfd   <- WfdPar$Lfd
  wbasis <- Wfd0$basis

  #  set up WLAMBDA

  wlambda <- WfdPar$lambda

  #  check landmark target values

  wrange <- wbasis$rangeval
  if (any(rangeval != wrange)) stop(
		"Ranges for FD and WFDPAR do not match.")

  #  set up analysis

  n   <- min(c(101,10*nbasis))
  x   <- seq(rangeval[1],rangeval[2],length=n)

  y       <- eval.fd(x, fdobj, returnMatrix=returnMatrix)
  yregmat <- y
  hfunmat <- matrix(0,n,ncurve)
  wlambda <- max(wlambda,1e-10)

  xval    <- c(rangeval[1],x0marks,rangeval[2])
  nwbasis <- wbasis$nbasis
  Wcoef   <- matrix(0,nwbasis,ncurve)
  nval    <- length(xval)
  wval    <- rep(1,nval)

  #  --------------------------------------------------------------------
  #                  Iterate through curves to register
  #  --------------------------------------------------------------------

  cat("Progress:  Each dot is a curve\n")

  for (icurve in 1:ncurve) {
    cat(".")
    #  set up landmark times for this curve
    yval   <- c(rangeval[1],ximarks[icurve,],rangeval[2])
    #  smooth relation between this curve"s values and target"s values
    if (monwrd) {
       #  use monotone smoother
       Wfds      <- smooth.morph(xval, yval, WfdPar)
       Wfd       <- Wfds$Wfdobj
       h         <- monfn(x, Wfd, returnMatrix=returnMatrix)
       b         <- (rangeval[2]-rangeval[1])/(h[n]-h[1])
       a         <- rangeval[1] - b*h[1]
       h         <- a + b*h
       h[c(1,n)] <- rangeval
       wcoefi    <- Wfd$coef
       Wcoef[,icurve] <- wcoefi
    } else {
       #  use unconstrained smoother
       warpfd <- smooth.basis(xval, yval, WfdPar, wval)$fd
       #  set up warping function by evaluating at sampling values
       h         <- as.vector(eval.fd(x, warpfd, returnMatrix=returnMatrix))
       b         <- (rangeval[2]-rangeval[1])/(h[n]-h[1])
       a         <- rangeval[1] - b*h[1]
       h         <- a + b*h
       h[c(1,n)] <- rangeval
       #  check for monotonicity
       deltah <- diff(h)
       if (any(deltah <= 0)) stop(
           paste("Non-increasing warping function estimated for curve",icurve,
                 " Try setting MONWRD to TRUE."))
       wcoefi    <- warpfd$coef
       Wcoef[,icurve] <- wcoefi
    }
    hfunmat[,icurve] <- h

    #  compute h-inverse  in order to register curves

    if (monwrd) {
       wcoef        <- Wfd$coefs
       Wfdinv       <- fd(-wcoef,wbasis)
       WfdParinv    <- fdPar(Wfdinv, wLfd, wlambda)
       Wfdinv       <- smooth.morph(h, x, WfdParinv)$Wfdobj
       hinv         <- monfn(x, Wfdinv, returnMatrix=returnMatrix)
       b            <- (rangeval[2]-rangeval[1])/(hinv[n]-hinv[1])
       a            <- rangeval[1] - b*hinv[1]
       hinv         <- a + b*hinv
       hinv[c(1,n)] <- rangeval
   } else {
       hinvfd       <- smooth.basis(h, x, WfdPar)$fd
       hinv         <- as.vector(eval.fd(x, hinvfd, returnMatrix=returnMatrix))
       b            <- (rangeval[2]-rangeval[1])/(hinv[n]-hinv[1])
       a            <- rangeval[1] - b*hinv[1]
       hinv         <- a + b*hinv
       hinv[c(1,n)] <- rangeval
       deltahinv <- diff(hinv)
       if (any(deltahinv <= 0)) stop(
           paste("Non-increasing warping function estimated for curve",icurve))
    }

    #  compute registered curves

    if (length(dim(coef)) == 2) {
        #  single variable case
        yregfd <- smooth.basis(hinv, y[,icurve], fdParobj)$fd
        yregmat[,icurve] <- eval.fd(x, yregfd, 0, returnMatrix=returnMatrix)
    }
    if (length(dim(coef)) == 3) {
        #  multiple variable case
        for (ivar in 1:nvar) {
            # evaluate curve as a function of h at sampling points
            yregfd <- smooth.basis(hinv, y[,icurve,ivar], fdParobj)$fd
            yregmat[,icurve,ivar] <- eval.fd(x, yregfd,
                                             returnMatrix=returnMatrix)
        }
     }
  }

  cat("\n")

  #  create functional data objects for the registered curves

  fdParobj    <- fdPar(basisobj, 2, ylambda)
  regfdobj    <- smooth.basis(x, yregmat, fdParobj)$fd
  regnames    <- fdobj$fdnames
  names(regnames)[3] <- paste("Registered",names(regnames)[3])
  regfdobj$fdnames <- regnames

  #  create functional data objects for the warping functions

  warpfdobj             <- smooth.basis(x, hfunmat, fdParobj)$fd
  warpfdnames           <- fdobj$fdnames
  names(warpfdnames)[3] <- paste("Warped",names(regnames)[1])
  warpfdobj$fdnames     <- warpfdnames

  Wfd <- fd(Wcoef, wbasis)

  return( list("regfd" = regfdobj, "warpfd" = warpfdobj, "Wfd" = Wfd) )
}
