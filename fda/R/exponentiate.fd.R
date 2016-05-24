#  -------------------------------------------------------------------
#  power method for "fd"
#  -------------------------------------------------------------------

"^.fd" <- function(e1, e2){
#  A positive integer pointwise power of a functional data object with
#  a B-splinebasis.  powerfd = fdobj^a.
#  Generic arguments e1 = fdobj and e2 = a.
#
#  The basis is tested for being a B-spline basis.  The function then
#  sets up a new spline basis with the same knots but with an order
#  that is M-1 higher than the basis for FDOBJ, where M = ceiling(a),
#  so that the order of differentiability for the new basis is
#  appropriate in the event that a is a positive integer, and also
#  to accommodate the additional curvature arising from taking a power.
#  The power of the values of the function over a fine mesh are computed,
#  and these are fit using the new basis.
#
#  Powers should be requested with caution, however, and especially if
#  a < 1, because, if there is strong local curvature in FDOBJ,
#  and if its basis is just barely adequate to capture this curvature,
#  then the power of the function may have considerable error
#  over this local area.  fdobj^a where a is close to zero is just
#  such a situation.
#
#  If a power of a functional data object is required for which the
#  basis is not a spline, it is better to either re-represent the
#  function in a spline basis, or, perhaps even better, to do the
#  math required to get the right basis and interpolate function
#  values over a suitable mesh.  This is especially true for fourier
#  bases.

#  Last modified 2012.07.01 by Spencer Graves
#  previously modified 3 November 2009

#  check first two arguments

  fdobj = e1
  a     = e2
  tol   = 1e-4

  if ((!(inherits(fdobj, "fd"))))
        stop("First argument for ^ is not a functional data object.")
  if ((!(is.numeric(a))))
        stop("Second argument for ^ is not numeric.")

#  extract basis

  basisobj = fdobj$basis

#  test the basis for being of B-spline type

  if (basisobj$type != "bspline"){
    e12 <- exponentiate.fd(e1, e2,
       tolint=.Machine$double.eps^0.75, basisobj=e1$basis,
       tolfd=sqrt(.Machine$double.eps)*
          sqrt(sum(e1$coefs^2)+.Machine$double.eps)^abs(e2),
       maxbasis=NULL, npoints=NULL)
    return(e12)
#    stop("FDOBJ does not have a spline basis.")
#    a1 <- round(a)
#    if(abs(a-a1)>.Machine$double.eps^.75)
#        stop('Fractional power not allowed of an fd object ',
#             'without a spline basis.')
#    if(a1<0)
#        stop('Negative powers not allowed of an fd object ',
#             'without a spline basis.')
#    if(a1==0){
#      if(a==0){
#          rng <- basisobj$rangeval
#          fdNames <- list(fdobj$fdnames$args, NULL, NULL)
#          fdout <- fd(1, const, fdNames)
#          return(fdout)
#      } else {
#          stop('very small nonzero powers not allowed of an fd object ',
#               'without a spline basis;  requested power = ', a)
#      }
#    }
#    fdout <- fdobj
#    for(i in seq(length=a1-1)) fdout <- fdout*fdobj
#    return(fdout)
  }

  nbasis        = basisobj$nbasis
  rangeval      = basisobj$rangeval
  interiorknots = basisobj$params
  norder        = nbasis - length(interiorknots)

#  Number of points at which to evaluate the power.  Even low
#  order bases can generate steep slopes and sharp curvatures,
#  especially if powers less than 1 are involved.

  nmesh = max(10*nbasis+1,501)

#  determine number curves and variables

  coefmat = fdobj$coef
  coefd   = dim(coefmat)
  ncurve  = coefd[2]
  if (length(coefd) == 2) {
    nvar = 1
  } else {
    nvar = coefd[3]
  }

#  evaluate function over this mesh

  tval = seq(rangeval[1],rangeval[2],len=nmesh)
  fmat = eval.fd(tval, fdobj)
  fmatNeg <- (fmat<0)
# eliminate negatives from roundoff only
  if(any(fmatNeg) && all(fmat[fmatNeg]>(-.Machine$double.eps)))
      fmat[fmatNeg] <- 0

#  find the minimum value over this mesh.  If the power is less than
#  one, return an error message.

  fmin = min(c(fmat))

#  a == 0:  set up a constant basis and return the unit function(s)

  if (a == 0) {
    newbasis = create.constant.basis(rangeval)
    if (nvar == 1) {
      powerfd = fd(matrix(1,1,ncurve), newbasis)
    } else {
      powerfd = fd(array(1,c(1,ncurve,nvar)), newbasis)
    }
    return(powerfd)
  }

#  a == 1:  return the function

  if (a == 1) {
    powerfd = fdobj
    return(powerfd)
  }

#  Otherwise:

  m = ceiling(a)

#  Check the size of the power.  If greater than one, estimating the
#  functional data object is relatively safe since the curvatures
#  involved are mild.  If not, then taking the power is a dangerous
#  business.

  if (m == a && m > 1) {

    #  a is an integer greater than one

    newnorder = (norder-1)*m + 1
    if (length(interiorknots) < 9) {
        newbreaks = seq(rangeval[1], rangeval[2], len=11)
    } else {
        newbreaks = c(rangeval[1], interiorknots, rangeval[2])
    }
    nbreaks   = length(newbreaks)
    newnbasis = newnorder + nbreaks - 2
    newbasis  = create.bspline.basis(rangeval, newnbasis, newnorder,
                                     newbreaks)
    ymat    = fmat^a
    ytol    = max(abs(c(ymat)))*tol
    powerfd = smooth.basis(tval, ymat, newbasis)$fd
    ymathat = eval.fd(tval,powerfd)
    ymatres = ymat - ymathat
    maxerr  = max(abs(c(ymatres)))
    while  (maxerr > ytol && nbreaks < nmesh) {
        newnbasis = newnorder + nbreaks - 2
        newbasis  = create.bspline.basis(rangeval, newnbasis, newnorder,
                                         newbreaks)
        newfdPar  = fdPar(newbasis, 2, 1e-20)
        powerfd   = smooth.basis(tval, ymat, newfdPar)$fd
        ymathat   = eval.fd(tval,powerfd)
        ymatres   = ymat - ymathat
        maxerr    = max(abs(c(ymatres)))
        if (nbreaks*2 <= nmesh) {
        newbreaks = sort(c(newbreaks,
            (newbreaks[1:(nbreaks-1)]+newbreaks[2:nbreaks])/2))
        } else {
            newbreaks = tval
        }
        nbreaks = length(newbreaks)
    }

    if (maxerr > ytol)
        warning("The maximum error exceeds the tolerance level.")

    return(powerfd)

  } else {

    #  a is fractional or negative

    #  check for negative values and a fractional power

    if (a > 0 && fmin < 0) stop(
         paste("There are negative values",
               "and the power is a positive fraction."))

    #  check for zero or negative values and a negative power

    if (a < 0 && fmin <= 0) stop(
        paste("There are zero or negative values",
               "and the power is negative."))

    if (length(interiorknots) < 9) {
        newbreaks = seq(rangeval[1], rangeval[2], n=11)
    } else {
        newbreaks = c(rangeval[1], interiorknots, rangeval[2])
    }
    nbreaks   = length(newbreaks)
    newnorder = max(4, norder+m-1)
    newnbasis = newnorder + nbreaks - 2
    newbasis  = create.bspline.basis(rangeval, newnbasis, newnorder,
                                     newbreaks)
    nmesh     = max(10*nbasis+1,101)
    tval      = seq(rangeval[1],rangeval[2],len=nmesh)
    fmat      = eval.fd(tval, fdobj)
    fmatNeg <- (fmat<0)
    if(any(fmatNeg) && all(fmat[fmatNeg]>(-.Machine$double.eps))){
        fmat[fmatNeg] <- 0
    }
#
    ymat      = fmat^a
    ytol      = max(abs(c(ymat)))*tol
    newfdPar  = fdPar(newbasis, 2, 1e-20)
    powerfd   = smooth.basis(tval, ymat, newfdPar)$fd
    ymathat   = eval.fd(tval,powerfd)
    ymatres   = ymat - ymathat
    maxerr    = max(abs(c(ymatres)))
    while (maxerr > ytol && nbreaks < nmesh) {
        newnbasis = newnorder + nbreaks - 2
        newbasis  = create.bspline.basis(rangeval, newnbasis,
                                         newnorder, newbreaks)
        newfdPar  = fdPar(newbasis, 2, 1e-20)
        powerfd   = smooth.basis(tval, ymat, newfdPar)$fd
        ymathat   = eval.fd(tval,powerfd)
        ymatres   = ymat - ymathat
        maxerr    = max(abs(ymatres))*tol
        if (nbreaks*2 <= nmesh) {
            newbreaks = sort(c(newbreaks,
            (newbreaks[1:(nbreaks-1)]+newbreaks[2:nbreaks])/2))
        } else {
            newbreaks = tval
        }
        nbreaks = length(newbreaks)
    }

    if (maxerr > ytol) {
        warning("The maximum error exceeds the tolerance level.")
    }

    return(powerfd)

  }

}

#  -------------------------------------------------------------------
#  sqrt method for "fd"
#  -------------------------------------------------------------------

sqrt.fd <- function(x)
{
#  Arguments:
#  x ...  A functional data object
#  Returns:
#  FDAROOT  ...  A functional data object that is the square root of x
#  Last modified:  27 october 2009

    if ((!(inherits(x, "fd")))) stop(
      "Argument for ^ is not a functional data object.")

    fdaroot <- x^0.5

    return(fdaroot)
}

exponentiate.fd <- function(e1, e2,
       tolint=.Machine$double.eps^0.75, basisobj=e1$basis,
       tolfd=sqrt(.Machine$double.eps)*
          sqrt(sum(e1$coefs^2)+.Machine$double.eps)^abs(e2),
       maxbasis=NULL, npoints=NULL){
##
## e2=0?
##
  e1basis <- e1$basis
  rng <- e1basis$rangeval
  coefmat <- e1$coef
  coefd <- dim(coefmat)
  ncurve <- coefd[2]
  if(length(coefd)==2){
      nvar <- 1
  } else nvar <- coefd[3]
#
  if(e2==0){
    const <- create.constant.basis(rng)
    fdNames <- list(e1$fdnames$args, NULL, NULL)
    if(nvar==1){
        fdout <- fd(matrix(1, 1, ncurve), const, fdNames)
    } else {
        fdout <- fd(array(1, c(1, ncurve, nvar)), const, fdNames)
    }
    return(fdout)
  }
##
## e2<0?
##
  if(e2<0){
    if(e1$type=='bspline')return(e1^e2)
    stop('Negative powers not allowed of fd objects without ',
         'a spline basis;  power = ', e2)
  }
##
## Fourier basis
##
  if(e1basis$type=='fourier'){
    e2. <- floor(e2)
    if(e2.<1){
      const <- create.constant.basis(rng)
      fdNames <- list(e1$fdnames$args, NULL, NULL)
      if(nvar==1){
        e120 <- fd(matrix(1, 1, ncurve), const, fdNames)
      } else {
        e120 <- fd(array(1, c(1, ncurve, nvar)), const, fdNames)
      }
      outbasis <- e1basis
    } else {
      e120 <- e1
      for(i in seq(length(e2.-1))) e120 <- e120*e1
      outbasis <- e120$basis
    }
    if((e2-e2.)==0)return(e120)
#
    rng <- outbasis$rangeval
    if(is.null(maxbasis)) maxbasis <- 2*outbasis$nbasis+1
    if(is.null(npoints))
        npoints <- max(10*maxbasis+1, 501)
#
    Time <- seq(rng[1], rng[2], length=npoints)
    e1. <- predict(e1, Time)
    e1.2 <- e1.^e2
    n.na <- sum(is.na(e1.2))
    if(n.na>0){
      stop('NAs generated in computing e1^e2 at ', n.na,
           ' of ', npoints, ' sample points.')
    }
#    fd1.2 <- Data2fd(Time, e1.2, outbasis)
    fd1.2 <- smooth.basis(Time, e1.2, outbasis)
    d1.2 <- (e1.2 - predict(fd1.2, Time))
    if(all(abs(d1.2)<tolfd))return(fd1.2$fd)
#
    morebases <- (maxbasis-outbasis$nbasis+1)
    for(i in 1:morebases){
      nbasisi <- outbasis$nbasis+2
      if (nbasisi == (nbasisi/2)*2) {
        nbasisi = nbasisi + 1;
      }
      outbasis <- create.fourier.basis(rng,
          nbasis=nbasisi, period=diff(rng) )
#      fd1.2 <- Data2fd(Time, e1.2, npoints)
      fd1.2 <- smooth.basis(Time, e1.2, outbasis)
      d1.2 <- (e1.2-predict(fd1.2, Time))
      maxd1.2 <- max(abs(d1.2))
      if(maxd1.2<tolfd) return(fd1.2$fd)
    }
    e1name <- deparse(substitute(e1))
    e2name <- deparse(substitute(e2))
    warning('Lack of precision in ', e1name, '^', e2name,
            ': max error at ', npoints, ' sample points = ',
            maxd1.2)
    return(fd1.2$fd)
  }
##
## positive integer
##
  if(abs(e2%%1)<=tolint){
    e2. <- round(e2)
    if(e2.==0)
      stop('powers near zero of fd objects without ',
           'a spline basis not allowed;  power = ', e2)
    fdout <- e1
    for(i in seq(length=e2.-1)) fdout <- fdout*e1
    return(fdout)
  }
}
