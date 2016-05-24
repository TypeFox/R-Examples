smooth.morph <- function(x, y, WfdPar, wt=rep(1,nobs),
                         conv=.0001, iterlim=20, dbglev=0)
{
#SMOOTH_MORPH smooths the relationship of Y to X
#  by fitting a monotone fn.  f(argvals) = b_0 + b_1 D^{-1} exp W(t)
#     where  W  is a function defined over the same range as ARGVALS,
#  W + ln b_1 = log Df and w = D W = D^2f/Df.
#  b_0 and b_1 are chosen so that f(t_1) = y_1 and f(t_n) = y_n.
#  The fitting criterion is penalized mean squared error:
#    PENSSE(lambda) = \sum [y_i - f(t_i)]^2 +
#                     \lambda * \int [L W]^2
#  W(x) is expanded by the basis in functional data object Wfd.

#  Arguments:
#  X       ...  vector of argument values
#  Y       ...  vector of function values to be fit
#  WFDPAR  ...  functional parameter object for W(x).  The coefficient array
#          for WFDPAROBJ$FD has a single column, and these are the starting
#          values for the iterative minimization of mean squared error.
#  WT      ...  a vector of weights
#  CONV    ...  convergence criterion, 0.0001 by default
#  ITERLIM ...  maximum number of iterations, 20 by default
#  DBGLEV  ...  Controls the level of output on each iteration.  If 0,
#               no output, if 1, output at each iteration, if higher, output
#               at each line search iteration. 1 by default.

#  Returns a list containing:
#  WFDOBJ    ...  functional data object for W(x).  Its coefficients are
#                   those that optimize fit.
#  FLIST     ...  List containing final criterion value, gradient, and
#                 gradient norm.
#  ITERNUM   ...  number of iterations
#  ITERHIST  ...  ITERNUM+1 by 5 array containing iteration history

# last modified 13 may 2012 by Spencer Graves
# previously modified 2 May 2012 by Jim Ramsay

#  initialize some arrays

  x      <- as.vector(x)
  nobs   <- length(x)         #  number of observations

#  check WfdPar

  if (!(inherits(WfdPar, "fdPar"))) stop(
		"Argument WFDPAROBJ is not a functional data object.")

#  extract information from WfdPar

  Wfdobj   <- WfdPar$fd
  Lfdobj   <- WfdPar$Lfd
  lambda   <- WfdPar$lambda

#  check LFDOBJ

  Lfdobj   <- int2Lfd(Lfdobj)

#  extract further information

  basisobj <- Wfdobj$basis     #  basis for W(x)
  nbasis   <- basisobj$nbasis  #  number of basis functions
  type     <- basisobj$type
  rangeval <- basisobj$rangeval
  cvec     <- Wfdobj$coefs
  ncvec    <- length(cvec)

  if (type == "bspline" || type == "fourier") {
      active <- c(FALSE, rep(TRUE, nbasis-1))
  } else {
      active <- rep(TRUE, nbasis)
  }

#  check some arguments

  if (any(diff(x) <= 0)) stop("Arguments are not strictly increasing.")
  if (x[1] < rangeval[1] | x[nobs] > rangeval[2]) stop(
		"Argument values are out of range.")
  if (any(wt < 0))  stop("One or more weights are negative.")
  if (all(wt == 0)) stop("All weights are zero.")

#  set up some variables

  wtroot <- sqrt(wt)
  climit <- c(-100*rep(1,nbasis), 100*rep(1,nbasis))
  inact  <- !active   #  indices of inactive coefficients

#  set up cell for storing basis function values

  JMAX <- 15
  basislist <- vector("list", JMAX)

#  initialize matrix Kmat defining penalty term

  if (lambda > 0) {
      Kmat <- lambda*getbasispenalty(basisobj, Lfdobj)
  } else {
      Kmat <- NULL
  }

#  Compute initial function and gradient values

  fngradlist <- fngrad.smooth.morph(y, x, wt, Wfdobj, lambda,
                                    Kmat, inact, basislist)
  Flist <- fngradlist[[1]]
  Dyhat <- fngradlist[[2]]

#  compute the initial expected Hessian

  hessmat <- hesscal.smooth.morph(Dyhat, wtroot, lambda, Kmat, inact)

#  evaluate the initial update vector for correcting the initial cvec

  result   <- linesearch.smooth.morph(Flist, hessmat, dbglev)
  deltac   <- result[[1]]
  cosangle <- result[[2]]

#  initialize iteration status arrays

  iternum <- 0
  status  <- c(iternum, Flist$f, Flist$norm)
  if (dbglev >= 1) {
    cat("\nIter.   PENSSE   Grad Length")
    cat("\n")
    cat(iternum)
    cat("        ")
    cat(round(status[2],4))
    cat("      ")
    cat(round(status[3],4))
  }
  # if (dbglev == 0 && iterlim > 1)
  #       cat("Progress:  Each dot is an iteration\n")

  iterhist <- matrix(0,iterlim+1,length(status))
  iterhist[1,]  <- status
  if (iterlim == 0)
    return ( list( "Wfdobj" = Wfdobj, "Flist" = Flist,
                 "iternum" = iternum, "iterhist" = iterhist ) )

#  -------  Begin iterations  -----------

  MAXSTEPITER <- 10
  MAXSTEP     <- 100
  trial       <- 1
  reset       <- FALSE
  linemat     <- matrix(0,3,5)
  cvecold     <- cvec
  Foldlist    <- Flist
  dbgwrd      <- dbglev >= 2

  for (iter in 1:iterlim)
  {

      # if (dbglev == 0 && iterlim > 1) cat(".")
      iternum <- iternum + 1
      #  initialize logical variables controlling line search
      dblwrd <- c(0,0)
      limwrd <- c(0,0)
      stpwrd <- 0
      ind    <- 0
      ips    <- 0
      #  compute slope at 0 for line search
      linemat[2,1] <- sum(deltac*Flist$grad)
      #  normalize search direction vector
      sdg     <- sqrt(sum(deltac^2))
      deltac  <- deltac/sdg
      linemat[2,1] <- linemat[2,1]/sdg
      # initialize line search vectors
      linemat[,1:4] <- outer(c(0, linemat[2,1], Flist$f),rep(1,4))
      stepiter <- 0
      if (dbglev >= 2) {
          cat("\n")
          cat(paste("                 ", stepiter, "  "))
          cat(format(round(t(linemat[,1]),6)))
          cat("\n")
      }
      #  return with error condition if initial slope is nonnegative
      if (linemat[2,1] >= 0) {
        if (dbgwrd >= 2) print("Initial slope nonnegative.")
        ind <- 3
        break
      }
      #  return successfully if initial slope is very small
      if (linemat[2,1] >= -1e-7) {
        if (dbglev >= 2) print("Initial slope too small")
        ind <- 0
        break
      }
      #  first step set to trial
      linemat[1,5]  <- trial
      #  Main iteration loop for linesearch
      cvecnew <- cvec
      Wfdnew  <- Wfdobj
      for (stepiter in 1:MAXSTEPITER)
      {
      #  ensure that step does not go beyond limits on parameters
        limflg  <- FALSE
        #  check the step size
        result <-
              stepchk(linemat[1,5], cvec, deltac, limwrd, ind,
                      climit, active, dbgwrd)
        linemat[1,5] <- result[[1]]
        ind          <- result[[2]]
        limwrd       <- result[[3]]
        if (linemat[1,5] <= 1e-7)
        {
          #  Current step size too small ... terminate
          if (dbglev >= 2) {
            print("Stepsize too small")
#            print(avec[5])
            print(linemat[1, 5])
          }
          if (limflg) ind <- 1 else ind <- 4
          break
        }
        #  compute new function value and gradient
        cvecnew <- cvec + linemat[1,5]*deltac
        Wfdnew$coefs <- as.matrix(cvecnew)
        fngradlist   <- fngrad.smooth.morph(y, x, wt, Wfdnew, lambda,
                                          Kmat, inact, basislist)
        Flist <- fngradlist[[1]]
        Dyhat <- fngradlist[[2]]
        linemat[3,5] <- Flist$f
        #  compute new directional derivative
        linemat[2,5] <- sum(deltac*Flist$grad)
        if (dbglev >= 2) {
          cat(paste("                 ", stepiter, "  "))
          cat(format(round(t(linemat[,5]),6)))
          cat("\n")
        }
        #  compute next line search step, also test for convergence
        result  <- stepit(linemat, ips, dblwrd, MAXSTEP)
        linemat <- result[[1]]
        ips     <- result[[2]]
        ind     <- result[[3]]
        dblwrd  <- result[[4]]
        trial   <- linemat[1,5]
        #  ind == 0  mean convergence
        if (ind == 0 | ind == 5) break
     }
     #  end iteration loop
     cvec    <- cvecnew
     Wfdobj  <- Wfdnew
     #  check that function value has not increased
     if (Flist$f > Foldlist$f) {
        # if it has, terminate iterations with a message
        if (dbglev >= 2) {
          cat("Criterion increased: ")
          cat(format(round(c(Foldlist$f, Flist$f),4)))
          cat("\n")
        }
        #  reset parameters and fit
        cvec         <- cvecold
        Wfdobj$coefs <- cvec
        Flist        <- Foldlist
        deltac       <- -Flist$grad
        if (reset) {
          # This is the second time in a row that
          #  this has happened ... quit
          if (dbglev >= 2) cat("Reset twice, terminating.\n")
          return ( list( "Wfdobj" = Wfdobj, "Flist" = Flist,
                         "iternum" = iternum, "iterhist" = iterhist ) )
        } else {
          reset <- TRUE
        }
     } else {
       if (abs(Foldlist$f - Flist$f) < conv) {
	       # cat("\n")
	       break
       }
       cvecold  <- cvec
       Foldlist <- Flist
       hessmat  <- hesscal.smooth.morph(Dyhat, wtroot, lambda, Kmat, inact)
       #  udate the line search direction
       result   <- linesearch.smooth.morph(Flist, hessmat, dbglev)
       deltac   <- result[[1]]
       cosangle <- result[[2]]
       reset    <- FALSE
     }
     #  store iteration status
     status <- c(iternum, Flist$f, Flist$norm)
     iterhist[iter+1,] <- status
     if (dbglev >= 1) {
        cat("\n")
        cat(iternum)
        cat("        ")
        cat(round(status[2],4))
        cat("      ")
        cat(round(status[3],4))
     }
  }
  return ( list( "Wfdobj" = Wfdobj, "Flist" = Flist,
                 "iternum" = iternum, "iterhist" = iterhist ) )
}

#  ----------------------------------------------------------------

linesearch.smooth.morph <- function(Flist, hessmat, dbglev)
{
  deltac   <- -symsolve(hessmat,Flist$grad)
  cosangle <- -sum(Flist$grad*deltac)/sqrt(sum(Flist$grad^2)*sum(deltac^2))
  if (dbglev >= 2) {
    cat(paste("\nCos(angle) =",format(round(cosangle,4))))
    if (cosangle < 1e-7) {
      if (dbglev >=2)  cat("\nCosine of angle too small\n")
      deltac <- -Flist$grad
    }
  }
  return(list(deltac, cosangle))
}

#  ----------------------------------------------------------------

fngrad.smooth.morph <- function(y, x, wt, Wfdobj, lambda,
                                   Kmat, inact, basislist)
{
  nobs   <- length(x)
  width  <- y[nobs] - y[1]
  cvec   <- Wfdobj$coefs
  nbasis <- length(cvec)
  h      <- monfn(x, Wfdobj)
  Dyhat  <- mongrad(x, Wfdobj)
  #  adjust h and Dyhat for normalization
  hmax   <- h[nobs]
  Dymax  <- Dyhat[nobs,]
  Dyhat  <- width*(hmax*Dyhat - outer(h, Dymax))/hmax^2
  h      <- y[1] + width*h/hmax
  #  update residuals and function values
  res    <- y - h
  f      <- mean(res^2*wt)
  temp   <- Dyhat*outer(wt,rep(1,nbasis))
  grad   <- -2*crossprod(temp,res)/nobs
  if (lambda > 0) {
    grad <- grad +         2 * Kmat %*% cvec
    f    <- f    + t(cvec) %*% Kmat %*% cvec
  }
  if (any(inact)) grad[inact] <- 0
  norm <- sqrt(sum(grad^2)) #  gradient norm
  Flist <- list("f"=f,"grad"=grad,"norm"=norm)
  return(list(Flist, Dyhat))
}

#  ----------------------------------------------------------------

hesscal.smooth.morph <- function(Dyhat, wtroot, lambda, Kmat, inact)
{
  Dydim   <- dim(Dyhat)
  nobs    <- Dydim[1]
  nbasis  <- Dydim[2]
  temp    <- Dyhat*(wtroot %*% matrix(1,1,nbasis))
  hessmat <- 2*crossprod(temp)/nobs
  #  adjust for penalty
  if (lambda > 0) hessmat <- hessmat + 2*Kmat
  #  adjust for inactive coefficients
  if (any(inact)) {
    eyemat               <- diag(rep(1,nbasis))
    hessmat[inact,     ] <- 0
    hessmat[     ,inact] <- 0
    hessmat[inact,inact] <- eyemat[inact,inact]
  }
  return(hessmat)
}
