smooth.monotone <- function(argvals, y, WfdParobj, wtvec=rep(1,n),
                            zmat=NULL, conv=.0001, iterlim=50,
                            active=rep(TRUE,nbasis), dbglev=1,
                            returnMatrix=FALSE)
{
#  Smooths the relationship of Y to ARGVALS using weights in WTVEC by
#  fitting a monotone function of the form
#                   f(x) = b_0 + b_1 D^{-1} exp W(x)
#     where  W  is a function defined over the same range as ARGVALS,
#                 W + ln b_1 = log Df and w = D W = D^2f/Df.
#  The constant term b_0 in turn can be a linear combinations of
#  covariates:
#                         b_0 = zmat * c.
#  The fitting criterion is penalized mean squared error:
#    PENSSE(lambda) = \sum w_i[y_i - f(x_i)]^2 +
#                     \lambda * \int [L W(x)]^2 dx
#  where L is a linear differential operator defined in argument Lfdobj,
#  and w_i is a positive weight applied to the observation.
#  The function W(x) is expanded by the basis in functional data object
#  Wfdobj.   The coefficients of this expansion are called
#  "coefficients" in the comments, while the b's are called "regression
#  coefficients"

#  Arguments:
#  ARGVALS ...  Argument value array of length N, where N is the number
#               of observed curve values for each curve.  It is assumed
#               that these argument values are common to all observed
#               curves.  If this is not the case, you will need to
#               run this function inside one or more loops, smoothing
#               each curve separately.
#  Y       ...  Function value array (the values to be fit).
#               If the functional data are univariate, this array will
#               be an N by NCURVE matrix, where N is the number of
#               observed curve values for each curve and NCURVE is the
#               number of curves observed.
#               If the functional data are muliivariate, this array will
#               be an N by NCURVE by NVAR matrix, where NVAR the number
#               of functions observed per case.  For example, for the
#               gait data, NVAR = 2, since we observe knee and hip
#               angles.
#  WFDPAROBJ... A functional parameter or fdPar object.  This object
#               contains the specifications for the functional data
#               object to be estimated by smoothing the data.  See
#               comment lines in function fdPar for details.
#               The functional data object WFD in WFDPAROBJ is used
#               to initialize the optimization process.
#               Its coefficient array contains the starting values for
#               the iterative minimization of mean squared error.
#  ZMAT    ...  An N by NCOV matrix of covariate values for the constant
#               term.  It defaults to NULL, in this case the constant
#               term is the value of BETA[1] for all values of a given
#               curve.
#  WTVEC   ...  A vector of weights, a vector of N one's by default.
#  CONV    ...  Convergence criterion, 0.0001 by default
#  ITERLIM ...  maximum number of iterations, 50 by default.
#  ACTIVE  ...  indices among 1:NBASIS of parameters to optimize.
#               Defaults to 1:NBASIS.
#  DBGLEV  ...  Controls the level of output on each iteration.  If 0,
#               no output, if 1, output at each iteration, if higher,
#               output at each line search iteration. 1 by default.
#  RETURNMATRIX ... If False, a matrix in sparse storage model can be returned
#               from a call to function BsplineS.  See this function for
#               enabling this option.

#  Returns are:
#  WFD     ...  Functional data object for W.
#               Its coefficient matrix an N by NCURVE (by NVAR) matrix
#               (or array), depending on whether the functional
#               observations are univariate or multivariate.
#  BETA    ...  The regression coefficients b_0 and b_1 for each
#               smoothed curve.
#               If the curves are univariate and
#                  ... ZMAT is NULL,       BETA is 2   by NCURVE.
#                  ... ZMAT has P columns, BETA is P+1 by NCURVE.
#               If the curves are multivariate and
#                  ... ZMAT is NULL,       BETA is 2   by NCURVE by NVAR.
#                  ... ZMAT has P columns, BETA is P+1 by NCURVE by NVAR.
#  YHATFD ...   A functional data object for the monotone curves that
#               smooth the data
#  FLIST  ...   A list object or a vector of list objects, one for
#               each curve (and each variable if functions are
#               multivariate).
#               Each list object has slots:
#                 f    ... The sum of squared errors
#                 grad ... The gradient
#                 norm ... The norm of the gradient
#  Y2CMAP ...   For each estimated curve (and variable if functions are
#               multivariate, this is an N by NBASIS matrix containing
#               a linear mappping from data to coefficients that can be
#               used for computing point-wise confidence intervals.
#               If NCURVE = NVAR = 1, a matrix is returned.  Otherwise
#               an NCURVE by NVAR list is returned, with each
#               slot containing this mapping.
#  When multiple curves and variables are analyzed, the lists containing
#  FLIST and Y2CMAP objects are indexed linear with curves varying
#  inside variables.

# last modified 11 May 2012 by Spencer Graves

#  check ARGVALS

argvals <- argcheck(argvals)
n       <- length(argvals)
onesobs <- matrix(1,n,1)

#  at least three points are necessary for monotone smoothing

if (n < 3) stop('ARGVALS does not contain at least three values.')

#  check Y

ychk   <- ycheck(y, n)
y      <- ychk$y
ncurve <- ychk$ncurve
nvar   <- ychk$nvar
ndim   <- ychk$ndim

#  check WfdParobj and get LAMBDA

WfdParobj <- fdParcheck(WfdParobj)
lambda    <- WfdParobj$lambda

#  the starting values for the coefficients are in FD object WFDOBJ

Wfdobj   <- WfdParobj$fd
Lfdobj   <- WfdParobj$Lfd
basisobj <- Wfdobj$basis     #  basis for W(argvals)
nbasis   <- basisobj$nbasis  #  number of basis functions

#  set up initial coefficient array

coef0    <- Wfdobj$coefs

#  check WTVEC

wtvec <- wtcheck(n, wtvec)$wtvec

#  check ZMAT

if (!is.null(zmat)) {
  zdim <- dim(zmat)
  if (zdim[1] != n) stop("First dimension of ZMAT not correct.")
  ncov   <- zdim[2]   #  number of covariates
} else {
  ncov <- 1
}

#  set up some variables

ncovp1 <- ncov + 1  #  index for regression coef. for monotone fn.
wtroot <- sqrt(wtvec)
wtrtmt <- wtroot %*% matrix(1,1,ncovp1)
yroot  <- y*as.numeric(wtroot)
climit <- c(-100*rep(1,nbasis), 100*rep(1,nbasis))
inact  <- !active   #  indices of inactive coefficients

#  set up list for storing basis function values

JMAX <- 15
basislist <- vector("list", JMAX)

#  initialize matrix Kmat defining penalty term

if (lambda > 0) {
  Kmat <- lambda*eval.penalty(basisobj, Lfdobj)
} else {
  Kmat <- matrix(0,nbasis,nbasis)
}

#  --------------------------------------------------------------------
#              loop through variables and curves
#  --------------------------------------------------------------------

#  set up arrays and lists to contain returned information

if (ndim == 2) {
    coef <- matrix(0,nbasis,ncurve)
    beta <- matrix(0,ncovp1,ncurve)
} else {
    coef <- array(0,c(nbasis,ncurve,nvar))
    beta <- array(0,c(ncovp1,ncurve,nvar))
}

if (ncurve > 1 || nvar > 1 ) {
    Flist <- vector("list",ncurve*nvar)
} else {
    Flist <- NULL
}

if (ncurve > 1 || nvar > 1)  {
    y2cMap <- vector("list",ncurve*nvar)
} else {
    y2cMap <- NULL
}

if (dbglev == 0 && ncurve > 1) cat("Progress:  Each dot is a curve\n")

for (ivar in 1:nvar) {
  for (icurve in 1:ncurve) {
    if (ndim == 2) {
        yi    <- y[,icurve]
        cveci <- coef0[,icurve]
    } else {
        yi    <- y[,icurve,ivar]
        cveci <- coef0[,icurve,ivar]
    }

  #  Compute initial function and gradient values

  result <- fngrad.smooth.monotone(yi, argvals, zmat, wtvec, cveci, lambda,
                                   basisobj, Kmat, inact, basislist, returnMatrix)
  Flisti <- result[[1]]
  betai  <- result[[2]]
  Dyhat  <- result[[3]]

  #  compute the initial expected Hessian

  hessmat <- hesscal.smooth.monotone(betai, Dyhat, wtroot,
                                     lambda, Kmat, inact)

  #  evaluate the initial update vector for correcting the initial cveci

  result   <- linesearch.smooth.monotone(Flisti, hessmat, dbglev)
  deltac   <- result[[1]]
  cosangle <- result[[2]]

  #  initialize iteration status arrays

  iternum <- 0
  status  <- c(iternum, Flisti$f, Flisti$norm, betai)
  if (dbglev >= 1) {
    if (ncurve > 1 || nvar > 1) {
          if (ncurve > 1 && nvar > 1) {
            cat("\n")
            curvetitle <- paste('Results for curve',icurve,'and variable',ivar)
          }
          if (ncurve > 1 && nvar == 1) {
            cat("\n")
            curvetitle <- paste('Results for curve',icurve)
          }
          if (ncurve == 1 && nvar > 1) {
            cat("\n")
            curvetitle <- paste('Results for variable',ivar)
          }
    }
    else curvetitle <- 'Results'

    cat("\n",curvetitle,"\n")
    cat("\nIter.   PENSSE   Grad Length Intercept   Slope\n")
    cat(iternum)
    cat("        ")
    cat(round(status[2],4))
    cat("      ")
    cat(round(status[3],4))
    cat("      ")
    cat(round(betai[1],4))
    cat("      ")
    cat(round(betai[ncovp1],4))
  } else {
    cat(".")
  }

#  -------  Begin iterations  -----------

  MAXSTEPITER <- 10
  MAXSTEP     <- 100
  trial       <- 1
  reset       <- FALSE
  linemat     <- matrix(0,3,5)
  betaold     <- betai
  cvecold     <- cveci
  Foldlist    <- Flisti
  dbgwrd      <- dbglev >= 2

  if (iterlim > 0) {
  for (iter in 1:iterlim) {
      iternum <- iternum + 1
      #  initialize logical variables controlling line search
      dblwrd <- c(FALSE,FALSE)
      limwrd <- FALSE
      stpwrd <- FALSE
      ind    <- 0
      ips    <- 0
      #  compute slope at 0 for line search
      linemat[2,1] <- sum(deltac*Flisti$grad)
      #  normalize search direction vector
      sdg     <- sqrt(sum(deltac^2))
      deltac  <- deltac/sdg
      dgsum   <- sum(deltac)
      linemat[2,1] <- linemat[2,1]/sdg
      # initialize line search vectors
      linemat[,1:4] <- outer(c(0, linemat[2,1], Flisti$f),rep(1,4))
      stepiter <- 0
      if (dbglev >= 2) {
          cat("\n")
          cat(paste("                 ", stepiter, "  "))
          cat(format(round(t(linemat[,1]),6)))
      }
      #  break with error condition if initial slope is nonnegative
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
      for (stepiter in 1:MAXSTEPITER) {
        #  ensure that step does not go beyond limits on parameters
        limflg  <- FALSE
        #  check the step size
        result <-
              stepchk(linemat[1,5], cveci, deltac, limwrd, ind,
                      climit, active, dbgwrd)
        linemat[1,5] <- result[[1]]
        ind          <- result[[2]]
        limwrd       <- result[[3]]
        if (linemat[1,5] <= 1e-7)
        {
          #  Current step size too small ... terminate
          Flisti  <- Foldlist
          cvecnew <- cveci
          gvecnew <- Flisti$grad
          if (dbglev >= 2) {
            print("Stepsize too small")
            print(linemat[1,5])
          }
          if (limflg) ind <- 1 else ind <- 4
          break
        }
        #  compute new function value and gradient
        cvecnew <- cveci + linemat[1,5]*deltac
        result  <- fngrad.smooth.monotone(yi, argvals, zmat, wtvec, cvecnew, lambda,
                                          basisobj, Kmat, inact, basislist,
                                          returnMatrix)
        Flisti  <- result[[1]]
        betai   <- result[[2]]
        Dyhat   <- result[[3]]
        linemat[3,5] <- Flisti$f
        #  compute new directional derivative
        linemat[2,5] <- sum(deltac*Flisti$grad)
        if (dbglev >= 2) {
          cat("\n")
          cat(paste("                 ", stepiter, "  "))
          cat(format(round(t(linemat[,5]),6)))
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
        #  end of line search loop
     }
     cveci  <- cvecnew
     #  check that function value has not increased
     if (Flisti$f > Foldlist$f) {
        # if it has, terminate iterations with a message
        if (dbglev >= 2) {
          cat("Criterion increased: ")
          cat(format(round(c(Foldlist$f, Flisti$f),4)))
          cat("\n")
        }
        #  reset parameters and fit
        betai        <- betaold
        cveci        <- cvecold
        Wfdobj$coefs <- cveci
        Flisti       <- Foldlist
        deltac       <- -Flisti$grad
        if (reset) {
          # This is the second time in a row that
          #  this has happened ... quit
          if (dbglev >= 2) cat("Reset twice, terminating.\n")
          break
        } else {
          reset <- TRUE
        }
     } else {
       if (abs(Foldlist$f - Flisti$f) < conv) {
	       if (dbglev >= 1) cat("\n")
	       break
       }
       cvecold  <- cveci
       betaold  <- betai
       Foldlist <- Flisti
       hessmat  <- hesscal.smooth.monotone(betai, Dyhat, wtroot,
                                           lambda, Kmat, inact)
       #  update the line search direction
       result   <- linesearch.smooth.monotone(Flisti, hessmat, dbglev)
       deltac   <- result[[1]]
       cosangle <- result[[2]]
       reset    <- FALSE
     }
     #  display iteration status
     status <- c(iternum, Flisti$f, Flisti$norm, betai)
     if (dbglev >= 1) {
        cat("\n")
        cat(iternum)
        cat("        ")
        cat(round(status[2],4))
        cat("      ")
        cat(round(status[3],4))
        cat("      ")
        cat(round(betai[1],4))
        cat("      ")
        cat(round(betai[ncovp1],4))
     }
    }
    }

    #  save coefficients in arrays COEF and BETA

    if (ndim == 2) {
        coef[,icurve] <- cveci
        beta[,icurve] <- betai
    } else {
        coef[,icurve,ivar] <- cveci
        beta[,icurve,ivar] <- betai
    }

    #  save Flisti if required in a list,
    #      indexed with curves varying inside variables.

    if (ncurve == 1 && nvar == 1) {
        Flist <- Flisti
    } else {
        Flist[[(ivar-1)*ncurve+icurve]] <- Flisti
    }

    #  save y2cMap if required in a list,
    #      indexed with curves varying inside variables.

    y2cMapij <- solve(crossprod(Dyhat) + lambda*Kmat) %*%
                    t(Dyhat)/sqrt(n)
    if (ncurve == 1 && nvar == 1) {
        y2cMap <- y2cMapij
    } else {
        y2cMap[[(ivar-1)*ncurve+icurve]] <- y2cMapij
    }
  }
}

Wfdobj <- fd(coef, basisobj)

#  Set up yhatfd, a functional data object for the monotone curves
#  fitting the data.
#  This can only be done if the covariate matrix ZMAT is NULL, meaning that
#  the same constant term is used for all curve values.

if (is.null(zmat)) {

  rangeval <- basisobj$rangeval
  narg     <- 10*nbasis+1
  evalarg  <- seq(rangeval[1], rangeval[2], len=narg)
  hmat     <- eval.monfd(evalarg, Wfdobj, 0, returnMatrix)
  if (ndim == 2) {
    yhatmat <- matrix(0,narg,ncurve)
    for (icurve in 1:ncurve) {
      yhatmat[,icurve] <- beta[1,icurve] +
                          beta[2,icurve]*hmat[,icurve]
    }
    yhatcoef <- project.basis(yhatmat, evalarg, basisobj)
    yhatfd   <- fd(yhatcoef, basisobj)
  } else {
    yhatcoef <- array(0,c(nbasis,ncurve,nvar))
    yhatmati <- matrix(0,narg,ncurve)
    for (ivar in 1:nvar) {
      for (icurve in 1:ncurve) {
        yhatmati[,icurve] <- beta[1,icurve,ivar] +
                             beta[2,icurve,ivar]*hmat[,icurve,ivar]
      }
      yhatcoef[,,ivar] <- project.basis(yhatmati, evalarg, basisobj)
    }
    yhatfd <- fd(yhatcoef, basisobj)
  }
} else {
  yhatfd <- NULL
}

monFd <- list( "Wfdobj"  = Wfdobj,  "beta"   = beta, "yhatfd" = yhatfd,
               "Flist"   = Flist,   "y2cMap" = y2cMap,
               "argvals" = argvals, "y"      = y )
class(monFd) <- 'monfd'
monFd
}

#  ----------------------------------------------------------------

linesearch.smooth.monotone <- function(Flisti, hessmat, dbglev)
{
  deltac   <- -symsolve(hessmat,Flisti$grad)
  cosangle <- -sum(Flisti$grad*deltac)/sqrt(sum(Flisti$grad^2)*sum(deltac^2))
  if (dbglev >= 2) {
    cat(paste("\nCos(angle) =",format(round(cosangle,4))))
    if (cosangle < 1e-7) {
      if (dbglev >=2)  cat("\nCosine of angle too small\n")
      deltac <- -Flisti$grad
    }
  }
  return(list(deltac, cosangle))
}

#  ----------------------------------------------------------------

fngrad.smooth.monotone <- function(yi, argvals, zmat, wtvec, cveci, lambda,
                                   basisobj, Kmat, inact, basislist,
                                   returnMatrix=FALSE)
{
  if (!is.null(zmat)) {
    ncov   <- ncol(zmat)
    ncovp1 <- ncov + 1
  } else {
    ncov   <- 1
    ncovp1 <- 2
  }
  n      <- length(argvals)
  nbasis <- basisobj$nbasis
  Wfdobj <- fd(cveci, basisobj)
  h      <- monfn(argvals, Wfdobj, basislist, returnMatrix)
  Dyhat  <- mongrad(argvals, Wfdobj, basislist, returnMatrix)
  if (!is.null(zmat)) {
    xmat <- cbind(zmat,h)
  } else {
    xmat <- cbind(matrix(1,n,1),h)
  }
  Dxmat  <- array(0,c(n,ncovp1,nbasis))
  Dxmat[,ncovp1,] <- Dyhat
  wtroot <- sqrt(wtvec)
  wtrtmt <- wtroot %*% matrix(1,1,ncovp1)
  yroot  <- yi*as.numeric(wtroot)
  xroot  <- xmat*wtrtmt
  #  compute regression coefs.
  betai  <- lsfit(xmat, yi, wt=as.vector(wtvec), intercept=FALSE)$coef
  #  update fitted values
  yhat   <- xmat %*% betai
  #  update residuals and function values
  res    <- yi - yhat
  f      <- mean(res^2*wtvec)
  grad   <- matrix(0,nbasis,1)
  #print(betai)
  for (j in 1:nbasis) {
    Dxroot <- Dxmat[,,j]*wtrtmt
    yDx <- crossprod(yroot,Dxroot) %*% betai
    xDx <- crossprod(xroot,Dxroot)
    #print(crossprod(betai,(xDx+t(xDx))))
    #print(2*yDx)
    grad[j] <- crossprod(betai,(xDx+t(xDx))) %*% betai - 2*yDx
  }
  grad <- grad/n
  if (lambda > 0) {
    grad <- grad +          2 * Kmat %*% cveci
    f    <- f    + t(cveci) %*% Kmat %*% cveci
  }
  if (any(inact)) grad[inact] <- 0
  norm <- sqrt(sum(grad^2)) #  gradient norm
  Flisti <- list("f"=f,"grad"=grad,"norm"=norm)
  return(list(Flisti, betai, Dyhat))
}

#  ----------------------------------------------------------------

hesscal.smooth.monotone <- function(betai, Dyhat, wtroot, lambda,
                                    Kmat, inact)
{
  nbet    <- length(betai)
  Dydim   <- dim(Dyhat)
  n       <- Dydim[1]
  nbasis  <- Dydim[2]
  temp    <- betai[nbet]*Dyhat
  temp    <- temp*(wtroot %*% matrix(1,1,nbasis))
  hessmat <- 2*crossprod(temp)/n
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
