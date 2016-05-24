register.fd <- function(y0fd=NULL, yfd=NULL, WfdParobj=NULL,
                    conv=1e-4, iterlim=20, dbglev=1, periodic=FALSE, crit=2,
                    returnMatrix=FALSE)
{
#REGISTERFD registers a set of curves YFD to a target function Y0FD.
#  Arguments are:
#  Y0FD      ... Functional data object for target function.  It may be either
#                a single curve, or have the same dimensions as YFD.
#  YFD       ... Functional data object for functions to be registered
#  WFDPAROBJ ... Functional parameter object for function W defining warping
#                functions. The basis that is defined in WFDPAROBJ must be a
#                B-spline basis, and the number of basis functions must be at
#                least 2.
#                The coefficients, which supplied either in a functional
#                data object or as defaults to zero if a basis object is
#                supplied in the definition of WFDPAROBJ are used as the
#                initial values in the iterative computation of the final
#                warping functions.
#                NB:  The value of the first coefficient is NOT used.
#                This is because a warping function is normalized, and when
#                this happens, the impact of the first coefficient, if used,
#                would be eliminated.  This first position is used, however, to
#                contain the shift parameter in case the data are to be
#                treated as periodic.  At the end of the calculations,
#                the shift parameter is returned separately.
#                If WFDPAROBJ is not supplied, it defaults to a bspline
#                basis of order 2 with 2 basis functions.  This is equivalent
#                to using a linear function for W.
#  CONV    ... Convergence criterion
#  ITERLIM ... iteration limit for scoring iterations
#  DBGLEV  ... Level of output of computation history
#  PERIODIC... If one, curves are periodic and a shift parameter is fit.
#              Initial value for shift parameter is taken to be 0.
#              The periodic option should ONLY be used with a Fourier
#              basis for the target function Y0FD, the functions to be
#              registered, YFD.
#  CRIT    ... if 1 least squares,
#              if 2 log eigenvalue ratio,
#              if 3 least squares using residual y0 - yi*dhdt
#              Default is 2.
#  RETURNMATRIX ... If False, a matrix in sparse storage model can be returned
#               from a call to function BsplineS.  See this function for
#               enabling this option.
#                Default:  0
#  Returns:

#  REGLIST ...  A list with fields:
#    REGLIST$REGFD  ... A functional data object for the registered curves
#    REGLIST$WARPFD ... A Functional data object for warping functions h
#    REGLIST$WFD    ... A Functional data object for functions W defining
#                         warping fns
#    REGLIST$SHIFT  ... Shift parameter value if curves are periodic
#    REGLIST$Y0FD   ... Argument Y0FD
#    REGLIST$YFD    ... Argument YFD

#  Last modified 14 November 2012 by Jim Ramsay

##
## 1.  Check y0fd and yfd
##

#  if YFD is not supplied, and therefore defaults to NULL,  it is
#  assumed that the first argument is to be taken as YFD, and that
#  the target function Y0FD defaults to NULL.  In this event,  Y0FD
#  is set up as the mean of the functions in YFD.

  if(is.null(yfd)){
    yfd  <- y0fd
    y0fd <- NULL
  }

#  Check that YFD is a functional data object

  if (!(inherits(yfd, "fd")))
      stop("'yfd' must be a functional data object.  ",
           "Instead, class(yfd) = ", class(yfd))

#  If target Y0FD is not supplied, and therefore defaults to NULL, replace
#  it be the mean of the functions in YFD

  if(is.null(y0fd)) {
    y0fd <- mean.fd(yfd)
  } else {
    if (!(inherits(y0fd, "fd")))
        stop("First argument is not a functional data object.",
             "Instead, class(y0fd) = ", class(y0fd))
  }

#  get dimensions of functions in YFD to be registered

  ycoefs <- yfd$coefs
  if (is.vector(ycoefs)) ycoefs = as.matrix(ycoefs)
  ydim   <- dim(ycoefs)
  ncurve <- ydim[2]
  ndimy  <- length(ydim)
  if (ndimy == 3) {
      nvar <- ydim[3]
  } else {
      nvar <- 1
  }
  if (ndimy > 3) stop("'yfd' is more than 3-dimensional.")

#  Extract basis information from YFD

  ybasis  <- yfd$basis
  ynbasis <- ybasis$nbasis
  yrange  <- ybasis$rangeval
  if (periodic && !(ybasis$type == "fourier"))
    stop("'periodic' is TRUE but 'type' is not 'fourier'; ",
         "periodic B-splines are not currently part of 'fda'")

#  Get dimensions of target function object Y0FD

  y0coefs0 <- y0fd$coefs
  if (is.vector(y0coefs0)) y0coefs0 = as.matrix(y0coefs0)
  y0dim0   <- dim(y0coefs0)
  ndimy00  <- length(y0dim0)
  if (ndimy00 > ndimy) stop("Y0FD has more dimensions than YFD")
  #  Determine whether the target function is full or not
  if (y0dim0[2] == 1) {
      fulltarg <- FALSE
  } else {
      if (y0dim0[2] == ydim[2]) {
          fulltarg <- TRUE
      } else {
          stop(paste("Second dimension of coefficient matrix for Y0FD",
                     "is neither 1 nor equal to the number of functions",
                     "to be registered."))
      }
  }
  if (ndimy00 == 3 && ydim[3] != y0dim0[3]) stop(
      "Third dimension of YOFD does not match that of YFD.")

#  Extract basis information from Y0FD

  y0basis  <- y0fd$basis
  y0nbasis <- y0basis$nbasis
  y0range  <- y0basis$rangeval
  #  check that target range matches function range
  if (!all(y0range == yrange)) stop(
      "Range for Y0FD does not match range for YFD.")

##
## 2.  Check WfdParobj
##

  if (is.null(WfdParobj)) {
      #  default WfdParobj to a B-spline basis of order 2 with 2 basis functions
      wbasis    <- create.bspline.basis(yrange,2,2)
      Wfd0      <- fd(matrix(0,2,ncurve),wbasis)
      WfdParobj <- fdPar(Wfd0)
  }
  WfdParobj <- fdParcheck(WfdParobj)
  Wfd0   <- WfdParobj$fd
  wcoef  <- Wfd0$coefs
  if (is.vector(wcoef)) wcoef <- as.matrix(wcoef)
  wbasis <- Wfd0$basis
  wtype  <- wbasis$type
  if (wtype != "bspline") stop("Basis for Wfd is not a B-spline basis.")
  wnbasis <- wbasis$nbasis
  norder  <- wnbasis - length(wbasis$params)
  if (wnbasis < 2) stop(
      "At least two basis functions for W are required.")
  if (norder < 2) stop(
      "The order of the basis functions for W must be at least 2.")
  wtype  <- wbasis$type
  rangex <- wbasis$rangeval
  wdim   <- dim(wcoef)
  if (length(wdim) > 2) stop("WFDPAROBJ contains a multivariate function.")
  if (wdim[2] == 1) {
      wcoef      <- wcoef %*% matrix(1,1,ncurve)
      Wfd0$coefs <- wcoef
  } else {
      if (wdim[2] != ncurve) stop(
          "WFDPAROBJ and YFD containing differing numbers of functions.")
  }


##
## 3.  Do the work
##

#  set up a fine mesh of argument values

NFINEMIN <- 201
nfine <- 10*ynbasis + 1
if (nfine < NFINEMIN) nfine <- NFINEMIN
xlo   <- rangex[1]
xhi   <- rangex[2]
width <- xhi - xlo
xfine <- seq(xlo, xhi, len=nfine)

#  set up indices of coefficients that will be modified in ACTIVE

wcoef1   <- wcoef[1,]
if (periodic) {
   active   <- 1:wnbasis
   wcoef[1] <- 0
   shift    <- 0
} else {
   active <- 2:wnbasis
}

#  initialize matrix Kmat defining penalty term

lambda <- WfdParobj$lambda
if (lambda > 0) {
   Lfdobj <- WfdParobj$Lfd
   Kmat <- getbasispenalty(wbasis, Lfdobj)
   ind  <- 2:wnbasis
   Kmat <- lambda*Kmat[ind,ind]
} else {
   Kmat <- NULL
}

#  set up limits on coefficient sizes

climit <- 50*c(-rep(1,wnbasis), rep(1,wnbasis))

#  set up cell for storing basis function values

JMAX <- 15
basislist <- vector("list", JMAX)

yregcoef <- yfd$coefs

#  loop through the curves

wcoefnew <- wcoef
if (dbglev == 0 && ncurve > 1) cat("Progress:  Each dot is a curve\n")

for (icurve in 1:ncurve) {
  if (dbglev == 0 && ncurve > 1) cat(".")
  if (dbglev >= 1 && ncurve > 1)
      cat(paste("\n\n-------  Curve ",icurve,"  --------\n"))
  if (ncurve == 1) {
    yfdi  <- yfd
    y0fdi <- y0fd
    Wfdi  <- Wfd0
    cvec  <- wcoef
  } else {
    Wfdi <- Wfd0[icurve]
    cvec <- wcoef[,icurve]
    if (nvar == 1) {
      yfdi <- yfd[icurve]
    } else {
      yfdi <- yfd[icurve]
      yfdi$coef <- array(yfdi$coef,c(ynbasis,1,nvar))
    }
    if (fulltarg) {
      if (nvar == 1) {
        y0fdi <- y0fd[icurve]
      } else {
        y0fdi <- y0fd[icurve,]
      }
    } else {
      y0fdi <- y0fd
    }
  }

  #  evaluate curve to be registered at fine mesh

  yfine  <- matrix(eval.fd(xfine, yfdi, 0, returnMatrix),nfine,nvar)

  #  evaluate target curve at fine mesh

  y0fine <- matrix(eval.fd(xfine, y0fdi, 0, returnMatrix),nfine,nvar)

  #  evaluate objective function for starting coefficients

  #  first evaluate warping function and its derivative at fine mesh

  ffine  <-   monfn(xfine, Wfdi, basislist, returnMatrix)
  Dffine <- mongrad(xfine, Wfdi, basislist, returnMatrix)
  fmax   <- ffine[nfine]
  Dfmax  <- Dffine[nfine,]
  hfine  <- xlo + width*ffine/fmax
  Dhfine <- width*(fmax*Dffine - outer(ffine,Dfmax))/fmax^2
  hfine[1]     <- xlo
  hfine[nfine] <- xhi

  #  register curves given current Wfdi

  yregfdi <- regyfn(xfine, yfine, hfine, yfdi, Wfdi, periodic)

  #  compute initial criterion value and gradient

  Flist <- regfngrad(xfine, y0fine, Dhfine, yregfdi, Wfdi,
                     Kmat, periodic, crit, returnMatrix)

  #  compute the initial expected Hessian

  if (crit == 2) {
     D2hwrtc <- monhess(xfine, Wfdi, basislist, returnMatrix)
     D2fmax  <- D2hwrtc[nfine,]
     fmax2 <- fmax*fmax
     fmax3 <- fmax*fmax2
     m <- 1
     if (wnbasis > 1) {
        for (j in 2:wnbasis) {
           m <- m + 1
           for (k in 2:j) {
              m <- m + 1
              D2hwrtc[,m] <- width*(2*ffine*Dfmax[j]*Dfmax[k]
                   - fmax*(Dffine[,j]*Dfmax[k] + Dffine[,k]*Dfmax[j])
                   + fmax2*D2hwrtc[,m] - ffine*fmax*D2fmax[m])/fmax3
           }
        }
     }
  } else {
     D2hwrtc <- NULL
  }

  hessmat <- reghess(xfine, y0fine, Dhfine, D2hwrtc, yregfdi,
                     Kmat, periodic, crit, returnMatrix)

  #  evaluate the initial update vector for correcting the initial cvec

  result   <- linesearch(Flist, hessmat, dbglev)
  deltac   <- result[[1]]
  cosangle <- result[[2]]
  #  initialize iteration status arrays

  iternum <- 0
  status <- c(iternum, Flist$f, Flist$norm)
  if (dbglev >= 1) {
        cat("\nIter.    Criterion   Grad Length")
        cat("\n")
        cat(iternum)
        cat("        ")
        cat(round(status[2],4))
        cat("      ")
        cat(round(status[3],4))
  }
  iterhist <- matrix(0,iterlim+1,length(status))
  iterhist[1,]  <- status
  if (iterlim == 0) break

  #  -------  Begin main iterations  -----------

  MAXSTEPITER <- 5
  MAXSTEP <- 100
  trial   <- 1
  reset   <- 0
  linemat <- matrix(0,3,5)
  cvecold <- cvec
  Foldlist <- Flist
  dbgwrd  <- dbglev >= 2
  #  ---------------  beginning of optimization loop  -----------
  for (iter in 1:iterlim) {
      iternum <- iternum + 1
      #  set logical parameters
      dblwrd <- c(FALSE,FALSE)
      limwrd <- c(FALSE,FALSE)
      ind <- 0
      ips <- 0
      #  compute slope
      linemat[2,1] <- sum(deltac*Foldlist$grad)
      #  normalize search direction vector
      sdg          <- sqrt(sum(deltac^2))
      deltac       <- deltac/sdg
      linemat[2,1] <- linemat[2,1]/sdg
      # initialize line search vectors
      linemat[,1:4] <- outer(c(0, linemat[2,1], Flist$f),rep(1,4))
      stepiter  <- 0
      if (dbglev >= 2) {
          cat("\n")
          cat(paste("                 ", stepiter, "  "))
          cat(format(round(t(linemat[,1]),4)))
      }
      #  return with stop condition if initial slope is nonnegative
      if (linemat[2,1] >= 0) {
        if (dbglev >= 2) cat("\nInitial slope nonnegative.")
        ind <- 3
        break
      }
      #  return successfully if initial slope is very small
      if (linemat[2,1] >= -min(c(1e-3,conv))) {
        if (dbglev >= 2) cat("\nInitial slope too small")
        ind <- 0
        break
      }
      #  first step set to trial
      linemat[1,5]  <- trial
      #  ------------  begin line search iteration loop  ----------
      cvecnew <- cvec
      Wfdnewi <- Wfdi
      for (stepiter in 1:MAXSTEPITER) {
        #  check the step size and modify if limits exceeded
        result <- stepchk(linemat[1,5], cvec, deltac, limwrd, ind,
                   climit, active, dbgwrd)
        linemat[1,5] <- result[[1]]
        ind          <- result[[2]]
        limwrd       <- result[[3]]
        if (ind == 1) break    # break of limit hit twice in a row
        if (linemat[1,5] <= 1e-7) {
           #  Current step size too small  terminate
           if (dbglev >= 2)
             cat("\nStepsize too small: ", round(linemat[1,5],4))
           break
        }
        #  update parameter vector
        cvecnew <- cvec + linemat[1,5]*deltac
        #  compute new function value and gradient
        Wfdnewi[[1]] <- cvecnew
        #  first evaluate warping function and its derivative at fine mesh
        cvectmp <- cvecnew
        cvectmp[1] <- 0
        Wfdtmpi <- Wfdnewi
        Wfdtmpi[[1]] <- cvectmp
        ffine  <-   monfn(xfine, Wfdtmpi, basislist, returnMatrix)
        Dffine <- mongrad(xfine, Wfdtmpi, basislist, returnMatrix)
        fmax   <- ffine[nfine]
        Dfmax  <- Dffine[nfine,]
        hfine  <- xlo + width*ffine/fmax
        Dhfine <- width*(fmax*Dffine - outer(ffine,Dfmax))/fmax^2
        hfine[1]     <- xlo
        hfine[nfine] <- xhi
        #  register curves given current Wfdi
        yregfdi <- regyfn(xfine, yfine, hfine, yfdi, Wfdnewi, periodic)
        Flist    <- regfngrad(xfine, y0fine, Dhfine, yregfdi, Wfdnewi,
                             Kmat, periodic, crit, returnMatrix)
        linemat[3,5] <- Flist$f
        #  compute new directional derivative
        linemat[2,5] <- sum(deltac*Flist$grad)
        if (dbglev >= 2) {
          cat("\n")
          cat(paste("                 ", stepiter, "  "))
          cat(format(round(t(linemat[,5]),4)))
        }
        #  compute next line search step, also testing for convergence
        result  <- stepit(linemat, ips, dblwrd, MAXSTEP)
        linemat <- result[[1]]
        ips     <- result[[2]]
        ind     <- result[[3]]
        dblwrd  <- result[[4]]
        trial   <- linemat[1,5]
        #  ind == 0 implies convergence
        if (ind == 0 || ind == 5) break
     }
     #  ------------  end line search iteration loop  ----------
     cvec   <- cvecnew
     Wfdi   <- Wfdnewi
     #  test for function value made worse
     if (Flist$f > Foldlist$f) {
        #  Function value worse  warn and terminate
        ier <- 1
        if (dbglev >= 2) {
          cat("Criterion increased, terminating iterations.\n")
          cat(paste("\n",round(c(Foldlist$f, Flist$f),4)))
        }
        #  reset parameters and fit
        cvec   <- cvecold
        Wfdi[[1]] <- cvecold
        Flist   <- Foldlist
        deltac <- -Flist$grad
        if (dbglev > 2) {
          for (i in 1:wnbasis) cat(cvec[i])
          cat("\n")
        }
        if (reset == 1) {
           #  This is the second time in a row that this
           #     has happened   quit
           if (dbglev >= 2) cat("Reset twice, terminating.\n")
            break
        } else {
           reset <- 1
        }
     } else {
        #  function value has not increased,  check for convergence
        if (abs(Foldlist$f-Flist$f) < conv) {
           wcoef[,icurve]    <- cvec
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
           break
        }
        #  update old parameter vectors and fit list
        cvecold <- cvec
        Foldlist <- Flist
        #  update the expected Hessian
        if (crit == 2) {
           cvectmp <- cvec
           cvectmp[1] <- 0
           Wfdtmpi[[1]] <- cvectmp
           D2hwrtc <- monhess(xfine, Wfdtmpi, basislist, returnMatrix)
           D2fmax  <- D2hwrtc[nfine,]
           #  normalize 2nd derivative
           fmax2 <- fmax*fmax
           fmax3 <- fmax*fmax2
           m <- 1
           if (wnbasis > 1) {
              for (j in 2:wnbasis) {
                 m <- m + 1
                 for (k in 2:j) {
                    m <- m + 1
                    D2hwrtc[,m] <- width*(2*ffine*Dfmax[j]*Dfmax[k]
                   - fmax*(Dffine[,j]*Dfmax[k] + Dffine[,k]*Dfmax[j])
                   + fmax2*D2hwrtc[,m] - ffine*fmax*D2fmax[m])/fmax3
                 }
              }
           }
        } else {
           D2hwrtc <- NULL
        }
        hessmat <- reghess(xfine, y0fine, Dhfine, D2hwrtc, yregfdi,
                           Kmat, periodic, crit, returnMatrix)
        #  update the line search direction vector
        result   <- linesearch(Flist, hessmat, dbglev)
        deltac   <- result[[1]]
        cosangle <- result[[2]]
        reset <- 0
     }
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
  #  ---------------  end of optimization loop  -----------
  wcoef[,icurve] <- cvec
  if (nvar == 1) {
     yregcoef[,icurve]  <- yregfdi$coefs
  } else {
     yregcoef[,icurve,] <- yregfdi$coefs
  }
}

cat("\n")

#  --------------------   end of variable loop  -----------

#  create functional data objects for the registered curves

regfdnames <- yfd$fdnames
regfdnames[[3]] <- paste("Registered ",regfdnames[[3]])
ybasis  <- yfd$basis
regfd   <- fd(yregcoef, ybasis, regfdnames)

#  set up vector of time shifts

if (periodic) {
  shift <- c(wcoef[1,])
  wcoef[1,] <- wcoef1
} else {
  shift <- rep(0,ncurve)
}

#  functional data object for functions W(t)

Wfd <- fd(wcoef, wbasis)

#  functional data object for warping functions

warpmat = eval.monfd(xfine, Wfd)
warpmat = rangex[1] + (rangex[2]-rangex[1])*
           warpmat/outer(rep(1,nfine),warpmat[nfine,]) +
           outer(rep(1,nfine),shift)
if (wnbasis > 1) {
   warpfdobj  = smooth.basis(xfine, warpmat, wbasis)$fd
} else {
   wbasis    = create.monomial.basis(rangex, 2)
   warpfdobj = smooth.basis(xfine, warpmat, wbasis)$fd
}
warpfdnames       <- yfd$fdnames
warpfdnames[[3]]  <- paste("Warped",warpfdnames[[1]])
warpfdobj$fdnames <- warpfdnames

reglist <- list("regfd"=regfd, "warpfd"=warpfdobj, "Wfd"=Wfd,
                "shift"=shift, "y0fd"  =y0fd,      "yfd"=yfd)

return(reglist)
}

#  ----------------------------------------------------------------

regfngrad <- function(xfine, y0fine, Dhwrtc, yregfd, Wfd,
                      Kmat, periodic, crit, returnMatrix=FALSE)
{
  y0dim <- dim(y0fine)
  if (length(y0dim) == 3) nvar <- y0dim[3] else nvar <- 1
  nfine <- length(xfine)
  cvec  <- Wfd$coefs
  ncvec <- length(cvec)
  onecoef <- matrix(1,1,ncvec)

  if (periodic) {
     Dhwrtc[,1] <- 1
  } else {
     Dhwrtc[,1] <- 0
  }
  yregmat  <- eval.fd(xfine, yregfd, 0, returnMatrix)
  Dyregmat <- eval.fd(xfine, yregfd, 1, returnMatrix)

  #  loop through variables computing function and gradient values

  Fval <- 0
  gvec <- matrix(0,ncvec,1)
  for (ivar in 1:nvar) {
    y0ivar  <-   y0fine[,ivar]
    ywrthi  <-  yregmat[,ivar]
    Dywrthi <- Dyregmat[,ivar]
    aa      <- mean(y0ivar^2)
    bb      <- mean(y0ivar*ywrthi)
    cc      <- mean(ywrthi^2)
    Dywrtc  <- (Dywrthi %*% onecoef)*Dhwrtc
    if (crit == 1) {
      res  <- y0ivar - ywrthi
      Fval <- Fval + aa - 2*bb + cc
      gvec <- gvec - 2*crossprod(Dywrtc, res)/nfine
    } else if (crit == 2) {
      ee   <- aa + cc
      ff   <- aa - cc
      dd   <- sqrt(ff^2 + 4*bb^2)
      Fval <- Fval + ee - dd
      Dbb  <- crossprod(Dywrtc, y0ivar)/nfine
      Dcc  <- 2.0 * crossprod(Dywrtc, ywrthi)/nfine
      Ddd  <- (4*bb*Dbb - ff*Dcc)/dd
      gvec <- gvec + (Dcc - Ddd)
      # } else if (crit == 3) {
        #  least squares with root Dh weighting criterion
        #  dhdt not defined here ... fix later
      # res  = y0ivar - ywrthi*sqrt(dhdt)
      # Fval = Fval + mean(res^2)
      # gvec = gvec - 2*crossprod(Dywrtc, res)/nfine
    } else {
      stop("Invalid value for fitting criterion CRIT.")
    }
  }
  if (!is.null(Kmat)) {
     if (ncvec > 1) {
        ind   <- 2:ncvec
        ctemp <- cvec[ind,1]
        Kctmp <- Kmat%*%ctemp
        Fval  <- Fval + t(ctemp)%*%Kctmp
        gvec[ind] <- gvec[ind] + 2*Kctmp
     }
  }

#  set up FLIST list containing function value and gradient

  Flist      <- list(f=0, grad=rep(0,ncvec), norm=0)
  Flist$f    <- Fval
  Flist$grad <- gvec
  #  do not modify initial coefficient for B-spline and Fourier bases
  if (!periodic)  Flist$grad[1] <- 0
  Flist$norm <- sqrt(sum(Flist$grad^2))
  return(Flist)
}

#  ---------------------------------------------------------------

reghess <- function(xfine, y0fine, Dhfine, D2hwrtc, yregfd,
                    Kmat, periodic, crit, returnMatrix=FALSE)
{
	#cat("\nreghess")
  y0dim <- dim(y0fine)
  if (length(y0dim) == 3) nvar <- y0dim[3] else nvar <- 1
  nfine   <- length(xfine)
  wnbasis   <- dim(Dhfine)[2]
  onecoef <- matrix(1,1,wnbasis)
  npair   <- wnbasis*(wnbasis+1)/2

  if (periodic) {
     Dhfine[,1] <- 1
  } else {
     Dhfine[,1] <- 0
  }
  yregmat  <- eval.fd(yregfd, xfine, 0, returnMatrix)
  Dyregmat <- eval.fd(yregfd, xfine, 1, returnMatrix)
  if (nvar > 1) {
	   y0fine   <- y0fine[,1,]
	   yregmat  <- yregmat[,1,]
	   Dyregmat <- Dyregmat[,1,]
  }

  if (crit == 2) {
     D2yregmat <- eval.fd(yregfd, xfine, 2, returnMatrix)
     if (nvar > 1) D2yregmat <- D2yregmat[,1,]
     if (periodic) {
        D2hwrtc[,1] <- 0
        if (wnbasis > 1) {
           for (j in 2:wnbasis) {
              m <- j*(j-1)/2 + 1
              D2hwrtc[,m] <- Dhfine[,j]
           }
        }
     } else {
        D2hwrtc[,1] <- 1
        if (wnbasis > 1) {
           for (j in 2:wnbasis) {
              m <- j*(j-1)/2 + 1
              D2hwrtc[,m] <- 0
           }
        }
     }
  }

  hessvec <- matrix(0,npair,1)
  for (ivar in 1:nvar) {
    y0i        <-   y0fine[,ivar]
    yregmati   <-  yregmat[,ivar]
    Dyregmati  <- Dyregmat[,ivar]
    Dywrtc <- ((Dyregmati %*% onecoef)*Dhfine)
    if (crit == 1) {
      hessmat <-  2*crossprod(Dywrtc, Dywrtc)/nfine
      m <- 0
       for (j in 1:wnbasis) {
          for (k in 1:j) {
             m <- m + 1
             hessvec[m] <- hessvec[m] + hessmat[j,k]
          }
      }
    } else {
      D2yregmati <- D2yregmat[,ivar]
      aa     <- mean(y0i^2)
      bb     <- mean(y0i*yregmati)
      cc     <- mean(    yregmati^2)
      Dbb    <- crossprod(Dywrtc, y0i)/nfine
      Dcc    <- 2.0 * crossprod(Dywrtc, yregmati)/nfine
      D2bb   <- matrix(0,npair,1)
      D2cc   <- matrix(0,npair,1)
      crossprodmat <- matrix(0,nfine,npair)
      DyD2hmat     <- matrix(0,nfine,npair)
      m <- 0
      for (j in 1:wnbasis) {
        for (k in 1:j) {
          m <- m + 1
          crossprodmat[,m] <- Dhfine[,j]*Dhfine[,k]*D2yregmati
          DyD2hmat[,m] <- Dyregmati*D2hwrtc[,m]
          temp <- crossprodmat[,m] + DyD2hmat[,m]
          D2bb[m] <- mean(y0i*temp)
          D2cc[m] <- 2*mean(yregmati*temp +
                     Dyregmati^2*Dhfine[,j]*Dhfine[,k])
        }
      }
      ee     <- aa + cc
      ff     <- aa - cc
      ffsq   <- ff*ff
      dd     <- sqrt(ffsq + 4*bb*bb)
      ddsq   <- dd*dd
      ddcu   <- ddsq*dd
      m <- 0
      for (j in 1:wnbasis) {
        for (k in 1:j) {
          m <- m + 1
          hessvec[m] <- hessvec[m] + D2cc[m] -
            (4*Dbb[j]*Dbb[k] + 4*bb*D2bb[m] + Dcc[j]*Dcc[k] -
                   ff* D2cc[m])/dd +
            (4*bb*Dbb[j] - ff*Dcc[j])*(4*bb*Dbb[k] - ff*Dcc[k])/ddcu
        }
      }
    }
  }
  hessmat <- matrix(0,wnbasis,wnbasis)
  m <- 0
  for (j in 1:wnbasis) {
    for (k in 1:j) {
      m <- m + 1
      hessmat[j,k] <- hessvec[m]
      hessmat[k,j] <- hessvec[m]
    }
  }
  if (!is.null(Kmat)) {
     if (wnbasis > 1) {
        ind <- 2:wnbasis
        hessmat[ind,ind] <- hessmat[ind,ind] + 2*Kmat
     }
  }
  if (!periodic) {
     hessmat[1,]  <- 0
     hessmat[,1]  <- 0
     hessmat[1,1] <- 1
  }
  return(hessmat)
}

#  ----------------------------------------------------------------

regyfn <- function(xfine, yfine, hfine, yfd, Wfd, periodic)
{
	#cat("\nregyfn")
coef  <- Wfd$coefs
shift <- coef[1]
coef[1] <- 0
Wfd[[1]] <- coef

if (all(coef == 0)) {
   if (periodic) {
      if (shift == 0) {
         yregfd <- yfd
         return(yregfd)
      }
   } else {
      yregfd <- yfd
      return(yregfd)
   }
}

#  Estimate inverse of warping function at fine mesh of values
#  28 dec 000
#  It makes no real difference which
#     interpolation method is used here.
#  Linear is faster and sure to be monotone.
#  Using WARPSMTH added nothing useful, and was abandoned.
nfine       <- length(xfine)
hinv        <- approx(hfine, xfine, xfine)$y
hinv[1]     <- xfine[1]
hinv[nfine] <- xfine[nfine]

#  carry out shift if period and shift != 0
basis  <- yfd$basis
rangex <- basis$rangeval
ydim <- dim(yfine)
#if (length(ydim) == 3) yfine <- yfine[,1,]
if (periodic & shift != 0) yfine <- shifty(xfine, yfine, shift)
#  make FD object out of Y
ycoef  <- project.basis(yfine, hinv, basis, 1)
yregfd <- fd(ycoef, basis)
return(yregfd)
}

#  ----------------------------------------------------------------

linesearch <- function(Flist, hessmat, dbglev)
{
deltac   <- -solve(hessmat,Flist$grad)
cosangle <- -sum(Flist$grad*deltac)/sqrt(sum(Flist$grad^2)*sum(deltac^2))
if (dbglev >= 2) cat(paste("\nCos(angle) = ",round(cosangle,2)))
if (cosangle < 1e-7) {
   if (dbglev >=2) cat("\nangle negative")
   deltac <- -Flist$grad
}
return(list(deltac, cosangle))
}

#  ---------------------------------------------------------------------

shifty <- function(x, y, shift)
{
#SHIFTY estimates value of Y for periodic data for
#       X shifted by amount SHIFT.
#  It is assumed that X spans interval over which functionis periodic.
#  Last modified 6 February 2001

ydim <- dim(y)
if (is.null(ydim)) ydim <- 1
if (length(ydim) > 3) stop("Y has more than three dimensions")

if (shift == 0) {
   yshift <- y
   return(yshift)
}

n   <- ydim[1]
xlo <- min(x)
xhi <- max(x)
wid <- xhi - xlo
if (shift > 0) {
   while (shift > xhi)  shift <- shift - wid
   ind <- 2:n
   x2  <- c(x, x[ind]+wid)
   xshift <- x + shift
   if (length(ydim) == 1) {
	  y2 <- c(y, y[ind])
      yshift <- approx(x2, y2, xshift)$y
   }
   if (length(ydim) == 2) {
	   nvar <- ydim[2]
	   yshift <- matrix(0,n,nvar)
      for (ivar in 1:nvar) {
         y2 <- c(y[,ivar], y[ind,ivar])
         yshift[,ivar] <- approx(x2, y2, xshift)$y
      }
   }
   if (length(ydim) == 3) {
	   nrep <- ydim[2]
	   nvar <- ydim[3]
      yshift <- array(0,c(n,nrep,nvar))
      for (irep in 1:nrep) for (ivar in 1:nvar) {
         y2 <- c(y[,irep,ivar], y[ind,irep,ivar])
         yshift[,irep,ivar] <- approx(x2, y2, xshift)$y
      }
   }
} else {
   while (shift < xlo - wid) shift <- shift + wid
   ind <- 1:(n-1)
   x2 <- c(x[ind]-wid, x)
   xshift <- x + shift
   if (length(ydim) == 1) {
      y2 <- c(y[ind], y)
      yshift <- approx(x2, y2, xshift)$y
   }
   if (length(ydim) == 2) {
	   nvar <- ydim[2]
	   yshift <- matrix(0,n,nvar)
	   for (ivar in 1:nvar) {
		   y2 <- c(y[ind,ivar],y[,ivar])
		   yshift[,ivar] <- approx(x2, y2, xshift)$y
	   }
   }
   if (length(ydim) == 3) {
	   nrep <- ydim[2]
	   nvar <- ydim[3]
      yshift <- array(0, c(n,nrep,nvar))
      for (irep in 1:nrep) for (ivar in 1:nvar) {
         y2 <- c(y[ind,irep,ivar], y[,irep,ivar])
         yshift[,irep,ivar] <- approx(x2, y2, xshift)$y
      }
   }
}
return(yshift)
}

