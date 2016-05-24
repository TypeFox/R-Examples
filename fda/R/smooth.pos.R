smooth.pos <- function(argvals, y, WfdParobj, wtvec=rep(1,n), conv=1e-4,
                       iterlim=50, dbglev=1, returnMatrix=FALSE) {
#  Smooths the relationship of Y to ARGVALS using weights in WTVEC by fitting a
#     positive function of the form
#                      f(x) = exp W(x)
#     where  W  is a function defined over the same range as ARGVALS,
#                         W = log Df.
#  The fitting criterion is penalized mean squared error:
#    PENSSE(lambda) = \sum w_i[y_i - f(x_i)]^2 +
#                     \lambda * \int [L W(x)]^2 dx
#  where L is a linear differential operator defined in argument Lfdobj,
#  and w_i is a positive weight applied to the observation.
#  The function W(x) is expanded by the basis in functional data object
#    Wfdobj.

#  Arguments:
#  ARGVALS ...  Argument value array of length N, where N is the number of
#               observed curve values for each curve.  It is assumed that
#               that these argument values are common to all observed
#               curves.  If this is not the case, you will need to
#               run this function inside one or more loops, smoothing
#               each curve separately.
#  Y       ...  Function value array (the values to be fit).
#               If the functional data are univariate, this array will be
#               an N by NCURVE matrix, where N is the number of observed
#               curve values for each curve and NCURVE is the number of
#               curves observed.
#               If the functional data are muliivariate, this array will be
#               an N by NCURVE by NVAR matrix, where NVAR the number of
#               functions observed per case.  For example, for the gait
#               data, NVAR = 2, since we observe knee and hip angles.
#  WFDPAROBJ... A functional parameter or fdPar object.  This object
#               contains the specifications for the functional data
#               object to be estimated by smoothing the data.  See
#               comment lines in function fdPar for details.
#               The functional data object WFD in WFDPAROBJ is used
#               to initialize the optimization process.
#               Its coefficient array contains the starting values for
#               the iterative minimization of mean squared error.
#               This means that it must be set up with the coefficients
#               for each replication in Y and each variable.  Do not
#               supply only a basis object for the FDOBJ argument in
#               setting up the FDPAROBJ argument.
#  WTVEC   ...  a vector of weights, a vector of N one's by default.
#  CONV    ...  convergence criterion, 0.0001 by default
#  ITERLIM ...  maximum number of iterations, 50 by default.
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
#  FLIST ... A list object or a vector of list objects, one for
#            each curve (and each variable if functions are multivariate).
#            Each list object has slots:
#                 f    ... The sum of squared errors
#                 grad ... The gradient
#                 norm ... The norm of the gradient
#  When multiple curves and variables are analyzed, the lists containing
#  FLIST objects are indexed linear with curves varying inside
#  variables.

#  Last modified 16 April 2014  by Jim Ramsay


#  check ARGVALS

argvals <- argcheck(argvals)
n       <- length(argvals)
onesobs <- matrix(1,n,1)

#  at least two points are necessary for monotone smoothing

if (n < 2) stop('ARGVALS does not contain at least two values.')

#  check Y

y      = as.matrix(y)
ychk   = ycheck(y, n)
y      = ychk$y
ncurve = ychk$ncurve
nvar   = ychk$nvar
ndim   = ychk$ndim

#  check WfdParobj and get LAMBDA

WfdParobj = fdParcheck(WfdParobj)
lambda    = WfdParobj$lambda

#  the starting values for the coefficients are in FD object WFDOBJ

Wfdobj   <- WfdParobj$fd
Lfdobj   <- WfdParobj$Lfd
basisobj <- Wfdobj$basis     #  basis for W(argvals)
nbasis   <- basisobj$nbasis  #  number of basis functions

#  set up initial coefficient array

coef0 <- Wfdobj$coefs

#  check WTVEC

wtvec = wtcheck(n, wtvec)$wtvec

#  set up some arrays

climit  <- c(rep(-400,nbasis),rep(400,nbasis))
active  <- 1:nbasis

#  set up initial coefficient: use that of Wfd0 if its dimensions
#  are correct, otherwise set to a zero array.

coef0 = Wfdobj$coefs  #  initial coefficients
if (ndim == 2) {
    if (dim(coef0)[2] != ncurve || length(dim(coef0)) != 2) {
        coef0 = matrix(0,nbasis,ncurve)
        warning(paste("Dimensions of coefficient array inconsistent",
                      "with data Y, initial coefficients set to 0."))
    }
} else {
    if (length(dim(coef0)) != 3) {
        coef0 = array(0,c(nbasis,ncurve,nvar))
        warning(paste("Dimensions of coefficient array inconsistent",
                      "with data Y, initial coefficients set to 0."))
    }
    if  (dim(coef0)[2] != ncurve || dim(coef0)[3] != nvar) {
        coef0 = array(0,c(nbasis,ncurve,nvar))
        warning(paste("Dimensions of coefficient array inconsistent",
                      "with data Y, initial coefficients set to 0."))
    }
}
coef <- coef0

#  initialize matrix Kmat defining penalty term

if (lambda > 0) Kmat <- lambda*eval.penalty(basisobj, Lfdobj)
else            Kmat <- matrix(0,nbasis,nbasis)

#  --------------------------------------------------------------------
#              loop through variables and curves
#  --------------------------------------------------------------------

#  set up arrays and lists to contain returned information

if (ncurve > 1 || nvar > 1 ) Flist = vector("list",ncurve*nvar)
else                         Flist = NULL


for (ivar in 1:nvar) {
   if (ndim == 2) {
      sclfac = mean(c(y)^2)
   } else {
      sclfac = mean(c(y[,,ivar])^2)
   }
   for (icurve in 1:ncurve) {
      if (ndim == 2) {
         yi    = y[,icurve]
         cveci = coef0[,icurve]
      } else {
         yi    = y[,icurve,ivar]
         cveci = coef0[,icurve,ivar]
      }

	#  evaluate log likelihood
	#    and its derivatives with respect to these coefficients

	result <- PENSSEfun(argvals, yi, basisobj, cveci, Kmat, wtvec)
	PENSSE   <- result[[1]]
	DPENSSE  <- result[[2]]

	#  compute initial badness of fit measures

	f0    <- PENSSE
	gvec0 <- DPENSSE
	Foldlist <- list(f = f0, grad = gvec0, norm = sqrt(mean(gvec0^2)))

	#  compute the initial expected Hessian

	hmat0 <- PENSSEhess(argvals, yi, basisobj, cveci, Kmat, wtvec)

	#  evaluate the initial update vector for correcting the initial bmat

	deltac   <- -solve(hmat0,gvec0)
	cosangle <- -sum(gvec0*deltac)/sqrt(sum(gvec0^2)*sum(deltac^2))

	#  initialize iteration status arrays

	iternum <- 0
	status <- c(iternum, Foldlist$f, -PENSSE, Foldlist$norm)
      if (dbglev >= 1) {
        if (ncurve > 1 || nvar > 1) {
          if (ncurve > 1 && nvar > 1) {
            cat("\n")
            curvetitle = paste('Results for curve',icurve,'and variable',ivar)
          }
          if (ncurve > 1 && nvar == 1) {
            cat("\n")
            curvetitle = paste('Results for curve',icurve)
          }
          if (ncurve == 1 && nvar > 1) {
            cat("\n")
            curvetitle = paste('Results for variable',ivar)
          }
          cat("\n")
          cat(curvetitle)
        }
        cat("\n")
        cat("\nIter.   PENSSE   Grad Length")
        cat("\n")
        cat(iternum)
        cat("        ")
        cat(round(status[2],4))
        cat("      ")
        cat(round(status[3],4))
     }
	#  -------  Begin iterations  -----------

	MAXSTEPITER <- 10
	MAXSTEP     <- 400
	trial       <- 1
	linemat     <- matrix(0,3,5)
      Flisti      <- Foldlist
      gvec        <- gvec0
      dbgwrd      <- dbglev > 1

	if (iterlim == 0) {
            cat("\n")
      } else {
	for (iter in 1:iterlim) {
   		iternum <- iternum + 1
	   	#  take optimal stepsize
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
      	#  break with error condition if (initial slope is nonnegative
      	if (linemat[2,1] >= 0) {
        	  print("Initial slope nonnegative.")
        	  ind <- 3
        	  break
      	}
      	#  return successfully if (initial slope is very small
      	if (linemat[2,1] >= -1e-5) {
        	  if (dbglev > 1) print("Initial slope too small")
        	  break
      	}
      	#  first step set to trial
      	linemat[1,5]  <- trial
      	#  Main iteration loop for linesrch
      	for (stepiter in 1:MAXSTEPITER) {
        	  #  ensure that step does not go beyond limits on parameters
        	  limflg  <- 0
        	  #  check the step size
        	  result <- stepchk(linemat[1,5], cveci, deltac, limwrd, ind,
                                climit, active, dbgwrd)
		  linemat[1,5] <- result[[1]]
		  ind          <- result[[2]]
		  limwrd       <- result[[3]]
       	  if (linemat[1,5] <= 1e-7) {
          		#  Current step size too small  terminate
          		Flisti  <- Foldlist
          		cvecnew <- cveci
          		gvecnew <- gvec
                  if (dbglev >= 2) {
                    print("Stepsize too small")
                    print(linemat[1,5])
                  }
          		if (limflg) ind <- 1 else ind <- 4
          		break
              }
        	  #  compute new function value and gradient
              cvecnew  <- cveci + linemat[1,5]*deltac
		  result   <- PENSSEfun(argvals, yi, basisobj, cvecnew, Kmat, wtvec)
		  PENSSE   <- result[[1]]
		  DPENSSE  <- result[[2]]
        	  Flisti$f <- PENSSE
        	  gvecnew  <- DPENSSE
              Flisti$grad <- gvecnew
        	  Flisti$norm <- sqrt(mean(gvecnew^2))
        	  linemat[3,5] <- Flisti$f
        	  #  compute new directional derivative
        	  linemat[2,5] <- sum(deltac*gvecnew)
              if (dbglev >= 2) {
                cat("\n")
                cat(paste("                 ", stepiter, "  "))
                cat(format(round(t(linemat[,5]),6)))
              }
        	  #  compute next step
		  result  <- stepit(linemat, ips, dblwrd, MAXSTEP)
		  linemat <- result[[1]]
		  ips     <- result[[2]]
		  ind     <- result[[3]]
		  dblwrd  <- result[[4]]
        	  trial   <- linemat[1,5]
        	  #  ind == 0 implies convergence
        	  if (ind == 0 | ind == 5) break
              #  end of line search loop
     	      }
     	      cveci <- cvecnew
     	      gvec  <- gvecnew
     	      #  test for convergence
     	      if (abs(Flisti$f - Foldlist$f) < sclfac*conv) {
 	        cat("\n")
	        break
     	      }
     	      if (Flisti$f >= Foldlist$f) break
     	      #  compute the Hessian
     	      hmat <- PENSSEhess(argvals, yi, basisobj, cveci, Kmat, wtvec)
     	      #  evaluate the update vector
     	      deltac <- -solve(hmat,gvec)
     	      cosangle  <- -sum(gvec*deltac)/sqrt(sum(gvec^2)*sum(deltac^2))
     	      if (cosangle < 0) {
       	  if (dbglev > 1) print("cos(angle) negative")
       	  deltac <- -gvec
     	      }
     	      Foldlist <- Flisti
            #  display iteration status
            status <- c(iternum, Flisti$f, Flisti$norm)
            if (dbglev >= 1) {
              cat("\n")
              cat(iternum)
              cat("        ")
              cat(round(status[2],4))
              cat("      ")
              cat(round(status[3],4))
            }
            #  end of iteration loop
      }
      }

      #  save coefficients in arrays COEF and BETA

      if (ndim == 2) {
        coef[,icurve] = cveci
      } else {
        coef[,icurve,ivar] = cveci
      }

      #  save Flisti

      if (ncurve == 1 && nvar == 1) {
        Flist = Flisti
      } else {
        Flist[[(ivar-1)*ncurve+icurve]] = Flisti
      }

    }
  }

  Wfdobj = fd(coef, basisobj)

  posFd <- list("Wfdobj"=Wfdobj, "Flist"=Flist,
                "argvals"=argvals, "y"=y)

  class(posFd) <- 'posfd'

  posFd
}

#  ---------------------------------------------------------------

PENSSEfun <- function(argvals, yi, basisobj, cveci, Kmat, wtvec,
                      returnMatrix=FALSE) {
	#  Computes the log likelihood and its derivative with
	#    respect to the coefficients in CVEC
	N       <- length(argvals)
	nbasis  <- basisobj$nbasis
	phimat  <- getbasismatrix(argvals, basisobj, 0, returnMatrix)
	Wvec    <- phimat %*% cveci
	EWvec   <- exp(Wvec)
	res     <- yi - EWvec
	PENSSE  <- mean(wtvec*res^2) + t(cveci) %*% Kmat %*% cveci
	DPENSSE <- -2*crossprod(phimat,wtvec*res*EWvec)/N + 2*Kmat %*% cveci
	return( list(PENSSE, DPENSSE) )
}

#  ---------------------------------------------------------------

PENSSEhess <- function(argvals, yi, basisobj, cveci, Kmat, wtvec,
                       returnMatrix=FALSE) {
	#  Computes the expected Hessian
   	n       <- length(argvals)
   	nbasis  <- basisobj$nbasis
   	phimat  <- getbasismatrix(argvals, basisobj, 0, returnMatrix)
	Wvec    <- phimat %*% cveci
	EWvec   <- exp(Wvec)
	res     <- yi - EWvec
      temp    <- wtvec*res*EWvec
	Dres    <- matrix(temp,n,nbasis) * phimat
  	D2PENSSE  <- 2*crossprod(Dres)/n + 2*Kmat
	return(D2PENSSE)
}


