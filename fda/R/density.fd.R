density <- function(x, ...)UseMethod('density')

density.fd <- function(x, WfdParobj, conv=0.0001, iterlim=20,
                      active=1:nbasis, dbglev=1, returnMatrix=FALSE, ...) {
# DENSITYFD estimates the density of a sample of scalar observations.

#  These observations may be one of two forms:
#   1.  a vector of observatons x_i
#   2.  a two-column matrix, with the observations x_i in the
#       first column, and frequencies f_i in the second.
#   Option 1. corresponds to all f_i = 1.

#  Arguments are:
#  X         ... data value array, either a vector or a two-column
#                matrix.
#  WFDPAROBJ ... functional parameter object specifying the initial log
#              density, the linear differential operator used to smooth
#              smooth it, and the smoothing parameter.
#  CONV      ... convergence criterion
#  ITERLIM   ... iteration limit for scoring iterations
#  ACTIVE    ... indices among 1:NBASIS of parameters to optimize
#  DBGLEV    ... level of output of computation history
#  RETURNMATRIX ... If False, a matrix in sparse storage model can be returned
#               from a call to function BsplineS.  See this function for
#               enabling this option.

#  Returns:
#  A list containing
#  WFDOBJ ...   functional data basis object defining final density
#  C      ...   normalizing constant for density p = exp(WFDOBJ)/C
#  FLIST  ...   Struct object containing
#               FSTR$f     final log likelihood
#               FSTR$norm  final norm of gradient
#  ITERNUM   Number of iterations
#  ITERHIST  History of iterations

#  To plot the density function or to evaluate it, evaluate WFDOBJ,
#  exponentiate the resulting vector, and then divide by the normalizing
#  constant C.

# last modified 3 March 2014 by Jim Ramsay

#  check WfdParobj

if (!inherits(WfdParobj, "fdPar")) {
	if (inherits(WfdParobj, "fd") || inherits(WfdParobj, "basisfd")) {
	    WfdParobj <- fdPar(WfdParobj)
	} else {
          stop("WFDPAROBJ is not a fdPar object")
      }
}

#  set up WFDOBJ

Wfdobj   <- WfdParobj$fd

#  set up LFDOBJ

Lfdobj <- WfdParobj$Lfd
Lfdobj <- int2Lfd(Lfdobj)

#  set up BASIS

basisobj <- Wfdobj$basis
nbasis   <- basisobj$nbasis
rangex   <- basisobj$rangeval

x    <- as.matrix(x)
xdim <- dim(x)
N    <- xdim[1]
m    <- xdim[2]

if (m > 2 && N > 2)
    	stop("Argument X must have either one or two columns.")

if ((N == 1 | N == 2) & m > 1) {
    x <- t(x)
    n <- N
    N <- m
    m <- n
}

if (m == 1) {
    f <- rep(1,N)
} else {
    f    <- x[,2]
    fsum <- sum(f)
    f    <- f/fsum
    x    <- x[,1]
}
f = as.matrix(f)

#  check for values outside of the range of WFD0

inrng <- (1:N)[x >= rangex[1] & x <= rangex[2]]
if (length(inrng) != N) {
    print(c(length(inrng), N))
    print(c(rangex[1], rangex[2], min(x), max(x)))
    warning("Some values in X out of range and not used.")
}
x     <- x[inrng]
f     <- f[inrng]
nobs  <- length(x)

#  set up some arrays

climit <- c(rep(-50,nbasis),rep(400,nbasis))
cvec0  <- Wfdobj$coefs
dbgwrd <- dbglev > 1

zeromat <- zerobasis(nbasis)

#  initialize matrix Kmat defining penalty term

lambda <- WfdParobj$lambda
if (lambda > 0) Kmat <- lambda*getbasispenalty(basisobj, Lfdobj)

#  evaluate log likelihood
#    and its derivatives with respect to these coefficients

result <- loglfnden(x, f, basisobj, cvec0, returnMatrix)
logl   <- result[[1]]
Dlogl  <- result[[2]]

#  compute initial badness of fit measures

fun  <- -logl
gvec <- -Dlogl
if (lambda > 0) {
   gvec <- gvec + 2*(Kmat %*% cvec0)
   fun    <- fun + t(cvec0) %*% Kmat %*% cvec0
}
Foldstr <- list(f = fun, norm = sqrt(mean(gvec^2)))
gvec0 <- t(zeromat) %*% as.matrix(gvec)

#  compute the initial expected Hessian

hmat <- Varfnden(x, basisobj, cvec0, returnMatrix)
if (lambda > 0) hmat <- hmat + 2*Kmat
hmat0 = t(zeromat) %*% hmat %*% zeromat

#  evaluate the initial update vector for correcting the initial bmat

deltac0  <- -solve(hmat0,gvec0)
deltac   <- zeromat %*% as.matrix(deltac0)
cosangle <- -sum(gvec0*deltac0)/sqrt(sum(gvec0^2)*sum(deltac0^2))

#  initialize iteration status arrays

iternum <- 0
status <- c(iternum, Foldstr$f, -logl, Foldstr$norm)
if (dbglev > 0) {
  cat("Iteration  Criterion  Neg. Log L  Grad. Norm\n")
  cat("      ")
  cat(format(iternum))
  cat("    ")
  cat(format(status[2:4]))
  cat("\n")
}
iterhist <- matrix(0,iterlim+1,length(status))
iterhist[1,]  <- status

#  quit if ITERLIM == 0

if (iterlim == 0) {
    Flist     <- Foldstr
    iterhist <- iterhist[1,]
    C        <- normden.phi(basisobj, cvec0, returnMatrix=returnMatrix)
    return( list(Wfdobj=Wfdobj, C=C, Flist=Flist, iternum=iternum,
                   iterhist=iterhist) )
}

#  -------  Begin iterations  -----------

STEPMAX <- 5
MAXSTEP <- 400
trial   <- 1
cvec    <- cvec0
linemat <- matrix(0,3,5)

for (iter in 1:iterlim) {
   	iternum <- iternum + 1
  	#  take optimal stepsize
	  dblwrd <- c(0,0)
	  limwrd <- c(0,0)
  	stpwrd <- 0
	  ind    <- 0
	  #  compute slope
      Flist <- Foldstr
      linemat[2,1] <- sum(deltac*gvec)
      #  normalize search direction vector
      sdg     <- sqrt(sum(deltac^2))
      deltac  <- deltac/sdg
      dgsum   <- sum(deltac)
      linemat[2,1] <- linemat[2,1]/sdg
      #  return with stop condition if (initial slope is nonnegative
      if (linemat[2,1] >= 0) {
        print("Initial slope nonnegative.")
        ind <- 3
        iterhist <- iterhist[1:(iternum+1),]
        break
      }
      #  return successfully if (initial slope is very small
      if (linemat[2,1] >= -1e-5) {
        if (dbglev>1) print("Initial slope too small")
        iterhist <- iterhist[1:(iternum+1),]
        break
      }
    	#  load up initial search matrix
      linemat[1,1:4] <- 0
      linemat[2,1:4] <- linemat[2,1]
      linemat[3,1:4] <- Foldstr$f
     	#  output initial results for stepsize 0
      stepiter  <- 0
      if (dbglev > 1) {
	      cat("              ")
	      cat(format(stepiter))
	      cat(format(linemat[,1]))
	      cat("\n")
	    }
      ips <- 0
      #  first step set to trial
      linemat[1,5]  <- trial
      #  Main iteration loop for linesrch
      for (stepiter in 1:STEPMAX) {
        #  ensure that step does not go beyond limits on parameters
        limflg  <- 0
        #  check the step size
        result <- stepchk(linemat[1,5], cvec, deltac, limwrd, ind,
                            climit, active, dbgwrd)
	      linemat[1,5] <- result[[1]]
	      ind          <- result[[2]]
	      limwrd       <- result[[3]]
        if (linemat[1,5] <= 1e-9) {
          	#  Current step size too small  terminate
          	Flist   <- Foldstr
          	cvecnew <- cvec
          	gvecnew <- gvec
          	if (dbglev > 1) print(paste("Stepsize too small:", linemat[1,5]))
          	if (limflg) ind <- 1 else ind <- 4
          	break
        }
        cvecnew <- cvec + linemat[1,5]*deltac
        #  compute new function value and gradient
	      result  <- loglfnden(x, f, basisobj, cvecnew, returnMatrix)
	      logl    <- result[[1]]
	      Dlogl   <- result[[2]]
        Flist$f <- -logl
        gvecnew <- -as.matrix(Dlogl)
        if (lambda > 0) {
            gvecnew <- gvecnew + 2*Kmat %*% cvecnew
            Flist$f <- Flist$f + t(cvecnew) %*% Kmat %*% cvecnew
        }
        gvecnew0 <- t(zeromat) %*% gvecnew
        Flist$norm <- sqrt(mean(gvecnew0^2))
        #  compute new directional derivative
        linemat[2,5] <- sum(deltac*gvecnew)
        linemat[3,5] <- Flist$f
        if (dbglev > 1) {
           cat("              ")
           cat(format(stepiter))
	         cat(format(linemat[,5]))
	         cat("\n")
	      }
        #  compute next step
	      result  <- stepit(linemat, ips, dblwrd, MAXSTEP)
	      linemat <- result$linemat
	      ips     <- result$ips
	      ind     <- result$ind
	      dblwrd  <- result$dblwrd
        trial   <- linemat[1,5]
        #  ind == 0 implies convergence
        if (ind == 0 | ind == 5) break
        #  end of line search loop
     	}

    	#  update current parameter vectors

     	cvec  <- cvecnew
     	gvec  <- gvecnew
      gvec0 <- t(zeromat) %*% as.matrix(gvec)
	    Wfdobj$coefs <- cvec
     	status <- c(iternum, Flist$f, -logl, Flist$norm)
     	iterhist[iter+1,] <- status
	    cat("      ")
	    cat(format(iternum))
	    cat("    ")
	    cat(format(status[2:4]))
	    cat("\n")

     	#  test for convergence

     	if (abs(Flist$f-Foldstr$f) < conv) {
          iterhist <- iterhist[1:(iternum+1),]
  	      C <- normden.phi(basisobj, cvec, returnMatrix=returnMatrix)
	        denslist <- list("Wfdobj" = Wfdobj, "C" = C, "Flist" = Flist,
			                     "iternum" = iternum, "iterhist" = iterhist)
	        return( denslist )
     	}
     	if (Flist$f >= Foldstr$f) break
     	#  compute the Hessian
     	hmat <- Varfnden(x, basisobj, cvec, returnMatrix)
     	if (lambda > 0) hmat <- hmat + 2*Kmat
      hmat0 <- t(zeromat) %*% hmat %*% zeromat
     	#  evaluate the update vector
     	deltac0 <- -solve(hmat0,gvec0)
     	cosangle  <- -sum(gvec0*deltac0)/sqrt(sum(gvec0^2)*sum(deltac0^2))
     	if (cosangle < 0) {
       	if (dbglev > 1) print("cos(angle) negative")
       	deltac0 <- -gvec0
     	}
      deltac <- zeromat %*% as.matrix(deltac0)
     	Foldstr <- Flist
		#  end of iterations
  	}
      #  compute final normalizing constant
 	C <- normden.phi(basisobj, cvec, returnMatrix=returnMatrix)
	denslist <- list("Wfdobj" = Wfdobj, "C" = C, "Flist" = Flist,
			             "iternum" = iternum, "iterhist" = iterhist)
 	return( denslist )
}

#  -----------------------------------------------------------------------------

loglfnden <- function(x, f, basisobj, cvec, returnMatrix=FALSE) {
	#  Computes the log likelihood and its derivative with
	#    respect to the coefficients in CVEC
   	N       <- length(x)
   	nbasis  <- basisobj$nbasis
   	fmat    <- outer(f, rep(1,nbasis))
   	fsum    <- sum(f)
   	nobs    <- length(x)
   	phimat  <- getbasismatrix(x, basisobj)
   	Cval    <- normden.phi(basisobj, cvec, , returnMatrix=returnMatrix)
   	logl    <- sum((phimat %*% cvec) * f - fsum*log(Cval)/N)
    EDw     <- expectden.phi(basisobj, cvec, Cval, returnMatrix=returnMatrix)
   	Dlogl   <- apply((phimat - outer(rep(1,nobs),EDw))*fmat,2,sum)
	return( list(logl, Dlogl) )
}

#  -----------------------------------------------------------------------------

Varfnden <- function(x, basisobj, cvec, returnMatrix=FALSE) {
	#  Computes the expected Hessian
   	nbasis  <- basisobj$nbasis
   	nobs    <- length(x)
   	Cval    <- normden.phi(basisobj, cvec, returnMatrix=returnMatrix)
   	EDw     <- expectden.phi(basisobj, cvec, Cval, returnMatrix=returnMatrix)
   	EDwDwt  <- expectden.phiphit(basisobj, cvec, Cval, returnMatrix=returnMatrix)
   	Varphi  <- nobs*(EDwDwt - outer(EDw,EDw))
	return(Varphi)
}

#  -----------------------------------------------------------------------------

normden.phi <- function(basisobj, cvec, JMAX=15, EPS=1e-7, returnMatrix=FALSE) {

#  Computes integrals of
#      p(x) = exp phi"(x) %*% cvec
#  by numerical integration using Romberg integration

  	#  check arguments, and convert basis objects to functional data objects

  	if (!inherits(basisobj, "basisfd") )
    	stop("First argument must be a basis function object.")

	  nbasis <- basisobj$nbasis
  	oneb   <- matrix(1,1,nbasis)
  	rng    <- basisobj$rangeval

  	#  set up first iteration

  	width <- rng[2] - rng[1]
  	JMAXP <- JMAX + 1
  	h <- matrix(1,JMAXP,1)
  	h[2] <- 0.25
  	#  matrix SMAT contains history of discrete approximations to the integral
  	smat <- matrix(0,JMAXP,1)
  	#  the first iteration uses just the }points
  	x  <- rng
  	nx <- length(x)
  	ox <- matrix(1,nx,1)
  	fx <- getbasismatrix(x, basisobj, returnMatrix=returnMatrix)
  	wx <- fx %*% cvec
  	wx[wx < -50] <- -50
  	px <- exp(wx)
  	smat[1]  <- width*sum(px)/2
  	tnm <- 0.5
  	j   <- 1

  	#  now iterate to convergence
  	for (j in 2:JMAX) {
    	tnm  <- tnm*2
    	del  <- width/tnm
    	if (j == 2) {
      		x <- (rng[1] + rng[2])/2
    	} else {
      		x <- seq(rng[1]+del/2, rng[2], del)
    	}
    	fx <- getbasismatrix(x, basisobj, returnMatrix=returnMatrix)
    	wx <- fx %*% cvec
    	wx[wx < -50] <- -50
    	px <- exp(wx)
    	smat[j] <- (smat[j-1] + width*sum(px)/tnm)/2
    	if (j >= 5) {
      		ind <- (j-4):j
			result <- polintarray(h[ind],smat[ind],0)
			ss  <- result[[1]]
			dss <- result[[2]]
      		if (!any(abs(dss) >= EPS*max(abs(ss)))) {
        		#  successful convergence
        		return(ss)
      		}
    	}
    	smat[j+1] <- smat[j]
    	h[j+1]    <- 0.25*h[j]
 	}
  	warning(paste("No convergence after ",JMAX," steps in NORMDEN.PHI"))
	return(ss)
}

#  -----------------------------------------------------------------------------

expectden.phi <- function(basisobj, cvec, Cval=1, nderiv=0,
                     JMAX=15, EPS=1e-7, returnMatrix=FALSE) {
    #  Computes expectations of basis functions with respect to density
    #      p(x) <- Cval^{-1} exp t(c)*phi(x)
    #  by numerical integration using Romberg integration

    #  check arguments, and convert basis objects to functional data objects

    if (!inherits(basisobj, "basisfd"))
    	stop("First argument must be a basis function object.")

    nbasis <- basisobj$nbasis
    rng    <- basisobj$rangeval
    oneb   <- matrix(1,1,nbasis)

    #  set up first iteration

    width <- rng[2] - rng[1]
    JMAXP <- JMAX + 1
    h <- matrix(1,JMAXP,1)
    h[2] <- 0.25
    #  matrix SMAT contains the history of discrete approximations to the integral
    smat <- matrix(0,JMAXP,nbasis)
    sumj <- matrix(0,1,nbasis)
    #  the first iteration uses just the }points
    x  <- rng
    nx <- length(x)
    ox <- matrix(1,nx,nx)
    fx <- getbasismatrix(x, basisobj, 0, returnMatrix)
    wx <- fx %*% cvec
    wx[wx < -50] <- -50
    px <- exp(wx)/Cval
    if (nderiv == 0) {
    	Dfx <- fx
    } else {
    	Dfx <- getbasismatrix(x, basisobj, 1, returnMatrix)
    }
    sumj <- t(Dfx) %*% px
    smat[1,]  <- width*as.vector(sumj)/2
    tnm <- 0.5
    j   <- 1

    #  now iterate to convergence

    for (j in 2:JMAX) {
      tnm  <- tnm*2
    	del  <- width/tnm
    	if (j == 2) {
        x <- (rng[1] + rng[2])/2
    	} else {
        x <- seq(rng[1]+del/2, rng[2], del)
    	}
    	nx <- length(x)
    	fx <- getbasismatrix(x, basisobj, 0, returnMatrix)
    	wx <- fx %*% cvec
    	wx[wx < -50] <- -50
    	px <- exp(wx)/Cval
    	if (nderiv == 0) {
        Dfx <- fx
    	} else {
        Dfx <- getbasismatrix(x, basisobj, 1, returnMatrix)
    	}
    	sumj <- t(Dfx) %*% px
    	smat[j,] <- (smat[j-1,] + width*as.vector(sumj)/tnm)/2
    	if (j >= 5) {
        ind <- (j-4):j
        temp <- smat[ind,]
	  result <- polintarray(h[ind],temp,0)
	  ss  <- result[[1]]
	  dss <- result[[2]]
        if (!any(abs(dss) > EPS*max(abs(ss)))) {
          #  successful convergence
          return(ss)
        }
    	}
    	smat[j+1,] <- smat[j,]
    	h[j+1]     <- 0.25*h[j]
    }
    warning(paste("No convergence after ",JMAX," steps in EXPECTDEN.PHI"))
    return(ss)
}

#  -----------------------------------------------------------------------------

expectden.phiphit <- function(basisobj, cvec, Cval=1, nderiv1=0, nderiv2=0,
                              JMAX=15, EPS=1e-7,
                              returnMatrix=FALSE) {

#  Computes expectations of cross product of basis functions with
#  respect to density
#      p(x) = Cval^{-1} int [exp t(c) %*% phi(x)] phi(x) t(phi(x)) dx
#  by numerical integration using Romberg integration

  	#  check arguments, and convert basis objects to functional data objects

  	if (!inherits(basisobj, "basisfd"))
    	stop("First argument must be a basis function object.")

  	nbasis <- basisobj$nbasis
  	rng    <- basisobj$rangeval
  	oneb   <- matrix(1,1,nbasis)

  	#  set up first iteration

  	width <- rng[2] - rng[1]
  	JMAXP <- JMAX + 1
  	h <- matrix(1,JMAXP,1)
  	h[2] <- 0.25
  	#  matrix SMAT contains history of discrete approximations to the integral
  	smat <- array(0,c(JMAXP,nbasis,nbasis))
  	#  the first iteration uses just the }points
  	x  <- rng
  	nx <- length(x)
  	fx <- getbasismatrix(x, basisobj, 0, returnMatrix)
  	wx <- fx %*% cvec
  	wx[wx < -50] <- -50
  	px <- exp(wx)/Cval
  	if (nderiv1 == 0) {
    	Dfx1 <- fx
  	} else {
    	Dfx1 <- getbasismatrix(x, basisobj, 1, returnMatrix)
  	}
  	if (nderiv2 == 0) {
    	Dfx2 <- fx
  	} else {
    	Dfx2 <- getbasismatrix(x, basisobj, 2, returnMatrix)
  	}
  	oneb <- matrix(1,1,nbasis)
  	sumj <- t(Dfx1) %*% ((px %*% oneb) * Dfx2)
  	smat[1,,]  <- width*as.matrix(sumj)/2
  	tnm <- 0.5
  	j   <- 1

  	#  now iterate to convergence
  	for (j in 2:JMAX) {
    	tnm  <- tnm*2
    	del  <- width/tnm
    	if (j == 2) {
      		x <- (rng[1] + rng[2])/2
    	} else {
      		x <- seq(rng[1]+del/2, rng[2], del)
    	}
    	nx <- length(x)
    	fx <- getbasismatrix(x, basisobj, 0, returnMatrix)
    	wx <- fx %*% cvec
    	wx[wx < -50] <- -50
    	px <- exp(wx)/Cval
    	if (nderiv1 == 0) {
      		Dfx1 <- fx
    	} else {
      		Dfx1 <- getbasismatrix(x, basisobj, 1, returnMatrix)
    	}
    	if (nderiv2 == 0) {
      		Dfx2 <- fx
    	} else {
      		Dfx2 <- getbasismatrix(x, basisobj, 2, returnMatrix)
    	}
    	sumj <- t(Dfx1) %*% ((px %*% oneb) * Dfx2)
    	smat[j,,] <- (smat[j-1,,] + width*as.matrix(sumj)/tnm)/2
    	if (j >= 5) {
      		ind <- (j-4):j
      		temp <- smat[ind,,]
	   		result <- polintarray(h[ind],temp,0)
	   		ss  <- result[[1]]
	   		dss <- result[[2]]
      		if (!any(abs(dss) > EPS*max(max(abs(ss))))) {
        		#  successful convergence
        		return(ss)
      		}
    	}
    	smat[j+1,,] <- smat[j,,]
    	h[j+1] <- 0.25*h[j]
  	}
  	warning(paste("No convergence after ",JMAX," steps in EXPECTDEN.PHIPHIT"))
	return(ss)
}

#  -----------------------------------------------------------------------------

polintarray <- function(xa, ya, x0) {
  	#  YA is an array with up to 4 dimensions
  	#     with 1st dim the same length same as the vector XA
  	n     <- length(xa)
  	yadim <- dim(ya)
  	if (is.null(yadim)) {
		yadim <- n
		nydim <- 1
  	} else {
    	nydim <- length(yadim)
  	}
  	if (yadim[1] != n) stop("First dimension of YA must match XA")
  	difx <- xa - x0
  	absxmxa <- abs(difx)
  	ns <- min((1:n)[absxmxa == min(absxmxa)])
  	cs <- ya
  	ds <- ya
  	if (nydim == 1) y <- ya[ns]
  	if (nydim == 2) y <- ya[ns,]
  	if (nydim == 3) y <- ya[ns,,]
  	if (nydim == 4) y <- ya[ns,,,]
  	ns <- ns - 1
  	for (m in 1:(n-1)) {
    	if (nydim == 1) {
      		for (i in 1:(n-m)) {
        		ho    <- difx[i]
        		hp    <- difx[i+m]
        		w     <- (cs[i+1] - ds[i])/(ho - hp)
        		ds[i] <- hp*w
        		cs[i] <- ho*w
      		}
      		if (2*ns < n-m) {
        		dy <- cs[ns+1]
      		} else {
        		dy <- ds[ns]
        		ns <- ns - 1
      		}
  		}
  		if (nydim == 2) {
      		for (i in 1:(n-m)) {
        		ho     <- difx[i]
        		hp     <- difx[i+m]
        		w      <- (cs[i+1,] - ds[i,])/(ho - hp)
        		ds[i,] <- hp*w
        		cs[i,] <- ho*w
      		}
      		if (2*ns < n-m) {
        		dy <- cs[ns+1,]
      		} else {
        		dy <- ds[ns,]
        		ns <- ns - 1
      		}
  		}
   		if (nydim == 3) {
      		for (i in 1:(n-m)) {
        		ho       <- difx[i]
        		hp       <- difx[i+m]
        		w        <- (cs[i+1,,] - ds[i,,])/(ho - hp)
        		ds[i,,] <- hp*w
        		cs[i,,] <- ho*w
      		}
      		if (2*ns < n-m) {
        		dy <- cs[ns+1,,]
      		} else {
        		dy <- ds[ns,,]
        		ns <- ns - 1
      		}
  		}
   		if (nydim == 4) {
      		for (i in 1:(n-m)) {
        		ho      <- difx[i]
        		hp      <- difx[i+m]
        		w       <- (cs[i+1,,,] - ds[i,,,])/(ho - hp)
        		ds[i,,,] <- hp*w
        		cs[i,,,] <- ho*w
      		}
      		if (2*ns < n-m) {
        		dy <- cs[ns+1,,,]

      		} else {
        		dy <- ds[ns,,,]
        		ns <- ns - 1
      		}
  		}
   		y <- y + dy
	}
   	return( list(y, dy) )
}
