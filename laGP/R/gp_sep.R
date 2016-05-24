#*******************************************************************************
#
# Local Approximate Gaussian Process Regression
# Copyright (C) 2013, The University of Chicago
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
#
# Questions? Contact Robert B. Gramacy (rbgramacy@chicagobooth.edu)
#
#*******************************************************************************


## newGPsep:
##
## build an initial separable GP representation on the C-side
## using the X-Z data and d/g paramterization.  

newGPsep <- function(X, Z, d, g, dK=FALSE)
  {
    n <- nrow(X)
    m <- ncol(X)
    if(is.null(n)) stop("X must be a matrix")
    if(length(Z) != n) stop("must have nrow(X) = length(Z)")
    if(length(d) == 1) d <- rep(d, m)
    else if(length(d) != m) stop("must have length(d) = ncol(X)")
    
    out <- .C("newGPsep_R",
              m = as.integer(m),
              n = as.integer(n),
              X = as.double(t(X)),
              Z = as.double(Z),
              d = as.double(d),
              g = as.double(g),
              dK = as.integer(dK),
              gpsepi = integer(1),
              package = "laGP")

    ## return C-side GP index
    return(out$gpsepi)
  }


## buildkGPsep:
##
## allocates/calculates the C-side derivative info (only) for particular 
## separable GP

buildKGPsep <- function(gpsepi)
  {
    .C("buildKGPsep_R",
       gpisep = as.integer(gpsepi),
       PACKAGE = "laGP")
    invisible(NULL)
  }


## deletedkGPsep:
##
## deletes the C-side derivative info (only) for particular separable GP

deletedkGPsep <- function(gpsepi)
  {
    .C("deletedKGPsep_R",
       gpi = as.integer(gpsepi),
       PACKAGE = "laGP")
    invisible(NULL)
  }

## deleteGPsep:
##
## deletes the C-side of a particular separable GP

deleteGPsep <- function(gpsepi)
  {
    .C("deleteGPsep_R",
       gpsepi = as.integer(gpsepi), package="laGP")
    invisible(NULL)
  }


## deleteGPseps:
##
## deletes all gpseps on the C side

deleteGPseps <- function()
  {
    .C("deleteGPseps_R", package="laGP")
    invisible(NULL)
  }


## llikGPsep:
##
## calculate the log likelihood of the GP

llikGPsep <- function(gpsepi, dab=c(0,0), gab=c(0,0))
  {
    r <- .C("llikGPsep_R",
            gpsepi = as.integer(gpsepi),
            dab = as.double(dab),
            gab = as.double(gab),
            llik = double(1),
            package = "laGP")

    return(r$llik)
  }


## getmGPsep
##
## acces the input dimension of a separable GP
##
## totall new to GPsep

getmGPsep <- function(gpsepi)
  {
    .C("getmGPsep_R", gpsepi = as.integer(gpsepi), m = integer(1), package="laGP")$m
  }


## getdGPsep
##
## acces the separable lengthscale of a separable gp
##
## totally new to GPsep

getdGPsep <- function(gpsepi)
  {
    m <- getmGPsep(gpsepi) 
    .C("getdGPsep_R", gpsepi = as.integer(gpsepi), d = double(m), package="laGP")$d
  }


## getgGPsep
##
## acces the input dimension of a separable GP
##
## totally new to GPsep

getgGPsep <- function(gpsepi)
  { 
    .C("getgGPsep_R", gpsepi = as.integer(gpsepi), g = double(1), package="laGP")$g
  }


## dllikGPsep:
##
## calculate the first and second derivative of the
## log likelihood of the GP with respect to d, the
## lengthscale parameter
##
## SIMILAR to dllikGP except with vector d and gpsep
## isntead of gp

dllikGPsep <- function(gpsepi, ab=c(0,0), param=c("d", "g"), d2nug=FALSE)
  {
    param <- match.arg(param)
    if(param == "d") {
      dim <- getmGPsep(gpsepi)
      r <- .C("dllikGPsep_R",
            gpsepi = as.integer(gpsepi),
            ab = as.double(ab),
            d = double(dim),
            package = "laGP")
      return(r$d)
    } else {
      if(d2nug) d2 <- 1
      else d2 <- 0
      r <- .C("dllikGPsep_nug_R",
            gpsepi = as.integer(gpsepi),
            ab = as.double(ab),
            d = double(1),
            d2 = as.double(d2),
            package = "laGP")
      if(d2nug) return(list(d=r$d, d2=r$r2))
      else return(r$d)
    }
  }



## newparamsGPsep:
##
## change the separable GP lengthscale and nugget parameerization
## (without destroying the object and creating a new one)

newparamsGPsep <- function(gpsepi, d, g=-1)
  {
    if(all(d <= 0) & g < 0) stop("one of d or g must be new")
    m <- getmGPsep(gpsepi)
    if(length(d) != m) stop("length(d) !=", m)

    r <- .C("newparamsGPsep_R",
            gpi = as.integer(gpsepi),
            d = as.double(d),
            g = as.double(g),
            package = "laGP")

    invisible(NULL)
  }


## mleGPsep.R:
##
## updates the separable GP to use its MLE lengthscale
## parameterization using the current data;
## 
## differs substantially from mleGP in that L-BFGS-B from
## optim is used to optimize over the separable lengthscale;
## an option is also provided to include the nugget in that
## optimization, or do to a mleGP style profile optimization
## for the the nugget instead
##
## in this ".R" version the optim command is used; in mleGPsep
## an internal C-side call to lbfgsb is used

mleGPsep.R <- function(gpsepi, param=c("d", "g"), 
                  tmin=sqrt(.Machine$double.eps), 
                  tmax=-1, ab=c(0,0), maxit=100, verb=0)
  {
    param <- match.arg(param)

    if(param == "d") { ## lengthscale with L-BFGS-B given nugget

      theta <- getdGPsep(gpsepi)
      if(length(ab) != 2 || any(ab < 0)) stop("ab should be a positive 2-vector")   
 
      ## objective
      f <- function(theta, gpsepi, dab) 
        {
          newparamsGPsep(gpsepi, d=theta)
          -llikGPsep(gpsepi, dab=dab)
        }
      ## gradient of objective
      g <- function(theta, gpsepi, dab)
        {
          newparamsGPsep(gpsepi, d=theta)
          -dllikGPsep(gpsepi, param="d", ab=dab)
        }

      ## for compatibility with mleGP
      tmax[tmax < 0] <- Inf

      ## possibly print progress meter
      if(verb > 0) {
        cat("(d=[", paste(signif(theta, 3), collapse=","), "], llik=", 
          llikGPsep(gpsepi, dab=ab), ") ", sep="")
      }

      ## call R's optim function
      out <- optim(theta, fn=f, gr=g, method="L-BFGS-B",
        control=list(trace=max(verb-1,0), maxit=maxit), lower=tmin, upper=tmax, 
        gpsepi=gpsepi, dab=ab)

      ## sanity check completion of scheme
      if(sqrt(mean((out$par - getdGPsep(gpsepi))^2)) > sqrt(.Machine$double.eps))
        warning("stored d not same as d-hat")

      ## check that we moved somewhere
      if(sqrt(mean((out$par - theta)^2)) < sqrt(.Machine$double.eps)) {
        out$convergence <- 0
        out$counts <- c(0,0)
        out$message <- "optim initialized at minima"
      }

      ## possibly print progress meter
      if(verb > 0) {
        cat("-> ", out$counts[1], " lbfgs.R its -> (d=[", 
          paste(signif(theta, 3), collapse=","), "], llik=", 
          llikGPsep(gpsepi, dab=ab), ")\n", sep="")
      }

    }
    else { ## nugget conditionally on lengthscale

      ## sanity check
      if(length(ab) != 2 || any(ab < 0)) stop("ab should be a positive 2-vector");

      r <- .C("mleGPsep_nug_R",
            gpsepi = as.integer(gpsepi),
            verb = as.integer(verb),
            tmin = as.double(tmin),
            tmax = as.double(tmax),
            ab = as.double(ab),
            g = double(1),
            its = integer(1),
            package = "laGP")
    }

    ## build object for returning
    if(param == "d") return(list(d=out$par, its=max(out$counts), msg=out$message, conv=out$convergence))
    else return(list(g=r$g, its=r$its))
  }


## mleGPsep:
##
## updates the separable GP to use its MLE lengthscale
## parameterization using the current data;
## 
## differs substantially from mleGP in that lbfgsb is used 
## to optimize over the separable lengthscale;
## an option is also provided to include the nugget in that
## optimization, or do to a mleGP style profile optimization
## for the the nugget instead
##
## this is a mostly C verision

mleGPsep <- function(gpsepi, param=c("d", "g"), 
                  tmin=sqrt(.Machine$double.eps), 
                  tmax=-1, ab=c(0,0), maxit=100, verb=0)
  {
    param <- match.arg(param)

    if(param == "d") { ## lengthscale with L-BFGS-B given nugget

      ## sanity checking
      m <- getmGPsep(gpsepi)
      if(length(tmax) == 1) tmax <- rep(tmax, m)
      else if(length(tmax) != m) stop("length(tmax) should be 1 or m")
      if(length(tmin) == 1) tmin <- rep(tmin, m)
      else if(length(tmin) != m) stop("length(tmin) should be 1 or m")
      if(length(ab) != 2 || any(ab < 0)) stop("ab should be a positive 2-vector")   

      out <- .C("mleGPsep_R",
                gpsepi = as.integer(gpsepi),
                maxit = as.integer(maxit),
                verb = as.integer(verb),
                dmin = as.double(tmin),
                dmax = as.double(tmax),
                ab = as.double(ab),
                par = double(m),
                counts = integer(2),
                msg = character(1),
                convergence = integer(1),
                package = "laGP")

      ## sanity check completion of scheme
      if(sqrt(mean((out$par - getdGPsep(gpsepi))^2)) > sqrt(.Machine$double.eps))
        warning("stored d not same as theta-hat")
    }
    else { ## nugget conditionally on lengthscale

      ## sanity check
      if(length(ab) != 2 || any(ab < 0)) stop("ab should be a positive 2-vector");

      r <- .C("mleGPsep_nug_R",
            gpsepi = as.integer(gpsepi),
            verb = as.integer(verb),
            tmin = as.double(tmin),
            tmax = as.double(tmax),
            ab = as.double(ab),
            g = double(1),
            its = integer(1),
            package = "laGP")
    }

    ## build object for returning
    if(param == "d") return(list(d=out$par, its=max(out$counts), msg=out$msg, conv=out$convergence))
    else return(list(g=r$g, its=r$its))
  }


## jmleGPsep.R:
##
## joint MLE for lengthscale (d) and nugget (g) parameters;
## updates the internal GP parameterization (since mleGP does);
## R-only version

jmleGPsep.R <- function(gpsepi, N=100, drange=c(sqrt(.Machine$double.eps), 10), 
  grange=c(sqrt(.Machine$double.eps), 1), dab=c(0,0), gab=c(0,0), maxit=100, 
  mleGPsep=mleGPsep.R, verb=0)
  {
    ## sanity check N
    if(length(N) != 1 && N > 0) 
      stop("N should be a positive scalar integer")
    m <- getmGPsep(gpsepi)
    dmle <- matrix(NA, nrow=N, ncol=m)
    gmle <- dits <- dconv <- gits <- rep(NA, N)

    ## sanity check tmin and tmax
    if(length(drange) != 2) stop("drange should be a 2-vector for c(min,max)")
    if(length(grange) != 2) stop("grange should be a 2-vector for c(min,max)")

    ## loop over outer interations
    for(i in 1:N) {
      d <- mleGPsep(gpsepi, param="d", tmin=drange[1], tmax=drange[2],
                    ab=dab, maxit=maxit, verb=verb)
      dmle[i,] <- d$d; dits[i] <- d$its; dconv[i] <- d$conv
      g <- mleGPsep(gpsepi, param="g", tmin=grange[1], tmax=grange[2],
                    ab=gab, verb=verb)
      gmle[i] <- g$g; gits[i] <- g$its
      if((gits[i] <= 2 && (dits[i] <= m+1 && dconv[i] == 0)) || dconv[i] > 1) break;
    }

    ## check if not converged
    if(i == N) warning("max outer its (N=", N, ") reached", sep="")
    else {
      dmle <- dmle[1:i,]; dits <- dits[1:i]; dconv <- dconv[1:i]
      gmle <- gmle[1:i]; gits <- gits[1:i]
    }

    ## total iteration count
    totits <- sum(c(dits, gits), na.rm=TRUE)

    ## assemble return objects
    return(list(mle=data.frame(d=dmle[i,,drop=FALSE], g=gmle[i], tot.its=totits, 
      conv=dconv[i]), prog=data.frame(dmle=dmle, dits=dits, dconv=dconv, gmle=gmle, 
      gits=gits)))
  }


## jmleGPsep
##
## interface to C-version for jmleGPsep; 
## right now doesn't take an N argument -- the C-side hard-codes
## N=100

jmleGPsep <- function(gpsepi, drange=c(sqrt(.Machine$double.eps), 10),
   grange=c(sqrt(.Machine$double.eps), 1), dab=c(0,0), gab=c(0,0), maxit=100,
   verb=0)
  {
    ## sanity check tmin and tmax
    m <- getmGPsep(gpsepi)
    if(length(drange) != 2) stop("drange should be a two vector for c(dmin, dmax)")
    dmin <- rep(drange[1], m)
    dmax <- rep(drange[2], m)
    if(length(grange) != 2) stop("grange should be a 2-vector for c(gmin, gmax)")

    ## sanity check ab
    if(length(dab) != 2 || any(dab < 0)) stop("dab should be a positive 2-vector")
    if(length(gab) != 2 || any(gab < 0)) stop("gab should be a positive 2-vector")

    ## call the C-side function
    r <- .C("jmleGPsep_R",
            gpsepi = as.integer(gpsepi),
            maxit = as.integer(maxit),
            verb = as.integer(verb),
            dmin = as.double(dmin),
            dmax = as.double(dmax),
            grange = as.double(grange),
            dab = as.double(dab),
            gab = as.double(gab),
            d = double(m),
            g = double(1),
            dits = integer(1),
            gits = integer(1),
            dconv = integer(1))

    return(data.frame(d=t(r$d), g=r$g, tot.its=r$dits+r$gits,
                      dits=r$dits, gits=r$gits, dconv=r$dconv))
  }


## predGPsep
##
## obtain the parameters to a multivariate-t
## distribution describing the predictive surface
## of the fitted GP model

predGPsep <- function(gpsepi, XX, lite=FALSE)
  {
    nn <- nrow(XX)

    if(lite) {  ## lite means does not compute full Sigma, only diag
      out <- .C("predGPsep_R",
                gpsepi = as.integer(gpsepi),
                m = as.integer(ncol(XX)),
                nn = as.integer(nn),
                XX = as.double(t(XX)),
                lite = as.integer(TRUE),
                mean = double(nn),
                s2 = double(nn),
                df = double(1),
                llik = double(1))
      
      ## coerce matrix output
      return(list(mean=out$mean, s2=out$s2, df=out$df, llik=out$llik))

    } else { ## compute full predictive covariance matrix

      out <- .C("predGPsep_R",
                gpsepi = as.integer(gpsepi),
                m = as.integer(ncol(XX)),
                nn = as.integer(nn),
                XX = as.double(t(XX)),
                lite = as.integer(FALSE),
                mean = double(nn),
                Sigma = double(nn*nn),
                df = double(1),
                llik = double(1))
      
      ## coerce matrix output
      Sigma <- matrix(out$Sigma, ncol=nn)
      
      ## return parameterization
      return(list(mean=out$mean, Sigma=Sigma, df=out$df, llik=out$llik))
    }
  }


## updateGPsep:
##
## add X-Z pairs to the C-side GPsep represnetation
## using only O(n^2) for each pair

updateGPsep <- function(gpsepi, X, Z, verb=0)
  {
    if(length(Z) != nrow(X))
      stop("bad dims")

    out <- .C("updateGPsep_R",
              gpsepi = as.integer(gpsepi),
              m = as.integer(ncol(X)),
              n = as.integer(nrow(X)),
              X = as.double(t(X)),
              Z = as.double(Z),
              verb = as.integer(verb),
              PACKAGE = "laGP")

    invisible(NULL)
  }


## alGPsep:
##
## calculate the E(Y) and EI(Y) for an augmented Lagrangian 
## composite objective function with linear objective (in X), or
## estimate objective (fhat) and constraint separable GP (gpsepi) 
## predictive surfaces

alGPsep <- function(XX, fgpsepi, fnorm, Cgpsepis, Cnorms, lambda, alpha, ymin, 
                 nomax=FALSE, N=100, fn=NULL, Bscale=1)
  {
    ## dims
    m <- ncol(XX)
    nn <- nrow(XX)
    nCgpseps <- length(Cgpsepis)

    ## checking lengths for number of constraint gps
    if(length(Cnorms) != nCgpseps) stop("length(Cgpsepis) != length(Cnorms)")
    if(length(lambda) != nCgpseps) stop("length(Cgpsepis) != length(lambda)")
    if(length(alpha) != 1) stop("length(alpha) != 1")

    ## checking scalars
    if(length(nomax) != 1) stop("nomax should be a scalar logical or negative number")
    if(length(N) != 1 || N <= 0) stop("N should be a positive integer scalar")
    if(length(ymin) != 1) stop("ymin should be a scalar")
    if(length(fnorm) != 1) stop("fnorm should be a scalar")

    ## run fn to get cheap objectives and constraints
    if(fgpsepi < 0 || any(Cgpsepis < 0)) {
      if(is.null(fn)) stop("fn must be provided when fgpsepi or Cgpsepis < -1")
      out <- fn(XX*Bscale, known.only=TRUE)
      if(fgpsepi < 0) {
        if(is.null(out$obj)) stop("fgpsepi < 0 but out$obj from fn() is NULL")
        obj <- out$obj
      } else obj <- NULL
      if(any(Cgpsepis < 0)) {
        C <- out$c
        if(ncol(C) != sum(Cgpsepis < 0)) stop("ncol(C) != sum(Cgpsepis < 0)")
      } else C <- NULL
    } else { obj <- C <- NULL }

    ## call the C-side
    out <- .C("alGPsep_R",
      m = as.integer(m),
      XX = as.double(t(XX)),
      nn = as.integer(nn),
      fgpsepi = as.integer(fgpsepi),
      ff = as.double(obj),
      fnorm = as.double(fnorm),
      nCgpseps = as.integer(nCgpseps),
      Cgpsepis = as.integer(Cgpsepis),    
      CC = as.double(C),  
      Cnorms = as.double(Cnorms),
      lambda = as.double(lambda),
      alpha = as.double(alpha),
      ymin = as.double(ymin),
      nomax = as.integer(nomax),
      N = as.integer(N),
      eys = double(nn),
      eis = double(nn),
      PACKAGE = "laGP")
    
    ## done
    return(data.frame(ey=out$eys, ei=out$eis))
  }


## eicGPsep:
##
## calculate EI(f) and p(Y(c) <= 0) for known linear or esitmated
## objective f and vectorized constraints C via separable GP (gpsepi)
## predictive surfaces; returns log probabilities (lplex) and 
## and log EIs

eicGPsep <- function(XX, fgpsepi, fnorm, Cgpsepis, Cnorms, fmin, fn=NULL, Bscale=1)
  {
    ## doms
    m <- ncol(XX)
    nn <- nrow(XX)
    nCgpseps <- length(Cgpsepis)

    ## checking lengths for number of constraint gps
    if(length(Cnorms) != nCgpseps) stop("length(Cgpsepis) != length(Cnorms)")
    ## checking scalars
    if(length(fmin) != 1) stop("ymin should be a scalar")
    if(length(fnorm) != 1) stop("fnorm should be a scalar")

    ## run fn to get cheap objectives and constraints
    if(fgpsepi < 0 || any(Cgpsepis < 0)) {
      if(is.null(fn)) stop("fn must be provided when fgpsepi or Cgpsepis < -1")
      out <- fn(XX*Bscale, known.only=TRUE)
      if(fgpsepi < 0) {
        if(is.null(out$obj)) stop("fgpsepi < 0 but out$obj from fn() is NULL")
        obj <- out$obj
      } else obj <- NULL
      if(any(Cgpsepis < 0)) {
        C <- out$c
        if(ncol(C) != sum(Cgpsepis < 0)) stop("ncol(C) != sum(Cgpsepis < 0)")
      } else C <- NULL
    }

    ## calculate expected improvement part
    if(fgpsepi < 0) {
      obj <- rowSums(XX) * fnorm
      if(!is.finite(fmin)) fmin <- quantile(obj, p=0.9)
      I <- fmin - obj
      ei <- pmax(I, 0)
    } else {
      p <- predGPsep(fgpsepi, XX=XX, lite=TRUE)
      pm <- p$mean * fnorm
      ps <- sqrt(p$s2) * fnorm
      if(!is.finite(fmin)) fmin <- quantile(pm, p=0.9)
      u <- (fmin  - pm)/ps
      ei <- ps*dnorm(u) + (fmin-pm)*pnorm(u)
    }

    ## calculate constraint part
    lplez <- matrix(NA, nrow=nrow(XX), nCgpseps)
    ik <- 1
    for(j in 1:nCgpseps) {
      if(Cgpsepis[j] < 0) {
        lplez[,j] <- log(C[,ik] <= 0)
        ik <- ik + 1
      } else {
        pc <- predGPsep(Cgpsepis[j], XX=XX, lite=TRUE)
        lplez[,j] <- pnorm(0, pc$mean, sqrt(pc$s2), log.p=TRUE)
      }
    }
    
    ## done
    return(data.frame(lei=log(ei), lplez=lplez))
  }


## alcGPsep:
##
## wrapper used to calculate the ALCs in C using
## the pre-stored separable GP representation.  
## Note that this only calculates the s2' component 
## of ds2 = s2 - s2'

alcGPsep <- function(gpsepi, Xcand, Xref=Xcand, 
  parallel=c("none", "omp", "gpu"), verb=0)
  {
    m <- ncol(Xcand)
    if(ncol(Xref) != m) stop("Xcand and Xref have mismatched cols")
    ncand <- nrow(Xcand)

    parallel <- match.arg(parallel)
    if(parallel == "omp") cf <- "alcGPsep_omp_R"
    else if(parallel == "gpu") stop("not implemented") ## cf <- "alcsepGP_gpu_R"
    else cf <- "alcGPsep_R"

    out <- .C(cf,
              gpsepi = as.integer(gpsepi),
              m = as.integer(m),
              Xcand = as.double(t(Xcand)),
              ncand = as.integer(ncand),
              Xref = as.double(t(Xref)),
              nref = as.integer(nrow(Xref)),
              verb = as.integer(verb),
              alcs = double(ncand),
              PACKAGE = "laGP")
    
    return(out$alcs)
  }


## mleGPsep.switch:
## 
## switch function for mle calculaitons by laGPsep.R

mleGPsep.switch <- function(gpsepi, method, d, g, verb) 
  { 
    ## do nothing if no MLE required
    if(!(d$mle || g$mle)) return(NULL)

    ## calculate derivatives
    if(d$mle && method != "mspe" && method != "efi") buildKGPsep(gpsepi)

    if(d$mle && g$mle) { ## joint lengthscale and nugget
      return(jmleGPsep(gpsepi, drange=c(d$min,d$max), grange=c(g$min, g$max), 
                    dab=d$ab, gab=g$ab))
    } else { ## maybe one or the other
      if(d$mle) { ## lengthscale only
        dmle <- mleGPsep(gpsepi, param="d", d$min, d$max, d$ab, verb=verb)
        return(data.frame(d=matrix(dmle$d, nrow=1), dits=dmle$its))
      } 
      if(g$mle) { ## nugget only
        gmle <- mleGPsep(gpsepi, param="g", g$min, g$max, g$ab, verb=verb)
        return(data.frame(g=gmle$g, gits=gmle$its))
      } 
    }
  }



## alcrayGPsep:
##
## wrapper used to optimize AIC via a ray search using
## the pre-stored separable GP representation.  Return 
## the convex combination s in (0,1) between Xstart and Xend;

alcrayGPsep <- function(gpsepi, Xref, Xstart, Xend, verb=0)
  {
    ## coerse to matrices
    if(is.null(ncol(Xref))) Xref <- matrix(Xref, nrow=1)
    if(is.null(ncol(Xstart))) Xstart <- matrix(Xstart, nrow=1)
    if(is.null(ncol(Xend))) Xend <- matrix(Xend, nrow=1)
   
    ## check dimensions of matrices
    m <- ncol(Xstart)
    if(ncol(Xref) != m) stop("Xstart and Xref have mismatched cols")
    if(ncol(Xend) != m) stop("Xend and Xref have mismatched cols")
    if(nrow(Xref) != 1) stop("only one reference location allowed for ray search")
    numrays <- nrow(Xstart)
    if(nrow(Xend) != numrays) stop("must have same number of starting and ending locations")

    ## call the C routine
    out <- .C("alcrayGPsep_R",
              gpsepi = as.integer(gpsepi),
              m = as.integer(m),
              Xref = as.double(t(Xref)),
              numrays = as.integer(numrays),
              Xstart = as.double(t(Xstart)),
              Xend = as.double(t(Xend)),
              verb = as.integer(verb),
              s = double(1),
              r = integer(1))
    
    ## return the convex combination
    return(list(r=out$r, s=out$s))
  }


## lalcrayGPsep.R:
##
## calculates a ray emiating from a random nearest (of start)
## neighbor(s) to Xref in Xcand.  The ending point of the ray
## is 10 times the (opposite) distance from Xstart to Xref,
## then alcrayGP (either C or R version) is called to optimize
## over the ray.  The candidate in Xcand which is closest
## to the solution is returned.  This works differently than
## lalcrayGPsep, since the starts of the rays are random from 
## 1:offset -- otherwise this is nearly identical to lalcrayGP.R

lalcrayGPsep.R <- function(gpsepi, Xref, Xcand, rect, offset=1, 
  numrays=ncol(Xref), verb=0) 
  {
    ## sanity checks
    m <- ncol(Xref)
    if(nrow(Xref) != 1) stop("alcray only applies for one Xref")
    if(m != ncol(Xcand)) stop("ncol(Xref) != ncol(Xcand)")
    if(ncol(rect) != m) stop("ncol(rect) != ncol(Xref)")
    if(length(offset) != 1 || offset < 1 || offset > nrow(Xcand))
      stop("offset should be a scalar integer >= 1 and <= nrow(Xcand)") 
    if(length(numrays) != 1 || numrays < 1)
      stop("numrays should be an integer scalar >= 1")

    ## adjust numrays
    if(numrays > nrow(Xcand)) numrays <- nrow(Xcand)

    ## get starting and ending point of ray
    Xstart <- Xcand[sample(1:offset, numrays),,drop=FALSE]
    Xend <- ray.end(numrays, Xref, Xstart, rect)

    ## solve for the best convex combination of Xstart and Xend
    so <- alcrayGPsep(gpsepi, Xref, Xstart, Xend, verb)
    Xstar <- matrix((1-so$s)*Xstart[so$r,] + so$s*Xend[so$r,], nrow=1)

    ## return the index of the closest Xcand to Xstar
    w <- which.min(distance(Xstar, Xcand)[1,])
    return(w)
  }


## lalcrayGPsep:
##
## wrapper to a C-side function used to calculate a ray emiating 
## from a random nearest (of start) neighbor(s) to Xref in Xcand.  
## The ending point of the ray is 10 times the (opposite) distance 
## from Xstart to Xref, then alcrayGPsep (on the C-side) is called to 
## optimize over the ray.  The candidate in Xcand which is closest
## to the solution is returned -- nearly identical to lalcrayGPsep

lalcrayGPsep <- function(gpsepi, Xref, Xcand, rect, offset=1, 
  numrays=ncol(Xref), verb=0)
  {
    ## sanity checks
    m <- ncol(Xref)
    ncand <- nrow(Xcand)
    if(nrow(Xref) != 1) stop("alcray only applies for one Xref")
    if(m != ncol(Xcand)) stop("ncol(Xref) != ncol(Xcand)")
    if(ncol(rect) != m) stop("ncol(rect) != ncol(Xref)")
    if(length(offset) != 1 || offset < 1 || offset > ncand)
      stop("offset should be a scalar integer >= 1 and <= nrow(Xcand)") 
    if(length(numrays) != 1 || numrays < 1)
      stop("numrays should be an integer scalar >= 1")

    out <- .C("lalcrayGPsep_R",
              gpsepi = as.integer(gpsepi),
              m = as.integer(m),
              Xcand = as.double(t(Xcand)),
              ncand = as.integer(ncand),
              Xref = as.double(t(Xref)),
              offset = as.integer(offset-1),
              numrays = as.integer(numrays),
              rect = as.double(t(rect)),
              verb = as.integer(verb),
              w = integer(1),
              PACKAGE = "laGP")

    return(out$w+1)
  }