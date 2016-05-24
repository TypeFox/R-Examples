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


## newGP:
##
## build an initial GP representation on the C-side
## using the X-Z data and d/g paramterization.  Calling
## this function writes over the previous GP representation

newGP <- function(X, Z, d, g, dK=FALSE)
  {
    n <- nrow(X)
    if(is.null(n)) stop("X must be a matrix")
    if(length(Z) != n) stop("must have nrow(X) = length(Z)")
    
    out <- .C("newGP_R",
              m = as.integer(ncol(X)),
              n = as.integer(n),
              X = as.double(t(X)),
              Z = as.double(Z),
              d = as.double(d),
              g = as.double(g),
              dK = as.integer(dK),
              gpi = integer(1),
              PACKAGE = "laGP")

    ## return C-side GP index
    return(out$gpi)
  }


## buildkGP:
##
## allocates/calculates the C-side derivative info (only) for particular GP

buildKGP <- function(gpi)
  {
    .C("buildKGP_R",
       gpi = as.integer(gpi),
       PACKAGE = "laGP")
    invisible(NULL)
  }


## deletedkGP:
##
## deletes the C-side derivative info (only) for particular GP

deletedkGP <- function(gpi)
  {
    .C("deletedKGP_R",
       gpi = as.integer(gpi),
       PACKAGE = "laGP")
    invisible(NULL)
  }


## deleteGP:
##
## deletes the C-side of a particular GP

deleteGP <- function(gpi)
  {
    .C("deleteGP_R",
       gpi = as.integer(gpi),
       PACKAGE = "laGP")
    invisible(NULL)
  }


## deleteGPs:
##
## deletes all gps on the C side

deleteGPs <- function()
  {
    .C("deleteGPs_R", PACKAGE="laGP")
    invisible(NULL)
  }

## copyGP:
##
## allocate a new GP with a copy of the contents of an
## old one

copyGP <- function(gpi)
  {
    r <- .C("copyGP_R",
            gpi = as.integer(gpi),
            newgpi = integer(1),
            PACKAGE = "laGP")

    return(r$newgpi)
  }


## newparamsGP:
##
## change the GP lengthscale and nugget parameerization
## (without destroying the object and creating a new one)

newparamsGP <- function(gpi, d=-1, g=-1)
  {
    if(d <= 0 & g < 0) stop("one of d or g must be new")
    r <- .C("newparamsGP_R",
            gpi = as.integer(gpi),
            d = as.double(d),
            g = as.double(g),
            PACKAGE = "laGP")

    invisible(NULL)
  }


## llikGP:
##
## calculate the log likelihood of the GP
llikGP <- function(gpi, dab=c(0,0), gab=c(0,0))
  {
    r <- .C("llikGP_R",
            gpi = as.integer(gpi),
            dab = as.double(dab),
            gab = as.double(gab),
            llik = double(1),
            PACKAGE = "laGP")

    return(r$llik)
  }


## jmleGP.R:
##
## joint MLE for lengthscale (d) and nugget (g) parameters;
## updates the internal GP parameterization (since mleGP does);
## R-only version

jmleGP.R <- function(gpi, N=100, drange=c(sqrt(.Machine$double.eps), 10), 
  grange=c(sqrt(.Machine$double.eps), 1), dab=c(0,0), gab=c(0,0), verb=0)
  {
    ## sanity check N
    if(length(N) != 1 && N > 0) 
      stop("N should be a positive scalar integer")
    dmle <- gmle <- rep(NA, N)
    dits <- gits <- rep(NA, N)

    ## sanity check tmin and tmax
    if(length(drange) != 2) stop("drange should be a 2-vector for c(min,max)")
    if(length(grange) != 2) stop("grange should be a 2-vector for c(min,max)")

    ## loop over outer interations
    for(i in 1:N) {
      d <- mleGP(gpi, param="d", tmin=drange[1], tmax=drange[2],
                 ab=dab, verb=verb)
      dmle[i] <- d$d; dits[i] <- d$its
      g <- mleGP(gpi, param="g", tmin=grange[1], tmax=grange[2],
                 ab=gab, verb=verb)
      gmle[i] <- g$g; gits[i] <- g$its
      if(gits[i] <= 1 && dits[i] <= 1) break;
    }

    ## check if not convergedf
    if(i == N) warning("max outer its (N=", N, ") reached", sep="")
    else {
      dmle <- dmle[1:i]; dits <- dits[1:i]
      gmle <- gmle[1:i]; gits <- gits[1:i]
    }

    ## total iteration count
    totits <- sum(c(dits, gits), na.rm=TRUE)

    ## assemble return objects
    return(list(mle=data.frame(d=dmle[i], g=gmle[i], tot.its=totits), 
      prog=data.frame(dmle=dmle, dits=dits, gmle=gmle, gits=gits)))
  }


## jmleGP
##
## interface to C-version for jmleGP; 
## right now doesn't take an N argument -- the C-side hard-codes
## N=100

jmleGP <- function(gpi, drange=c(sqrt(.Machine$double.eps), 10), 
  grange=c(sqrt(.Machine$double.eps), 1), dab=c(0,0), gab=c(0,0), verb=0)
  {
    ## sanity check tmin and tmax
    if(length(drange) != 2) stop("drange should be a 2-vector for c(d,g)")
    if(length(grange) != 2) stop("grange should be a 2-vector for c(d,g)")

    ## sanity check ab
    if(length(dab) != 2 || any(dab < 0)) stop("dab should be a positive 2-vector")
    if(length(gab) != 2 || any(gab < 0)) stop("gab should be a positive 2-vector")

    ## call the C-side function
    r <- .C("jmleGP_R",
            gpi = as.integer(gpi),
            verb = as.integer(verb),
            drange = as.double(drange),
            grange = as.double(grange),
            dab = as.double(dab),
            gab = as.double(gab),
            d = double(1),
            g = double(1),
            dits = integer(1),
            gits = integer(1),
            package = "laGP")

    return(data.frame(d=r$d, g=r$g, tot.its=r$dits+r$gits,
                      dits=r$dits, gits=r$gits))
  }



## mleGP:
##
## updates the GP to use its MLE lengthscale
## parameterization using the current data

mleGP <- function(gpi, param=c("d", "g"), 
                  tmin=sqrt(.Machine$double.eps), 
                  tmax=-1, ab=c(0,0), verb=0)
  {
    param <- match.arg(param)
    if(param == "d") param <- 1
    else param <- 2

    ## sanity check
    if(length(ab) != 2 || any(ab < 0)) stop("ab should be a positive 2-vector")

    r <- .C("mleGP_R",
            gpi = as.integer(gpi),
            param = as.integer(param),
            verb = as.integer(verb),
            tmin = as.double(tmin),
            tmax = as.double(tmax),
            ab = as.double(ab),
            theta = double(1),
            its = integer(1),
            package = "laGP")

    if(param == 1) return(list(d=r$theta, its=r$its))
    else return(list(g=r$theta, its=r$its))
  }


## dllikGP:
##
## calculate the first and second derivative of the
## log likelihood of the GP with respect to d, the
## lengthscale parameter

dllikGP <- function(gpi, ab=c(0,0), param=c("d", "g"))
  {
    param <- match.arg(param)
    if(param == "d") cf <- "dllikGP_R"
    else cf <- "dllikGP_nug_R"

    r <- .C(cf,
            gpi = as.integer(gpi),
            ab = as.double(ab),
            d = double(1),
            d2 = double(1),
            package = "laGP")

    return(data.frame(d=r$d, d2=r$d2))
  }


## mleGP.switch:
## 
## switch function for mle calculaitons by laGP.R

mleGP.switch <- function(gpi, method, d, g, verb) 
  { 
    ## do nothing if no MLE required
    if(!(d$mle || g$mle)) return(NULL)

    ## calculate derivatives
    if(d$mle && method != "mspe" && method != "efi") buildKGP(gpi)

    ## switch
    if(d$mle && g$mle) { 
      ## joint lengthscale and nugget
      return(jmleGP(gpi, drange=c(d$min,d$max), grange=c(g$min, g$max), 
                    dab=d$ab, gab=g$ab))
    } else { ## maybe one or the other
      if(d$mle) { ## lengthscale only
        dmle <- mleGP(gpi, param="d", d$min, d$max, d$ab, verb=verb)
        return(data.frame(d=dmle$d, dits=dmle$its))
      } 
      if(g$mle) { ## nugget only
        gmle <- mleGP(gpi, param="g", g$min, g$max, g$ab, verb=verb)
        return(data.frame(g=gmle$g, gits=gmle$its))
      } 
    }
  }


## updateGP:
##
## add X-Z pairs to the C-side GP represnetation
## using only O(n^2) for each pair

updateGP <- function(gpi, X, Z, verb=0)
  {
    if(length(Z) != nrow(X))
      stop("bad dims")

    out <- .C("updateGP_R",
              gpi = as.integer(gpi),
              m = as.integer(ncol(X)),
              n = as.integer(nrow(X)),
              X = as.double(t(X)),
              Z = as.double(Z),
              verb = as.integer(verb),
              PACKAGE = "laGP")

    invisible(NULL)
  }


## predGP
##
## obtain the parameters to a multivariate-t
## distribution describing the predictive surface
## of the fitted GP model

predGP <- function(gpi, XX, lite=FALSE)
  {
    nn <- nrow(XX)

    if(lite) {
      out <- .C("predGP_R",
                gpi = as.integer(gpi),
                m = as.integer(ncol(XX)),
                nn = as.integer(nn),
                XX = as.double(t(XX)),
                lite = as.integer(TRUE),
                mean = double(nn),
                s2 = double(nn),
                df = double(1),
                llik = double(1),
                PACKAGE = "laGP")
      
      ## coerce matrix output
      return(list(mean=out$mean, s2=out$s2, df=out$df, llik=out$llik))

    } else {

      out <- .C("predGP_R",
                gpi = as.integer(gpi),
                m = as.integer(ncol(XX)),
                nn = as.integer(nn),
                XX = as.double(t(XX)),
                lite = as.integer(FALSE),
                mean = double(nn),
                Sigma = double(nn*nn),
                df = double(1),
                llik = double(1),
                PACKAGE = "laGP")
      
      ## coerce matrix output
      Sigma <- matrix(out$Sigma, ncol=nn)
      
      ## return parameterization
      return(list(mean=out$mean, Sigma=Sigma, df=out$df, llik=out$llik))
    }
  }


## alcGP:
##
## wrapper used to calculate the ALCs in C using
## the pre-stored GP representation.  Note that this only
## calculates the s2' component of ds2 = s2 - s2'

alcGP <- function(gpi, Xcand, Xref=Xcand, parallel=c("none", "omp", "gpu"), 
                  verb=0)
  {
    m <- ncol(Xcand)
    if(ncol(Xref) != m) stop("Xcand and Xref have mismatched cols")
    ncand <- nrow(Xcand)

    parallel <- match.arg(parallel)
    if(parallel == "omp") cf <- "alcGP_omp_R"
    else if(parallel == "gpu") cf <- "alcGP_gpu_R"
    else cf <- "alcGP_R"

    out <- .C(cf,
              gpi = as.integer(gpi),
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


## ray.bounds:
##
## get the end of the ray eminating from Xstart that goes
## away (in the direction) from Xref which is 10 times the
## distance but not beyond the bounding rectangle

ray.end <- function(numrays, Xref, Xstart, rect)
  {
    for(r in 1:numrays) {
      Xend[r,] <- as.numeric(10*(Xstart[r,] - Xref) + Xstart[r,])
      while(any(Xend[r,] < rect[1,]) | any(Xend[r,] > rect[2,])) {
        w <- which(Xend[r,] < rect[1,])
        if(length(w) > 0) { col <- 1
        } else { w <- which(Xend[r,] > rect[2,]); col <- 2 }
        w <- w[1]
        sc <- (rect[col,w] - Xstart[r,w])/(Xend[r,w] - Xstart[r,w])
        Xend[r,] <- (Xend[r,] - Xstart[r,])*sc + Xstart[r,]
      }
    }

    return(Xend)
  }


## lalcrayGP.R:
##
## calculates a ray emiating from a random nearest (of start)
## neighbor(s) to Xref in Xcand.  The ending point of the ray
## is 10 times the (opposite) distance from Xstart to Xref,
## then alcrayGP (either C or R version) is called to optimize
## over the ray.  The candidate in Xcand which is closest
## to the solution is returned.  This works differently than
## lalcrayGP, since the starts of the rays are random from 
## 1:offset

lalcrayGP.R <- function(gpi, Xref, Xcand, rect, offset=1, numrays=ncol(Xref), verb=0) 
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
    so <- alcrayGP(gpi, Xref, Xstart, Xend, verb)
    Xstar <- matrix((1-so$s)*Xstart[so$r,] + so$s*Xend[so$r,], nrow=1)

    ## return the index of the closest Xcand to Xstar
    w <- which.min(distance(Xstar, Xcand)[1,])
    return(w)
  }


## lalcrayGP:
##
## wrapper to a C-side function used to calculate a ray emiating 
## from a random nearest (of start) neighbor(s) to Xref in Xcand.  
## The ending point of the ray is 10 times the (opposite) distance 
## from Xstart to Xref, then alcrayGP (on the C-side) is called to 
## optimize over the ray.  The candidate in Xcand which is closest
## to the solution is returned

lalcrayGP <- function(gpi, Xref, Xcand, rect, offset=1, numrays=ncol(Xref), verb=0)
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

    out <- .C("lalcrayGP_R",
              gpi = as.integer(gpi),
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


## alcrayGP:
##
## wrapper used to optimize AIC via a ray search using
## the pre-stored GP representation.  Return the convex
## combination s in (0,1) between Xstart and Xend

alcrayGP <- function(gpi, Xref, Xstart, Xend, verb=0)
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
    out <- .C("alcrayGP_R",
              gpi = as.integer(gpi),
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


## alGP:
##
## calculate the E(Y) and EI(Y) for an augmented Lagrangian 
## composite objective function with linear objective (in X)
## and constraint GP (gpi) predictive surfaces

alGP <- function(XX, fgpi, fnorm, Cgpis, Cnorms, lambda, alpha, ymin, 
                 nomax=FALSE, N=100, fn=NULL, Bscale=1)
  {
    ## dimensions
    m <- ncol(XX)
    nn <- nrow(XX)
    nCgps <- length(Cgpis)

    ## checking lengths for number of gps
    if(length(Cnorms) != nCgps) stop("length(Cgpis) != length(Cnorms)")
    if(length(lambda) != nCgps) stop("length(Cgpis) != length(lambda)")
    if(length(alpha) != 1) stop("length(alpha) != 1")

    ## checking scalars
    if(length(nomax) != 1) stop("nomax should be a scalar logical or negative nuumber")
    if(length(N) != 1 || N <= 0) stop("N should be a positive integer scalar")
    if(length(ymin) != 1) stop("ymin should be a scalar")
    if(length(fnorm) != 1) stop("fnorm should be a scalar")

    ## run fn to get cheap objectives and constraints
    if(fgpi < 0 || any(Cgpis < 0)) {
      if(is.null(fn)) stop("fn must be provided when fgpi or Cgpis < -1")
      out <- fn(XX*Bscale, known.only=TRUE)
      if(fgpi < 0) {
        if(is.null(out$obj)) stop("fgpsepi < 0 but out$obj from fn() is NULL")
        obj <- out$obj
      } else obj <- NULL
      if(any(Cgpis < 0)) {
        C <- out$c
        if(ncol(C) != sum(Cgpis < 0)) stop("ncol(C) != sum(Cgpis < 0)")
      } else C <- NULL
    } else { obj <- C <- NULL }

    ## call the C-side
    out <- .C("alGP_R",
      m = as.integer(m),
      XX = as.double(t(XX)),
      nn = as.integer(nn),
      fgpi = as.integer(fgpi),
      fnorm = as.double(fnorm),
      nCgps = as.integer(nCgps),
      Cgpis = as.integer(Cgpis),
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


## eicGP:
##
## calculate EI(f) and p(Y(c) <= 0) for known linear or esitmated
## objective f and vectorized constraints C via isotropic GP (gpsi)
## predictive surfaces; returns log probabilities (lplex) and 
## EIs on the original scale

eicGP <- function(XX, fgpi, fnorm, Cgpis, Cnorms, fmin, fn=NULL, Bscale=1)
  {
    ## doms
    m <- ncol(XX)
    nn <- nrow(XX)
    nCgps <- length(Cgpis)

    ## checking lengths for number of constraint gps
    if(length(Cnorms) != nCgps) stop("length(Cgpis) != length(Cnorms)")
    ## checking scalars
    if(length(fmin) != 1) stop("ymin should be a scalar")
    if(length(fnorm) != 1) stop("fnorm should be a scalar")

    ## run fn to get cheap objectives and constraints
    if(fgpi < 0 || any(Cgpis < 0)) {
      if(is.null(fn)) stop("fn must be provided when fgpi or Cgpis < -1")
      out <- fn(XX*Bscale, known.only=TRUE)
      if(fgpi < 0) {
        if(is.null(out$obj)) stop("fgpi < 0 but out$obj from fn() is NULL")
        obj <- out$obj
      } else obj <- NULL
      if(any(Cgpis < 0)) {
        C <- out$c
        if(ncol(C) != sum(Cgpis < 0)) stop("ncol(C) != sum(Cgpis < 0)")
      } else C <- NULL
    }

    ## calculate expected improvement part
    if(fgpi <= 0) {
      obj <- rowSums(XX) * fnorm
      if(!is.finite(fmin)) fmin <- quantile(obj, p=0.9)
      I <- fmin - obj
      ei <- pmax(I, 0)
    } else {
      p <- predGP(fgpi, XX=XX, lite=TRUE)
      pm <- p$mean * fnorm
      ps <- sqrt(p$s2) * fnorm
      if(!is.finite(fmin)) fmin <- quantile(pm, p=0.9)
      u <- (fmin  - pm)/ps
      ei <- ps*dnorm(u) + (fmin-pm)*pnorm(u)
    }

    ## calculate constraint part
    lplez <- matrix(NA, nrow=nrow(XX), nCgps)
    for(j in 1:nCgps) {
      pc <- predGP(Cgpis[j], XX=XX, lite=TRUE)
      lplez[,j] <- pnorm(0, pc$mean, sqrt(pc$s2), log.p=TRUE)
    }
    
    ## done
    return(data.frame(ei=ei, lplez=lplez))
  }


## mspeGP:
##
## wrapper used to calculate the MSPEs in C using
## the pre-stored GP representation.  

mspeGP <- function(gpi, Xcand, Xref=Xcand, fi=TRUE, verb=0)
  {
    m <- ncol(Xcand)
    if(ncol(Xref) != m) stop("Xcand and Xref have mismatched cols")
    ncand <- nrow(Xcand)

    out <- .C("mspeGP_R",
              gpi = as.integer(gpi),
              m = as.integer(m),
              Xcand = as.double(t(Xcand)),
              ncand = as.integer(ncand),
              Xref = as.double(t(Xref)),
              nref = as.integer(nrow(Xref)),
              fi = as.integer(fi),
              verb = as.integer(verb),
              mspes = double(ncand),
              PACKAGE = "laGP")
    
    return(out$mspes)
  }


## dmus2GP:
##
## obtain the derivative of the predictive scale
## of the fitted GP model

dmus2GP <- function(gpi, XX)
  {
    nn <- nrow(XX)

    out <- .C("dmus2GP_R",
              gpi = as.integer(gpi),
              m = as.integer(ncol(XX)),
              nn = as.integer(nn),
              XX = as.double(t(XX)),
              mu = double(nn),
              dmu = double(nn),
              d2mu = double(nn),
              s2 = double(nn),
              ds2 = double(nn),
              d2s2 = double(nn),
              PACKAGE = "laGP")
      
    return(data.frame(mu=out$mu, dmu=out$dmu, d2mu=out$d2mu,
                      s2=out$s2, ds2=out$ds2, d2s2=out$d2s2))
  }


## efiGP:
##
## obtain the expected (approx) Fisher information for
## the fitted GP model; returns the absolute value (i.e.,
## determinant)

efiGP <- function(gpi, Xcand)
  {
    nn <- nrow(Xcand)

    out <- .C("efiGP_R",
              gpi = as.integer(gpi),
              m = as.integer(ncol(Xcand)),
              nn = as.integer(nn),
              Xcand = as.double(t(Xcand)),
              efi = double(nn),
              PACKAGE = "laGP")

    ## remove silly values
    return(out$efi)
  }
