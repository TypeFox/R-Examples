#******************************************************************************* 
#
# Particle Learning of Gaussian Processes
# Copyright (C) 2010, University of Cambridge
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
# Questions? Contact Robert B. Gramacy (bobby@statslab.cam.ac.uk)
#
#*******************************************************************************


## lpredprob.GP:
##
## for the PL resample step -- evaluate the (log) predictive
## density of the new z=(x,y) point given suff stats
## and params Zt up to time t

lpredprob.GP <- function(z, Zt, prior)
  {
    ## get the parameters to the predictive equations
    tp <- pred.GP(z$x, Zt, prior, PL.env$pall$Y)
    
    ## calculate the predictive probability
    return(dt((z$y - tp$m)/sqrt(tp$s2), tp$df, log=TRUE))
  }


## propagate.GP:
##
## for the PL propagate step -- add z to Zt and
## calculate the relevant updates to the sufficient
## statistics

propagate.GP <- function(z, Zt, prior)
  {
    ## increment t
    Zt$t <- Zt$t + 1
    
    ## call the update function on new data
    Zt <- updat.GP(Zt, prior, PL.env$pall$Y)
    
    ## propose changes to the correlation parameters
    Zt <- draw.GP(Zt, prior, l=3, h=4, thin=1)

    ## return the propagated particle
    return(Zt)
  }


## prior.GP:
##
## default prior specification for the GP model

prior.GP <- function(m, cov=c("isotropic", "separable", "sim"))
  {
    cov <- match.arg(cov)
    prior <- list(bZero=FALSE, s2p=c(0,0), grate=20, cov=cov)
    if(cov == "isotropic") prior$drate <- 5
    else {
      if(m == 1) stop("use isotropic when m=1")
      prior$drate <- rep(5, m)
      if(cov == "sim") {
        prior$bZero <- TRUE
        prior$drate <- sqrt(1/prior$drate)
      }
    }
    return(prior)
  }


## unif.propose.pos:
##
## proposals from uniform sliding winows over the
## positive reals

unif.propose.pos <- function(x, l=3, h=4)
{
  x.new <- runif(length(x), l/h*x, h/l*x)
  lfwd <- dunif(x.new, l/h*x, h/l*x, log=TRUE)
  lbak <- dunif(x, l/h*x.new, h/l*x.new, log=TRUE)
  return(list(x=x.new, lfwd=sum(lfwd), lbak=sum(lbak)))
}


## mvnorm.propose.rw
##
## random walk multivariate normal proposal for the
## sim index parameter for the sim covariance function

mvnorm.propose.rw <- function(x, cov=diag(0.2, length(x)))
  {
    
    x.new <- rmvnorm(1, mean=x, sigma=cov)
    return(list(x=x.new, lfwd=0, lbak=0))
  }


## draw.GP:
##
## MH-style draw for the range (d) and nugget (g) parameters
## do the correlation function (K) -- THIS FUNCTION IS VERY
## INEFFICIENT via RE-INIT and UPDATE

draw.GP <- function(Zt, prior, l=3, h=4, thin=10, Y=NULL)
  {
    ## perhaps initialize instead of draw
    if(is.null(Zt)) return(init.GP(prior=prior))
    if(any(prior$drate < 0) && prior$grate < 0) return(Zt)
    
    ## determine which Y to use
    if(is.null(Y)) Y <- PL.env$pall$Y
    else Zt <- updat.GP(Zt, prior, Y)

    ## increase thin in sim case
    if(prior$cov == "sim") thin <- thin * length(Zt$d)
    
    ## take thin draws from the posterior
    for(i in 1:thin) {

      ## check if sampling d
      if(all(prior$drate > 0)) {
      ## propose a change to d
        if(prior$cov == "sim") {
          d.new <- mvnorm.propose.rw(Zt$d)
          d.new$x <- d.new$x * sample(c(-1,1), length(d.new$x), replace=TRUE)
        } else d.new <- unif.propose.pos(Zt$d, l, h)
        
        ## accept or reject
        if(all(d.new$x != 0)) {
          Zt.prop <- init.GP(prior, d.new$x, Zt$g, Y)
          
          ## accept or reject d-change via MH
          if(runif(1) < exp(Zt.prop$lpost + d.new$lbak - Zt$lpost - d.new$lfwd))
            Zt <- Zt.prop
        }
      }

      ## propose a change to g
      if(prior$grate < 0) next; ## skip if freezing nugget
      g.new <- unif.propose.pos(Zt$g, l, h)
      if(all(g.new$x > 0)) {
        Zt.prop <- init.GP(prior, Zt$d, g.new$x, Y)
        
        ## accept or reject g-change via MH
        if(runif(1) < exp(Zt.prop$lpost + g.new$lbak - Zt$lpost - g.new$lfwd)) 
          Zt <- Zt.prop
      } 
    }

    ## return the possibly updated Zt
    return(Zt)
  }


## lpost.GP:
##
## calculate the log posterior probability of the
## GP model in the Zt particle

lpost.GP <- function(n, m, Vb, phi, ldetK, d, g, prior)
  {
    ## sanity check
    if(any(prior$s2p < 0)) stop("bad s2p")
    
    ## adjust m if not using linear model
    if(prior$bZero) m <- 0
    s2p <- prior$s2p
    
    ## deal with the prior (uniform in the range)
    lp <- 0
    if(all(prior$drate > 0)) lp <- lp + sum(dexp(abs(d), rate=prior$drate, log=TRUE))
    if(prior$grate > 0) lp <- lp + dexp(abs(g), rate=prior$grate, log=TRUE)
    
    ## deal with the (integrated) likelihood part
    lp <- lp + 0.5* (det(Vb, log=TRUE) - (n-m)*log(2*pi))
    lp <- lp - 0.5*ldetK
    lp <- lp + lgamma(0.5*(s2p[1]+n-m))
    lp <- lp - 0.5*(s2p[1]+n-m)*log(0.5*(phi+s2p[2]))
    if(all(s2p > 0))
      lp <- lp + 0.5*s2p[1]*log(0.5*s2p[2]) - lgamma(0.5*s2p[1]);
    
    ## sanity check
    if(!is.finite(lp)) stop("bad log posterior")
    
    ## return the log posterior proabability
    return(lp)
  }


## updat.GP:
##
## recalculate all of the suffficient statistics
## (etc) for the X & Y information in Zt

updat.GP <- function(Zt, prior, Y)
  {
    ## re-calculate the covariance matrix
    util <- util.GP(Zt, prior, Y, ldetK=TRUE)

    ## update the calculation of the posterior probability
    Zt$lpost <- lpost.GP(Zt$t, ncol(PL.env$pall$X), util$Vb, util$phi,
                         util$ldetK, Zt$d, Zt$g, prior)

    return(Zt)
  }


## init.GP:
##
## create a new particle for data X & Y, which
## comprises the kitchen-sink (suff stats) variable

init.GP <- function(prior, d=NULL, g=NULL, Y=NULL)
  {
    ## determine which Y to use
    if(is.null(Y)) Y <- PL.env$pall$Y

    ## defaults from prior
    if(is.null(d)) {
      d <- abs(1/prior$drate)
      ## perhaps negate some for sim cov
      if(prior$cov == "sim")
        d <- d * sample(c(-1,1), length(d), replace=TRUE)
    }
    if(is.null(g)) g <- abs(1/prior$grate)

    ## initialize particle parameters
    Zt <- list(d=d, g=g)
    if(!is.null(PL.env$pall$Y)) Zt$t <- sum(!is.na(PL.env$pall$Y))
    else Zt$t <- nrow(PL.env$pall$X)

    ## calculate the sufficient statistics
    Zt <- updat.GP(Zt, prior, Y)

    ## return the newly allocated particle
    return(Zt)
  }


## pred.GP:
##
## obtain the parameters to the (Student-t) predictive
## distribution based upon a particle at the XX predictive
## input locations

pred.GP <- function(XX, Zt, prior, Y=NULL, quants=FALSE,
                       Sigma=FALSE, sub=1:Zt$t)
  {
    ## get the right X and Y
    if(is.null(Y)) Y <- PL.env$pall$Y
    
    ## coerse the XX input
    XX <- matrix(XX, ncol=ncol(PL.env$pall$X))

    ## allocate space for the Student-t parameters
    ## from each particle
    I <- nrow(XX)
    tm <- ts2 <- rep(NA, I)
    
    ## predictive degrees of freedom
    if(prior$bZero) m <- 0
    else m <- ncol(PL.env$pall$X)
    tdf <- Zt$t - m - 1

    ## utility for calculations below
    if(prior$cov == "isotropic") { ## isotropic
      k <- covar(X1=PL.env$pall$X[sub,], X2=XX, d=Zt$d, g=0)
      kk <- drop(covar(X1=XX, d=Zt$d, g=Zt$g))
    } else if(prior$cov == "separable") { ## separable
      k <- covar.sep(X1=PL.env$pall$X[sub,], X2=XX, d=Zt$d, g=0)
      kk <- drop(covar.sep(X1=XX, d=Zt$d, g=Zt$g))
    } else { ## sim (rank 1)
      k <- covar.sim(X1=PL.env$pall$X[sub,], X2=XX, d=Zt$d, g=0)
      kk <- drop(covar.sim(X1=XX, d=Zt$d, g=Zt$g))
    }

    ## utility for calculations below
    util <- util.GP(Zt, prior, Y, sub, retKi=TRUE)

    ## calculate predictive quantities using Ki
    ktKi <- t(k) %*% util$Ki
    ktKik <- ktKi %*% k
    ktKiYmFbmu <- drop(ktKi %*% util$YmFbmu)
    if(any(!is.finite(ktKiYmFbmu))) stop("bad ktKiYmFbmu")

    ## predictive mean
    FF <- cbind(1,XX)
    tm <- drop(FF %*% util$bmu + ktKiYmFbmu)
    if(any(!is.finite(tm))) stop("bad tm")

    ## predictive variance
    ## REALLY NEED TO CHECK WHICH ONES OF THESE IS CORRECT
    badj <- diag(1,nrow(XX)) + FF %*% util$Vb %*% t(FF)
    ## badj <- kk + FF %*% util$Vb %*% t(FF)
    zphi <- (prior$s2p[2] + util$phi)*(kk - ktKik)
    tSigma <- badj * zphi / (prior$s2p[1] + tdf)
    ## stop("blah")
    
    ## assemble the data frame
    if(Sigma) pred.df <- list(m=tm, Sigma=tSigma, df=tdf)
    else {
      ts2 <- diag(tSigma) ## extract variance
      ts2[ts2 < 0] <- 0   ## correct for roundoff error
      if(any(!is.finite(ts2))) stop("bad ts2")
      pred.df <- data.frame(m=tm, s2=ts2, df=tdf)
    }

    ## maybe append condnum
    if(!is.null(util$condnum)) {
      if(Sigma) pred.df$condnum <- rep(util$condnum, nrow(XX))
      else pred.df <- cbind(pred.df, data.frame(condnum=rep(util$condnum, nrow(pred.df))))
    }
    
    ## get quantiles
    if(quants) {
      if(Sigma) stop("cannot get quants when Sigma=TRUE")
      pred.df <- cbind(pred.df, tquants(pred.df))
    }

    ## return the df
    return(pred.df)
  }


## ieci.GP:
##
## use the 1-step-ahead predictive equations to calculate
## the IECI statistic under a GP model.  

ieci.GP <- function(Xcand, Xref, Zt, prior, Y=NULL, w=NULL, verb=1)
  {
    ## get the right Y
    if(is.null(Y)) Y <- PL.env$pall$Y
    
    ## coerse the Xcand & Xref input
    Xcand <- matrix(Xcand, ncol=ncol(PL.env$pall$X))
    Xref <- matrix(Xref, ncol=ncol(PL.env$pall$X))

    ## predictive degrees of freedom
    if(prior$bZero) m <- 0
    else m <- ncol(PL.env$pall$X)
    tdf <- Zt$t - m - 1
    util <- util.GP(Zt, prior, Y, retKi=TRUE)

    ## utility for calculations below
    if(prior$cov == "isotropic") { ## isotropic
      k <- covar(X1=PL.env$pall$X, X2=Xref, d=Zt$d, g=0)#g=Zt$g)
    } else if(prior$cov == "separable") { ## separable
      k <- covar.sep(X1=PL.env$pall$X, X2=Xref, d=Zt$d, g=0)#, g=Zt$g)
    } else { ## sim (rank 1)
      k <- covar.sim(X1=PL.env$pall$X, X2=Xref, d=Zt$d, g=0)#, g=Zt$g)
    }

    ## build up the K quantities
    ktKi <- t(k) %*% util$Ki
    ktKik <- diag(ktKi %*% k)
    ktKiYmFbmu <- drop(ktKi %*% util$YmFbmu)
    if(any(!is.finite(ktKiYmFbmu))) stop("bad ktKiYmFbmu")

    ## predictive mean
    FF <- cbind(1,Xref)
    tm <- drop(FF %*% util$bmu + ktKiYmFbmu)
    if(any(!is.finite(tm))) stop("bad tm")

    ## this assumes that Xref is representative of PL.env$pall$X
    fmin <- min(tm)

    ## initial steps in the predictive variance calculation
    badj <- 1 + diag(FF %*% util$Vb %*% t(FF))  ## could be simplified further
    
    ## IECI calculation for each entry in Xcand
    ieci <- calc.iecis(ktKik, k, Xcand, PL.env$pall$X, util$Ki, Xref, Zt$d, Zt$g,
                       prior$s2p, util$phi, badj, tm, tdf, fmin, w, verb)

    ## sanity checks
    if(any(is.nan(ieci))) stop("NaN in ieci")
    
    ## return the df
    return(ieci)
  }


## alc.GP:
##
## use the 1-step-ahead predictive equations to calculate
## the alc statistic under a GP model.  

alc.GP <- function(Xcand, Xref, Zt, prior, Y=NULL, w=NULL, verb=1)
  {
    ## get the right Y
    if(is.null(Y)) Y <- PL.env$pall$Y
    
    ## coerse the Xcand & Xref input
    Xcand <- matrix(Xcand, ncol=ncol(PL.env$pall$X))
    Xref <- matrix(Xref, ncol=ncol(PL.env$pall$X))

    ## predictive degrees of freedom
    if(prior$bZero) m <- 0
    else m <- ncol(PL.env$pall$X)
    tdf <- Zt$t - m - 1
    util <- util.GP(Zt, prior, Y, retKi=TRUE)

    ## utility for calculations below
    if(prior$cov == "isotropic") { ## isotropic
      k <- covar(X1=PL.env$pall$X, X2=Xref, d=Zt$d, g=0)#g=Zt$g)
    } else if(prior$cov == "separable") { ## separable
      k <- covar.sep(X1=PL.env$pall$X, X2=Xref, d=Zt$d, g=0)#, g=Zt$g)
    } else { ## sim (rank 1)
      k <- covar.sim(X1=PL.env$pall$X, X2=Xref, d=Zt$d, g=0)#, g=Zt$g)
    }

    ## initial steps in the predictive variance calculation
    FF <- cbind(1,Xref)
    badj <- 1 + diag(FF %*% util$Vb %*% t(FF))  ## could be simplified further
    
    ## ALC calculation for each entry in Xcand
    alc <- calc.alcs(k, Xcand, PL.env$pall$X, util$Ki, Xref, Zt$d, Zt$g,
                     prior$s2p, util$phi, badj, tdf, w, verb)

    ## sanity checks
    if(any(is.nan(alc))) stop("NaN in alc")
    
    ## return the df
    return(alc)
  }


## tquants
##
## turns the (0.05, 0.95) quantiles of a Student-t
## distribution parameterized in tp

tquants <- function(tp)
  {
    q1 <- sqrt(tp$s2)*qt(0.05, tp$df) + tp$m
    q2 <- sqrt(tp$s2)*qt(0.95, tp$df) + tp$m
    return(data.frame(q1=q1, q2=q2))
  }


## params.GP:
##
## collects the GP parameters, range (d) and nugged (d)
## in a data frame for the purposes of making a histogram
## to track particle depletion

params.GP <- function()
  {
    ## allocate a data frame for the params
    P <- length(PL.env$peach)
    g <- lpost <- rep(NA, P)
    d <- matrix(NA, ncol=length(PL.env$peach[[1]]$d), nrow=P)
    
    ## collect the parameters from the particles
    for(p in 1:P) {
      d[p,] <- PL.env$peach[[p]]$d
      g[p] <- PL.env$peach[[p]]$g
      lpost[p] <- PL.env$peach[[p]]$lpost
    }

    ## return the particles
    return(data.frame(d=d, g=g, lpost=lpost))
  }


## rectscale:
##
## scale the input matrix (X) to lie in [0,1]^d by
## translating it according to the rectangle provided

rectscale <- function(X, rect)
  {
    X <- as.matrix(X)
    m <- ncol(X)
    if(!is.null(rect)) {
      rect <- matrix(rect, nrow=m)
      for(j in 1:m)
        X[,j] <- (X[,j]-rect[j,1])/(rect[j,2]-rect[j,1])
    }
    return(X)
  }


## rectunscale:
##
## un-scale the input matrix (X) to lie in [0,1]^d by
## translating it according to the rectangle provided

rectunscale <- function(X, rect)
  {
    X <- as.matrix(X)
    m <- ncol(X)
    if(!is.null(rect)) {
      rect <- matrix(rect, nrow=m)
      for(j in 1:m)
        X[,j] <- X[,j]*(rect[j,2]-rect[j,1]) + rect[j,1]
    }
    return(X)
  }


## data.GP:
##
## extract the appropriate columns from the X matrix
## and Y vector -- designed to be generic for other cases
## where we would want to get the next observation (end=NULL)
## or a range of observations from begin to end

data.GP <- function(begin, end=NULL, X, Y)
  {
    if(is.null(end) || begin == end)
      return(list(x=X[begin,], y=Y[begin]))
    else if(begin > end) stop("must have begin <= end")
    else return(list(x=as.matrix(X[begin:end,]), y=Y[begin:end]))
  }


## addpall.GP:
##
## add data to the pall data structure used as utility
## by all particles

addpall.GP <- function(Z)
  {
    PL.env$pall$X <- rbind(PL.env$pall$X, Z$x)
    PL.env$pall$Y <- c(PL.env$pall$Y, Z$y)
    PL.env$pall$D <- NULL
  }


## util.GP:
##
## calculate the full set of quantities needed to evaluate
## the likelihood of the GP model

util.GP <- function(Zt, prior, Y=NULL, sub=1:Zt$t, ldetK=FALSE, retKi=TRUE)
  {
    ## sanity check
    if(!is.null(PL.env$pall$Y) && Zt$t != sum(!is.na(PL.env$pall$Y)))
      stop("Zt$t and sum(!is.na(PL.env$pall$Y)) mismatch")
    else if(Zt$t != nrow(PL.env$pall$X)) stop("Zt$t and nrow(PL.env$pall$X) mismatch")
    m <- ncol(PL.env$pall$X)

    ## get the right X and Y
    if(is.null(Y)) Y <- PL.env$pall$Y
    
    ## calculate the covariance matrix
    if(prior$cov == "isotropic") { ## isotropic
      if(is.null(PL.env$pall$D)) PL.env$pall$D <- distance(PL.env$pall$X)
      K <- dist2covar.symm(D=PL.env$pall$D[sub,sub], d=Zt$d, g=Zt$g)
    } else if(prior$cov == "separable") { ## separable
      K <- covar.sep(X1=PL.env$pall$X[sub,], d=Zt$d, g=Zt$g)
    } else { ## sim (rank 1)
      K <- covar.sim(X1=PL.env$pall$X[sub,], d=Zt$d, g=Zt$g)      
    }

    ## Ki <- K^{-1}
    Ki <- solve(K)

    ## calculate the weighted least squares part
    Y <- Y[sub]
    if(! prior$bZero) {
      F <- cbind(1, PL.env$pall$X[sub,])
      FtKi <- t(F) %*% Ki 
      Vb <- solve(FtKi %*% F)
      bmu <- Vb %*% (FtKi %*% Y)
    } else { bmu <- rep(0, m+1); Vb <- matrix(0, m+1, m+1) }

    ## re-calculate the phi (variance) suff stats
    phi <- drop(t(Y) %*% Ki %*% Y)
    if(!prior$bZero) phi <- drop(phi - t(bmu) %*% solve(Vb) %*% bmu)

    ## possibly calculate the log determinant of K
    if(ldetK) ldetK <- det(K, log=TRUE)
    else ldetK <- NULL

    ## first work out the regular predictive equations
    if(! prior$bZero) { YmFbmu <- drop(Y - F %*% bmu) }
    else { YmFbmu <- Y }

    ## decide whether or not to return Ki
    if(!retKi) { Ki <- FALSE; condnum <- NULL }
    else condnum <- kappa(Ki)
    
    ## return the list of quantities calculated
    return(list(Vb=Vb, bmu=bmu, phi=phi, Ki=Ki, condnum=condnum, ldetK=ldetK, YmFbmu=YmFbmu))
  }


## EI:
##
## calculate the expected improvement

EI <- function(tp, fmin)
  {
    tp <- as.data.frame(t(tp))
    diff <- fmin - tp$m
    tsd <- sqrt(tp$s2)
    diffs <- diff/tsd
    scale <- (tp$df*tsd + diff^2/tsd)/(tp$df-1)
    ei <- diff*pt(diffs, tp$df)
    ei <- ei + scale*dt(diffs, tp$df)
    return(ei)
  }


## ei.adapt:
##
## return the index into Xcand that has the most potential to
## improve the estimate of the minimum

ei.adapt <- function(Xcand, rect, prior, verb)
  {
    ## calculate the average maximum entropy point
    if(verb > 0)
      cat("taking design point ", nrow(PL.env$pall$X)+1, " by EI\n", sep="")

    ## stale candidates
    Xcands <- rectscale(Xcand, rect)

    ## get predictive distribution information
    outp <- papply(XX=Xcands, fun=pred.GP, prior=prior, verb=verb)

    ## gather the entropy info
    ei <- rep(0, nrow(as.matrix(Xcand)))
    for(p in 1:length(outp)) {
      fmin <- min(outp[[p]]$m)
      ei <- ei + calc.eis(outp[[p]][,1:3], fmin)
      ## ei <- ei + apply(outp[[p]][,1:3], 1, EI, fmin=fmin)
    }

    ## return the candidate with the most potential
    return(ei/length(outp))
  }


## var.adapt:
##
## return the index into Xcand that has the largest variance

var.adapt <- function(Xcand, rect, prior, verb)
  {
    ## calculate the average maximum entropy point
    if(verb > 0)
      cat("taking design point ", nrow(PL.env$pall$X)+1, " by EI\n", sep="")

    ## stale candidates
    Xcands <- rectscale(Xcand, rect)
    
    ## get predictive distribution information
    outp <- papply(XX=Xcands, fun=pred.GP, prior=prior, verb=verb)

    ## gather the predictive variance info
    condnum <- m <- m2 <- v <- rep(0, nrow(as.matrix(Xcand)))
    for(p in 1:length(outp)) {
      m2 <- m2 + outp[[p]]$m^2
      m <- m + outp[[p]]$m
      v <- v + outp[[p]]$df*outp[[p]]$s2/(outp[[p]]$df - 2)
      condnum <- condnum + outp[[p]]$condnum
    }

    ## return the candidate with the most potential
    n <- length(outp)
    return(list(m=m/n, v=v/n + m2/n - (m/n)^2, condnum=condnum/n))
  }


## icei.adapt
##
## return the index into Xcand that has the most potential to
## improve the estimate of the minimum via integrated expected conitional
## improvement

ieci.adapt <- function(Xcand, rect, prior, verb, Xref=NULL, fun=ieci.GP)
  {
    ## calculate the average maximum entropy point
    if(verb > 0)
      cat("taking design point ", nrow(PL.env$pall$X)+1, " by IECI\n", sep="")

    ## adjust the candidates (X) and reference locations (Xref)
    Xcands <- rectscale(Xcand, rect)
    if(is.null(Xref)) Xrefs <- Xcands
    else Xrefs <- rectscale(Xref, rect) 

    ## get predictive distribution information
    iecis <- papply(Xcand=Xcands, Xref=Xrefs, fun=fun, prior=prior, verb=verb)
    
    ## gather the entropy info for each x averaged over
    ieci <- rep(0, nrow(Xcands))
    for(p in 1:length(iecis)) ieci <- ieci + iecis[[p]]

    ## calculate mean EIs to in order to make a better progress meter
    mei <- mean(ei.adapt(Xref, rect, prior, verb=0))

    ## return the candidate with the most potential
    return(mei-ieci/length(iecis))
    ##return(-ieci/length(iecis))
  }


## alc.adapt:
##
## return the index into Xcand that has the most potential to
## improve predictive variance at reference locations

alc.adapt <- function(Xcand, rect, prior, verb, Xref=NULL, fun=alc.GP)
  {
    ## calculate the average maximum entropy point
    if(verb > 0)
      cat("taking design point ", nrow(PL.env$pall$X)+1, " by ALC\n", sep="")

    ## adjust the candidates (X) and reference locations (Xref)
    Xcands <- rectscale(Xcand, rect)
    if(is.null(Xref)) Xrefs <- Xcands
    else Xrefs <- rectscale(Xref, rect)

    ## get predictive distribution information
    alcs <- papply(Xcand=Xcands, Xref=Xrefs, fun=fun, prior=prior, verb=verb)

    ## gather the alc info for each x averaged over
    alc <- rep(0, nrow(Xcands))
    for(p in 1:length(alcs)) alc <- alc + alcs[[p]]

    ## calculate mean vars to in order to make a better progress meter
    ## mvar <- mean(var.adapt(Xref, rect, prior, verb=0)$v)

    ## return the candidate with the most potential
    return(alc/length(alcs))
    ## return(mvar-alc/length(alcs))
    ##return(-alc/length(alcs))
  }


## mindist.adapt:
##
## return the index into Xcand that is the minimum distance
## from Xref

mindist.adapt <- function(Xcand, rect, prior, verb, Xref)
  {
    ## calculate the average maximum entropy point
    if(verb > 0)
      cat("taking design point ", nrow(PL.env$pall$X)+1, " by mindist\n", sep="")

    ## check Xref should be length 1
    if(nrow(Xref) != 1) stop("Xref should have one vector")
    
    ## adjust the candidates (X) and reference locations (Xref)
    Xcands <- rectscale(Xcand, rect)
    if(is.null(Xref)) Xrefs <- Xcands
    else Xrefs <- rectscale(Xref, rect)

    ## calculate distances
    D <- 1/distance(Xrefs, Xcands)

    ## return the candidate with the most potential
    return(D)
  }


## data.GP.improv:
##
## use the current state of the particles to calculate
## the next adaptive sample from the posterior predictive
## distribution based on the expected improvement
##
## THIS FUNCTION SHOULD BE BROKEN IN TWO -- ONE FOR OPTIMIZATION
## AND ONE FOR LOCAL DESIGN WITH APLEY

data.GP.improv <- function(begin, end=NULL, f, rect, prior, adapt=ei.adapt,
                           cands=40, save=TRUE, oracle=TRUE, verb=2, interp=interp.loess)
  {
    if(!is.null(end) && begin > end) stop("must have begin <= end")
    else if(is.null(end) || begin == end) { ## adaptive sample

      ## choose some adaptive sampling candidates
      if(is.na(cands)) Xc <- PL.env$Xcand
      else Xc <- lhs(cands, rect)

      if(oracle) {  ## add a cleverly chosen candidate
        xstars <- findmin.GP(PL.env$pall$X[nrow(PL.env$pall$X),], prior)
        xstar <- drop(rectunscale(rbind(xstars), rect))
        Xc <- rbind(Xc, xstar)
      } else if(!is.null(formals(adapt)$Xref)) {
        ## calculate the predictive variance of the reference location
        mv <- var.adapt(formals(adapt)$Xref, rect, prior, verb=0)
      } else mv <- NULL

      ## calculate the index with the best entropy reduction potential
      as <- adapt(Xc, rect, prior, verb)
      indx <- which.max(as)
      
      ## return the new adaptive sample
      x <- matrix(Xc[indx,], nrow=1)
      xs <- rectscale(x, rect)

      ## possibly remove the candidate from a fixed set
      if(is.na(cands)) PL.env$Xcand <- PL.env$Xcand[-indx,]
      
      ## maybe plot something
      if(verb > 1) {
        par(mfrow=c(1,1))
        if(ncol(Xc) > 1 && nrow(PL.env$pall$X) > 10) { ## 2-d+ data
          image(interp(Xc[,1], Xc[,2], as), main="AS surface")
          points(rectunscale(PL.env$pall$X, rect))
          points(Xc, pch=18)
          if(oracle) points(xstar[1], xstar[2], pch=17, col="blue")
          points(x[,1], x[,2], pch=18, col="green")
        } else if(ncol(Xc) == 1){ ## 1-d data
          o <- order(drop(Xc))
          plot(drop(Xc[o,]), as[o], type="l", lwd=2,
               xlab="x", ylab="AS stat", main="AS surface")
          points(drop(rectunscale(PL.env$pall$X, rect)), rep(min(as), nrow(PL.env$pall$X)))
          points(x, min(as), pch=18, col="green")
          if(oracle) points(xstar, min(as), pch=17, col="blue")
          legend("topright", c("chosen point", "oracle candidate"),
                 pch=c(18,17), col=c("green", "blue"), bty="n")
        }
      }

      ## maybe save the max log EI or IECI
      if(save) {
        if(oracle) PL.env$psave$xstar <- rbind(PL.env$psave$xstar, xstar)
        else if(!is.null(mv)) {
          PL.env$psave$vref <- rbind(PL.env$psave$vref, mv$v)
          PL.env$psave$mref <- rbind(PL.env$psave$mref, mv$m)
          PL.env$psave$condnum <- rbind(PL.env$psave$condnum, mv$condnum)
        }
        PL.env$psave$max.as <- c(PL.env$psave$max.as, max(as))
      }

      ## return the adaptively chosen location
      return(list(x=xs, y=f(x)))

    } else {  ## create an initial design

      ## MED from a subset when cands is NULL, uses global Xcand
      if(is.na(cands)) {
        if(verb > 0) cat("initializing with size", end-begin+1, "MED\n")
        out <- dopt.gp(end-begin+1, X=NULL, Xcand=PL.env$Xcand)
        X <- out$XX
        PL.env$Xcand <- PL.env$Xcand[-out$fi,]
      }  else {  ## calculate a LHS 
        if(verb > 0) cat("initializing with size", end-begin+1, "LHS\n")
        X <- lhs(end-begin+1, rect)
      }
      
      ## get the class labels
      Y <- f(X)
      return(list(x=rectscale(X, rect), y=Y))
    }
  }


## pred.mean.GP:
##
## calculate the predictive mean given the util
## structure for the particle and range and nuggets

pred.mean.GP <- function(x, util, cov, dparam, gparam)
  {
    ## coerse the XX input
    XX <- matrix(x, ncol=ncol(PL.env$pall$X))

    ## calculate covariance vectors and matrices
    if(cov == "isotropic") { ## isotropic
      k <- covar(X1=PL.env$pall$X, X2=XX, d=dparam, g=0)
      kk <- drop(covar(X1=XX, d=dparam, g=gparam))
    } else if(cov == "separable") { ## separable
      k <- covar.sep(X1=PL.env$pall$X, X2=XX, d=dparam, g=0)
      kk <- drop(covar.sep(X1=XX, d=dparam, g=gparam))
    } else { ## sim (rank 1)
      k <- covar.sim(X1=PL.env$pall$X, X2=XX, d=dparam, g=0)
      kk <- drop(covar.sim(X1=XX, d=dparam, g=gparam))
    }
    
    ## calculate predictive quantities using Ki
    ktKi <- t(k) %*% util$Ki
    ktKik <- ktKi %*% k
    ktKiYmFbmu <- drop(ktKi %*% util$YmFbmu)
    if(any(!is.finite(ktKiYmFbmu))) stop("bad ktKiYmFbmu")
    
    ## predictive mean
    FF <- cbind(1,XX)
    tm <- drop(FF %*% util$bmu + ktKiYmFbmu)
    if(any(!is.finite(tm))) stop("bad tm")
    
    ## print(c(x,tm))
    ## return predictive mean
    return(drop(tm))
  }


## getmap.GP:
##
## return MAP particle

getmap.GP <- function()
  {
    ## calculate the MAP particle
    mi <- 1
    if(length(PL.env$peach) > 1) {
      for(p in 2:length(PL.env$peach))
        if(PL.env$peach[[p]]$lpost > PL.env$peach[[mi]]$lpost) mi <- p
    }
    return(PL.env$peach[[mi]])
  }


## findmin.GP:
##
## find the minimum of the predictive surface for the MAP particle

findmin.GP <- function(xstart, prior)
  {
    ## calculate the MAP particle
    Zt <- getmap.GP()
    
    ## utility for calculations below
    util <- util.GP(Zt, prior, PL.env$pall$Y, retKi=TRUE)

    ## call the optim function
    m <- ncol(PL.env$pall$X)
    xstar <- optim(xstart, pred.mean.GP, method="L-BFGS-B",
                   lower=rep(0,m), upper=rep(1,m), util=util,
                   cov=prior$cov, dparam=Zt$d, gparam=Zt$g)$par
    return(xstar)
  }
