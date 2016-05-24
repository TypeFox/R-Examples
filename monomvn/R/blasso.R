#******************************************************************************* 
#
# Estimation for Multivariate Normal Data with Monotone Missingness
# Copyright (C) 2007, University of Cambridge
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


## bridge:
##
## Bayesian Ridge regression.  Simply calls the blasso function with
## the ridge argument modified to be TRUE and the rd argment of c(0,0)
## for a Jeffrey's prior on the squared ridge penalty lambda2

'bridge' <-
function(X, y, T=1000, thin=NULL, RJ=TRUE, M=NULL, beta=NULL,
         lambda2=1, s2=var(y-mean(y)), mprior=0, rd=NULL,
         ab=NULL, theta=0, rao.s2=TRUE, icept=TRUE, normalize=TRUE,
         verb=1)
  {
    blasso(X=X, y=y, T=T, thin=thin, RJ=RJ, M=M, beta=beta,
           lambda2=lambda2, s2=s2, case="ridge", mprior=mprior,
           rd=rd, ab=ab, theta=theta, rao.s2=rao.s2, icept=icept,
           normalize=normalize, verb=verb)

    ## need to rename the outputs and re-work the S3 functions
    ## appropriately (to make ridge-specific)
  }


## bhs:
##
## Bayesian Horseshoe regression.  Simply calls the blasso function with
## the ridge argument modified to be TRUE and the rd argment of c(0,0)
## for a Jeffrey's prior on the squared ridge penalty lambda2

'bhs' <-
function(X, y, T=1000, thin=NULL, RJ=TRUE, M=NULL, beta=NULL,
         lambda2=1, s2=var(y-mean(y)), mprior=0, ab=NULL,
         theta=0, rao.s2=TRUE, icept=TRUE, normalize=TRUE, verb=1)
  {
    blasso(X=X, y=y, T=T, thin=thin, RJ=RJ, M=M, beta=beta,
           lambda2=lambda2, s2=s2, case="hs", mprior=mprior,
           rd=c(0,0), ab=ab, theta=theta, rao.s2=rao.s2, icept=icept,
           normalize=normalize, verb=verb)

    ## need to rename the outputs and re-work the S3 function
    ## appropriately (to make horseshoe-specific)
  }



## blasso:
##
## function for sampling from the posterior distribution of the
## regression coefficients (beta) and variance (s2) under the
## Bayesian Lasso linear model of Park & Casella

'blasso' <-
function(X, y, T=1000, thin=NULL, RJ=TRUE, M=NULL, beta=NULL,
         lambda2=1, s2=var(y-mean(y)), case=c("default", "ridge", "hs", "ng"),
         mprior=0, rd=NULL, ab=NULL, theta=0, rao.s2=TRUE, icept=TRUE,
         normalize=TRUE, verb=1)
  {
    ## (quitely) double-check that blasso is clean before-hand
    blasso.cleanup()
    
    ## what to do if fatally interrupted?
    on.exit(blasso.cleanup())

    ## dimensions of the inputs
    X <- as.matrix(X)
    m <- ncol(X)
    n <- nrow(X)

    ## save the call
    cl <- match.call()
    
    ## check/conform arguments
    X <- as.matrix(X)
    y <- as.numeric(y)
    if(length(y) != nrow(X))
      stop("must have nrow(X) == length(y)")

    ## check T
    if(length(T) != 1 || T <= 1)
      stop("T must be a scalar integer > 1")

    ## check RJ (reversible jump)
    if(length(RJ) != 1 || !is.logical(RJ))
      stop("RJ must be a scalar logical\n")

    ## check lambda2
    if(length(lambda2) != 1 || lambda2 < 0)
      stop("lambda2 must be a non-negative scalar")

    ## check case
    case <- match.arg(case)
    if(case == "ridge" && lambda2 == 0)
      stop("specifying case=\"ridge\" and lambda2=0 doesn't make any sense")
    if(case == "ng") gamma <- rep(2, T)
    else gamma <- double(0)
    
    ## check M or default
    if(is.null(M)) {
      if(RJ) M <- as.integer(min(m, n-1))
      else M <- m
    }
    M <- as.integer(M)
    if(length(M) != 1 || M < 0 || M > m)
      stop("M must be a positive integer 0 <= M <= ncol(X)")
    if(!RJ && M != m) {
      M <- m
      warning("must have M=", M, " == ncol(X)=", m, " when RJ=FALSE",
              immediate.=TRUE)
    }
    
    ## check beta or default
    if(is.null(beta)) beta <- rep(!RJ, m)
    if(length(beta) != m)
      stop("must have length(beta) == ncol(X)")

    ## general big-p small-n problem -- must have lasso on
    if(lambda2 == 0 && m >= n)
      stop("big p small n problem; must have lambda2 > 0")
    else if(lambda2 == 0) lambda2 <- double(0)
    else lambda2 <- as.double(rep(lambda2, T))
    
    ## check for a valid regression
    if(!RJ) {  ## when not doing reversible jump (RJ)
      if(length(lambda2) > 0 && any(beta == 0)) {
        warning("must start with non-zero beta when RJ=FALSE, using beta=1\n",
                immediate.=TRUE)
        beta  <- rep(!RJ, m)
      }
    } else if(sum(beta != 0) > M) { ## when doing RJ
      beta <- rep(0, m)
      warning("initial beta must have M or fewer non-zero entries",
              immediate.=TRUE)
    }

    ## check thin or default
    if(is.null(thin)) {
      ## thin a bit for the Lasso latent variables
      if(RJ || (length(lambda2) > 0 && case!="ridge")) thin <- M
      else if(length(lambda2) > 0) thin <- 2
      else thin <- 1
      ## thin a bit for the Student-t latent variables
      if(theta > 0) thin <- thin + n
    }
    if(length(thin) != 1 || thin < 1)
      stop("thin must be a scalar integer >= 1")

    ## check s2
    if(length(s2) != 1 || s2 <= 0)
      stop("s2 must be a positive scalar")

    ## check tau2i or default
    if(case=="ridge" || length(lambda2) == 0) tau2i <- double(0)
    else {
      tau2i <- rep(1, m)
      tau2i[beta == 0] <- -1
      tau2i <- as.double(rep(tau2i, T))
    }

    ## check mprior and allocate pi appropriately
    if(any(mprior < 0)) stop("must have all(0 <= mprior)");
    if(length(mprior) == 1) {
      if(mprior != 0 && RJ == FALSE)
        warning(paste("setting mprior=", mprior,
                      " ignored since RJ=FALSE", sep=""))
      if(mprior > 1) stop("must have scalar 0 <= mprior < 1")
      mprior <- c(mprior, 0)
    } else if(length(mprior) != 2)
      stop("mprior should be a scalar or 2-vector in [0,1]")
    if(mprior[2] != 0) pi <- as.double(rep(mprior[1]/(mprior[1]+mprior[2]), T))
    else pi <- double(0)
    
    ## check r and delta (rd)
    if(is.null(rd)) {
      if(case=="ridge") { ## if using ridge regression IG prior
        rd <- c(0,0)
        ## big-p small-n setting for ridge
        if(m >= n) rd <- c(5, 10) 
      } else rd <- c(2, 0.1) ## otherwise lasso G prior
    } 
    ## double-check rd
    if(length(rd) != 2 || (length(tau2i) > 0 && any(rd <= 0))) {
      if(length(rd) == 1 && rd == FALSE) rd <- c(-1,-1) ## fixed lambda2
      else { ## actual error handling
        if(case != "hs") stop("rd must be a positive 2-vector")
        if(case == "ng" && rd[1] != 2) stop("must have rd[1] = 2 for NG prior")
      }
    }

    ## check ab or default
    if(is.null(ab) || all(ab == -1)) {
      if(all(ab != -1)) ab <- c(0,0)
      if(all(ab == -1) || (!RJ && lambda2 > 0 && m >= n)) {
        ab[1] <- 3/2
        ab[2] <- Igamma.inv(ab[1], 0.95*gamma(ab[1]), lower=FALSE)*sum(y^2)
      }
    }

    ## double check ab
    if(length(ab) != 2 || any(ab < 0))
      stop("ab must be a non-negative 2-vector")
    if(case!="ridge" && !RJ && m >= n && any(ab <= 0))
      stop("must have ab > c(0,0) when case!=\"ridge\", !RJ, and ncol(X) >= length(y)")

    ## check theta and possibly allocate omega2
    if(length(theta) != 1 || theta < 0)
      stop("theta must be a non-negative scalar")
    if(theta > 0) {
      omega2 <- as.double(rep(rep(1,n), T))
      nu <- as.double(rep(1.0/theta, T))
    } else nu <- omega2 <- double(0)

    ## check rao.s2
    if(length(rao.s2) != 1 || !is.logical(rao.s2))
      stop("rao.s2 must be a scalar logical")
    
    ## check normalize
    if(length(normalize) != 1 || !is.logical(normalize))
      stop("normalize must be a scalar logical")

    ## check icept
    if(length(icept) != 1 || !is.logical(icept))
      stop("icept must be a scalar logical")

    ## check verb
    if(length(verb) != 1 || verb < 0)
      stop("verb must be non-negative a scalar integer")

    ## call the C routine
    r <- .C("blasso_R",
            T = as.integer(T),
            thin = as.integer(thin),
            cols = as.integer(m),
            n = as.integer(n),
            X = as.double(t(X)),
            y = as.double(y),
            lambda2.len = as.integer(length(lambda2)),
            lambda2 = lambda2,
            gamma.len = as.integer(length(gamma)),
            gamma = gamma,
            mu = double(T),
            RJ = as.integer(RJ),
            M = as.integer(M),
            beta = as.double(rep(beta, T)),
            m = as.integer(rep(sum(beta!=0), T)),
            s2 = as.double(rep(s2, T)),
            tau2i.len = as.integer(length(tau2i)),
            tau2i = tau2i,
            hs = as.integer(case == "hs"),
            omega2.len = as.integer(length(omega2)),
            omega2 = omega2,
            nu.len = as.integer(length(nu)),
            nu = nu,
            pi.len = as.integer(length(pi)),
            pi = pi,
            lpost = double(T),
            llik = double(T),
            llik.norm.len = as.integer(T * (theta != 0)),
            llik.norm = double(T * (theta != 0)),
            mprior = as.double(mprior),
            r = as.double(rd[1]),
            delta = as.double(rd[2]),
            a = as.double(ab[1]),
            b = as.double(ab[2]),
            theta = as.double(theta),
            rao.s2 = as.integer(rao.s2),
            icept = as.integer(icept),
            normalize = as.integer(normalize),
            verb = as.integer(verb),
            PACKAGE = "monomvn")

    ## copy the inputs back into the returned R-object
    r$X <- X
    r$y <- y

    ## turn the beta vector of samples into a matrix
    r$beta <- matrix(r$beta, nrow=T, ncol=m, byrow=TRUE,
                     dimnames=list(NULL,paste("b.", 1:m, sep="")))

    ## turn the tau2i vector of samples into a matrix
    if(r$lambda2[1] != 0 && length(r$tau2i) > 0) {
      r$tau2i <- matrix(r$tau2i, nrow=T, ncol=m, byrow=TRUE,
                        dimnames=list(NULL,paste("tau2i.", 1:m, sep="")))
      ## put NAs where tau2i has -1
      r$tau2i[r$tau2i == -1] <- NA
    } else if(length(r$tau2i) > 0) {
      r$lambda <- r$tau2i <- NULL 
    } else r$tau2i <- NULL

    ## turn the omega2 vector of samples into a matrix
    if(r$theta > 0) 
      r$omega2 <- matrix(r$omega2, nrow=T, ncol=n, byrow=TRUE,
                     dimnames=list(NULL,paste("omega2.", 1:n, sep="")))
    else { r$theta <- NULL; r$omega2 <- NULL; r$nu <- NULL }
        
    ## first llik and lpost not available
    if(theta == 0) r$llik.norm <- NULL
    else r$llik.norm[1] <- NA
    r$llik[1] <- r$lpost[1] <- NA

    ## make logicals again
    r$normalize = as.logical(r$normalize)
    r$icept = as.logical(r$icept)
    r$RJ <- as.logical(r$RJ)
    r$rao.s2 <- as.logical(r$rao.s2)
    if(case != "ridge") r$hs <- as.logical(r$hs)
    else r$hs <- NULL

    ## transform mprior and pi
    if(r$mprior[2] == 0) {
      r$mprior <- r$mprior[-2]
      r$pi <- NULL
    }

    ## null-out redundancies
    r$pi.len <- r$col <- r$n <- r$cols <- r$verb <- NULL
    r$lambda2.len <- r$gamma.len <- r$tau2i.len <- r$omega2.len <- r$nu.len <- NULL
    if(length(r$lambda2) == 0) r$lambda2 <- NULL
    if(length(r$gamma) == 0) r$gamma <- NULL

    ## assign call and class
    r$call <- cl
    class(r) <- "blasso"
    
    return(r)
  }


## blasso.cleanup
##
## gets called when the C-side is aborted by the R-side and enables
## the R-side to clean up the memory still allocaed to the C-side,
## as well as whatever files were left open on the C-side

"blasso.cleanup" <-  function()
{
  .C("blasso_cleanup", PACKAGE="monomvn")
}
