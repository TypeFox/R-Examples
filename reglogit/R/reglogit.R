#*******************************************************************************
#
# Regularized logistic regression
# Copyright (C) 2011, The University of Chicago
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


## functions for Gibbs sampling the parameters to the
## regularized logit regression model

## rmultnorm
##
## faster version of rmvnorm from the mvtnorm package

rmultnorm <- function(n,mu,sigma)
  {
    p <- length(mu)
    z <- matrix(rnorm(n * p),nrow=n)
    ch <- chol(sigma,pivot=T)
    piv <- attr(ch,"pivot")
    zz <- (z%*%ch)
    zzz <- 0*zz
    zzz[,piv] <- zz
    zzz + matrix(mu,nrow=n,ncol=p,byrow=T)
}


## draw.beta:
##
## draw the regression coefficients conditional on the
## latent zs, lambdas, and omegas, and the penalty
## parameter nu, and thermo parameter kappa

draw.beta <- function(X, z, lambda, omega, nu, sigma,
                      kappa, kp)
  {
    ## negative nu and no sigma indicates L2 prior if sigma not null
    if(nu < 0 && !is.null(sigma)) n <- 1
    
    ## covariance of the conditional
    tXLi <- t(X * (1/lambda))
    Bi <- tXLi %*% X
    if(nu > 0) {
      bdi <- (kp/nu)^2 * 1/(sigma^2 * omega)
      if(ncol(X) == length(sigma) + 1) bdi <- c(0, bdi)
      if(class(Bi)[1] == "dgCMatrix")
        Bi <- Bi + .sparseDiagonal(nrow(Bi)) * bdi
      else diag(Bi) <- diag(Bi) + bdi
    } ## else Bi <- 0
    B <- base::solve(Bi) # + tXLi %*% X)

    ## mean of the conditional
    if(is.null(z)) { ## z=0 PDF representation

      ## choices between different (a,b) parameterizations
      a  <- 0.5; b <- kappa - 0.5
      kappa <- (a - 0.5*(a-b))

      ## special handling in Binomial case
      if(length(kappa) == 1) uk <- colSums(X*kappa)
      else uk <- drop(t(kappa) %*% X)
      b <- B %*% uk

    } else ## CDF representation
      b <- B %*% tXLi %*% (z - 0.5*(1-kappa)*lambda)

    ## draw from the MVN
    ## return(rmvnorm(1, b, B))
    return(rmultnorm(1, b, B))
  }


## draw.nu:
##
## draw from the full conditional distribution of the
## penalty parameter nu

draw.nu <- function(beta, sigma, kp, nup=list(a=0, b=0))
  {
    ## intercept?
    p <- length(sigma)
    if(p == length(beta) - 1) beta <- beta[-1]

    ## temper the prior
    a <- kp*(nup$a+1)-1
    b <- kp*nup$b

    ## update using beta
    a <- a + kp*p
    ## a <- a + p
    b <- b + kp*sum(abs(beta)/sigma)

    ## return a sample from an inverse gamma
    nu <- 1.0/rgamma(1, shape=a, scale=1/b)

    ## return sampled nu
    return(nu)
  }


## draw.z:
##
## draw the latent z-values given beta, y, and the
## latent lambdas

draw.z <- function(X, beta, lambda, kappa, C=NULL)
{
  n <- nrow(X)
  xbeta <- X %*% beta
  if(!is.null(C)) xbeta <- xbeta - log(C)
  
  return(.C("draw_z_R",
            n = as.integer(n),
            xbeta = as.double(xbeta),
            beta = as.double(beta),
            lambda = as.double(lambda),
            kappa = as.double(kappa),
            kl = as.integer(length(kappa)),
            z = double(n),
            PACKAGE = "reglogit")$z)
}


## draw.omega:
##
## draw the latent omega variables given the betas
## and sigmas

draw.omega <- function(beta, sigma, nu, kp)
  {
    if(any(nu < 0)) stop("nu must be positive")

    ## deal wth intercept in prior?
    if(length(sigma) == length(beta) - 1) beta <- beta[-1]
    m <- length(beta)
    
    return(.C("draw_omega_R",
              m = as.integer(m),
              beta = as.double(beta),
              sigma = as.double(sigma),
              nu = as.double(nu/kp),
              omega = double(m),
              PACKAGE = "reglogit")$omega)
  }



## draw.lambda:
##
## draw the latent lambda variables that help us
## implement the logit link, via Metropolis-Hastings

draw.lambda <- function(X, beta, lambda, kappa, kmax,
                        zzero=TRUE, thin=kappa, C=NULL)
  {
    ## calculate X * beta
    xbeta <- X %*% beta
    if(!is.null(C)) xbeta <- xbeta - log(C)
    n <- length(xbeta)

    ## call C code
    return(.C("draw_lambda_R",
              n = as.integer(n),
              xbeta = as.double(xbeta),
              kappa = as.double(kappa),
              kl = as.integer(length(kappa)),
              kmax = as.integer(kmax),
              zzero = as.integer(zzero),
              thin = as.integer(thin),
              tl = as.integer(length(thin)),
              lambda = as.double(lambda),
              PACKAGE = "reglogit")$lambda)
  }



## my.rinvgauss:
##
## draw from an Inverse Gaussian distribution following
## Gentle p193

my.rinvgauss <- function(n, mu, lambda)
  {
    return(.C("rinvgauss_R",
              n = as.integer(n),
              mu = as.double(mu),
              lambda = as.double(lambda),
              x = double(n),
              PACKAGE = "reglogit")$x)
  }


## calc.lpost:
##
## calculate the log posterior probability of the
## parameters in argument

calc.lpost <- function(yX, beta, nu, kappa, kp, sigma,
                       nup=list(a=0, b=0))
  {
    ## likelihood part
    llik <- - sum(kappa*log(1 + exp(-drop(yX %*% beta))))

    ## prior part
    
    ## deal with intercept in prior
    p <- length(sigma)
    if(p == length(beta) - 1) beta <- beta[-1]

    ## for beta
    if(nu > 0) lprior <- - kp*p*log(nu) - (kp/nu) * sum(abs(beta/sigma))
    ## if(nu > 0) lprior <- - p*log(nu) - (kp/nu) * sum(abs(beta/sigma))
    ## else if(!is.null(sigma)) lprior <- - sum((kp*beta/sigma)^2)
    else lprior <- 0

    ## inverse gamma prior for nu
    if(!is.null(nup))
      ## lprior <- lprior - kp*(2*(nup$a+1)*log(nu) + nup$b/(nu^2))
      lprior <- lprior - kp*((nup$a+1)*log(nu) + nup$b/nu)

    ## return log posterior
    return(llik + lprior)
  }


## calc.mlpost:
##
## calculate the log posterior probability of the
## parameters in argument, for multinomial reglogit

calc.mlpost <- function(yX, beta, nu, kappa, kp, sigma,
                       nup=list(a=0, b=0))
  {
    ## deal with vector beta
    if(is.null(dim(beta))) beta <- matrix(beta, ncol=dim(yX)[3])
    
    ## likelihood part
    if(!is.null(dim(kappa))) kappa <- rowSums(kappa)
    Cs <- rep(1, nrow(yX))
    Qm1 <- ncol(beta)
    for(q in 1:Qm1) Cs <- Cs + exp(-drop(yX[,,q] %*% beta[,q]))
    llik <- - sum(kappa*log(Cs))

    ## prior part
    
    ## deal with intercept in prior
    p <- length(sigma)
    if(p == nrow(beta) - 1) beta <- beta[-1,,drop=FALSE]

    ## for beta
    lprior <- 0.0
    if(any(nu > 0))
      for(q in 1:Qm1) {
        lprior <- lprior - kp*p*log(nu[q]) - (kp/nu[q]) * sum(abs(beta[,q]/sigma))
      }
    ## if(nu > 0) lprior <- - p*log(nu) - (kp/nu) * sum(abs(beta/sigma))
    ## else if(!is.null(sigma)) lprior <- - sum((kp*beta/sigma)^2)
    else lprior <- 0

    ## inverse gamma prior for nu
    if(!is.null(nup))
      for(q in 1:Qm1) {
      ## lprior <- lprior - kp*(2*(nup$a+1)*log(nu) + nup$b/(nu^2))
        lprior <- lprior - kp*((nup$a+1)*log(nu[q]) + nup$b/nu[q])
      }

    ## return log posterior
    return(llik + lprior)
  }


## calc.Cs
##
## for calculating the sum of the likelihood contributions
## for the Q-2 other reglogits in the multinomial implementation
## of Holmes and Held

## maybe make C version
calc.Cs <- function(yX, beta, lambda, kappa)
  {
    Qm1 <- dim(yX)[3]
    if(ncol(kappa) != Qm1) stop("wups")
    Cs <- matrix(1, nrow=nrow(lambda), ncol=Qm1)
    for(q in 1:Qm1) {
      for(nq in setdiff(1:Qm1, q)) {
        Cs[,q] <- Cs[,q] + exp(yX[,,nq] %*% beta[,nq] -
          0.5*(1-drop(kappa[,nq]))*lambda[,nq])
      }
    }
    return(Cs)
  }


## preprocess:
##
## pre-process the X and y and kappa data depending on
## how binomial responses might be processed

preprocess <- function(X, y, N, flatten, kappa)
  {
    ## process y's and set up design matrices
    if(is.null(N) || flatten == TRUE) { 
      
      if(flatten) { ## flattened binomial case
        if(is.null(N)) stop("flatten only applies for N != NULL")
        if(any(y > N)) stop("must have y <= N")
        ye <- NULL
        Xe <- NULL
        for(i in 1:length(N)) {
          Xe <- rbind(Xe, matrix(rep(X[i,], N[i]), nrow=N[i], byrow=TRUE))
          ye <- c(ye, rep(1, y[i]), rep(0, N[i]-y[i]))
        }
        y <- ye; X <- Xe
      } else if(any(y > 1)) stop("must have y in {0,1}")

      ## create y in {-1,+1}
      ypm1 <- y
      ypm1[y == 0] <- -1
      
      ## create y *. X
      yX <- X * ypm1 
      
    } else {  ## unflattened binomial case
      
      ## sanity checks
      if(length(N) != length(y))
        stop("should have length(N) = length(y)")
      if(any(N < y)) stop("should have all N >= y")
      
      ## construct yX and re-use kappa
      n <- length(y)
      kappa <- kappa*c(y, N-y)
      knz <- kappa != 0
      yX <- rbind(X, -X)[knz,]
      y <- c(rep(1,n),rep(0,n))[knz]
      kappa <- kappa[knz]
    }

    ## return the data we're going to work with
    return(list(yX=yX, y=y, kappa=kappa))
  }


## mpreprocess:
##
## pre-process the X and y and kappa data depending on
## how binomial responses might be processed

mpreprocess <- function(X, y, flatten, kappa)
  {
    ## process y's and set up design matrices

    ## sanity check
    if(is.null(dim(y))) stop("should have matrix y")

    ## dumb initialization
    yX <- ynew <- NULL
    Q <- ncol(y)

    ## check if there is a need to flatten
    N <- rowSums(y)
    
    if(max(N) == 1 || flatten) { ## flattened binomial case

      for(q in 2:Q) { ## for each label but the last

        if(flatten) {
          ye <- Xe <- NULL
          for(i in 1:nrow(y)) { ## for each obs
            Xe <- rbind(Xe, matrix(rep(X[i,], N[i], nrow=N[i], byrow=TRUE)))
            ye <- c(ye, rep(1, y[i,q]), rep(0, N[i]-y[i,q]))
          }
          y[,q] <- ye; X <- Xe
        } 
        
        ## create matrix y in {-1,+1}
        ypm1 <- y[,q]
        ypm1[y[,q] == 0] <- -1
        
        ## create y *. X
        if(is.null(yX)) {
          yX <- array(NA, dim=c(nrow(X), ncol(X), Q-1))
        } else if(nrow(X) != nrow(yX[,,q-1])) stop("dim problem")

        ## populate yX
        yX[,,q-1] <- X * ypm1
      }
      knew <- matrix(kappa, nrow=1, ncol=Q-1)
      
    } else {  ## unflattened multinomial case

      for(q in 2:Q) { ## for each label but the last
      
        ## construct yX and re-use kappa
        n <- nrow(y)
        kk <- kappa*c(y[,q], N-y[,q])
        knz <- kk != 0

        ## create y *. X
        if(is.null(yX)) {
          ne <- sum(knz)
          yX <- array(NA, dim=c(ne, ncol(X), Q-1))
          ynew <- matrix(NA, nrow=ne, ncol=Q-1)
          knew <- matrix(NA, nrow=ne, ncol=Q-1)
        } else if(sum(knz) != nrow(yX[,,q])) stop("dim problem")

        ## populate yX
        yX[,,q-1] <- rbind(X, -X)[knz,]
        ynew[,q-1] <- c(rep(1,n), rep(0,n))[knz]
        knew[,q-1] <- kk[knz]
      }
      y <- ynew
    }

    ## return the data we're going to work with
    return(list(yX=yX, y=y, kappa=knew))
  }


## reglogit:
##
## function for Gibbs sampling from the a logistic
## regression model, sigma=NULL and nu < 0 indicates
## no prior

reglogit <- function(T, y, X, N=NULL, flatten=FALSE, sigma=1, nu=1,
                  kappa=1, icept=TRUE, normalize=TRUE,
                  zzero=TRUE, powerprior=TRUE, kmax=442,
                  bstart=NULL, lt=NULL, nup=list(a=2, b=0.1),
                  save.latents=FALSE, verb=100)
{
  ## checking T
  if(length(T) != 1 || !is.numeric(T) || T < 1)
    stop("T should be a positive integer scalar")
  
  ## check ys
  if(any(y < 0)) stop("ys must be non-negative integers")
  
  ## getting and checking data size
  m <- ncol(X)
  n <- length(y)
  if(n != nrow(X)) stop("dimension mismatch")

  ## possibly deal with sparse X
  if(class(X)[1] == "dgCMatrix") {
    sparse <- TRUE
    X <- as.matrix(X)
  } else sparse <- FALSE

  ## design matrix processing
  X <- as.matrix(X)
  if(normalize) {
    one <- rep(1, n)
    normx <- sqrt(drop(one %*% (X^2)))
    if(any(normx == 0)) stop("degenerate X cols")
    X <- as.matrix(as.data.frame(scale(X, FALSE, normx)))
  } else normx <- rep(1, ncol(X))

  ## init sigma
  if(length(sigma) == 1) sigma <- rep(sigma, ncol(X))

  ## add on intecept term?
  if(icept) {
    X <- cbind(1, X)
    normx <- c(1, normx)
  }

  ## put the starting beta value on the right scale
  if(!is.null(bstart) && normalize) bstart <- bstart*normx
  else if(is.null(bstart)) bstart <- rnorm(m+icept)

  ## allocate beta, and set starting position
  beta <- matrix(NA, nrow=T, ncol=m+icept)
  beta[1,] <- bstart
  map <- list(beta=bstart)
  
  ## check if we are inferring nu
  if(!is.null(nup)) {
    nus <- rep(NA, T)
    if(nu <= 0) stop("starting nu must be positive")
    nus[1] <- nu
  } else nus <- NULL
  
  ## check for agreement between nu and sigma
  if(is.null(sigma) && nu > 0)
    stop("must define sigma for regularization")
  else if(nu == 0) sigma <- NULL

  ## check save.latents
  if(length(save.latents) != 1 || !is.logical(save.latents))
    stop("save.latents should be a scalar logical")

  ## initial values of the regularization latent variables
  if(nu > 0) {  ## omega lasso/L1 latents
    ot <- 1
    if(save.latents) {
      omega <- matrix(NA, nrow=T, ncol=m)
      omega[1,] <- ot
    }
  } else { omega <- NULL; ot <- rep(1, m) }

  ## save the original kappa
  if(powerprior) kp <- kappa
  else kp <- 1

  ## process y's and set up design matrices
  nd <- preprocess(X, y, N, flatten, kappa)
  yX <- nd$yX; y <- nd$y; kappa <- nd$kappa
  n <- nrow(yX)
  if(sparse) yX.sparse <- Matrix(yX, sparse=TRUE)
 
  ## initialize lambda latent variables
  if(is.null(lt)) lt <- rep(1, n)
  map$lambda <- lt
  if(save.latents) {
    lambda <- matrix(NA, nrow=T, ncol=n)
    lambda[1,] <- lt
  }

  ## initial values for the logit latents
  if(!zzero) {
    zt <- y
    if(save.latents) {
      z <- matrix(NA, nrow=T, ncol=n)
      z[1,] <- zt
    }
  } else { z <- zt <- NULL }

  ## allocate and initial log posterior calculation
  lpost <- rep(NA, T)
  map$lpost <- lpost[1] <-
    calc.lpost(yX, map$beta, nu, kappa, kp, sigma, nup)
  if(!is.null(nus)) map$nu <- nu
  
  ## initialize internal RNG
  .C("newRNGstates")

  ## Gibbs sampling rounds
  for(t in 2:T) {

    ## progress meter
    if(t %% verb == 0) cat("round", t, "\n")
    
    ## if regularizing, draw the latent omega values
    if(nu > 0) {
      ot <- draw.omega(beta[t-1,], sigma, nu, kp)
      if(save.latents) omega[t,] <- ot
    }

    ## if logistic, draw the latent lambda values
    lt <- draw.lambda(yX, beta[t-1,], lt, kappa, kmax, zzero)
    if(save.latents) lambda[t,] <- lt

    ## draw the latent z values
    if(!zzero) {
      zt <- draw.z(yX, beta[t-1,], lt, kappa)
      if(save.latents) z[t,] <- zt
    }

    ## draw the regression coefficients
    if(sparse) beta[t,] <- draw.beta(yX.sparse, zt, lt, ot, nu, sigma, kappa, kp)
    else beta[t,] <- draw.beta(yX, zt, lt, ot, nu, sigma, kappa, kp)

    ## maybe draw samples from nu
    if(!is.null(nus)) nu <- nus[t] <- draw.nu(beta[t,], sigma, kp, nup)
    
    ## calculate the posterior probability to find the map
    lpost[t] <- calc.lpost(yX, beta[t,], nu, kappa, kp, sigma, nup)

    ## update the map
    if(lpost[t] > map$lpost) {
      map <- list(beta=beta[t,], lpost=lpost[t], lambda=lt)
      if(!is.null(nus)) map$nu <- nus[t]
    }
      
  }

  ## destroy internal RNG
  .C("deleteRNGstates")

  ## un-normalize
  if(normalize) {
    beta <- as.matrix(as.data.frame(scale(beta, FALSE, scale=normx)))
    map$beta <- map$beta/normx
  }     

  ## construct the return object
  r <- list(X=X, y=y, beta=beta, lpost=lpost,
            map=map, kappa=kappa)
  if(save.latents) r$lambda <- lambda
  if(nu > 0 && save.latents) r$omega <- omega
  if(normalize) r$normx <- normx
  if(!zzero && save.latents) r$z <- z
  if(!is.null(nus)) r$nu <- nus
  r$icept <- icept

  ## assign a class to the object
  class(r) <- "reglogit"
  
  return(r)
}


## regmlogit:
##
## function for Gibbs sampling from the a logistic
## (multinomial) regression model, sigma=NULL and
## nu < 0 indicates no prior

regmlogit <- function(T, y, X, flatten=FALSE, sigma=1, nu=1,
                      kappa=1, icept=TRUE, normalize=TRUE,
                      zzero=TRUE, powerprior=TRUE, kmax=442,
                      bstart=NULL, lt=NULL, nup=list(a=2, b=0.1),
                      save.latents=FALSE, verb=100)
{
  ## checking T
  if(length(T) != 1 || !is.numeric(T) || T < 1)
    stop("T should be a positive integer scalar")
  
  ## getting and checking data size 
  m <- ncol(X)
  if(is.null(dim(y))) { ## convert y into a matrix
    if(any(y <= 0)) stop("vector ys must be positive integer class labels")
    ymat <- matrix(0, nrow=length(y), ncol=max(y))
    for(i in 1:nrow(ymat)) ymat[i,y[i]] <- 1  ## COULD DO THIS FASTER
    y <- ymat
  }

  ## check matrix form of y
  if(any(apply(y, 1, function(x) { any(x < 0) || sum(x) <= 0 })))
    stop("matrix ys must have at least one non-zero in each row")
  Q <- ncol(y)
  n <- nrow(y)
  if(n != nrow(X)) stop("dimension mismatch")
  
  ## possibly deal with sparse X
  if(class(X)[1] == "dgCMatrix") {
    sparse <- TRUE
    X <- as.matrix(X)
  } else sparse <- FALSE

  ## design matrix processing
  X <- as.matrix(X)
  if(normalize) {
    one <- rep(1, n)
    normx <- sqrt(drop(one %*% (X^2)))
    if(any(normx == 0)) stop("degenerate X cols")
    X <- as.matrix(as.data.frame(scale(X, FALSE, normx)))
  } else normx <- rep(1, ncol(X))
  
  ## init sigma
  if(length(sigma) == 1) sigma <- rep(sigma, ncol(X))
  
  ## add on intecept term?
  if(icept) {
    X <- cbind(1, X)
    normx <- c(1, normx)
  }
  
  ## put the starting beta value on the right scale
  if(!is.null(bstart)) {
    if(is.null(dim(bstart))) {
      bstart <- matrix(bstart, ncol=Q-1)
      if(nrow(bstart) != m+icept) stop("bstart should be (m+icept)x(Q-1) matrix")
    }
    if(normalize) bstart <- apply(bstart, 2, function(x, n) x*n, n=normx)
  } else if(is.null(bstart)) bstart <- matrix(rnorm((Q-1)*(m+icept)), ncol=Q-1)
  
  ## allocate beta, and set starting position
  beta <- array(NA, dim=c(T, m+icept, Q-1))
  beta[1,,] <- bstart
  map <- list(beta=beta[1,,])

  ## check nu
  if(length(nu) == 1) nu <- rep(nu, Q-1)
  if(length(nu) != Q-1)
    stop("nu must be scalar or (Q-1)-vector")
  
  ## check if we are inferring nu
  if(!is.null(nup)) {
    nus <- matrix(NA, nrow=T, ncol=Q-1)
    ## must infer nu for all q in 1:(Q-1), could change later
    if(any(nu <= 0)) stop("all starting nu must be positive")
    nus[1,] <- nu
  } else nus <- NULL
  
  ## check for agreement between nu and sigma
  if(is.null(sigma) && any(nu > 0))
    stop("must define sigma for regularization")
  else if(all(nu == 0)) sigma <- NULL

  ## check save.latents
  if(length(save.latents) != 1 || !is.logical(save.latents))
    stop("save.latents should be a scalar logical")

  ## initial values of the regularization latent variables
  if(all(nu > 0)) {  ## omega lasso/L1 latents
    ot <- matrix(1, nrow=m, ncol=Q-1)
    if(save.latents) {
      omega <- array(NA, dim=c(T, m, Q-1))
      omega[1,,] <- ot; 
    }
  } else { omega <- NULL; ot <- matrix(1, nrow=m, ncol=Q-1) }
  
  ## save the original kappa
  if(powerprior) kp <- kappa
  else kp <- 1
  
  ## process y's and set up design matrices
  nd <- mpreprocess(X, y, flatten, kappa)
  yX <- nd$yX; y <- nd$y; kappa <- nd$kappa
  n <- dim(yX)[1]
  if(sparse) {
    yX.sparse <- array(NA, dim=dim(yX))
    for(q in 1:(Q-1)) yX.sparse[,,q] <- Matrix(yX[,,q], sparse=TRUE)
  }
  
  ## initialize lambda latent variables
  if(is.null(lt)) lt <- matrix(1, nrow=n, ncol=Q-1)
  map$lambda <- lt
  if(save.latents) {
    lambda <- array(NA, dim=c(T, n, Q-1))
    lambda[1,,] <- lt
  }
  
  ## initial values for the logit latents
  if(!zzero) {
    zt <- y
    if(save.latents) {
      z <- array(NA, dim=c(T, n, Q-1))
      z[1,,] <- zt
    }
  } else { z <- zt <- NULL }
  
  ## allocate and initial log posterior calculation
  lpost <- matrix(NA, nrow=T, ncol=Q-1)
  map$lpost <- lpost[1,] <- 
    calc.mlpost(yX, map$beta, nu, kappa, kp, sigma, nup)
  if(!is.null(nus)) map$nu <- nu

  ## Gibbs sampling rounds
  for(t in 2:T) {
    
    ## progress meter
    if(t %% verb == 0) cat("round", t, "\n")

    ## loop over class labels
    for(q in 1:(Q-1)) {

      ## calculate Cs
      Cs <- calc.Cs(yX, beta[t-1,,], lt, kappa)
      
      ## if regularizing, draw the latent omega values
      if(all(nu > 0)) {
        ot[,q] <- draw.omega(beta[t-1,,q], sigma, nu, kp)
        if(save.latents) omega[t,,q] <- ot[,q]
      }
      
      ## if logistic, draw the latent lambda values
      lt[,q] <- draw.lambda(yX[,,q], beta[t-1,,q], lt[,q], kappa[,q], kmax, 
        zzero, C=Cs[,q])
      if(save.latents) lambda[t,,q] <- lt[,q]

      ## draw the latent z values
      if(!zzero) {
        zt[,q] <- draw.z(yX[,,q], beta[t-1,,q], lt[,q], kappa, C=Cs[,q])
        if(save.latents) z[t,,q] <- zt[,q]
      }

      ## draw the regression coefficients
      if(sparse) beta[t,,q] <- draw.beta(yX.sparse[,,q], zt[,q], lt[,q], ot[,q], nu[q],
                              sigma, kappa[,q], kp)
      else beta[t,,q] <- draw.beta(yX[,,q], zt[,q], lt[,q], ot[,q], nu[q],
                              sigma, kappa[,q], kp)
      
      ## maybe draw samples from nu
      if(!is.null(nus)) nu[q] <- nus[t,q] <- draw.nu(beta[t,,q], sigma, kp, nup)
    }

    ## calculate the posterior probability to find the map
    lpost[t] <- calc.mlpost(yX, beta[t,,], nu, kappa, kp, sigma, nup)
    
    ## update the map
    if(lpost[t] > map$lpost) {
      map <- list(beta=beta[t,,], lpost=lpost[t], lambda=lt)
      if(!is.null(nus)) map$nu <- nus[t,]
    }
      
  }
  
  ## un-normalize
  if(normalize) {
    for(q in 1:(Q-1)) {
      beta[,,q] <- as.matrix(as.data.frame(scale(beta[,,q], FALSE, scale=normx)))
      map$beta[,q] <- map$beta[,q]/normx
    }
  }     
  
  ## construct the return object
  r <- list(X=X, y=y, beta=beta, lpost=lpost,
            map=map, kappa=kappa)
  if(save.latents) r$lambda <- lambda
  if(all(nu > 0) & save.latents) r$omega <- omega
  if(normalize) r$normx <- normx
  if(!zzero && save.latents) r$z <- z
  if(!is.null(nus)) r$nu <- nus
  r$icept <- icept
  
  ## assign a class to the object
  class(r) <- "regmlogit"
  
  return(r)
}

## predict.reglogit:
##
## prediction method for binary regularized logistic
## regression

predict.reglogit <- function(object, XX,
                             burnin=round(0.1*nrow(object$beta)),
                             ...)
  {
    ## check XX
    if(object$icept) XX <- cbind(1, XX)
    if(ncol(XX) != ncol(object$X)) stop("XX dims don't match")

    ## check burnin
    if(length(burnin) != 1 || burnin < 0)
      stop("burnin must be a positive scalar")
    if(burnin+1 >= nrow(object$beta))
      stop("burnin too big; >= T-1 samples")

    ## extract beta
    beta <- object$beta[-(1:burnin),]

    ## calculate 1/(1+exp(-X %*% beta))
    probs <- 1-plogis(0, XX %*% t(beta))

    ## normalize and average over MCMC samples
    mprobs <- rowMeans(probs)

    ## get point-predictions and variance estimates
    pc <- round(mprobs)
    ent <- - mprobs * log(mprobs) - (1-mprobs)*log(1-mprobs)

    return(list(p=probs, mp=mprobs, c=pc, ent=ent))
  }



## predict.regmlogit:
##
## prediction method for multinomial regularized logistic
## regression

predict.regmlogit <- function(object, XX,
                              burnin=round(0.1*dim(object$beta)[1]),
                              ...)
  {
    ## check XX
    if(object$icept) XX <- cbind(1, XX)
    if(ncol(XX) != ncol(object$X)) stop("XX dims don't match")

    ## check burnin
    if(length(burnin) != 1 || burnin < 0)
      stop("burnin must be a positive scalar")
    if(burnin+1 >= dim(object$beta)[1])
      stop("burnin too big; >= T-1 samples")

    ## extract beta
    beta <- object$beta[-(1:burnin),,,drop=FALSE]

    ## calculate exp(X %*% beta)
    Q <- dim(beta)[3]+1
    probs <- array(1, dim=c(nrow(XX), dim(beta)[1], Q))
    psum <- matrix(1, nrow=nrow(XX), ncol=dim(beta)[1])
    for(q in 1:(Q-1)) {
      probs[,,q+1] <- exp(XX %*% t(beta[,,q]))
      psum <- psum + probs[,,q+1]
    }

    ## normalize and average over MCMC samples
    mprobs <- matrix(NA, nrow=nrow(XX), ncol=Q)
    for(q in 1:Q) {
      probs[,,q] <- probs[,,q]/psum
      mprobs[,q] <- rowMeans(probs[,,q])
    }

    ## get point-predictions and variance estimates
    pc <- apply(mprobs, 1, which.max)
    ent <- apply(mprobs, 1, function(p) -sum(p * log(p)) )

    return(list(p=probs, mp=mprobs, c=pc, ent=ent))
  }