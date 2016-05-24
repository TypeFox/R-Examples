################################################################################
#
# Bayesian estimation of spatial autoregressive probit model (SAR probit)
# using MCMC sampling
#
# Stefan Wilhelm <Stefan.Wilhelm@financial.com>
# using code from Miguel Godinho de Matos <miguelgodinhomatos@cmu.edu>
#
################################################################################

if (FALSE) {
 library(tmvtnorm)
 library(mvtnorm)
 library(Matrix)    # sparseMatrix
 source("sar_base.r")
 source("matrix_operations.r")
 source("stats_distributions.r")
 source("utility_functions.r") 
}

# estimated tr(W^i) for i=1,...,100
# see LeSage (2009), chapter 4, p.97/98
#
# "The expectation of the quadratic form u'Au equals tr(A).
#  Since u_i^2 follows a chi^2 distribution with one degree of freedom."
# Pre-calculate traces tr(W^i) i=1,...,100 for the x-impacts calculations
#
# @param W spatial weight matrix (n x n)
# @param o highest order of W^i = W^o
# @param iiter number of MCMC iterations (we run this 50 times)
# @return (n x o) matrix with tr(W^i) in each column, for i=1..,o
tracesWi <- function(W, o=100, iiter=50) {
  n <- nrow(W)
  trW_i <- matrix( data=0, nrow=n, ncol=o )   # n x o
  u  <- matrix(rnorm(n * iiter), nrow = n, ncol = iiter)   # (n x iiter)
  xx <- u
  trW_i[,1] <- apply(u * as.matrix(xx), 1, sum)    # tr(W^0) = 1
  for(i in 2:o ){
    xx <- W %*% xx  # (n x iter)
    trW_i[,i] <- apply(u * as.matrix(xx), 1, sum)  # u'(W^i)u; sum across all iterations
  }
  trW_i <- trW_i / iiter
  return(trW_i)
}

# PURPOSE: draw rho from conditional distribution p(rho | beta, z, y)
# ---------------------------------------------------
#  USAGE: 
#
#  results.rmin 
#  results.rmax maximum eigen value
#  results.time execution time
# ---------------------------------------------------
#where epe0 <- t(e0) %*% e0
#      eped <- t(ed) %*% ed
#      epe0d<- t(ed) %*% e0
#
# detval1 umbenannt in => rho_grid = grid/vector of rhos
# detval2 umbenannt in => lndet = vector of log-determinants
# detval1sq = rvec^2
# lnbprior = vector of log prior densities for rho (default: Beta(1,1))
#
draw_rho <- function (rho_grid, lndet, rho_gridsq, yy, epe0, eped, epe0d, 
    rho, nmk, nrho, lnbprior, u) 
{
    # This is the scalar concentrated log-likelihood function value
    # lnL(rho). See eqn(3.7), p.48
    z <- epe0 - 2 * rho_grid * epe0d + rho_gridsq * eped
    z <- -nmk * log(z)
    den <- lndet + z + lnbprior          # vector of log posterior densities for rho vector
    # posterior density post(rho | data) \propto likelihood(data|rho,beta,z) * prior(rho)
    # log posterior density post(rho | data) \propto loglik(data|rho,beta,z) + log prior(rho)
    n <- nrho
    adj <- max(den)
    den <- den - adj                     # adjustieren/normieren der log density auf maximum 0; 
    x <- exp(den)                        # density von p(rho) --> pdf
    isum <- sum(yy * (x[2:n] - x[1:(n - 1)])/2)
    z <- abs(x/isum)
    den <- cumsum(z)
    rnd <- u * sum(z)
    ind <- which(den <= rnd)
    idraw <- max(ind)
    if (idraw > 0 && idraw < nrho) {
        results <- rho_grid[idraw]
    }
    else {
        results <- rho
    }
    return(results)
}

# Bayesian estimation of SAR probit model
#
# @param formula 
sarprobit <- function(formula, W, data, subset, ...) {
  cl <- match.call()                     # cl ist object of class "call"
  mf <- match.call(expand.dots = FALSE)  # mf ist object of class "call"
  m  <- match(c("formula", "data", "subset"), names(mf), 0L)        # m is index vector
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())         # from here mf is a data.frame
  mt <- attr(mf, "terms")                # mt is object of class "terms" and "formula"
  y <- model.response(mf, "numeric")
  if (!is.null(W) && !is.numeric(W) && !inherits(W, "sparseMatrix") && nrow(W) != NROW(y)) 
    stop(gettextf("'W' must be a numeric square matrix, dimension %d should equal %d (number of observations)",
      NROW(W), NROW(y)), domain = NA)
  
  X <- model.matrix(mt, mf, contrasts)
  sar_probit_mcmc(y, X, W, ...)    
}

# faster update of matrix S = (I - rho * W) for new values of rho
#
# @param S template matrix of (I - rho * W)
# @param ind indizes to replaced
# @param W spatial weights matrix W
# @return (I - rho * W)
update_I_rW <- function(S, ind, rho, W) {
  S@x[ind] <- (-rho*W)@x
  return(S)
}

# Estimate the spatial autoregressive probit model (SAR probit)
# z = rho * W * z + X \beta + epsilon
# where y = 1 if z >= 0 and y = 0 if z < 0 observable
#
# @param y
# @param X
# @param W spatial weight matrix
# @param ndraw number of MCMC iterations
# @param burn.in  number of MCMC burn-in to be discarded
# @param thinning MCMC thinning factor, defaults to 1
# @param m number of burn.in sampling in inner Gibbs sampling
# @param prior list of prior settings: 
#   prior$rho ~ Beta(a1,a2); 
#   prior$beta ~ N(c, T)
#   prior$lflag   
# lflag=0 --> default to 1997 Pace and Barry grid approach
# lflag=1 --> Pace and LeSage (2004) Chebyshev approximation
# lflag=2 --> Barry and Pace (1999) MC determinant approx
# @param start
# @param m
# @param computeMarginalEffects
# @param showProgress
sar_probit_mcmc <- function(y, X, W, ndraw=1000, burn.in=100, thinning=1, 
  prior=list(a1=1, a2=1, c=rep(0, ncol(X)), T=diag(ncol(X))*1e12, lflag = 0), 
  start=list(rho=0.75, beta=rep(0, ncol(X))),
  m=10, computeMarginalEffects=TRUE, showProgress=FALSE){  

  #start timer
  timet <- Sys.time()
  
  n  <- nrow( X )            # number of observations
  n1 <- nrow( X )          
  n2 <- nrow( W )
  k <- ncol( X )             # number of of parameters/exogenous variables
  I_n <- sparseMatrix(i=1:n, j=1:n, x=1) # sparse identity matrix
  if (is.null(colnames(X))) colnames(X) <- paste("x",1:k,sep="")
  
  #validate inputs
  if( length(c(which(y == 0 ),which(y == 1))) != length( y ) ){
    stop('sarprobit: not all y-values are 0 or 1')
  }
  if( n1 != n2 && n1 != n ){
    stop('sarprobit: wrong size of spatial weight matrix W')
  }
  # check if spatial weights matrix W does not contain zeros in the main diagonal
  if (!inherits(W, "sparseMatrix") || any(diag(W) != 0)) {
    stop('sarprobit: spatial weights matrix W must be a sparse matrix with zeros in the main diagonal')
  }
  # check if we have a constant term in X
  ind <- match( n, apply(X,2,sum))
  if( is.na(ind) ){
    cflag <- 0
    p     <- k
  }else if( ind == 1 ){
    cflag <- 1
    p     <- k - 1
  }else{
    stop('sarprobit: intercept term must be in first column of the X-matrix')
  }
  
  # MCMC start values
  rho  <- start$rho          # start value of rho
  beta <- start$beta         # start value of parameters, prior value, we could also sample from beta ~ N(c, T)
  
  # MCMC priors
  # conjugate prior beta ~ N(c, T)
  # parametrize, default to diffuse prior, for beta, e.g. T <- diag(k) * 1e12
  c <- rep(0, k)             # prior distribution of beta ~ N(c, T) : c = 0
  if (is.numeric(prior$c) && length(prior$c) == k) {
    c <- prior$c
  }
  if (is.matrix(prior$T) && ncol(prior$T) == k && isSymmetric(prior$T) && det(prior$T) > 0) {
    T <- prior$T               # prior distribution of beta ~ N(c, T) : T = I_n --> diffuse prior
  } else {
    T <- diag(k)*1e12
  }
  
  Tinv <- solve(T)           # T^{-1}
  
  # prepare computation of (I_n - rho * W)
  if (class(W) == "dgCMatrix") {
   I <- sparseMatrix(i=1:n,j=1:n,x=Inf)
   S <- (I - rho * W)
   ind  <- which(is.infinite(S@x))  # Stellen an denen wir 1 einsetzen müssen (I_n)
   ind2 <- which(!is.infinite(S@x))  # Stellen an denen wir -rho*W einsetzen müssen
   S@x[ind] <- 1
  } else {
   S <- I_n - rho * W
  }
  
  H <- t(S) %*% S            # precision matrix H for beta | rho, z, y
  QR <- qr(S)                # class "sparseQR"
  mu <- solve(QR, X %*% beta)
  
  # truncation points for z, depend only on y, can be precalculated
  lower <- ifelse(y > 0, 0,  -Inf)
  upper <- ifelse(y > 0, Inf,   0)

# prepare settings for drawing rho  
  rmin       <- -1   # use -1,1 rho interval as default
  rmax       <-  1
  
  lflag <- 0
  if (is.numeric(prior$lflag) && lflag %in% c(0, 1, 2)) lflag <- prior$lflag
  #lflag=0 --> default to 1997 Pace and Barry grid approach
  #lflag=1 --> Pace and LeSage (2004) Chebyshev approximation
  #lflag=2 --> Barry and Pace (1999) MC determinant approx
  tmp <- sar_lndet(lflag, W, rmin, rmax)
  detval <- tmp$detval
  
  # Some precalculated quantities for drawing rho
  # rho ~ Beta(a1, a2) prior
  a1         <-  1
  a2         <-  1
  if (is.numeric(prior$a1)) a1 <- prior$a1
  if (is.numeric(prior$a2)) a2 <- prior$a2
  
  lnbprior <- log(beta_prior(detval[,1],a1,a2))
  u        <- runif(thinning * ndraw + burn.in)   # u ~ U(0, 1)
  nrho     <- nrow(detval)  # do_ldet() liefert nur 2000 statt 2001 Gridpoints
  nmk      <- (n-k)/2
  rho_grid   <- detval[,1]  # rho grid values
  lndet    <- detval[,2]  # log-determinant grid values
  rho_gridsq <- rho_grid * rho_grid
  yy       <- (rho_grid[2:nrho] + rho_grid[1:(nrho-1)])
    
  # matrix to store the beta + rho parameters for each iteration/draw
  B <- matrix(NA, ndraw, k+1)
  
  # progress bar
  if (showProgress) {
    pb <- txtProgressBar(min=0, max=(thinning * ndraw + burn.in), initial=0, style=3)
  }
  
  # immutable matrices
  tX <- t(X)                       # X'               # k x n
  xpx  <- t(X) %*% X               # (X'X)            # k x k
  xpxI <- solve(xpx)               # (X'X)^{-1}       # k x k
  xxpxI <- X %*% xpxI              # X(X'X)^(-1)     # n x k (better, compromise)
  AA    <- solve(xpx + Tinv)       # (X'X + T^{-1})^{-1}
  
  # draw from multivariate normal beta ~ N(c, T). we can precalculate 
  # betadraws ~ N(0, T) befor running the chain and later just create beta as
  # beta = c + betadraws ~ N(c, T)
  betadraws <- rmvnorm(n=(burn.in + ndraw * thinning), mean=rep(0, k), sigma=AA)
  
  # matrices for direct and indirect impacts
  direct       <- matrix(NA, ndraw,p)    # n x p
  indirect     <- matrix(NA, ndraw,p)    # n x p
  total        <- matrix(NA, ndraw,p)    # n x p
  zmean        <- rep(0, n)
    
  # names of non-constant parameters
  if(cflag == 0) {
    namesNonConstantParams <- colnames(X)
  } else {
    namesNonConstantParams <- colnames(X)[-1]
  }
  colnames(total) <- namesNonConstantParams
  colnames(direct)   <- namesNonConstantParams
  colnames(indirect) <- namesNonConstantParams
    
  if (computeMarginalEffects) {
    # simulate Monte Carlo estimation of tr(W^i) for i = 1..o before MCMC iterations
    trW.i <- tracesWi(W, o=100, iiter=50)
  }
    
  # just to set a start value for z
  z <- rep(0, n)
  ones <- rep(1, n)
  
  for (i in (1 - burn.in):(ndraw * thinning)) {
  
  # 1. sample from z | rho, beta, y using precision matrix H
  # mu will be updated after drawing from rho
  
  # see LeSage (2009) for choice of burn-in size, often m=5 or m=10 is used!
  # we can also use m=1 together with start.value=z, see LeSage (2009), section 10.1.5
  if (m==1) {
    z <- as.double(rtmvnorm.sparseMatrix(n=1, mean=mu, H=H, 
      lower=lower, upper=upper, burn.in=m, start.value=z))
  } else {
    z <- as.double(rtmvnorm.sparseMatrix(n=1, mean=mu, H=H, 
      lower=lower, upper=upper, burn.in=m))
  }
    
  # 2. sample from beta | rho, z, y
  Sz <- as.double(S %*% z)               # (n x 1); dense
  c2 <- AA  %*% (tX %*% Sz + Tinv %*% c) # (n x 1); dense
  T <- AA   # no update basically on T, TODO: check this
  beta <- as.double(c2 + betadraws[i + burn.in, ])
  
  # 3. sample from rho | beta, z
  #---- DRAW RHO ----
  #see LeSage 2009 chapter 5 - page 132 for the explanation of the
  #code below which is used for numerical integration of the rho prior.
  #I changed from the original code to match the notation of the book
  #using c0 and cd below instead of b0 and bd ....
  xpz  <- tX %*% z           # X'z
  Wz   <- as.double(W %*% z) # Wz        # SW: coerce Wz to vector 
  # (from n x 1 sparse matrix! we do not need a sparse matrix here)
  xpWz <- tX %*% Wz          # X'Wz      # k x 1
  e0   <-  z - xxpxI %*% xpz  # z  - X(X'X)^-1X' z
  ed   <- Wz - xxpxI %*% xpWz # Wz - X(X'X)^(-1)X'Wz
  epe0 <- as.double(crossprod(e0))  # slightly faster than t(e0) %*% e0
  eped <- as.double(crossprod(ed))
  epe0d<- as.double(crossprod(ed, e0))
  rho  <- draw_rho(rho_grid, lndet, rho_gridsq, yy, epe0, eped, epe0d, rho, nmk=nmk, nrho=nrho, lnbprior, u=u[i + burn.in])
  
  ############################################################################## 
  
  # update S, H and QR decomposition of S and mu after each iteration; before effects
  S <- update_I_rW(S, ind=ind2, rho, W)  # update (I - rho * W)
  H <- t(S) %*% S      # H = S'S 
  QR <- qr(S)          # class "sparseQR"
  
  # solving equation 
  # (I_n - rho * W) mu = X beta 
  # instead of inverting S = I_n - rho * W as in mu = ( In -  rho W)^{-1} X beta.
  # QR-decomposition for sparse matrices
  mu <- solve(QR, X %*% beta)

  if (i > 0) {
    if (thinning == 1) {
      ind <- i
    }
    else if (i%%thinning == 0) {
      ind <- i%/%thinning
    } else {
      next
    }
    
    B[ind,] <- c(beta, rho)
    zmean   <- zmean + z
  
    # compute effects estimates (direct and indirect impacts) in each MCMC iteration
    if (computeMarginalEffects) {
      o <- 100
      rhovec <- rho^(0:(o-1)) # SW: (100 x 1)   mit [1, rho^1, rho^2 ..., rho^99], see LeSage(2009), eqn (4.145), p.115
      if( cflag == 1 ){ #has intercept
        beff <- beta[-1]      # beff is parameter vector without constant
      }else if(cflag == 0){
        beff <- beta          # no constant in model
      }
      # beff is parameter vector without constant!
      # See LeSage (2009), section 5.6.2., p.149/150 for spatial effects estimation in MCMC
      #   direct: M_r(D) = n^{-1} tr(S_r(W))           # SW: efficient approaches available, see chapter 4, pp.114/115
      #    total: M_r(T) = n^{-1} 1'_n S_r(W) 1_n      # SW: Problem: S_r(W) is dense, but can be solved via QR decomposition of S
      # indirect: M_r(I) = M_r(T) - M_r(D)
      # SW: See LeSage (2009), section 10.1.6, p.293 for Marginal effects in SAR probit
      pdfz <- dnorm(as.numeric(mu))                     # standard normal pdf phi(mu)
      dd   <- sparseMatrix(i=1:n, j=1:n, x=pdfz)       # dd is diagonal matrix with pdfz as diagonal (n x n)
  
      dir      <- as.double(t(pdfz) %*% trW.i %*% rhovec /n)  # (1 x n) * (n x o) * (o x 1)
      # direct impact : dy_i / d X_ir = phi((In -  rho W)^{-1} X beta_r) * beta_r
      avg_direct     <- dir * beff      # (p x 1)
         
      # We compute the average total effects without inverting S 
      # unlike in the LeSage Matlab Code, 
      # but using the QR decomposition of S which we already have!
      # average total effects = n^(-1) * 1_n' %*% (D %*% S^(-1) * b[r]) %*% 1_n
      #                       = n^(-1) * 1_n' %*% (D %*% x) * b[r]
      # where D=dd is the diagonal matrix containing phi(mu)
      # and x is the solution of S %*% x = 1_n, obtained from the QR-decompositon
      # of S. The average total effects is then the mean of (D %*% x) * b[r]
      # average total effects, which can be furthermore done for all b[r] in one operation.
      avg_total    <- mean(dd %*% qr.coef(QR, ones)) * beff
      avg_indirect <- avg_total - avg_direct    # (p x 1)
      
      total[ind, ]      <- avg_total    # an (ndraw-nomit x p) matrix
      direct[ind, ]     <- avg_direct   # an (ndraw-nomit x p) matrix
      indirect[ind, ]   <- avg_indirect # an (ndraw-nomit x p) matrix
    }
  ##############################################################################
    
  }
    
  if (showProgress) setTxtProgressBar(pb, i + burn.in) # update progress bar
  }
    
  if (showProgress)  close(pb) #close progress bar
  
  # fitted values for estimates (based on z rather than binary y like in fitted(glm.fit))
  # (on response scale y vs. linear predictor scale z...)
  beta  <- colMeans(B)[1:k]
  rho   <- colMeans(B)[k+1]
  S     <- (I_n - rho * W)
  fitted.values   <- solve(qr(S), X %*% beta)   # z = (I_n - rho * W)^{-1}(X * beta)
  fitted.response <- as.numeric(fitted.values >= 0) 
  # TODO: linear.predictors  vs. fitted.values
  
  # result
  results       <- NULL
  results$time  <- Sys.time() - timet
  results$nobs  <- n          # number of observations
  results$nvar  <- k          # number of explanatory variables
  results$y     <- y 
  results$zip   <- n - sum(y) # number of zero values in the y-vector
  results$beta  <- colMeans(B)[1:k]
  results$rho   <- colMeans(B)[k+1]
  results$coefficients <- colMeans(B)
  results$fitted.values <- fitted.values
  results$fitted.response <- fitted.response  # fitted values on response scale (binary y variable)
  results$ndraw <- ndraw
  results$nomit <- burn.in
  results$a1        <- a1
  results$a2        <- a2
  results$rmax      <- rmax 
  results$rmin      <- rmin
  results$tflag     <- 'plevel'
  results$lflag     <- lflag
  results$cflag     <- cflag
  results$lndet     <- detval
  results$names     <- c(colnames(X), 'rho')
  results$B         <- B        # (beta, rho) draws
  results$bdraw     <- B[,1:k]  # beta draws
  results$pdraw     <- B[,k+1]  # rho draws
  results$total     <- total
  results$direct    <- direct
  results$indirect  <- indirect
  results$W <- W
  results$X <- X
  #results$mlike     <- mlike    # log-likelihood based on posterior means

  #results$predicted <- # prediction required. The default is on the scale of the linear predictors
  class(results)    <- "sarprobit"
  return(results)
}

marginal.effects <- function (object, ...) 
 UseMethod("marginal.effects")

# compute marginal effects for every MCMC iteration of the 
# estimated SAR probit model
marginal.effects.sarprobit <- function(object, o=100, ...) {
  # check for class "sarprobit"
  if (!inherits(object, "sarprobit")) 
        stop("use only with \"sarprobit\" objects")
        
  nobs      <- object$nobs
  nvar      <- object$nvar   # number of explanatory variables
  p         <- ifelse(object$cflag == 0, nvar, nvar - 1) # number of non-constant variables
  ndraw     <- object$ndraw
  nomit     <- object$nomit
  betadraws <- object$bdraw
  rhodraws  <- object$pdraw
  X         <- object$X # data matrix
  W         <- object$W # spatial weight matrix
  I_n       <- sparseMatrix(i=1:nobs, j=1:nobs, x=1) # sparse identity matrix
  
  # Monte Carlo estimation of tr(W^i) for i = 1..o, tr(W^i) is (n x o) matrix
  trW.i <- tracesWi(W, o=o, iiter=50)
  
  # Matrices (n x k) for average direct, indirect and total effects for all k (non-constant) explanatory variables
  # TODO: check for constant variables
  D <- matrix(NA, ndraw, p)
  I <- matrix(NA, ndraw, p)
  T <- matrix(NA, ndraw, p)
  
  # names of non-constant parameters
  if(object$cflag == 0) {
    namesNonConstantParams <- colnames(X)
  } else {
    namesNonConstantParams <- colnames(X)[-1]
  }
  colnames(T) <- namesNonConstantParams
  colnames(D) <- namesNonConstantParams
  colnames(I) <- namesNonConstantParams
  
  ones <- rep(1, nobs)
  
  # loop all MCMC draws
  for(i in 1:ndraw) {
  
  # get parameters for this MCMC iteration
  # TODO: check for constant variables
  beta <- betadraws[i,]
  beff <- betadraws[i,-1]
  rho  <- rhodraws[i]
  
  # solving equation 
  # (I_n - rho * W) mu = X beta 
  # instead of inverting S = I_n - rho * W as in mu = (In -  rho W)^{-1} X beta.
  # QR-decomposition for sparse matrices
  S <- I_n - rho * W
  QR <- qr(S)  # class "sparseQR"
  mu <- qr.coef(QR, X %*% beta)      # n x 1
  
  # Marginal Effects for SAR Probit Model:
  # LeSage (2009), equation (10.10), p.294:
  # d E[y | x_r] / dx_r' = phi(S^{-1} I_n mean(x_r) beta_r) * S^{-1} I_n beta_r
  # This gives a (n x n) matrix
  # average direct effects = n^{-1} tr(S_r(W))
  
  # compute effects estimates (direct and indirect impacts) in each MCMC iteration
  rhovec <- rho^(0:(o-1)) # vector (o x 1) with [1, rho^1, rho^2 ..., rho^(o-1)], 
  # see LeSage(2009), eqn (4.145), p.115
    
  # See LeSage (2009), section 5.6.2., p.149/150 for spatial effects estimation in MCMC
  #   direct: M_r(D) = n^{-1} tr(S_r(W))           # efficient approaches available, see LeSage (2009), chapter 4, pp.114/115
  #    total: M_r(T) = n^{-1} 1'_n S_r(W) 1_n      # S_r(W) is dense n x n matrix!
  # indirect: M_r(I) = M_r(T) - M_r(D)             # Computation of dense n x n matrix M_r(T) for total effects required!
  # See LeSage (2009), section 10.1.6, p.293 for Marginal effects in SAR probit
  pdfz <- dnorm(as.numeric(mu))                     # standard normal pdf phi(mu) = phi( (In -  rho W)^{-1} X beta )  # (n x 1)
  dd   <- sparseMatrix(i=1:nobs, j=1:nobs, x=pdfz)  # dd is diagonal matrix with pdfz as diagonal (n x n)
  
  dir      <- as.double(t(pdfz) %*% trW.i %*% rhovec /nobs)  # (1 x n) * (n x o) * (o x 1) = (1 x 1)
  # direct impact : dy_i / d X_ir = phi((In -  rho W)^{-1} X beta_r) * (In -  rho W)^{-1} * beta_r
  avg_direct     <- dir * beff      # (p x 1)
    
  # We compute the average total effects without inverting S 
  # as in the LeSage Matlab Code, 
  # but using the QR decomposition of S which we already have!
  # average total effects = n^(-1) * 1_n' %*% (D %*% S^(-1) * b[r]) %*% 1_n
  #                       = n^(-1) * 1_n' %*% (D %*% x) * b[r]
  # where D=dd is the diagonal matrix containing phi(mu)
  # and x is the solution of S %*% x = 1_n, obtained from the QR-decompositon
  # of S. The average total effects is then the mean of (D %*% x) * b[r]
  # average total effects, which can be furthermore done for all b[r] in one operation.
  avg_total    <- mean(dd %*% qr.coef(QR, ones)) * beff
  avg_indirect       <- avg_total - avg_direct    # (p x 1)
  
  D[i,] <- avg_direct
  I[i,] <- avg_indirect
  T[i,] <- avg_total
  }
  
  summaryMarginalEffects <- function(x) {
    r <- cbind(
    apply(x, 2, mean),
    apply(x, 2, sd),
    apply(x, 2, mean)/apply(x, 2, sd))
    colnames(r) <- c("marginal.effect", "standard.error","z.ratio")
    return(r)
  }
  summary_direct <- summaryMarginalEffects(D)
  summary_indirect <- summaryMarginalEffects(I)
  summary_total <- summaryMarginalEffects(T)    
  
  return(list(direct=D, indirect=I, total=T,
   summary_direct=summary_direct,
   summary_indirect=summary_indirect,
   summary_total=summary_total)
  ) 
}

# summary method for class "sarprobit"
summary.sarprobit <- function(object, var_names=NULL, file=NULL, 
  digits = max(3, getOption("digits")-3), ...){
  # check for class "sarprobit"
  if (!inherits(object, "sarprobit")) 
        stop("use only with \"sarprobit\" objects")
        
  nobs      <- object$nobs
  nvar      <- object$nvar
  ndraw     <- object$ndraw
  nomit     <- object$nomit
  draws     <- object$B
  
  #bayesian estimation
  bout_mean <- object$coefficients                         #parameter mean column
  bout_sd   <- apply(draws, 2, sd)                         #parameter sd colum
  # build bayesian significance levels
  # for parameter > 0 count all draws > 0  and determine P(X <= 0)
  # for parameter <= 0 count all draws <0  and determine P(X >= 0)
  bout_sig <- 1 - apply(draws, 2, function(x) { ifelse (mean(x) > 0, sum(x > 0), sum(x < 0)) }) / ndraw
  #standard asymptotic measures
  bout_t    <- bout_mean / bout_sd             #t score b/se
  bout_tPval<- (1 - pt( abs(bout_t), nobs ))*2 #two tailed test = zero probability = z-prob
  #name definition
  if( is.null(var_names)){
    bout_names<- as.matrix(object$names)
  }else{
    bout_names<- as.matrix(var_names)
  }
  
  if(is.null(file)){file <- ""}#output to the console
  #HEADER
  write(sprintf("--------MCMC spatial autoregressive probit--------"), file, append=T)
  #sprintf("Dependent Variable")
  write(sprintf("Execution time  = %6.3f %s", object$time, attr(object$time, "units"))  , file, append=T)
  write(sprintf("N steps for TMVN= %6d"  , object$nsteps), file, append=T)
  write(sprintf("N draws         = %6d, N omit (burn-in)= %6d", ndraw, nomit), file, append=T)
  write(sprintf("N observations  = %6d, K covariates    = %6d", nobs, nvar)  , file, append=T)
  write(sprintf("# of 0 Y values = %6d, # of 1 Y values = %6d", object$zip, nobs - object$zip) , file, append=T)
  write(sprintf("Min rho         = % 6.3f, Max rho         = % 6.3f", object$rmin, object$rmax), file, append=T)
  write(sprintf("--------------------------------------------------"), file, append=T)
  write(sprintf(""), file, append=T)
  #ESTIMATION RESULTS
  coefficients <- cbind(bout_mean, bout_sd, bout_sig, bout_t, bout_tPval)
  dimnames(coefficients) <- list(bout_names, 
        c("Estimate", "Std. Dev", "p-level", "t-value", "Pr(>|z|)"))
  printCoefmat(coefficients, digits = digits,
    signif.stars = getOption("show.signif.stars"))      
  #if (getOption("show.signif.stars")) { 
    # The problem: The significance code legend printed by printCoefMat()
    # is too wide for two-column layout and it will not wrap lines...
    # The solution: using cat() instead of print() and use line breaks
    # cat(paste(strwrap(x, width = 70), collapse = "\\\\\n"), "\n")
    # http://r.789695.n4.nabble.com/Sweave-line-breaks-td2307755.html
    #  Signif <- symnum(1e-6, corr = FALSE, na = FALSE,
    #              cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
    #              symbols = c("***", "**", "*", ".", " "))
    #x <- paste("Signif. codes: ", attr(Signif, "legend"), "\n", sep="")
    #cat(paste(strwrap(x, width = getOption("width")), collapse = "\\\n"), "\n")
  #}
  return(invisible(coefficients))
}

# prints the marginal effects/impacts of the SAR probit fit
impacts.sarprobit <- function(obj, file=NULL, 
  digits = max(3, getOption("digits")-3), ...) {
  if (!inherits(obj, "sarprobit")) 
    stop("use only with \"sarprobit\" objects")
    
  if(is.null(file)){file <- ""}#output to the console
  write(sprintf("--------Marginal Effects--------"), file, append=T)  
  write(sprintf(""), file, append=T)      
  
  #(a) Direct effects
  write(sprintf("(a) Direct effects"), file, append=T)
  direct <- cbind(
   lower_005=apply(obj$direct, 2, quantile, prob=0.05),
   posterior_mean=colMeans(obj$direct),
   upper_095=apply(obj$direct, 2, quantile, prob=0.95)
  )
  printCoefmat(direct, digits = digits)
  
  # (b) Indirect effects
  write(sprintf(""), file, append=T)      
  write(sprintf("(b) Indirect effects"), file, append=T)
  indirect <- cbind(
   lower_005=apply(obj$indirect, 2, quantile, prob=0.05),
   posterior_mean=colMeans(obj$indirect),
   upper_095=apply(obj$indirect, 2, quantile, prob=0.95)
  )
  printCoefmat(indirect, digits = digits)
  
  # (c) Total effects
  write(sprintf(""), file, append=T)      
  write(sprintf("(c) Total effects"), file, append=T)
  total <- cbind(
   lower_005=apply(obj$total, 2, quantile, prob=0.05),
   posterior_mean=colMeans(obj$total),
   upper_095=apply(obj$total, 2, quantile, prob=0.95)
  )
  printCoefmat(total, digits = digits)
  
  return(invisible(list(direct=direct, indirect=indirect, total=total)))
}

# concatenate two objects of class "sarprobit"
# c.sarprobit works in the same way as boot:::c.boot().
c.sarprobit <- function(...) {
 args <- list(...)
 nm <- lapply(args, names)
 if (!all(sapply(nm, function(x) identical(x, nm[[1]]))))
   stop("arguments are not all the same type of \"sarprobit\" object")
 res <- args[[1]]
 res$time  <- as.difftime(max(as.numeric(sapply(args, "[[", "time"), units="secs")), units="secs") 
 res$ndraw <- sum(sapply(args, "[[", "ndraw")) 
 res$nomit <- sum(sapply(args, "[[", "nomit"))
 res$B     <- do.call(rbind, lapply(args, "[[", "B"))
 res$bdraw <- do.call(rbind, lapply(args, "[[", "bdraw"))
 res$pdraw <- do.call(rbind, lapply(args, "[[", "pdraw"))
 res$coefficients <- colMeans(res$B) 
 res$beta  <- res$coefficients[1:res$nvar]
 res$rho   <- res$coefficients[res$nvar+1]
 res
}


# extract the coefficients
coef.sarprobit <- function(object, ...) {
 if (!inherits(object, "sarprobit")) 
        stop("use only with \"sarprobit\" objects")
 return(object$coefficients)
}

# extract the coefficients
coefficients.sarprobit <- function(object, ...) {
 UseMethod("coef", object)
}

# plot MCMC results for class "sarprobit" (draft version);
# diagnostic plots for results (trace plots, ACF, posterior density function)
# method is very similar to plot.lm()
#
# @param x
# @param which
# @param ask
# @param trueparam a vector of "true" parameter values to be marked in posterior density plot
plot.sarprobit <- function(x, which=c(1, 2, 3), 
  ask = prod(par("mfcol")) < length(which) && dev.interactive(), ..., trueparam=NULL) {
 if (!inherits(x, "sarprobit")) 
        stop("use only with \"sarprobit\" objects")
 if (!is.numeric(which) || any(which < 1) || any(which > 3)) 
        stop("'which' must be in 1:3")
        
 names <- x$names
 B <- x$B
 k <- ncol(B)
  
 show <- rep(FALSE, 3)
 show[which] <- TRUE
 
 if (ask) {
   oask <- devAskNewPage(TRUE)
   on.exit(devAskNewPage(oask))
 }
 if (show[1L]) {
  # trace plots
  for (i in 1:k) {
    plot(1:nrow(B), B[,i], type="l", xlab="iteration", ylab=names[i], main=substitute("Trace plot of "*x, list(x=names[i])), ...)
    if (!is.null(trueparam)) abline(h=trueparam[i], col="red", lty=2)
  }
 }

 if (show[2L]) {
   # ACFs
   for (i in 1:k) {
     acf(B[,i], main=substitute("ACF of "*x, list(x=names[i])), ...)
   }
 }
 
 if (show[3L]) {
   # posterior distribution
   for (i in 1:k) {
     plot(density(B[,i]), main=substitute("Posterior distribution of "*x, list(x=names[i])), ...)
     if (!is.null(trueparam)) abline(v=trueparam[i], col="red", lty=2)
   }
 }
}

# return fitted values of SAR probit (on response scale vs. linear predictor scale)
# TODO: see fitted.lm() for comparison
fitted.sarprobit <- function(object, ...) {
  object$fitted.value
}

# Extract Log-Likelihood; see logLik.glm() for comparison
# Method returns object of class "logLik" with at least one attribute "df"
# giving the number of (estimated) parameters in the model.
# see Marsh (2000) equation (2.8), p.27 
logLik.sarprobit <- function(object, ...) {
  X <- object$X
  y <- object$y
  n <- nrow(X)
  k <- ncol(X)
  W <- object$W
  beta <- object$beta
  rho <- object$rho
  I_n <- sparseMatrix(i=1:n, j=1:n, x=1)
  S <- I_n - rho * W
  D <- diag(1/sqrt(diag(S %*% t(S))))  # D = diag(E[u u'])^{1/2}  (n x n)
  Xs <- D %*% solve(qr(S), X)          # X^{*} = D * (I_n - rho * W)^{-1} * X
  #Xs <- D %*% solve(S) %*% X
  F <- pnorm(as.double(Xs %*% beta))    # F(X^{*} beta)  # (n x 1)
  lnL <- sum(log(F[y == 1])) + sum(log((1 - F[y == 0]))) # see Marsh (2000), equation (2.8)
  #lnL <- sum(ifelse(y == 1, log(pnorm(xb)), log(1 - pnorm(xb))))
  out <- lnL
  class(out) <- "logLik"
  attr(out,"df") <- k+1                 # k parameters in beta, rho
 return(out)
}


