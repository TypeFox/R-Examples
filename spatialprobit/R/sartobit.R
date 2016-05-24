if (FALSE) {
 source("stats_distributions.r")
 source("rtnorm.R")
 source("SpatialProbit-MCMC.R")

}

# create a generic method
#sartobit <- function (y, ...)
# UseMethod("sartobit")
 
# Bayesian estimation of the SAR Tobit model
#
sartobit <- function(formula, W, data, ...) {
  cl <- match.call()                     # cl ist object of class "call"
  mf <- match.call(expand.dots = FALSE)  # mf ist object of class "call"
  m  <- match(c("formula", "data"), names(mf), 0L)        # m is index vector
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
  sar_tobit_mcmc(y, X, W, ...)    
} 

#Bayesian estimates of the spatial autoregressive tobit model
#          y = rho*W*y + XB + e, e = N(0,sige*Omega), Omega = inv[(I_n-rho*W)'*(I_n -rho*W)]
#          y is a nx1 vector of either 0 or > 0
#          B = N(c,T),
#          1/sige = Gamma(nu,d0),
#          rho = beta(a1,a2) prior
#
# Code based on James P. LeSage method sart_g.m; ported to R by Stefan Wilhelm
#
# --------------------------------------------------------------
# REFERENCES: LeSage and Pace (2009)
# Introduction to Spatial Econometrics, Chapter 10.
#----------------------------------------------------------------
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
#   prior$lflag = 0 for full lndet computation (default = 1, fastest)
#               = 1 for Chebyshev approx
#               = 2 for MC approx (fast for large problems)
# @param start
# @param m
# @param computeMarginalEffects
# @param showProgress
sar_tobit_mcmc <- function(y, X, W, ndraw=1000, burn.in=100, thinning=1,
  prior=list(a1=1, a2=1, c=rep(0, ncol(X)), T=diag(ncol(X))*1e12, lflag = 0),
  start=list(rho=0.75, beta=rep(0, ncol(X)), sige=1),
  m=10, computeMarginalEffects=FALSE, showProgress=FALSE){

  #start timer
  timet <- Sys.time()

  n  <- nrow( X )            # number of observations
  n1 <- nrow( X )
  n2 <- nrow( W )
  k <- ncol( X )             # number of of parameters/exogenous variables
  I_n <- sparseMatrix(i=1:n, j=1:n, x=1) # sparse identity matrix
  if (is.null(colnames(X))) colnames(X) <- paste("x",1:k,sep="")

  #validate inputs
  # TODO: uncomment
  #if( any(y < 0) ){
  #  stop('sartobit: not all y-values are greater or equal to 0')
  #}
  if( n1 != n2 && n1 != n ){
    stop('sartobit: wrong size of spatial weight matrix W')
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
    stop('sartobit: intercept term must be in first column of the X-matrix')
  }

  # MCMC sampling of beta
  rho  <- start$rho          # start value of row
  beta <- start$beta         # start value of parameters, prior value, we could also sample from beta ~ N(c, T)
  sige <- start$sige         # start value for sigma_e

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
  # prior for sige ~ IG(nu, d0)
  if (is.numeric(prior$nu)) {
    nu <- prior$nu
  } else {
    nu <- 0
  }
  if (is.numeric(prior$d0)) {
    d0 <- prior$d0
  } else {
    d0 <- 0
  }
  
  cc <- 0.1      # tuning factor for Metropolis-Hastings
  metflag <- 0   # 0 = Grid-based approach, 1 = Metropolis-Hastings
  nsample <- 5   # steps for updating z-values

  TI <- solve(T)             # T^{-1}
  TIc <- TI%*%c              # T^{-1}c
  S <- I_n - rho * W
  H <- t(S) %*% S            # precision matrix H for beta | rho, z, y
  QR <- qr(S)                # class "sparseQR"
  mu <- solve(QR, X %*% beta)
  
  rmin       <- -1   # use -1,1 rho interval as default
  rmax       <-  1

  lflag <- 0
  if (is.numeric(prior$lflag) && lflag %in% c(0, 1)) lflag <- prior$lflag
  #lflag=0 --> default to 1997 Pace and Barry grid approach
  #lflag=1 --> 2004 Pace and LeSage Chebyshev approx
  #lflag=2 --> 1999 Pace and Barry MC determinant approx
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
  nrho     <- nrow(detval)
  nmk      <- (n-k)/2
  detval1  <- detval[,1]  # rho grid values
  detval2    <- detval[,2]  # log-determinant grid values
  rho_gridsq <- detval1 * detval1
  yy       <- (detval1[2:nrho] + detval1[1:(nrho-1)])

  # matrix to store the beta + sige + rho parameters for each iteration/draw
  B <- matrix(NA, ndraw, k+2)
  
  ones <- rep(1, n)  # vector of ones

  # progress bar
  if (showProgress) {
    pb <- txtProgressBar(min=0, max=(thinning * ndraw + burn.in), initial=0, style=3)
  }
  
  # immutable matrices
  tX <- t(X)                       # X'               # k x n
  xpx  <- t(X) %*% X               # (X'X)            # k x k
  xpxI <- solve(xpx)               # (X'X)^{-1}       # k x k
  xxpxI    <- X %*% xpxI           # X(X'X)^(-1)      # n x k (better, compromise)

  # matrices for direct and indirect impacts
  direct       <- matrix(NA, ndraw,p)    # n x p
  indirect     <- matrix(NA, ndraw,p)    # n x p
  total        <- matrix(NA, ndraw,p)    # n x p

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
  
  yin <- y    # input values for y
  ind1 <- which(yin == 0)    # index vector for censored observations y=0
  nobs1 <- length(ind1)      # number of censored observations
  ind2 <- which(yin > 0)
  nobs2 <- length(ind2)

  for (i in (1 - burn.in):(ndraw * thinning)) {
  
    # update beta | y, rho, sige
    AI <- solve(xpx + sige*TI)
    Sy <- as.double(S%*%y)            # (n x 1)
    b  <- tX%*%Sy + sige*TIc          # (k x 1)
    b0 <- AI %*% b
    beta <- as.double(rmvnorm(n=1, mean=b0, sigma=sige*AI))
    Xb <- X%*%beta

    # update sige
    nu1 <- n + 2*nu
    e <- (Sy - Xb)              # (n x 1), residuals
    d1 <- 2*d0 + crossprod(e)
    chi <- rchisq(n=1,df=nu1)
    sige <- as.double(d1/chi)
        
    # draw rho
    if (metflag == 1) {
      # metropolis step to get rho update
      rhox <- c_sar(rho,y,Xb,sige,I_n,W,detval1,detval2)
      accept <- 0
      rho2 <- rho + cc*rnorm(n=1)
      while (accept == 0) {
        if ((rho2 > rmin) & (rho2 < rmax)) {
         accept <- 1
        } else {
         rho2 <- rho + cc*rnorm(n=1)
        }
      }
      rhoy <- c_sar(rho2,y,Xb,sige,I_n,W,detval1,detval2)
      ru <- runif(n=1,0,1)
      if ((rhoy - rhox) > exp(1)) {
        pp <- 1
      } else {
        ratio <- exp(rhoy-rhox)
        pp <- min(1,ratio)
      }
      if (ru < pp) {
        rho <- rho2
        acc <- acc + 1
      }
      iter <- i + burn.in
      acc_rate[iter] <- acc/iter;
      # update cc based on std of rho draws
      if (acc_rate[iter] < 0.4) {
       cc <- cc/1.1
      }
      if (acc_rate[iter] > 0.6) {
       cc <- cc*1.1
      }
    } # end draw rho Metropolis-Hastings
    
    if (metflag == 0) {
      # when metflag == 0,
      # we use numerical integration to perform rho-draw
      xpy  <- tX %*% y            # X'y
      Wy   <- as.double(W %*% y)  # Wy (n x 1) # SW: coerce Wz to vector
      # (from n x 1 sparse matrix! we do not need a sparse matrix here)
      xpWy <- tX %*% Wy           # X'Wy      # k x 1
      e0   <-  y - xxpxI %*% xpy  # y  - X(X'X)^-1X' y
      ed   <- Wy - xxpxI %*% xpWy # Wy - X(X'X)^(-1)X'Wy
      epe0 <- as.double(crossprod(e0))  # slightly faster than t(e0) %*% e0
      eped <- as.double(crossprod(ed))
      epe0d<- as.double(crossprod(ed, e0))
      rho  <- draw_rho(detval1, detval2, rho_gridsq, yy, epe0, eped, epe0d, rho, nmk=nmk, nrho=nrho, lnbprior, u=u[i + burn.in])
    }
    
    if (nobs1 > 0) {
    # update z-values for all zero-observations
    z <- rep(0, n)

    # loop over i
    S <- (I_n - rho*W)
    mu <- as.double(solve(qr(S), X %*% beta)) # (n x 1)
    
    # see LeSage(2009), chapter 10, 10.3 Spatial Tobit models, p.299-305
    # precision matrix H / Phi partioned as [ind1, ind2] where ind1 are censored obs
    # H = [ H_11  H_12 ]
    #     [ H_21  H_22 ]
    H <-(t(S)%*%S)/sige 
    H_11 <- H[ind1,ind1]   
    H_12 <- H[ind1,ind2]
    HI_11 <- 1/diag(H)[ind1]  
    
    # conditional mean and precision matrix for censored observations
    n1 <- length(ind1)
    mu_1 <- mu[ind1] -  HI_11 * H_12 %*% (y[ind2] - mu[ind2])  # (n1 x 1)
    m <- 5  
    z[ind1] <- rtmvnorm.sparseMatrix(n=1, mean=mu_1, H=H_11, 
      lower=rep(-Inf, n1), upper=rep(0, n1), burn.in=m, start.value=rep(0, n1))
    y[ind1] <- z[ind1]
    
    #browser()
    #plot(density(y[ind1]))
    #lines(density(ysave[ind1]), col="red")
    }
    
    if (i > 0) {
      if (thinning == 1) {
        ind <- i
      }
      else if (i%%thinning == 0) {
        ind <- i%/%thinning
      } else {
        next
      }
      B[ind,] <- c(beta, sige, rho)
      
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
        
        # Marginal effects in SAR Tobit
        probz <- pnorm(as.numeric(mu) / sige)     # Phi(Xß/sigma) = Prob(z <= Xß/sigma)
        dd <- sparseMatrix(i = 1:n, j = 1:n, x = probz)
        dir <- as.double(t(probz) %*% trW.i %*% rhovec/n)
        avg_direct <- dir * beff
        avg_total <- mean(dd %*% qr.coef(QR, ones)) * beff
        avg_indirect <- avg_total - avg_direct

        total[ind, ]      <- avg_total    # an (ndraw-nomit x p) matrix
        direct[ind, ]     <- avg_direct   # an (ndraw-nomit x p) matrix
        indirect[ind, ]   <- avg_indirect # an (ndraw-nomit x p) matrix
      }
    }

    if (showProgress) setTxtProgressBar(pb, i + burn.in) # update progress bar
  } # end loop

  if (showProgress)  close(pb) #close progress bar

  beta  <- colMeans(B)[1:k]
  sige  <- colMeans(B)[k+1]
  rho   <- colMeans(B)[k+2]
  
  S     <- (I_n - rho * W)
  fitted.values   <- solve(qr(S), X %*% beta) # E[z | beta] = (I_n - rho * W)^{-1}(X * beta)
  fitted.response <- pmax(as.numeric(fitted.values), 0)   # y = max(z, 0)
  
  # result
  results       <- NULL
  results$time  <- Sys.time() - timet
  results$nobs  <- n          # number of observations
  results$nvar  <- k          # number of explanatory variables
  results$y     <- y
  results$zip   <- nobs1      # number of zero values in the y-vector
  results$beta  <- colMeans(B)[1:k]
  results$sige  <- colMeans(B)[k+1]
  results$rho   <- colMeans(B)[k+2]
  results$coefficients <- colMeans(B)
  results$fitted.values   <- fitted.values    # E[z | beta]
  results$fitted.response <- fitted.response  # fitted values on response scale (censored y variable)
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
  results$names     <- c(colnames(X), "sige", 'rho')
  results$B         <- B        # (beta, rho) draws
  results$bdraw     <- B[,1:k]  # beta draws
  results$sdraw     <- B[,k+1]  # sige draws
  results$pdraw     <- B[,k+2]  # rho draws
  results$total     <- total
  results$direct    <- direct
  results$indirect  <- indirect
  results$W <- W
  results$X <- X

  class(results)    <- "sartobit"
  return(results)
}

# PURPOSE: evaluate the conditional distribution of rho given sige
#  spatial autoregressive model using sparse matrix algorithms
# ---------------------------------------------------
#  USAGE:cout = c_sar(rho,y,x,b,sige,I_n,W,detval,p,R)
#  where:  rho  = spatial autoregressive parameter
#          y    = dependent variable vector
#          W    = spatial weight matrix
#        detval = an (ngrid,2) matrix of values for det(I-rho*W)
#                 over a grid of rho values
#                 detval(:,1) = determinant values
#                 detval(:,2) = associated rho values
#          sige = sige value
#          c    = (optional) prior mean for rho
#          T    = (optional) prior variance for rho
# ---------------------------------------------------
#  RETURNS: a conditional used in Metropolis-Hastings sampling
#  NOTE: called only by sar_g
#  --------------------------------------------------
#  SEE ALSO: sar_g, c_far, c_sac, c_sem
# ---------------------------------------------------
c_sar <- function(rho,y,xb,sige,I_n,W,detval1,detval2,c=NULL,T=NULL) {

i <- findInterval(rho,detval1)
if (i == 0) index=1
else index=i
detm <- detval2[index]

if (is.null(c) && is.null(T)) {        # case of diffuse prior
 n <- length(y)
 z <- I_n - rho * W
 e <- z%*%y - xb
 epe <- crossprod(e)/(2*sige)
}
else if (c == NULL) {    # case of informative prior rho ~ N(c,T)
 z <- I_n - rho * W
 e <- z%*% - xb
 n <- length(y)
 T <- T*sige
 z <- (speye(n) - rho*W)*e
 epe <- ((t(z)%*%z)/2*sige) + 0.5*(((rho-c)^2)/T)
}
cout <- as.double(detm - epe)
return(cout)
}

# extract the coefficients
coef.sartobit <- function(object, ...) {
 if (!inherits(object, "sartobit")) 
        stop("use only with \"sartobit\" objects")
 return(object$coefficients)
}

# extract the coefficients
coefficients.sartobit <- function(object, ...) {
 UseMethod("coef", object)
}


# summary method for class "sartobit"
summary.sartobit <- function(object, var_names=NULL, file=NULL, 
  digits = max(3, getOption("digits")-3), ...){
  # check for class "sarprobit"
  if (!inherits(object, "sartobit")) 
        stop("use only with \"sartobit\" objects")
        
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
  write(sprintf("----MCMC spatial autoregressive Tobit model ----"), file, append=T)
  #sprintf("Dependent Variable")
  write(sprintf("Execution time  = %6.3f %s", object$time, attr(object$time, "units"))  , file, append=T)
  write(sprintf("N steps for TMVN= %6d"  , object$nsteps), file, append=T)
  write(sprintf("N draws         = %6d, N omit (burn-in)= %6d", ndraw, nomit), file, append=T)
  write(sprintf("N observations  = %6d, K covariates    = %6d", nobs, nvar)  , file, append=T)
  write(sprintf("# censored values = %6d, # observed values = %6d", object$zip, nobs - object$zip) , file, append=T)
  write(sprintf("Min rho         = % 6.3f, Max rho         = % 6.3f", object$rmin, object$rmax), file, append=T)
  write(sprintf("--------------------------------------------------"), file, append=T)
  write(sprintf(""), file, append=T)
  #ESTIMATION RESULTS
  coefficients <- cbind(bout_mean, bout_sd, bout_sig, bout_t, bout_tPval)
  dimnames(coefficients) <- list(bout_names, 
        c("Estimate", "Std. Dev", "p-level", "t-value", "Pr(>|z|)"))
  printCoefmat(coefficients, digits = digits,
    signif.stars = getOption("show.signif.stars"))      
  return(invisible(coefficients))
}

# compute marginal effects for every MCMC iteration of the 
# estimated SAR Tobit model
# @param object fitted model object
# @param o degree of approximating the tr(W^i)
marginal.effects.sartobit <- function (object, o = 100, ...) {
    if (!inherits(object, "sartobit")) 
        stop("use only with \"sartobit\" objects")
    nobs <- object$nobs
    nvar <- object$nvar
    p <- ifelse(object$cflag == 0, nvar, nvar - 1)
    ndraw <- object$ndraw
    nomit <- object$nomit
    betadraws <- object$bdraw
    sigedraws <- object$sdraw
    rhodraws <- object$pdraw
    
    X <- object$X
    W <- object$W
    I_n <- sparseMatrix(i = 1:nobs, j = 1:nobs, x = 1)
    trW.i <- tracesWi(W, o = o, iiter = 50)
    D <- matrix(NA, ndraw, p)
    I <- matrix(NA, ndraw, p)
    T <- matrix(NA, ndraw, p)
    if (object$cflag == 0) {
        namesNonConstantParams <- colnames(X)
    }
    else {
        namesNonConstantParams <- colnames(X)[-1]
    }
    colnames(T) <- namesNonConstantParams
    colnames(D) <- namesNonConstantParams
    colnames(I) <- namesNonConstantParams
    ones <- rep(1, nobs)
    for (i in 1:ndraw) {
        beta  <- betadraws[i, ]
        beff  <- betadraws[i, -1]    # TODO: check for constant
        rho   <- rhodraws[i]
        sigma <- sigedraws[i]
        S <- I_n - rho * W
        QR <- qr(S)
        mu <- qr.coef(QR, X %*% beta)
        rhovec <- rho^(0:(o - 1))
        probz <- pnorm(as.numeric(mu) / sigma)     # Phi(Xß/sigma) = Prob(z <= Xß/sigma)
        dd <- sparseMatrix(i = 1:nobs, j = 1:nobs, x = probz)
        dir <- as.double(t(probz) %*% trW.i %*% rhovec/nobs)
        avg_direct <- dir * beff
        avg_total <- mean(dd %*% qr.coef(QR, ones)) * beff
        avg_indirect <- avg_total - avg_direct
        D[i, ] <- avg_direct
        I[i, ] <- avg_indirect
        T[i, ] <- avg_total
    }
    summaryMarginalEffects <- function(x) {
        r <- cbind(apply(x, 2, mean), apply(x, 2, sd), apply(x, 
            2, mean)/apply(x, 2, sd))
        colnames(r) <- c("marginal.effect", "standard.error", 
            "z.ratio")
        return(r)
    }
    summary_direct <- summaryMarginalEffects(D)
    summary_indirect <- summaryMarginalEffects(I)
    summary_total <- summaryMarginalEffects(T)
    return(list(direct = D, indirect = I, total = T, summary_direct = summary_direct, 
        summary_indirect = summary_indirect, summary_total = summary_total))
}


# prints the marginal effects/impacts of the SAR Tobit fit
impacts.sartobit <- function(obj, file=NULL, 
  digits = max(3, getOption("digits")-3), ...) {
  if (!inherits(obj, "sartobit")) 
    stop("use only with \"sartobit\" objects")
    
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



# plot MCMC results for class "sartobit" (draft version);
# diagnostic plots for results (trace plots, ACF, posterior density function)
# method is very similar to plot.lm()
#
# @param x
# @param which
# @param ask
# @param trueparam a vector of "true" parameter values to be marked in posterior density plot
plot.sartobit <- function(x, which=c(1, 2, 3),
  ask = prod(par("mfcol")) < length(which) && dev.interactive(), ..., trueparam=NULL) {
 if (!inherits(x, "sartobit"))
        stop("use only with \"sartobit\" objects")
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

# return fitted values of SAR Tobit
fitted.sartobit <- function(object, ...) {
  object$fitted.value
}
