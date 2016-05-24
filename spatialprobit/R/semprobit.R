# Estimating a Probit Model with Spatial Errors (SEM Probit)
# autocorrelation in the error rather than on a lag variable
#
# References/Applications von Spatial Probit
#
# (1) Marsh (2000): 
#     "Probit with Spatial Correlation by Field Plot: Potato Leafroll Virus
#	    Net Necrosis in Potatoes"
# (2) Coughlin (2003): 
#     "Spatial probit and the geographic patterns of state lotteries"
# (3) Spatial Probit estimation of Freetrade agreements: http://ideas.repec.org/p/gii/giihei/heidwp07-2010.html
# (4) Probit with spatial correlation: http://www.jstor.org/stable/1400629
# (6) Jaimovich (2012)
#
# SAR Probit:
# (5) LeSage (2011)

# Comparison of SAR Probit vs. SEM Probit:
#
# (a) SAR Probit
# - model: 
#   z = pWz + xB + e, e ~ N(0, sige*I)
#   Sz = xB + e  where S = (I - pW)
# - DGP:
#   z = (I - pW)^(-1)xB + (I - pW)^(-1)e
# - observed:
#   y = 1, if (z >= 0)
#   y = 0, if (z < 0)
# - model parameters: p, B (sige=1 for identification)
# - conditional distributions:
#   (aa) p(z | B,p,y)  ~  TN((I - pW)^(-1) xB, [(I - pW)'(I - pW)]^(-1))
#   (bb) p(p | z, B)   ~  |I - pW| exp(-1/2*(Sz - xB)'(Sz - xB))
#        Sz - xB = e (e ~ N(i.i.d))
#   (cc) p(B | z,p,y)    ~ N(c*,T*)
#
#
# (b) SEM Probit:
# - autocorrelation in the error rather than on a lag variable:
# - Modell:
#   z = xB + u
#   u = pWu + e, e ~ N(0, sige*I)
# - observed:
#   y = 1, if (z >= 0)
#   y = 0, if (z < 0)
# - DGP: z = xB + (I - pW)^(-1)e 
# - model parameters: p, B, sige>0
# - conditional distributions:
#   (aa) p(z | B,p,y,sige) ~ TN(xB, sige*[(I - pW)'(I - pW)]^(-1))
#   (bb) p(p | z, B, sige) ~ |I - pW| exp(-1/(2*sige)*(Sz - SxB)'(Sz - SxB))
#                  z = xB + S^(-1)e
#                  Sz - SxB = e
#   (cc) Erwartungswert von von Normalverteilung p(beta | z, rho) ändert sich
# - No spatial spillover, marginal effects in SEM Probit

# Bayesian estimation of the probit model with spatial errors (SEM probit)
#
# @param formula 
semprobit <- function(formula, W, data, subset, ...) {
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
  sem_probit_mcmc(y, X, W, ...)    
}


# Estimate the probit model with spatial errors (SEM probit)
# z = xB + u
# u = pWu + e
#
# Also it requires the truncated multivariate normal for estimation as well:
#  z = xB + (I - pW)^(-1)e 
# where y = 1 if z >= 0 and y = 0 if z < 0 observable
#
# @param y
# @param X
# @param W spatial weight matrix
# @param ndraw number of MCMC iterations
# @param burn.in  number of MCMC burn-in to be discarded
# @param thinning MCMC thinning factor, defaults to 1
# @param m number of burn.in sampling in inner Gibbs sampling
# @param prior list of prior settings (a1, a2, c, T, nu, d0) for rho ~ Beta(a1,a2), beta ~ N(c, T),
#    1/sige ~ Gamma(nu,d0) prior on sigma, default: nu=0,d0=0 (diffuse prior)
#
# @param start  start value for the chain; list of three components rho and beta and sigma.
# @param m
# @param showProgress
sem_probit_mcmc <- function(y, X, W, ndraw=1000, burn.in=100, thinning=1, 
  prior=list(a1=1, a2=1, c=rep(0, ncol(X)), T=diag(ncol(X))*1e12,
  nu=0, d0=0, lflag = 0), 
  start=list(rho=0.75, beta=rep(0, ncol(X)), sige=1),
  m=10, showProgress=FALSE, univariateConditionals=TRUE){  

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
    stop("semprobit: not all y-values are 0 or 1")
  }
  if( n1 != n2 && n1 != n ){
    stop("semprobit: wrong size of spatial weight matrix W")
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
    stop("semprobit: intercept term must be in first column of the X-matrix")
  }
  
  # MCMC sampling of beta
  rho  <- start$rho          # start value of row
  beta <- start$beta         # start value of parameters, prior value, we could also sample from beta ~ N(c, T)
  sige <- start$sige         # start value for sigma_e
  
  # conjugate prior beta ~ N(c, T)
  # parametrize, default to diffuse prior, for beta, e.g. T <- diag(k) * 1e12
  c <- rep(0, k)             # prior distribution of beta ~ N(c, T) : c = 0
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
  
  
  TI <- solve(T)           # T^{-1}
  TIc <- TI%*%c            # T^{-1}c
  S <- I_n - rho * W
  H <- t(S) %*% S / sige   # precision matrix H for beta | rho, z, y, sige
  
  # truncation points for z, depend only on y, can be precalculated
  lower <- ifelse(y > 0, 0,  -Inf)
  upper <- ifelse(y > 0, Inf,   0)
  
  rmin       <- -1   # use -1,1 rho interval as default
  rmax       <-  1
  ldetflag   <-  0   # default to 1999 Pace and Barry MC determinant approx
  tmp <- sar_lndet(ldetflag, W, rmin, rmax)
  detval <- tmp$detval
  
  # Some precalculated quantities for drawing rho
  # rho ~ Beta(a1, a2) prior
  a1         <-  1.0
  a2         <-  1.0
  if(is.numeric(prior$a1)) a1 <- prior$a1
  if(is.numeric(prior$a2)) a2 <- prior$a2
  
  u        <- runif(thinning * ndraw + burn.in)   # u ~ U(0, 1)
  nrho     <- nrow(detval)
  nmk      <- (n-k)/2
  detval1  <- detval[,1]  # SW: avoid multiple accesses to detval[,1]
  detval2  <- detval[,2]
      
  # matrix to store the beta + sige + rho parameters for each iteration/draw
  B <- matrix(NA, ndraw, k+2)
  colnames(B) <- c(colnames(X), "sige", "rho")
  
  cc <- 0.2         # initial tuning parameter for M-H sampling
  acc <- 0          # number of accepted samples 
  acc_rate <- rep(NA, thinning * ndraw + burn.in)
  
  # progress bar
  if (showProgress) {
    pb <- txtProgressBar(min=0, max=(thinning * ndraw + burn.in), initial=0, style=3)
  }
  
  # names of non-constant parameters
  if(cflag == 0) {
    namesNonConstantParams <- colnames(X)
  } else {
    namesNonConstantParams <- colnames(X)[-1]
  }
    
  # just to set a start value for z
  #z <- rep(0, n)
  z <- y
  ones <- rep(1, n)
  W2diag <- diag(t(W)%*%W)
  
  ind0 <- which(y == 0)      # 0 obs.
  ind1 <- which(y == 1)      # 1 obs.
  
  for (i in (1 - burn.in):(ndraw * thinning)) {
    
  # update beta given rho, z  
  SX  <- S %*% X                          # n x k matrix; sparse
  tSX <- t(SX)                            # k x n matrix; sparse
  tSXSX <- as.matrix(tSX %*% SX)          # k x k matrix; dense
  AI <- solve(tSXSX + sige * TI)          # k x k matrix; dense
  Sz <- as.double(S %*% z)                # n x 1 vector
  b <- as.double(tSX %*% Sz + sige * TIc) # k x 1 vector
  b0 <- AI %*% b                          # k x 1 vector
  beta <- as.double(rmvnorm(n=1, mean=b0, sigma=sige*AI))
  
  # update sige given rho, beta, z
  nu1 <- n + 2*nu
  e <- as.double(S %*% (z - X %*% beta))  # n x 1 vector
  d1 <- 2*d0 + crossprod(e)
  chi <- rchisq(n=1,df=nu1)
  sige <- as.double(d1/chi)
  
  # Update H = 1/sige*t(S)%*%S after update of sige
  H <- t(S) %*% S / sige
  
  # update z-values given beta, sigma, rho and old z:
  #
  # univariate conditional distributions as univariate truncated normals:
  # p(z_{i}^(t) | z_{-i}^(t-1), beta, rho, sige)
  #
  # E[z_i^(t) | z_{-i}^(t-1)] = mu_i - H_{ii}^{-1} H_{i,-i} [ z_{-i}^(t-1) - mu_{-i}^(t-1) ]
  #                           = mu_i - H_{ii}^{-1} H_{i,.} [ z^(t-1) - mu^(t-1) ] + (z_{i}^(t-1) - mu_i^(t-1))
  # 
  # Vektorisiert für alle z_i | z_{-i} untereinandergeschrieben:
  # E[z^(t) | z_{-i}^(t-1) ]  = mu^(t-1)  - diag(H)^{-1} H [ z^(t-1) - mu^(t-1) ] + [ z^(t-1) - mu^(t-1) ]
  # 
  # H = 1/sige * (I_n - rho * W)'(I_n - rho * W) = 1/sige * (I_n - 2*rho*W + rho^2*W^2)
  # diag(H) = H_{ii} = 1/sige * diag(I_n - 2*rho*W + rho^2*W^2) 
  #                  = 1/sige * (diag(I_n) + rho^2*diag(W^2))
  # because diag(W) = 0!
  #
  #
  if (univariateConditionals) {
    # conditional variance  z_i | z_{-i}
    dsig <- 1/sige * (ones - rho * rho * W2diag)  # SW: Sollte es nicht ones + rho * rho * W2diag sein?
    zvar <- ones/dsig;            # conditional variances for each z_i | z = H_{ii}^{-1} = 1/diag(H) (n x 1)
                                  # TODO: check if sige is missing in zvar
                                  # TODO: Was passiert im Fall dsig < 0 --> zvar < 0 (negative variance) !!!
    #browser()
    # all.equal(zvar, 1 / diag(H))    # TRUE zvar = diag(H)^{-1} 
    
    # conditional mean  z_i | z_{-i}
    mu <- X %*% beta
    zmu <- z - mu                  
    A  <- (1/sige)* S %*% zmu     # a vector (n x 1)
    B2  <- t(S) %*% A             # B2 <- (1/sige) * t(S) %*% S %*% zmu
    Cz <- zmu - zvar*B2           # Cz = (z - mu) - diag(H)^{-1} * H * (z - mu)
    zm <- mu + Cz;                # mu + (z-mu) - zvar * t(S)[ (1/sige) * S * (z - mu) ]  =  mu + (z-mu)
    # zm2 <- mu - 1 / diag(H) * H %*% zmu  + zmu
    # all.equal(zm, zm2)          # TRUE
    
    z[zvar < 0] <- 0
    z[ind0] <- rtnorm(mu=zm[ind0], sd=sqrt(zvar[ind0]), a=-Inf, b=0)
    z[ind1] <- rtnorm(mu=zm[ind1], sd=sqrt(zvar[ind1]), a=0, b=Inf)
    z[is.infinite(z) | zvar < 0] <- 0    # some zvar become negative which causes sqrt(zvar) to be NA.
  }
  # TODO: check why sampling with rtmvnorm.sparseMatrix does not work. 
  # Chain is exploding in this case. Conditional variance is different from LeSage code.
  # multivariate truncated normal given beta, rho, sige
  if (!univariateConditionals) {
    mu <- X %*% beta
    H <- (1/sige)*t(S)%*%S
    if (m==1) {
     z <- as.double(rtmvnorm.sparseMatrix(n=1, mean=mu, H=H, 
          lower=lower, upper=upper, burn.in=m, start.value=z))
    } else {
     z <- as.double(rtmvnorm.sparseMatrix(n=1, mean=mu, H=H, 
       lower=lower, upper=upper, burn.in=m))
    }
  }
      
  # 3. sample from rho | beta, z, sige using Metropolis-Hastings with burn.in=20
  # update rho using metropolis-hastings
  # numerical integration is too slow here
  rhox <- c_sem(rho,z,X,beta,sige,I_n,W,detval1,detval2,ones,a1,a2)
  accept <- 0
  rho2 <- rho + cc * rnorm(1)
  while(accept == 0) {
    if ((rho2 > rmin) & (rho2 < rmax)) { 
      accept <- 1
    }  
    else {
      rho2 <- rho + cc * rnorm(1)
    } 
  }
  rhoy <- c_sem(rho2,z,X,beta,sige,I_n,W,detval1,detval2,ones,a1,a2)
  ru <- runif(1,0,1) # TODO: Precalculate ru
  if ((rhoy - rhox) > exp(1)) {
    p <- 1
  } else {
    ratio <- exp(rhoy-rhox)
    p <- min(1,ratio)
  }
  if (ru < p) {
    rho <- rho2
    acc <- acc + 1
  }
  iter <- i + burn.in
  acc_rate[iter] <- acc/iter
  # update cc based on std of rho draws
  if (acc_rate[iter] < 0.4) {
    cc <- cc/1.1;
  }
  if (acc_rate[iter] > 0.6) {
    cc <- cc*1.1;
  }
    
  ############################################################################## 
  
  # update S and H after update of rho
  S <- I_n - rho * W
  H <- t(S) %*% S / sige
  
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
  
  ##############################################################################
    
  }
    
  if (showProgress) setTxtProgressBar(pb, i + burn.in) # update progress bar
  }
    
  if (showProgress)  close(pb) #close progress bar
    
  # fitted values for estimates (based on z rather than binary y like in fitted(glm.fit))
  # (on response scale y vs. linear predictor scale z...)
  beta  <- colMeans(B)[1:k]
  sige  <- colMeans(B)[k+1]
  rho   <- colMeans(B)[k+2]
  S     <- (I_n - rho * W)
  fitted.values   <- X %*% beta                     # E[z | beta] = (X * beta)
  fitted.response <- as.numeric(fitted.values >= 0) # y = (z >= 0)
  # TODO: linear.predictors  vs. fitted.values
 
  # result
  results       <- NULL
  results$time  <- Sys.time() - timet
  results$nobs  <- n          # number of observations
  results$nvar  <- k          # number of explanatory variables
  results$y     <- y 
  results$zip   <- n - sum(y) # number of zero values in the y-vector
  results$beta  <- colMeans(B)[1:k]
  results$sige  <- sige
  results$rho   <- colMeans(B)[k+2]
  results$coefficients <- colMeans(B)
  results$fitted.values <- fitted.values    # fitted values
  results$fitted.response <- fitted.response  # fitted values on reponse scale (binary y variable)
  results$ndraw <- ndraw
  results$nomit <- burn.in
  results$a1        <- a1
  results$a2        <- a2
  results$nu        <- nu
  results$d0        <- d0
  results$rmax      <- rmax 
  results$rmin      <- rmin
  results$tflag     <- "plevel"
  results$lflag     <- ldetflag
  results$cflag     <- cflag
  results$lndet     <- detval
  results$names     <- c(colnames(X), "sige", "rho")
  results$B         <- B        # (beta, sige, rho) draws
  results$bdraw     <- B[,1:k]  # beta draws
  results$sdraw     <- B[,k+1]  # sige draws
  results$pdraw     <- B[,k+2]  # rho draws
  results$W <- W
  results$X <- X
  #results$predicted <- # prediction required. The default is on the scale of the linear predictors
  class(results)    <- "semprobit"
  return(results)
}

#% PURPOSE: evaluate the conditional distribution of rho given sige
#%  spatial autoregressive model using sparse matrix algorithms
#% ---------------------------------------------------
#%  USAGE:cout = c_sar(rho,y,x,b,sige,W,detval,a1,a2)
#%  where:  rho  = spatial autoregressive parameter
#%          y    = dependent variable vector
#%          W    = spatial weight matrix
#%        detval = an (ngrid,2) matrix of values for det(I-rho*W) 
#%                 over a grid of rho values 
#%                 detval(:,1) = determinant values
#%                 detval(:,2) = associated rho values
#%          sige = sige value
#%          a1    = (optional) prior parameter for rho
#%          a2    = (optional) prior parameter for rho
#% ---------------------------------------------------
#%  RETURNS: a conditional used in Metropolis-Hastings sampling
#%  NOTE: called only by sar_g
#%  --------------------------------------------------
#%  SEE ALSO: sar_g, c_far, c_sac, c_sem
#% ---------------------------------------------------
# Code from James P. LeSage; ported to R by Stefan Wilhelm
c_sem <- function(rho,y,X,b,sige,I_n,W,detval1,detval2,vi,a1,a2) {
 i <- findInterval(rho,detval1)
 if (i == 0) index=1
 else index=i
 detm = detval2[index] 
 z = I_n - rho*W;
 e = as.double(z %*% (y - X %*% b))
 ev = e * sqrt(vi)
 epe = (crossprod(ev))/(2*sige)
 cout =  as.double(detm - epe)  # log-density
 return(cout)
}

# extract the coefficients
coef.semprobit <- function(object, ...) {
 if (!inherits(object, "semprobit")) 
        stop("use only with \"semprobit\" objects")
 return(object$coefficients)
}

# extract the coefficients
coefficients.semprobit <- function(object, ...) {
 UseMethod("coef", object)
}


# define summary method
# summary method for class "semprobit"
summary.semprobit <- function(object, var_names=NULL, file=NULL, digits = max(3, getOption("digits")-3), ...){
  # check for class "sarprobit"
  if (!inherits(object, "semprobit")) 
        stop("use only with \"semprobit\" objects")
        
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
  write(sprintf("--------MCMC probit with spatial errors ---------"), file, append=T)
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
        c("Estimate", "Std. Dev", "Bayes p-level", "t-value", "Pr(>|z|)"))
  printCoefmat(coefficients, digits = digits,
    signif.stars = getOption("show.signif.stars"))      
#  if (getOption("show.signif.stars")) {               
#    # The solution: using cat() instead of print() and use line breaks
#    # cat(paste(strwrap(x, width = 70), collapse = "\\\\\n"), "\n")
#    # http://r.789695.n4.nabble.com/Sweave-line-breaks-td2307755.html
#    Signif <- symnum(1e-6, corr = FALSE, na = FALSE,
#                  cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
#                  symbols = c("***", "**", "*", ".", " "))
#    x <- paste("Signif. codes: ", attr(Signif, "legend"), "\n", sep="")
#    cat(paste(strwrap(x, width = getOption("width")), collapse = "\\\n"), "\n")
#  }
   return(invisible(coefficients))
} 

# plot MCMC results for class "semprobit" (draft version);
# diagnostic plots for results (trace plots, ACF, posterior density function)
# method is very similar to plot.lm()
#
# @param x
# @param which
# @param ask
# @param trueparam a vector of "true" parameter values to be marked in posterior density plot
plot.semprobit <- function(x, which=c(1, 2, 3), 
  ask = prod(par("mfcol")) < length(which) && dev.interactive(), ..., trueparam=NULL) {
 if (!inherits(x, "semprobit")) 
        stop("use only with \"semprobit\" objects")
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

# Extract Log-Likelihood; see logLik.glm() for comparison
# Method returns object of class "logLik" with at least one attribute "df"
# giving the number of (estimated) parameters in the model.
# see Marsh (2000) equation (2.8), p.27 
logLik.semprobit <- function(object, ...) {
  X <- object$X
  y <- object$y
  n <- nrow(X)
  k <- ncol(X)
  W <- object$W
  beta <- object$beta
  rho <- object$rho
  sige <- object$sige
  I_n <- sparseMatrix(i=1:n, j=1:n, x=1)
  S <- I_n - rho * W
  D <- diag(1/sqrt(sige*diag(S %*% t(S))))  # D = diag(E[u u'])^{1/2}  (n x n)
  Xs <- D %*% X                             # X^{*} = D %*% X
  F <- pnorm(as.double(Xs %*% beta))        # F(X^{*} beta)  # (n x 1)
  lnL <- sum(log(F[y == 1])) + sum(log((1 - F[y == 0]))) # see Marsh (2000), equation (2.8)
  out <- lnL
  class(out) <- "logLik"
  attr(out,"df") <- k+2                     # k parameters in beta, rho, sige
 return(out)
}

# return fitted values of SEM probit (on reponse scale vs. linear predictor scale)
fitted.semprobit <- function(object, ...) {
  object$fitted.value
}

