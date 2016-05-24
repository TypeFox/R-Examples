# SAR Ordered Probit / Ordered spatial probit model
#
# see LeSage (2009), section 10.2
# see Greene (2003), section 21.8 for non-spatial ordered probit
#
# model:
# (1)  z = rho * W  * z + X beta + eps
# (2) y_i can take J alternatives (j = 1,...,J) for
#     y_i = j, if phi_{j-1} <= z <= phi_j
# (3) phi is a vector
#
# SAR probit is special case with J=1 (2 alternatives) 
# and phi=c(-Inf, 0, Inf) is J+1 vector with phi[0] = -Inf, phi[1]=0 and phi[J]= +Inf
# Model paramaters to be estimated:
# beta, rho and vector phi (J-2 values: phi[2],...,phi[J-1])
#
# TO BE DONE
# 0. Is there an implementation of SAR Ordered Probit in LeSage Matlab Toolbox?
# 1. Does sigma_eps need to be estimated?
# 2. What is the difference for drawing rho if any?
# 3. Is there any change required for the marginal effects / impacts() method?
# 4. Summary method for class "sarprobit" anpassen? JA wegen phi und wegen #0 und #1 values (vs. table of y)
# 5. Guten Default Wert für phi überlegen...
# 6. Impacts() überlegen
#
# Bayesian estimation of the SAR Ordered Probit model
#
sarorderedprobit <- function(formula, W, data, subset, ...) {
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
  sar_ordered_probit_mcmc(y, X, W, ...)    
} 



################################################################################

sar_ordered_probit_mcmc <- function(y, X, W, ndraw=1000, burn.in=100, thinning=1, 
  prior=list(a1=1, a2=1, c=rep(0, ncol(X)), T=diag(ncol(X))*1e12, lflag = 0), 
  start=list(rho=0.75, beta=rep(0, ncol(X)), phi=c(-Inf, 0:(max(y)-1), Inf)),
  m=10, computeMarginalEffects=TRUE, showProgress=FALSE){  


# method parameters
#ndraw=1000
#burn.in=100
#thinning=1
#prior=list(a1=1, a2=1, c=rep(0, ncol(X)), T=diag(ncol(X))*1e12, lflag = 0)
#start=list(rho=0.75, beta=rep(0, ncol(X)), phi=c(-Inf, 0, 0.5, 2.5, Inf))
#m=10                            # Anzahl MCMC samples for truncated normal variates
#computeMarginalEffects=TRUE
#showProgress=TRUE

  #start timer
  timet <- Sys.time()
  
  n  <- nrow( X )            # number of observations
  n1 <- nrow( X )          
  n2 <- nrow( W )
  k <- ncol( X )             # number of of parameters/exogenous variables
  I_n <- sparseMatrix(i=1:n, j=1:n, x=1) # sparse identity matrix
  if (is.null(colnames(X))) colnames(X) <- paste("x",1:k,sep="")


# validate inputs
# (1) check y in 1...J
  # J <- 4
  J <- max(y)
  if( sum(y %in% 1:J) != length( y ) ){
    stop('sarorderedprobit: not all y-values are in 1...J')
  }
  if( n1 != n2 && n1 != n ){
    stop('sarorderedprobit: wrong size of spatial weight matrix W')
  }
  # check if spatial weights matrix W does not contain zeros in the main diagonal
  if (!inherits(W, "sparseMatrix") || any(diag(W) != 0)) {
    stop('sarorderedprobit: spatial weights matrix W must be a sparse matrix with zeros in the main diagonal')
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
    stop('sarorderedprobit: intercept term must be in first column of the X-matrix')
  }
  
  
# MCMC start values
rho  <- start$rho          # start value of rho
if (is.null(rho)) {
  rho <- 0.75
  warning("No start value set for rho. Setting rho=0.75")
}
beta <- start$beta         # start value of parameters, prior value, we could also sample from beta ~ N(c, T)
if (is.null(beta)) {
  beta <- rep(0, ncol(X))
  warning(sprintf("No start value set for beta. Setting beta=%s", paste(beta, collapse=",")))
}
phi  <- start$phi
if (is.null(phi)) {
  phi <- c(-Inf, 0:(J-1), +Inf)
}

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

S <- I_n - rho * W
H <- t(S) %*% S
QR <- qr(S)                # class "sparseQR"
mu <- solve(QR, X %*% beta)

# progress bar
if (showProgress) {
  pb <- txtProgressBar(min=0, max=(thinning * ndraw + burn.in), initial=0, style=3)
}

# immutable matrices
tX <- t(X)                    # X'               # k x n
xpx  <- t(X) %*% X            # (X'X)            # k x k
xpxI <- solve(xpx)            # (X'X)^{-1}       # k x k
xxpxI <- X %*% xpxI           # X(X'X)^(-1)     # n x k (better, compromise)
AA    <- solve(xpx + Tinv)    # (X'X + T^{-1})^{-1}

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


# MCMC parameters
params <- (k+1) + (J-1)            # parameters beta (k), rho (1), (J-1) cut parameters phi, aber nur (J-2) zu schätzen
# matrix to store the beta + rho parameters for each iteration/draw
B <- matrix(NA, ndraw, params)
colnames(B) <- c(paste("beta_", 1:k, sep=""), "rho", paste("y>=", 2:J, sep=""))

# just to set a start value for z
z <- rep(0, n)
ones <- rep(1, n)

# MCMC loop  
for (i in (1 - burn.in):(ndraw * thinning)) {
  
  # 1. sample from z | rho, beta, y using precision matrix H
  # mu will be updated after drawing from rho
  mu <- solve(qr(S), X %*% beta)
  
  # determine lower and upper bounds for z depending on the value of y and phi
  # eqn (10.14), p.299
  # besser: y = 1..J
  lower <- phi[y]
  upper <- phi[y+1]
  
  # see cbind(y, lower, upper)

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
  c <- AA  %*% (tX %*% S %*% z + Tinv %*% c)
  T <- AA   # no update basically on T, TODO: check this
  beta <- as.vector(rmvnorm(n=1, mean=c, sigma=T))

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
  # keep rho fixed for now
  #rho <- 0.75
  
  # 4. determine bounds/cut-points p(phi_j | phi_{-j}, z, y, beta) for j = 2,...,J-1
  # phi_j = 0 is set fixed!
  for (j in 2:(J-1)) {
    phi.lower <- max(max(z[y == j]),     phi[j-1+1])   # \bar{phi}_{j-1}, SW: +1 is needed as our vector index starts with 1
    phi.upper <- min(min(z[y == j + 1]), phi[j+1+1])   # \bar{phi}_{j+1}
    
    # Sample phi_{j | phi_{-j}, z, y, beta)
    phi[j + 1]   <- runif(n=1, min=phi.lower, max=phi.upper)
  }

  # update S and H
  S <- I_n - rho * W
  H <- t(S) %*% S      # H = S'S  / SW: crossprod(S) does not seem to work!
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
    # save parameters in this MCMC round
    B[ind,] <- c(beta, rho, phi[2:J])   # (k + 1) + (J - 1)
    #zmean   <- zmean + z
  }
  if (showProgress) setTxtProgressBar(pb, i + burn.in) # update progress bar
}
if (showProgress)  close(pb) #close progress bar
  
# fitted values for estimates (based on z rather than binary y like in fitted(glm.fit))
# (on response scale y vs. linear predictor scale z...)
beta  <- colMeans(B)[1:k]
rho   <- colMeans(B)[k+1]
#browser()
phi   <- c(-Inf, colMeans(B)[(k+2):((k+2) + (J-2))], Inf)

S     <- (I_n - rho * W)
fitted.values   <- solve(qr(S), X %*% beta) # z = (I_n - rho * W)^{-1}(X * beta)
fitted.response <- cut(as.double(fitted.values), breaks=phi, labels=FALSE, ordered_result = TRUE) # split according to phi values

#addMargins(table(y, fitted.response))
#   fitted.response
#y     1   2   3   4
#  1 153  34   9   0
#  2  19  20  16   4
#  3   4  16  27  13
#  4   5  10  23  47

# result
results       <- NULL
results$time  <- Sys.time() - timet
results$nobs  <- n          # number of observations
results$nvar  <- k          # number of explanatory variables
results$y     <- y 
results$beta  <- colMeans(B)[1:k]
results$rho   <- colMeans(B)[k+1]
results$phi   <- colMeans(B)[(k+2):((k+2) + (J-2))]
results$coefficients <- colMeans(B)
results$fitted.values <- fitted.values      # fitted latent values on linear predictor scale
results$fitted.response <- fitted.response  # fitted values on response scale (ordered y variable)
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
results$names     <- c(colnames(X), 'rho', paste("y>=", 2:J, sep=""))
results$B         <- B        # (beta, rho, phi) draws
results$bdraw     <- B[,1:k]  # beta draws
results$pdraw     <- B[,k+1]  # rho draws
results$phidraw   <- B[(k+2):((k+2) + (J-2))]
#results$total     <- total
#results$direct    <- direct
#results$indirect  <- indirect
results$W <- W
results$X <- X
#results$mlike     <- mlike    # log-likelihood based on posterior means

#results$predicted <- # prediction required. The default is on the scale of the linear predictors
class(results)    <- c("sarorderedprobit", "sarprobit")
return(results)
}

# summary method for class "sarorderedprobit"
summary.sarorderedprobit <- function(object, var_names=NULL, file=NULL, 
  digits = max(3, getOption("digits")-3), ...){
  # check for class "sarorderedprobit"
  if (!inherits(object, "sarorderedprobit")) 
        stop("use only with \"sarorderedprobit\" objects")
        
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
  write(sprintf("----MCMC spatial autoregressive ordered probit----"), file, append=T)
  #sprintf("Dependent Variable")
  write(sprintf("Execution time  = %6.3f %s", object$time, attr(object$time, "units"))  , file, append=T)
  write(sprintf("N steps for TMVN= %6d"  , object$nsteps), file, append=T)
  write(sprintf("N draws         = %6d, N omit (burn-in)= %6d", ndraw, nomit), file, append=T)
  write(sprintf("N observations  = %6d, K covariates    = %6d", nobs, nvar)  , file, append=T)
  write(sprintf("Min rho         = % 6.3f, Max rho         = % 6.3f", object$rmin, object$rmax), file, append=T)
  write(sprintf("--------------------------------------------------"), file, append=T)
  write(sprintf(""), file, append=T)
  #DISTRIBUTION OF y VALUES
  print(table(y=object$y))
  #ESTIMATION RESULTS
  coefficients <- cbind(bout_mean, bout_sd, bout_sig, bout_t, bout_tPval)
  dimnames(coefficients) <- list(bout_names, 
        c("Estimate", "Std. Dev", "p-level", "t-value", "Pr(>|z|)"))
  printCoefmat(coefficients, digits = digits,
    signif.stars = getOption("show.signif.stars"))      
  return(invisible(coefficients))
}

# return fitted values of SAR ordered probit (on response scale vs. linear predictor scale)
fitted.sarorderedprobit <- function(object, ...) {
  object$fitted.value
}




