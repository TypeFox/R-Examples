## functions.R
## David Ardia

## bayesGARCH
## __input__
## y : time series
## mu.alpha = c(0,0) : prior mean for alpha
## v.alpha = 1000 * diag(1,2) : prior covariance matrix for alpha
## mu.beta = 0 : prior mean for beta
## v.beta = 1000 : prior variance for beta
## c.nu = 0.01
## d.nu = 2 : prior mean of 102
## control = list() : a list with
##           n.chain = 1: number of chain
##           l.chain = 1000 : length of each chain
##           start.val = c(0.01, 0.1, 0.7, 20) : starting valuesn
##           addPriorConditions : a function
##           digits = 4 : digits for the output
##           refresh = 10 : refreshing frequency
## __output__
## r : a list of the class MCMC.list of length n.chain containing the MCMC outputs
"bayesGARCH" <- function(y, mu.alpha = c(0,0), Sigma.alpha = 1000 * diag(1,2),
                         mu.beta = 0, Sigma.beta = 1000,
                         lambda = 0.01, delta = 2, control = list()){
  if (missing(y))
    stop ("'y' is missing")
  if (!is.vector(y))
    stop ("'y' must be a vector")
  if (length(y) < 2)
    stop ("'y' must be a longer time series")
  if (any(is.na(y)))
    stop ("'y' contains 'NA' values")
  if (!is.vector(mu.alpha))
    stop ("'mu.alpha' must be a vector")
  if (length(mu.alpha) != 2) 
    stop ("'mu.alpha' is not of appropriate size")
  if (!is.matrix(Sigma.alpha))
    stop ("'Sigma.alpha' must be a matrix")
  if (!all(dim(Sigma.alpha) == c(2, 2))) 
    stop ("'Sigma.alpha' is not of appropriate size")
  if (Sigma.alpha[1,2] != Sigma.alpha[2,1])
    stop ("'Sigma.alpha' is not symmetric")
  if (any(eigen(Sigma.alpha)$values <= 0))
    stop ("'Sigma.alpha' is not positive definite")
  if (!is.vector(mu.beta) || length(mu.beta) != 1)
    stop ("'mu.beta' must be a scalar")
  if (!is.vector(Sigma.beta) || length(Sigma.beta) != 1)
    stop ("'Sigma.beta' must be a scalar")
  if (Sigma.beta <= 0)
    stop ("'Sigma.beta' must be positive")
  if (!is.vector(lambda) || length(lambda) != 1)
    stop ("'lambda' must be a scalar")
  if (lambda <= 0)
    stop ("'lambda' must be positive")
  if (!is.vector(delta) || length(delta) != 1)
    stop ("'delta' must be a scalar")
  if (delta < 2)
    stop ("'delta' cannot be lower than '2'")
  if (!is.list(control))
    stop ("'control' must be a list")
  if (!is.null(control$addPriorConditions))
    if (!is.function(control$addPriorConditions))
      stop ("'control$addPriorConditions' must be NULL or a function")
  
  fn.bayesGARCH(y, mu.alpha, solve(Sigma.alpha), mu.beta, 1/Sigma.beta, lambda, delta, control)
}

## sub function fn.bayesGARCH
"fn.bayesGARCH" <- function(y, mu.alpha, iv.alpha, mu.beta, iv.beta, c.nu, d.nu, control){

  con <- list(n.chain = NULL, l.chain = 10000, addPriorConditions = NULL, digits = 4, refresh = 10,
              start.val = matrix(c(0.01,0.1, 0.7, 100), 1, 4,
                dimnames = list("chain1", c("alpha0", "alpha1", "beta", "nu"))),
              hypers = list(mu.alpha = mu.alpha, iv.alpha = iv.alpha,
                mu.beta = mu.beta, iv.beta = iv.beta, c.nu = c.nu, d.nu = d.nu))
  
  con[names(control)] <- control
  if (is.vector(con$start.val)) 
    con$start.val <- matrix(con$start.val, nrow = 1, byrow = TRUE)
  if (is.null(con$n.chain)) 
    con$n.chain <- nrow(con$start.val)
  else con$start.val <- matrix(rep(con$start.val, con$n.chain), 
                               con$n.chain, byrow = TRUE)
  if (is.null(colnames(con$start.val))) 
    colnames(con$start.val) <- c("alpha0", "alpha1", "beta", "nu") 
  if (is.null(rownames(con$start.val))) 
    rownames(con$start.val) <- paste(rep("chain", con$n.chain), 1:con$n.chain, sep = "")
  if (is.null(con$addPriorConditions)) {con$addPriorConditions <- function(psi){TRUE}}
  r <- list()
  for (i in 1:con$n.chain) {
    k <- NULL
    chain <- fn.initializeChain(con$l.chain, con$start.val[i,])
    for (j in 2:con$l.chain) {
      k <- j - 1
      chain[j, ] <- fn.block(y, chain[k, 1:2], chain[k,3], chain[k, 4], con$hypers$mu.alpha, con$hypers$iv.alpha, 
                             con$hypers$mu.beta, con$hypers$iv.beta, con$hypers$c.nu, 
                             con$hypers$d.nu, con$addPriorConditions)
      if (con$refresh > 0 & j%%con$refresh == 0) 
        cat("chain: ", i, " iteration: ", j, " parameters: ", 
            round(chain[j, ], con$digits), "\n")
    }
    r[[i]] <- chain
  }
  names(r) <- paste(rep("chain", con$n.chain), 1:con$n.chain, sep = "")
  attr(r, "class") <- "mcmc.list"
  r
}

## cvGARCH
## __input__
## y : Tx1 time series
## theta : vector of scedastic's function parameter, theta := (alpha, beta)
## CSC : covariance statinarity condition
## __ouput__
## h : (T+1)x1 vector of conditional variances
"cvGARCH" <- function(y, theta, CSC = FALSE){
  if (missing(y))
    stop ("'y' is missing")
  if (!is.vector(y))
    stop ("'y' should be a vector")
  if (length(y) < 2)
    stop ("'y' sould be a longer time series")
  if (any(is.na(y)))
    stop ("'y' contains 'NA' values")
  if (missing(theta))
    stop ("'theta' is missing")
  if (!is.vector(theta))
    stop ("'theta' should be a vector")
  if (length(theta) != 3)
    stop ("not appropriate size for 'theta'")
  if (any(theta <= 0))
    stop ("some components of 'theta' are not positive")
  if (CSC && theta[2] + theta[3] >=1)
    warning ("the covariance stationarity condition 'alpha1 + beta < 1' is violated")
  
  fn.cvGARCH(y, theta)
}

"fn.cvGARCH" <- function(y, theta){
  n <- length(y)
  .C('fnGarchC',
     n = as.integer(n),
     order = as.integer(c(1,2)),
     alpha = as.double(theta[1:2]),
     delta = as.double(0), 
     beta = as.double(theta[3]),
     h = vector('double', n+1),
     u = as.double(y),
     sim = vector('integer',n+1),
     PACKAGE = 'bayesGARCH')$h
}

## fn.initializeChain
## __input__
## l.chain : length of the chain
## start.val : starting values
## __output__
## mcmc matrix 
"fn.initializeChain" <- function (l.chain, start.val){
  r <- matrix(NA, l.chain, length(start.val), dimnames = list(1:l.chain, names(start.val)))
  r[1,] <- start.val
  
  coda::mcmc(r, start = 1)
}

## fn.block
## __input__
## y : data
## alpha, beta, nu : parameters
## alpha0, iv.alpha0, beta0, iv.beta0, c.nu, d.nu : hyperparameters
## __output__
## block of parameters
"fn.block" <- function(y, alpha, beta, nu, alpha0, iv.alpha0, beta0, iv.beta0, c.nu, d.nu, addPriorConditions){
  
  w.new <- fn.w.full(y, alpha, beta, nu,
                     alpha0, iv.alpha0, beta0, iv.beta0, c.nu, d.nu)	
  
  alpha.new <- fn.alpha.full(y, alpha, beta, nu, w.new,
                             alpha0, iv.alpha0, beta0, iv.beta0, c.nu, d.nu)

  beta.new <- fn.beta.full(y, alpha.new, beta, nu, w.new,
                           alpha0, iv.alpha0, beta0, iv.beta0, c.nu, d.nu)

  nu.new <- fn.nu.full(y, alpha.new, beta.new, nu, w.new,
                       alpha0, iv.alpha0, beta0, iv.beta0, c.nu, d.nu)

  psi.new <- c(alpha.new, beta.new, nu.new)
  if (addPriorConditions(psi.new)) psi.new else c(alpha,beta,nu)
}

## fn.alpha.full
## __input__
## y : data
## alpha.old, beta, nu, w : parameters
## alpha0, iv.alpha0, beta0, iv.beta0, c.nu, d.nu : hyperparameters
## __output__
## MH draw for alpha
"fn.alpha.full" <- function(y, alpha.old, beta, nu, w,
                            alpha0, iv.alpha0, beta0, iv.beta0, c.nu, d.nu){
  r <- alpha.old
  n <- length(y)
  h <- fn.cvGARCH(y, c(alpha.old, beta))[1:n]
  
  post.old <- fn.post.garch(y, w*h, alpha.old, beta, nu,
                            alpha0, iv.alpha0, beta0, iv.beta0, c.nu, d.nu)

  tau <- 2 * h^2
  filter.alpha <- fn.filterAlpha(y, beta)

  Dd.D <- fn.Dd.D(y^2/w, filter.alpha, tau,
                  alpha0, iv.alpha0)
  
  alpha.new <- t( mvtnorm::rmvnorm(n = 1, mean = Dd.D$Dd, sigma = Dd.D$D) )

  if (all(alpha.new>0)){
    h.new <- fn.cvGARCH(y, c(alpha.new, beta))[1:n]

    tau.new <- 2 * h.new^2
		
    post.new <- fn.post.garch(y, w*h.new, alpha.new, beta, nu,
                              alpha0, iv.alpha0, beta0, iv.beta0, c.nu, d.nu)

    prop.new <- mvtnorm::dmvnorm(t(alpha.new), mean = Dd.D$Dd, sigma = Dd.D$D, log = TRUE)
		
    Dd.D <- fn.Dd.D(y^2/w, filter.alpha, tau.new,
                    alpha0, iv.alpha0)
		
    prop.old <- dmvnorm(alpha.old, mean = Dd.D$Dd, sigma = Dd.D$D, log = TRUE)

    ratio <- post.new - post.old + prop.old - prop.new

    r.u <- runif(1)
    if (r.u < min(1, exp(ratio)))
      r <- alpha.new
  }
  
  r
}

## fn.beta.full
## __input__
## y : data set
## alpha, beta.old, nu, w : parameters
## alpha0, iv.alpha0, beta0, iv.beta0, c.nu, d.nu : hyperparameters
## __output__
## MH draw for beta
"fn.beta.full" <- function(y, alpha, beta.old, nu, w,
                           alpha0, iv.alpha0, beta0, iv.beta0, c.nu, d.nu){
  n <- length(y)
  r <- beta.old
  h <- fn.cvGARCH(y, c(alpha, beta.old))[1:n]

  post.old <- fn.post.garch(y, w*h, alpha, beta.old, nu,
                            alpha0, iv.alpha0, beta0, iv.beta0, c.nu, d.nu)
  tau <- 2 * h^2
  beta.nls <- beta.old

  w1 <- fn.filterW(y, alpha, beta.nls)
  W <- fn.W(y^2-w1, -beta.nls)
  z <- w1 + W * beta.nls
  z <- z + y^2 * (1/w - 1)
  
  Dd.D <- fn.Dd.D(z, W, tau, beta0, iv.beta0)

  beta.new <- rnorm(n = 1, mean = Dd.D$Dd, sd = sqrt(Dd.D$D))
	
  if (beta.new>0){
    h.new <- fn.cvGARCH(y, c(alpha, beta.new))[1:n]
    tau.new <- 2 * h.new^2
    
    post.new <- fn.post.garch(y, w*h.new, alpha, beta.new, nu, 
                              alpha0, iv.alpha0, beta0, iv.beta0, c.nu, d.nu)
				
    prop.new <- dnorm(beta.new, mean = Dd.D$Dd, sd = sqrt(Dd.D$D), log = TRUE)
		
    w1 <- fn.filterW(y, alpha, beta.new)
    W <- fn.W(y^2-w1, -beta.new)
    z <- w1 + W * beta.new
    z <- z + y^2 * (1/w - 1)
    
    Dd.D <- fn.Dd.D(z, W, tau.new, beta0, iv.beta0)
		
    prop.old <- dnorm(beta.old, mean = Dd.D$Dd, sd = sqrt(Dd.D$D), log = TRUE)

    ratio <- post.new - post.old + prop.old - prop.new

    r.u <- runif(1)
    if (r.u < min(1, exp(ratio)))
      r <- beta.new  
  }
  
  r  	
}

## fn.w.full
## __input__
## u : data set
## alpha, beta, nu : parameters
## alpha0, iv.alpha0, beta0, iv.beta0, c.nu, d.nu : hyperparameters
## __output__
## draw for w
"fn.w.full" <- function(u, alpha, beta, nu,
                        alpha0, iv.alpha0, beta0, iv.beta0, c.nu, d.nu){
  n <- length(u)
  h <- fn.cvGARCH(u, c(alpha, beta))[1:n]

  r <- rgamma(n = n, shape = 0.5 * rep((nu+1),n),
              rate = 0.5 * ((u^2/h) + (nu-2)) )
  
  1/r
}

## fn.nu.full
## __input__
## u : data set
## alpha, beta, nu.old, w : parameters
## alpha0, iv.alpha0, beta0, iv.beta0, c.nu, d.nu : hyperparameters
## alpha.min, alpha.max, k.max : simulation parameters
## __output__
## simulated nu
"fn.nu.full" <- function(u, alpha, beta, nu.old, w,
                         alpha0, iv.alpha0, beta0, iv.beta0, c.nu, d.nu,
                         alpha.min = 0.00001, alpha.max = 500, k.max = 50000){
  n <- length(u)
  r <- nu.old
		
  phi <- 0.5 * sum(log(w) + 1/w) + c.nu
  
  optim.alpha <- uniroot(fn.neg.alpha,
                         interval = c(alpha.min, alpha.max),
                         tol = .Machine$double.eps,
                         maxiter = 1000,
                         n = n,
                         d = d.nu,
                         phi = phi)
  
  alpha <- optim.alpha$root
  value <- optim.alpha$f.root
  iter <- optim.alpha$iter
  if ( (alpha>alpha.min) && (alpha<alpha.max) ){
    k <- 0
    while (k<k.max){
      k <- k+1
      nu.new <- rexp(n = 1, rate = alpha)
      nu.new <- nu.new + d.nu
       
      p <- fn.accept.proba(n, alpha, d.nu, nu.new, phi)
    
      r.u <- runif(1)
      if(r.u < p){
        r = nu.new
        k = k.max
      }
    }
  }
  
  r
}

## fn.post.garch
## __input__
## u : data set
## h : conditional variance process
## alpha, beta, nu : parameters
## alpha0, iv.alpha0, beta0, iv.beta0, c.nu, d.nu : parameters
## __output__
## posterior density
"fn.post.garch" <- function(u, h, alpha, beta, nu,
                            alpha0, iv.alpha0, beta0, iv.beta0, c.nu, d.nu){

  r <- 0
  r <- -0.5 * sum(log(h))
  r <- r - 0.5 * sum(u^2/h)
  r <- r - 0.5 * t(alpha-alpha0) %*%
    iv.alpha0 %*% (alpha-alpha0)
  r <- r - 0.5 * (beta-beta0)^2 * iv.beta0
  r <- r - c.nu * nu
  
  r
}

## fn.neg.alpha
"fn.neg.alpha" <- function(alpha, n, d, phi){
  r <- 0
  r <- log( (1+alpha*(d-2))/(2*alpha) )
  r <- r + (1+alpha*d)/(1+alpha*(d-2))
  r <- r - digamma( (1+alpha*d)/(2*alpha) )
  r <- 0.5 * n * r + alpha - phi
  
  r
}

## fn.accept.proba
"fn.accept.proba" <- function(n, alpha, d, nu, phi){
  u <- (1+alpha*d) / (2*alpha)
  v <- nu/2
  w <- (1+alpha*(d-2)) / (2*alpha)
  
  r.ln <- 0
  r.ln <- n * (lgamma(u) - lgamma(v))
  r.ln <- r.ln + n * v * log((nu-2)/2)
  r.ln <- r.ln - n * u * log(w)
  r.ln <- r.ln + (nu - d) * (alpha - phi) + (phi/alpha) - 1
  
  exp(r.ln)
}

## fn.filterAlpha
"fn.filterAlpha" <- function(u, beta){
  n <- length(u)	
  filter.alpha <- .C('fnFilterAlphaC', 
                     n = as.integer(n),
                     u = as.double(u),
                     beta = as.double(beta),
                     vstar = vector('double', n),
                     lstar = vector('double', n),
                     PACKAGE = 'bayesGARCH')
  
  as.matrix(cbind(filter.alpha$lstar, 
                  c(0, filter.alpha$vstar[1:(n-1)])))
}

## fn.filterW
"fn.filterW" <- function(u, alpha, beta){
  n <- length(u)
  .C('fnFilterWC', 
     n = as.integer(n),
     u = as.double(u),
     alpha = as.double(alpha),
     beta = as.double(beta),
     w = vector('double', n),
     PACKAGE = 'bayesGARCH')$w
}

## fn.W
"fn.W" <- function(u, theta){
  n <- length(u)
  .C('fnQDiffC',
     n = as.integer(n),
     u = as.double(u),
     theta = as.double(theta),
     w = vector('double', n),
     w0 = as.double(0),
     PACKAGE = 'bayesGARCH')$w
}

## fn.Dd.D
## compute D and Dd
"fn.Dd.D" <- function(y, X, h, b0, iv.b0){
  r <- list()
  Xh <- X/h
  D.inv <- t(Xh) %*% X + iv.b0
  D <- solve(D.inv)
  d <- t(Xh) %*% y + iv.b0 %*% b0
  Dd <- D %*% d

  list(D = D, Dd = Dd)
}


