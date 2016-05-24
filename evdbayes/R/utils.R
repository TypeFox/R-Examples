## BETA AND GAMMA INFORMATION
"ibeta" <-
  function(mean, var, shape1, shape2)
{
  
  if(missing(shape1) && missing(shape2)) {
    if(min(mean) <= 0 || max(mean) >= 1)
      stop("`mean' must contain values in (0,1)")
    if(min(var) <= 0 || max(var) >= 0.25)
      stop("`var' must contain values in (0,0.25)")
    
    n <- max(length(mean), length(var))
    mean <- rep(mean, length.out=n)
    var <- rep(var, length.out=n)

    tmp <- var^(-1) * mean * (1-mean) - 1
    shape1 <- tmp * mean 
    shape2 <- tmp * (1 - mean)
    nind <- (tmp <= 0)
  }

  if(missing(mean) && missing(var)) {
    if(min(shape1) <= 0)
      stop("`shape1' must contain positive values")
    if(min(shape2) <= 0)
      stop("`shape2' must contain positive values")
    
    n <- max(length(shape1),length(shape2))
    shape1 <- rep(shape1, length.out=n)
    shape2 <- rep(shape2, length.out=n)
    
    mean <- shape1 / (shape1 + shape2)
    var <- mean * (1 - mean) / (shape1 + shape2 + 1)
    nind <- rep(FALSE, n)
  }
  
  mode <- ifelse(shape1 > 1 & shape2 > 1,
                 (shape1 - 1) / (shape1 + shape2 - 2), NA)
  rv <- matrix(c(shape1,shape2,mean,var,mode), ncol = 5,
               dimnames = list(1:n,c("shape1","shape2","mean","var","mode")))
  rv[nind, ] <- NA
  drop(rv)
}

"igamma" <-
  function(mean, var, shape, scale)
{
  if(missing(shape) && missing(scale)) {
    if(min(mean) <= 0)
      stop("`mean' must contain positive values")
    if(min(var) <= 0)
      stop("`var' must contain positive values")
    
    n <- max(length(mean),length(var))
    mean <- rep(mean, length.out=n)
    var <- rep(var, length.out=n)

    shape <- mean^2 / var
    scale <- var / mean
  }
  
  if(missing(mean) && missing(var)) {
    if(min(shape) <= 0)
      stop("`shape' must contain positive values")
    if(min(scale) <= 0)
      stop("`scale' must contain positive values")
    
    n <- max(length(shape),length(scale))
    shape <- rep(shape, length.out=n)
    scale <- rep(scale, length.out=n)
    
    mean <- shape * scale
    var <- shape * scale^2
  }
  mode <- ifelse(shape > 1, scale * (shape - 1), NA)
  rv <- matrix(c(shape,scale,mean,var,mode), ncol = 5,
               dimnames = list(1:n, c("shape","scale","mean","var","mode")))
  drop(rv)
}

## ACCEPT RATE COMPUTATIONS
ar.choice <- function(init, prior, lh = c("none","gev","gpd","pp","os"),
                      ..., psd, ar = rep(.4, npar), n = 1000,
                      tol = rep(.05, npar)){

  mc <- posterior(n = n, init = init, prior = prior, lh = lh,
                  ..., psd = psd)
  npar <- length(psd)
  ar.cur <- attributes(mc)$ar[1,1:npar]
  cat("Accept Rate values and proposal standard deviations at each iterations...\n")
  cat("Accept Rate", "\t", "Prop. Std\n")
  
  idx <- which(abs(ar.cur - ar) >= tol)
  cat(ar.cur, "\t", round(psd, 3), "\n")
  ##Initialization
  psd.inf <- rep(0, npar)
  psd.sup <- 5 * psd

  while( length(idx) > 0 ){
    
    psd.old <- psd
    
    for (i in idx){
      if ( (ar.cur - ar)[i] > 0 )
        psd[i] <- (psd.old[i] + psd.sup[i]) / 2
      
      else
        psd[i] <- (psd.inf[i] +  psd.old[i]) / 2
      
    }

    mc <- posterior(n = n, init = init, prior = prior, lh = lh,
                    ..., psd = psd)
    ar.cur <- attributes(mc)$ar[1,1:npar]
    cat(ar.cur, "\t", round(psd, 3), "\n" )
    
    for (i in idx){
      if ( (ar.cur - ar)[i] > 0 ){
        psd.inf[i] <- psd.old[i]
        psd.sup[i] <- psd.sup[i] * ar.cur[i] / ar[i]
      }
      
      else {
        psd.sup[i] <- psd.old[i]
        psd.inf[i] <- psd.inf[i] * ar.cur[i] / ar[i]
      }
    }

                                        #cat("psd Inf: ", psd.inf, "\n")
                                        #cat("psd Sup: ", psd.sup, "\n")
    idx <- which(abs(ar.cur - ar) >= tol)
  }
  
  return(list(psd = psd, ar = ar.cur))
}

## MAXIMIZE POSTERIOR DENSITY
"mposterior" <-
  function(init, prior, lh = c("none","gev","gpd","pp","os"),
           method = c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN"),
           lower = -Inf, upper = Inf, control = list(),
           hessian = FALSE, ...)
{
  if (!inherits(prior, "evprior")) 
    stop("`prior' must be a prior distribution")
  if(!prior$trendsd && (length(init) != 3 || mode(init) != "numeric"))
    stop("`init' must be a numeric vector of length three")
  if(prior$trendsd && (length(init) != 4 || mode(init) != "numeric"))
    stop("`init' must be a numeric vector of length four")
  if(init[2] <= 0)
    stop("initial value of scale parameter must be positive")

  lh <- match.arg(lh)
  ndpost <- function(par, prior, lh, ...) {
    if(par[2] <= 0) return(Inf)
    dprior <- do.call(prior$prior, c(list(par = par), prior[-1]))
    switch(lh,
           gev = -dprior - gevlik(par, ...),
           gpd = -dprior - gpdlik(par, ...),
           pp = -dprior - pplik(par, ...),
           os = -dprior - oslik(par, ...),
           none = -dprior)
  }

  inittest <- ndpost(par = init, prior = prior, lh = lh, ...) 
  if(is.infinite(inittest))
    stop("density is zero at initial parameter values")
  
  optim(init, ndpost, method = method, lower = lower,
        upper = upper, control = control, prior = prior,
        lh = lh, ...)
}
