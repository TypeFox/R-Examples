###########################################################
# check arguments 
###########################################################

check.args.cplm <- function(call,n.obs){
  ## checking arguments  
  if (!is.null(call$weights)){
    if (!is.numeric(call$weights)) 
        stop("'weights' must be a numeric vector")
    if (any(call$weights <= 0)) 
        stop("negative or zero weights not allowed")
  }
  if (!is.null(call$offset)) {
    if (length(call$offset) != n.obs) 
      stop(gettextf("number of 'offset' is %d should 
                    equal %d (number of observations)", 
                length(call$offset), n.obs), domain = NA)
    }
}

check.args.bcplm <- function(call, n.beta, n.chains){
  n.iter <- eval(call$n.iter)
  n.burnin <- eval(call$n.burnin)
  # check counts related inputs
  if (!is.null(call$n.chains) && (!is.numeric(call$n.chains) 
      || call$n.chains < 1))
    stop("'n.chains' must be greater than 1" )
  if (!is.null(n.burnin) && !is.null(n.iter) && 
    n.burnin >= n.iter)
  	stop("'n.burnin' should be less than 'n.iter'" )
  if (!is.null(call$prior.beta.mean) && 
  	length(call$prior.beta.mean) != n.beta)
  	stop(gettextf("'prior.beta.mean' should be of length %d"), n.beta)  
  if (!is.null(call$prior.beta.mean) && 
  	length(call$prior.beta.mean) != n.beta)
  	stop(gettextf("'prior.beta.mean' should be of length %d"), n.beta)
}

###########################################################
# Check initial values
###########################################################

# check initial values in cpglm
check.inits.cpglm <- function(inits, n.beta){
  
  if (any(is.na(match(c("beta", "phi", "p"), names(inits)))))
    stop("'inits' must contain 'beta', 'phi' and 'p'!")
  if (length(inits$beta) != n.beta)
    stop(gettextf("number of 'beta' in 'inits' is %d, but should 
                    equal %d (number of mean parameters)", 
                length(inits$beta), n.beta, domain = NA))
    
  if (length(inits$phi) > 1 || inits$phi <= 0) 
    stop("'phi' in 'inits' should be of length 1 and greater than 0")
  if (length(inits$p) > 1 || inits$p <= 1 || inits$p >= 2) 
    stop("'p' in 'inits' should be of length 1 and between 1 and 2")
}

# check initial values in cpglmm
check.inits.cpglmm <- function(inits, n.beta, n.term){
  check.inits.cpglm(inits, n.beta)
  if (!("Sigma" %in% names(inits)))
    stop("the 'Sigma' component in 'inits' is missing") 
  if (length(inits$Sigma) != n.term) 
    stop(gettextf("'Sigma' in 'inits' should be of length %d", n.term))
}

# check initial values in bcplm
check.inits.bcplm <- function(inits, n.beta, n.term, n.chains, is.cpglmm){
  if (length(inits) != n.chains)
    stop(gettextf("'inits' should be of length %d", n.chains))  
  lapply(inits, function(x) {
    if (any(is.na(match(c("beta", "phi", "p"), names(x)))))
      stop("elements of 'inits' should contain 'beta', 'phi' and 'p'")
    if (is.cpglmm && any(is.na(match(c("u", "Sigma"), names(inits)))))
      stop("elements of 'inits' should contain 'u', and 'Sigma'")  
    if (!is.cpglmm){
      check.inits.cpglm(x, n.beta)
    } else {
      check.inits.cpglmm(x, n.beta, n.term)
      if (!("u" %in% names(x)))
        stop("the 'u' component in 'inits' is missing")
    }
  })
}  


###########################################################
# default control options   
###########################################################

# set control parameters  
cplm.control <- function(max.iter = 300L,
                       max.fun = 2000L,               
                       bound.p = c(1.01, 1.99),
                       trace = 0,
                       PQL.init = TRUE){         
  if (!is.numeric(max.iter) || max.iter <= 0) 
        stop("value of 'max.iter' must be > 0")
  if (!is.numeric(max.fun) || max.fun <= 0) 
        stop("value of 'max.fun' must be > 0")
  if (!is.numeric(bound.p) || length(bound.p) != 2)
        stop("'bound.p' must be of length 2")
  if (min(bound.p) < 1 || max(bound.p) > 2)
        stop("invalid bounds in 'bound.p'")          
  if (!is.numeric(trace) && !is.logical(trace))
        stop("'trace' must be logical or numeric")
  
  list(max.iter = as.integer(max.iter),
       max.fun = as.integer(max.fun),
       bound.p = as.numeric(sort(bound.p)),
       trace = as.integer(trace),
       PQL.init = as.logical(PQL.init))
}


###########################################################
# numerical derivatives  
###########################################################

# function to compute gradient
grad <- function(parm, fun, ...){
  n <- length(parm)
  eps <- 0.001
  gd <- rep(NA, n)
  for (i in 1:n){
    parm[i] <- parm[i] - eps
    g1 <- fun(parm, ...)
    parm[i] <- parm[i] + 2 * eps
    g2 <- fun(parm, ...)
    gd[i] <- (g2 - g1) / (2 * eps)
    parm[i] <- parm[i] - eps
  }
  return(gd)
}

# function to compute hessian
hess <- function(parm, fun, ...){
  n <- length(parm)
  eps <- 0.001
  hn <- matrix(0, n, n)
  for (i in 1:n){
    parm[i] <- parm[i] - eps
    g1 <- grad(parm, fun, ...)
    parm[i] <- parm[i] + 2 * eps
    g2 <- grad(parm, fun, ...)
    hn[i,] <- (g2 - g1) / ( 2 * eps)
    parm[i] <- parm[i] - eps
  }
  return(hn)  
}


###########################################################
# glm related   
###########################################################

# construct model frame in cpglm   
cpglm.mf <- function(mf, contrasts){  
    m <- match(c("formula", "data", "subset", "weights",
                 "na.action", "offset"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame(2))
    mt <- attr(mf, "terms")
    Y <- model.response(mf, "any")
    X <- if (!is.empty.model(mt)) 
          model.matrix(mt, mf, contrasts)
    weights <- as.vector(model.weights(mf))
    offset <- as.vector(model.offset(mf))
    n.obs <- nrow(X)
    if (is.null(weights))
      weights <- rep(1, n.obs)    
    if (is.null(offset))
      offset <- rep(0, n.obs)
    return (list(mf = mf, wts = weights, off = offset,
                 Y = Y, X = X))
  }

# fit a Tweedie glm given a model frame 
cpglm.fit <- function(fr, p = 1.5, link.power = 0) {
  fm <- tweedie(var.power = p, link.power = link.power)
  int <- attr(attr(fr$mf,"terms"), "intercept") > 0L
  suppressWarnings(glm.fit(fr$X, fr$Y, weights = fr$wts, offset = fr$off,
                      family = fm, intercept = int))
}
  
# generate inital values for a Tweedie glm given a model frame
cpglm.init <- function(fr, link.power = 0){
  p <- 1.5
  fit <- cpglm.fit(fr, p, link.power)
  beta <- as.numeric(fit$coefficients)
  phi <- sum(fit$weights * fit$residuals^2) / fit$df.residual
  vbeta <- summary.glm(fit)$cov.scaled
  list(beta = beta, phi = phi, p = p, vcov = vbeta)
}

# generate inital values for bcplm
bcplm.init <- function(fr, link.power = 0, n.chains, bound.p, dm){
  init <- cpglm.init(fr, link.power)
  n.beta <- length(init$beta)
  bound.p[1] <- max(bound.p[1], 1.4)
  bound.p[2] <- min(bound.p[2], 1.6)
  init0 <- unname(unlist(init[1:3]))
  inits <- vector("list", n.chains)
  inits[[1]] <- init0
  if (n.chains > 1){
    for (i in 2:n.chains)
    inits[[i]] <- c(as.numeric(init$beta + rnorm(n.beta, 0, 0.5)),
                    runif(1, 0.8 * init$phi, 1.2 * init$phi),
                    runif(1, bound.p[1], bound.p[2]))
  }
  if (!is.null(dm)){
    s <- lapply(dm$ST, function(x) x %*% t(x))
    sv <- unlist(lapply(s, as.numeric))
    inits <- lapply(inits, function(x) c(x, rnorm(dm$dd[["q"]]), sv))
  }  
  return(list(inits = inits, vbeta = init$vcov))
}             
  
# compute fitted values of for bigglm  
fitted.bigglm <- function(object, data, ...){
  # get chunks of data
  tt <- terms(object)
  n <- object$n
  beta <- coef(object)
  cursor <- 0
  eta <- offset <- pwts <- c()
  datafun <- function(){
    if (cursor >= n)
        return(NULL)
    start <- cursor + 1
    cursor <<- cursor + min(object$call$chunksize, n - cursor)
    data[start:cursor, ]
  }
  # get stats for each chunk
  while(!is.null(chunk <- datafun())){
    mf <- model.frame(tt, chunk)
    mm <- model.matrix(tt, mf)
    if(is.null(off <- model.offset(mf))) 
      off <- rep(0, nrow(mm))  
    if (!is.null(object$weights))
      w <- model.frame(object$weights, chunk)[[1]] else 
      w <- rep(1, nrow(mm))
    eta <- c(eta, mm %*% beta + off)    
    offset <- c(offset, off)
    pwts <- c(pwts, w)
  }
  # compute stats to be returned
  mu <- object$family$linkinv(eta)
  dmu <- object$family$mu.eta(eta)
  wts <- pwts * dmu * dmu / (object$family$variance(mu))
  y <- eval(object$call$formula[[2]], data)
  res <- (y - mu) / dmu
  list(linear.predictors = eta, 
        fitted.values = mu,
        offset = offset,
        prior.weights = pwts,
        weights = wts,
        residuals = res )
}


###########################################################
# general utility functions  
###########################################################

# function to compute minus twice log density 
dtweedie.nlogl <- function(y, mu, phi, p) {
    -2 * sum(log(dtweedie(y = y, mu = mu, phi = phi, power = p)))
}

# function to take inverse of a matrix using svd 
svd.inv <- function(x){
  sx <- svd(x)
  return(sx$v %*% diag(1 / sx$d) %*% t(sx$u))	
}
    
# function to compute the link.power needed in tweedie
make.link.power <- function(link) {
  if (!is.character(link) && !is.numeric(link))
    stop("link.power must be either numeric or character.")
  if (is.character(link)){  
    okLinks <- c("log", "identity", "sqrt","inverse")
    if (link %in% okLinks) 
      switch(link, log = 0, identity = 1, sqrt = 0.5, inverse = -1) else
      stop("invalid link function!")
  } else 
    link  
}

# optimize an objective function using different optimizers
cplm_optim <- function(par, fn, gr = NULL, ..., 
                       lower = -Inf, upper = Inf, control = cplm.control(), 
                       optimizer = "nlminb"){
  optimizer <- match.arg(optimizer, c("nlminb", "L-BFGS-B", "bobyqa"))
  if (optimizer == "nlminb"){
    ans <- nlminb(par, fn, gradient = gr, ..., 
                  lower = lower, upper = upper,                
                  control = list(trace = control$trace,
                                 iter.max = control$max.iter,
                                 eval.max = control$max.fun))
    names(ans)[2] <- "value"
    return(ans[c("par", "value", "convergence", "message")])    
  } else if (optimizer == "L-BFGS-B"){
    ans <- optim(par, fn, gr = gr, ..., method = "L-BFGS-B",                     
                 lower = lower, upper = upper,
                 control = list(trace = control$trace,
                                maxit = control$max.iter))
    return(ans[c("par", "value", "convergence", "message")])
  } else if (optimizer == "bobyqa"){ 
    ans <- bobyqa(par, fn, lower = lower, upper = upper, 
                  control = list(iprint = control$trace,
                                 rhobe = 0.02, rhoend = 2e-7,
                                 maxfun = control$max.fun),
                  ...)
    names(ans)[c(2, 4, 5)] <- c("value", "convergence", "message")
    return(ans[c("par", "value", "convergence", "message")])    
  }
}


# 1-d random walk metropolis
metrop_rw <- function(n = 1, par = 0, sd = 1, fun, ..., 
                      lower = -Inf, upper = Inf){
  fn <- function(x) fun(x, ...)
  if (par < lower || par > upper)
    stop("starting value of x is invalid")
  # run the metropolis algorithm
  .Call("bcplm_metrop_rw", as.integer(n), as.double(par), as.double(sd), 
        as.double(lower), as.double(upper), fn, new.env())
}
