spGLM <- function(formula, family="binomial", weights, data = parent.frame(), coords, knots,
                  starting, tuning, priors, cov.model, 
                  amcmc, n.samples, verbose=TRUE, n.report=100, ...){
  
  ####################################################
  ##Check for unused args
  ####################################################
  formal.args <- names(formals(sys.function(sys.parent())))
  elip.args <- names(list(...))
  for(i in elip.args){
    if(! i %in% formal.args)
      warning("'",i, "' is not an argument")
  }

  ####################################################
  ##formula
  ####################################################
  if(missing(formula)){stop("error: formula must be specified")}
  
  if(class(formula) == "formula"){
    
    holder <- parseFormula(formula, data)
    Y <- holder[[1]]
    X <- as.matrix(holder[[2]])
    x.names <- holder[[3]]

  }else{
    stop("error: formula is misspecified")
  }

  p <- ncol(X)
  n <- nrow(X)
    
  ##make sure storage mode is correct
  storage.mode(Y) <- "double"
  storage.mode(X) <- "double"
  storage.mode(p) <- "integer"
  storage.mode(n) <- "integer"

  ####################################################
  ##family and weights
  ####################################################
  if(!family %in% c("binomial","poisson"))
    stop("error: family must be binomial or poisson")

  if(missing(weights)){weights <- rep(1, n)}
  if(length(weights) != n){stop("error: weights vector is misspecified")}
  
  storage.mode(weights) <- "integer"
  
  ####################################################
  ##sampling method
  ####################################################
  n.batch <- 0
  batch.length <- 0
  accept.rate <- 0
    
  if(missing(amcmc)){
    
    if(missing(n.samples)){stop("error: n.samples need to be specified")}

    n.batch <- n.samples
    batch.length <- 1
    is.amcmc <- FALSE
    
  }else{
    
    is.amcmc <- TRUE

    names(amcmc) <- tolower(names(amcmc))
    
    if(!"n.batch" %in% names(amcmc)){stop("error: n.batch must be specified in amcmc list")}
    n.batch <- amcmc[["n.batch"]]
    
    if(!"batch.length" %in% names(amcmc)){stop("error: batch.length must be specified in amcmc list")}
    batch.length <- amcmc[["batch.length"]]

    if(!"accept.rate" %in% names(amcmc)){
      warning("accept.rate was not specified in the amcmc list and was therefore set to the default 0.43")
      accept.rate <- 0.43
    }else{
      accept.rate <- amcmc[["accept.rate"]]
    }

    n.samples <- n.batch*batch.length
  }

  storage.mode(is.amcmc) <- "integer"
  storage.mode(n.samples) <- "integer"
  storage.mode(n.batch) <- "integer"
  storage.mode(batch.length) <- "integer"
  storage.mode(accept.rate) <- "double"

  ####################################################
  ##Fit non-spatial model if specified
  ####################################################
  if(missing(coords)){

    ##Starting
    beta.starting <- 0
    
    if(missing(starting)){stop("error: starting value list for the parameters must be specified")}
    
    names(starting) <- tolower(names(starting))   
    
    if("beta" %in% names(starting)){
      beta.starting <- starting[["beta"]]
      if(length(beta.starting) != p){stop(paste("error: starting values for beta must be of length ",p,sep=""))}
    }else{
      if(family=="poisson"){
        beta.starting <- coefficients(glm(Y~X-1, family="poisson"))
      }else{
        beta.starting <- coefficients(glm((Y/weights)~X-1, weights=weights, family="binomial"))
      }
    }
    
    storage.mode(beta.starting) <- "double"
    
    ##Priors
    beta.Norm <- 0
    beta.prior <- "flat"
    
    if(missing(priors)) {stop("error: prior list for the parameters must be specified")}
    
    names(priors) <- tolower(names(priors))
    
    if("beta.normal" %in% names(priors)){
      beta.Norm <- priors[["beta.normal"]]
      if(!is.list(beta.Norm) || length(beta.Norm) != 2){stop("error: beta.Norm must be a list of length 2")}
      if(length(beta.Norm[[1]]) != p ){stop(paste("error: beta.Norm[[1]] must be a vector of length, ",p, " with elements corresponding to betas' mean",sep=""))}
      if(length(beta.Norm[[2]]) != p ){stop(paste("error: beta.Norm[[2]] must be a vector of length, ",p, " with elements corresponding to betas' sd",sep=""))}
      beta.prior <- "normal"
    }
    
    ##Tuning values
    beta.tuning <- 0
    
    if(!missing(tuning)){
      
      names(tuning) <- tolower(names(tuning))
      
      if(!"beta" %in% names(tuning)){stop("error: beta must be specified in tuning value list")}
      beta.tuning <- tuning[["beta"]]
      
      if(is.matrix(beta.tuning)){
        if(nrow(beta.tuning) != p || ncol(beta.tuning) != p){
          stop(paste("error: if beta tuning is a matrix, it must be of dimension ",p,sep=""))
        }
        
        if(is.amcmc){
          beta.tuning <- diag(beta.tuning)
        }
        
      }else if(is.vector(beta.tuning)){
        if(length(beta.tuning) != p){
          stop(paste("error: if beta tuning is a vector, it must be of length ",p,sep=""))
        }
        
        if(!is.amcmc && p > 1){
          beta.tuning <- diag(beta.tuning)
        }
        
      }else{
        stop("error: beta tuning is misspecified")
      }
      
    }else{##no tuning provided
      
      if(!is.amcmc){
        stop("error: tuning value list must be specified")
      }
      
      beta.tuning <- rep(0.01,p)
    }
    
    storage.mode(beta.tuning) <- "double"
    
    ##Other stuff
    storage.mode(n.report) <- "integer"
    storage.mode(verbose) <- "integer"
    
    ##Send if off
    out <- .Call("nonSpGLM_AMCMC", Y, X, p, n, family, weights,
                 beta.prior, beta.Norm, beta.starting, beta.tuning,
                 n.batch, batch.length, accept.rate, verbose, n.report, is.amcmc)
    
    out$p.beta.samples <- mcmc(t(out$p.beta.samples))
    
    colnames(out$p.beta.samples) <- x.names
    
    out$weights <- weights
    out$family <- family
    out$Y <- Y
    out$X <- X
    
    class(out) <- "nonSpGLM"
    
    return(out)
  }
    
  ####################################################
  ##Distance matrices
  ####################################################
  
  ####################
  ##Coords
  ####################
  if(missing(coords)){stop("error: coords must be specified")}
  if(!is.matrix(coords)){stop("error: coords must n-by-2 matrix of xy-coordinate locations")}
  if(ncol(coords) != 2 || nrow(coords) != n){
    stop("error: either the coords have more than two columns or then number of rows is different than
          data used in the model formula")
  }
  
  coords.D <- iDist(coords)
  storage.mode(coords.D) <- "double"
  
  ####################
  ##Knots
  ####################
   is.pp <- FALSE
  
  if(!missing(knots)){
    
    if(is.vector(knots) && length(knots) %in% c(2,3)){
      
      ##allow single knot dim
      if(knots[1] > 1){
        x.knots <- seq(min(coords[,1]), max(coords[,1]), length.out=knots[1])
      }else{
        x.knots <- (max(coords[,1])-min(coords[,1]))/2
      }
      
      if(knots[2] > 1){
        y.knots <- seq(min(coords[,2]), max(coords[,2]), length.out=knots[2])
      }else{
        y.knots <- (max(coords[,2])-min(coords[,2]))/2
      }
      
      ##if not single knot then adjust out half distance on all sides
      if(length(knots) == 2){
        if(knots[1] > 1){
          x.int <- (x.knots[2]-x.knots[1])/2
          x.knots <- seq(min(x.knots)-x.int, max(x.knots)+x.int, length.out=knots[1])
        }
        
        if(knots[2] > 1){
          y.int <- (y.knots[2]-y.knots[1])/2
          y.knots <- seq(min(y.knots)-y.int, max(y.knots)+y.int, length.out=knots[2])
        }
        
        knot.coords <- as.matrix(expand.grid(x.knots, y.knots))
        is.pp <- TRUE
      }else{   
        if(knots[1] > 1){
          x.int <- knots[3]
          x.knots <- seq(min(x.knots)-x.int, max(x.knots)+x.int, length.out=knots[1])
        }
        
        if(knots[2] > 1){
          y.int <- knots[3]
          y.knots <- seq(min(y.knots)-y.int, max(y.knots)+y.int, length.out=knots[2])
        }
        
        knot.coords <- as.matrix(expand.grid(x.knots, y.knots))
        is.pp <- TRUE
      }
      
    }else if(is.matrix(knots) && ncol(knots) == 2){
      knot.coords <- knots
      is.pp <- TRUE
    }else{
      stop("error: knots is misspecified")
    }
  }
 
  m <- 0
  knots.D <- 0
  knots.coords.D <- 0
  
  if(is.pp){
    knots.D <- iDist(knot.coords)
    m <- nrow(knots.D)
    knots.coords.D <- iDist(knot.coords, coords)
  }

  storage.mode(m) <- "integer"
  storage.mode(knots.D) <- "double"
  storage.mode(knots.coords.D) <- "double"
  
  ####################################################
  ##Covariance model
  ####################################################
  if(missing(cov.model)){stop("error: cov.model must be specified")}
  if(!cov.model%in%c("gaussian","exponential","matern","spherical"))
    {stop("error: specified cov.model '",cov.model,"' is not a valid option; choose, from gaussian, exponential, matern, spherical.")}

  ####################################################
  ##Starting values
  ####################################################
  beta.starting <- 0
  sigma.sq.starting <- 0
  phi.starting <- 0
  nu.starting <- 0
  w.starting <- 0
 
  if(missing(starting)){stop("error: starting value list for the parameters must be specified")}
  
  names(starting) <- tolower(names(starting))   

  if("beta" %in% names(starting)){
    beta.starting <- starting[["beta"]]
    if(length(beta.starting) != p){stop(paste("error: starting values for beta must be of length ",p,sep=""))}
  }else{
    if(family=="poisson"){
      beta.starting <- coefficients(glm(Y~X-1, family="poisson"))
    }else{
      beta.starting <- coefficients(glm((Y/weights)~X-1, weights=weights, family="binomial"))
    }
  }
  
  if(!"sigma.sq" %in% names(starting)){stop("error: sigma.sq must be specified in starting value list")}
  sigma.sq.starting <- starting[["sigma.sq"]][1]
  
  if(!"phi" %in% names(starting)){stop("error: phi must be specified in starting value list")}
  phi.starting <- starting[["phi"]][1]
  
  if(cov.model == "matern"){
    if(!"nu" %in% names(starting)){stop("error: nu must be specified in starting value list")}
    nu.starting <- starting[["nu"]][1]
  }

  if(!"w" %in% names(starting)){
    warning("w is not specified in starting value list. Setting starting value to 0.")
  }else{
    w.starting <- starting[["w"]]
  }
    
  if(is.pp){
    if(length(w.starting) == 1){w.starting <- rep(w.starting, m)}  
    if(length(w.starting) != m){stop(paste("error: w in the starting value list must be a scalar of length 1 or vector of length ",m," (i.e., the number of predictive process knots)",sep=""))}
  }else{
    if(length(w.starting) == 1){w.starting <- rep(w.starting, n)}
    if(length(w.starting) != n){stop(paste("error: w in the starting value list must be a scalar of length 1 or vector of length ",n," (i.e., the number of predictive process knots)",sep=""))}    
  }

  storage.mode(beta.starting) <- "double"
  storage.mode(phi.starting) <- "double"
  storage.mode(sigma.sq.starting) <- "double"
  storage.mode(nu.starting) <- "double"
  storage.mode(w.starting) <- "double"

  ####################################################
  ##Priors
  ####################################################
  beta.Norm <- 0
  beta.prior <- "flat"
  sigma.sq.IG <- 0
  nu.Unif <- 0
  phi.Unif <- 0

  if(missing(priors)) {stop("error: prior list for the parameters must be specified")}
    
  names(priors) <- tolower(names(priors))

  if("beta.normal" %in% names(priors)){
    beta.Norm <- priors[["beta.normal"]]
    if(!is.list(beta.Norm) || length(beta.Norm) != 2){stop("error: beta.Norm must be a list of length 2")}
    if(length(beta.Norm[[1]]) != p ){stop(paste("error: beta.Norm[[1]] must be a vector of length, ",p, " with elements corresponding to betas' mean",sep=""))}
    if(length(beta.Norm[[2]]) != p ){stop(paste("error: beta.Norm[[2]] must be a vector of length, ",p, " with elements corresponding to betas' sd",sep=""))}
    beta.prior <- "normal"
  }

  if(!"sigma.sq.ig" %in% names(priors)){stop("error: sigma.sq.IG must be specified in priors value list")}
  sigma.sq.IG <- priors[["sigma.sq.ig"]]
  
  if(!is.vector(sigma.sq.IG) || length(sigma.sq.IG) != 2){stop("error: sigma.sq.IG must be a vector of length 2 in priors value list")}
  if(any(sigma.sq.IG <= 0)){stop("error: sigma.sq.IG must be a positive vector of length 2 in priors value list")}
 
  if(!"phi.unif" %in% names(priors)){stop("error: phi.Unif must be specified in priors value list")}
  phi.Unif <- priors[["phi.unif"]]
  
  if(!is.vector(phi.Unif) || length(phi.Unif) != 2){stop("error: phi.Unif must be a vector of length 2 in priors value list")}
  if(any(phi.Unif <= 0, phi.Unif[1] >= phi.Unif[2])){stop("error: phi.Unif must be a positive vector of length 2 with element 1 < element 2 in priors value list")}
  
  if(cov.model == "matern"){
    
    if(!"nu.unif" %in% names(priors)){stop("error: nu.Unif must be specified in priors")}
    nu.Unif <- priors[["nu.unif"]]
    
    if(!is.vector(nu.Unif) || length(nu.Unif) != 2){stop("error: nu.Unif must be a vector of length 2 in priors")}
    if(any(nu.Unif <= 0, nu.Unif[1] >= nu.Unif[2])){stop("error: nu.Unif must be a positive vector of length 2 with element 1 < element 2 in priors value list")}
  }

  storage.mode(sigma.sq.IG) <- "double"
  storage.mode(nu.Unif) <- "double"
  storage.mode(phi.Unif) <- "double"

  ####################################################
  ##Tuning values
  ####################################################
  beta.tuning <- 0
  phi.tuning <- 0
  sigma.sq.tuning <- 0
  nu.tuning <- 0
  w.tuning <- 0
   
  if(!missing(tuning)){
    
    names(tuning) <- tolower(names(tuning))
    
    if(!"beta" %in% names(tuning)){stop("error: beta must be specified in tuning value list")}
    beta.tuning <- tuning[["beta"]]
    
    if(is.matrix(beta.tuning)){
      if(nrow(beta.tuning) != p || ncol(beta.tuning) != p){
        stop(paste("error: if beta tuning is a matrix, it must be of dimension ",p,sep=""))
      }
      
      if(is.amcmc){
        beta.tuning <- diag(beta.tuning)
      }
      
    }else if(is.vector(beta.tuning)){
      if(length(beta.tuning) != p){
        stop(paste("error: if beta tuning is a vector, it must be of length ",p,sep=""))
      }
      
      if(!is.amcmc && p > 1){
        beta.tuning <- diag(beta.tuning)
      }
      
    }else{
      stop("error: beta tuning is misspecified")
    }
    
    if(!"sigma.sq" %in% names(tuning)){stop("error: sigma.sq must be specified in tuning value list")}
    sigma.sq.tuning <- tuning[["sigma.sq"]][1]
    
    if(!"phi" %in% names(tuning)){stop("error: phi must be specified in tuning value list")}
    phi.tuning <- tuning[["phi"]][1]
    
    if(cov.model == "matern"){
      if(!"nu" %in% names(tuning)){stop("error: nu must be specified in tuning value list")}
      nu.tuning <- tuning[["nu"]][1]
    }    
    
    if(!"w" %in% names(tuning)){stop("error: w must be specified in tuning value list")}
    w.tuning <- tuning[["w"]]
    
   
    if(!"w" %in% names(tuning)){stop("error: w must be specified in tuning value list")}
    w.tuning <- tuning[["w"]]

    if(is.pp){
      if(length(w.tuning) == 1){w.tuning <- rep(w.tuning, m)}  
      if(length(w.tuning) != m){stop(paste("error: w in the tuning value list must be a scalar of length 1 or vector of length ",m," (i.e., the number of predictive process knots)",sep=""))}
    }else{
      if(length(w.tuning) == 1){w.tuning <- rep(w.tuning, n)}
      if(length(w.tuning) != n){stop(paste("error: w in the tuning value list must be a scalar of length 1 or vector of length ",n," (i.e., the number of predictive process knots)",sep=""))}    
    }
   
  }else{##no tuning provided
    
    if(!is.amcmc){
      stop("error: tuning value list must be specified")
    }
    
    beta.tuning <- rep(0.01,p)
    phi.tuning <- 0.01
    sigma.sq.tuning <- 0.01
    nu.tuning <- 0.01
    
    if(is.pp){
      w.tuning <- rep(0.01,m)
    }else{
      w.tuning <- rep(0.01,n)
    }
    
  }
  
  storage.mode(beta.tuning) <- "double"
  storage.mode(phi.tuning) <- "double"
  storage.mode(sigma.sq.tuning) <- "double"
  storage.mode(nu.tuning) <- "double"
  storage.mode(w.tuning) <- "double"
  
  ####################################################
  ##Other stuff
  #################################################### 
  storage.mode(n.report) <- "integer"
  storage.mode(verbose) <- "integer"

  ####################################################
  ##Pack it up and off it goes
  ####################################################
  ptm <- proc.time()
  
  if(is.pp){
    
    if(!is.amcmc){
      out <- .Call("spPPGLM", Y, X, p, n, coords.D, family, weights,
                   m, knots.D, knots.coords.D,               
                   beta.prior, beta.Norm, sigma.sq.IG, nu.Unif, phi.Unif,
                   phi.starting, sigma.sq.starting, nu.starting, beta.starting, w.starting,
                   phi.tuning, sigma.sq.tuning, nu.tuning, beta.tuning, w.tuning,
                   cov.model, n.samples, verbose, n.report)  
    }else{
      out <- .Call("spPPGLM_AMCMC", Y, X, p, n, coords.D, family, weights,
                   m, knots.D, knots.coords.D,               
                   beta.prior, beta.Norm, sigma.sq.IG, nu.Unif, phi.Unif,
                   phi.starting, sigma.sq.starting, nu.starting, beta.starting, w.starting,
                   phi.tuning, sigma.sq.tuning, nu.tuning, beta.tuning, w.tuning,
                   cov.model, n.batch, batch.length, accept.rate, verbose, n.report)
    }
  }else{
    
    if(!is.amcmc){
      out <- .Call("spGLM", Y, X, p, n, coords.D, family, weights,
                   beta.prior, beta.Norm, sigma.sq.IG, nu.Unif, phi.Unif,
                   phi.starting, sigma.sq.starting, nu.starting, beta.starting, w.starting,
                   phi.tuning, sigma.sq.tuning, nu.tuning, beta.tuning, w.tuning,
                   cov.model, n.samples, verbose, n.report)
    }else{
      out <- .Call("spGLM_AMCMC", Y, X, p, n, coords.D, family, weights,
                   beta.prior, beta.Norm, sigma.sq.IG, nu.Unif, phi.Unif,
                   phi.starting, sigma.sq.starting, nu.starting, beta.starting, w.starting,
                   phi.tuning, sigma.sq.tuning, nu.tuning, beta.tuning, w.tuning,
                   cov.model, n.batch, batch.length, accept.rate, verbose, n.report)
    }
  }

  run.time <- proc.time() - ptm

  out$p.beta.theta.samples <- mcmc(t(out$p.beta.theta.samples))
  
  col.names <- rep("null",ncol(out$p.beta.theta.samples))
  
  col.names[1:p] <- x.names
  if(cov.model != "matern"){
    col.names[(p+1):(p+2)] <- c("sigma.sq", "phi")
  }else{
    col.names[(p+1):(p+3)] <- c("sigma.sq", "phi", "nu")
  }
    
  colnames(out$p.beta.theta.samples) <- col.names

  out$weights <- weights
  out$family <- family
  out$Y <- Y
  out$X <- X
  out$coords <- coords
  out$is.pp <- is.pp
  if(is.pp){
    out$knot.coords <- knot.coords
  }
  out$cov.model <- cov.model
  out$run.time <- run.time
  
  class(out) <- "spGLM"
  out  
}

