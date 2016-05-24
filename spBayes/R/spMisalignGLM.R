spMisalignGLM <- function(formula, family="binomial", weights, data = parent.frame(), coords,
                          starting, tuning, priors, cov.model,
                          amcmc, n.samples,
                          verbose=TRUE, n.report=100, ...){
  
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
  
  if(is.list(formula) && is.list(coords)){

    if(length(formula) != length(coords)){
      stop("error: formula and coords are misspecified")
    }
    
    mod.dat <- mkMisalignYX(formula, data)
    Y <- mod.dat[[1]]
    X <- mod.dat[[2]]
    misalign.n <- mod.dat[[3]]
    misalign.p <- mod.dat[[4]]
    x.names <- mod.dat[[5]]
    m <- length(formula)

    storage.mode(misalign.n) <- "integer"
    storage.mode(misalign.p) <- "integer"
    
  }else{
    stop("error: formula is misspecified")
  }
  
  p <- ncol(X)
  n <- nrow(X)
  n.ltr <- m*(m+1)/2

  ##make sure storage mode is correct
  storage.mode(Y) <- "double"
  storage.mode(X) <- "double"
  storage.mode(m) <- "integer"
  storage.mode(p) <- "integer"
  storage.mode(n) <- "integer"

  ####################################################
  ##family and weights
  ####################################################
  if(!family %in% c("binomial","poisson"))
    stop("error: family must be binomial or poisson")

  if(missing(weights)){
    weights <- rep(1, n)
  }else if(is.list(weights) && length(weights) == m && all(sapply(weights, length) == misalign.n)){
    weights<- do.call(c, weights)
  }else{
    stop(paste("error: weights must be a list length ", m," consisting of weight vectors of lengths ", paste(misalign.n, collapse=" "), sep=""))
  }
  
  storage.mode(weights) <- "integer"
  
  ####################################################
  ##sampling method
  ####################################################
  n.batch <- 0
  batch.length <- 0
  accept.rate <- 0
  is.amcmc <- TRUE

  ##require amcmc for now
  if(missing(amcmc)){
    
    stop("error: amcmc must be be specified")
    ## if(missing(n.samples)){stop("error: n.samples needs to be specified")}
    
    ## n.batch <- n.samples
    ## batch.length <- 1
    ## is.amcmc <- FALSE
        
  }else{
    
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

  }

  storage.mode(is.amcmc) <- "integer"
  storage.mode(n.batch) <- "integer"
  storage.mode(batch.length) <- "integer"
  storage.mode(accept.rate) <- "double"

  ####################################################
  ##Distance matrices
  ####################################################
  
  ####################
  ##Coords
  ####################
  if(missing(coords)){stop("error: coords must be specified")}
  
  coords <- as.matrix(do.call(rbind, coords))

  if(ncol(coords) != 2 || nrow(coords) != n){
    stop("error: either the coords have more than two columns or then number of rows is different than data used in the model formula")
  }
  
  coords.D <- iDist(coords)
  storage.mode(coords.D) <- "double"
  
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
  A.starting <- 0
  phi.starting <- 0
  nu.starting <- 0
  w.starting <- 0
  
  if(missing(starting)){stop("error: starting value list for the parameters must be specified")}
  
  names(starting) <- tolower(names(starting))   

  if("beta" %in% names(starting)){
    beta.starting <- starting[["beta"]]
    if(length(beta.starting) != p){stop(paste("error: starting values for beta must be of length ",p,sep=""))}
  }else{
    stop("error: beta must be specified in starting")
  }
   
  if(!"a" %in% names(starting)){stop("error: A must be specified in starting")}
  A.starting <- starting[["a"]]
  if(length(A.starting) != n.ltr){stop(paste("error: A must be of length ",n.ltr," in starting value list",sep=""))}
  
  if(!"phi" %in% names(starting)){stop("error: phi must be specified in starting")}
  phi.starting <- starting[["phi"]]
  if(length(phi.starting) != m){stop(paste("error: phi must be of length ",m," in starting value list",sep=""))}
  
  if(cov.model == "matern"){
    if(!"nu" %in% names(starting)){stop("error: nu must be specified in starting")}
    nu.starting <- starting[["nu"]]
    if(length(nu.starting) != m){stop(paste("error: nu must be of length ",m," in starting value list",sep=""))}
  }

  if(!"w" %in% names(starting)){stop("error: w must be specified in starting value list")}
  w.starting <- starting[["w"]]
  
  if(is.numeric(w.starting) && length(w.starting) == 1){
    w.starting <- rep(w.starting, n)
  }else if(is.list(w.starting) && length(w.starting) == m && all(sapply(w.starting, length) == misalign.n)){
       w.starting <- do.call(c, w.starting)
  }else{
    stop(paste("error: w in the starting value list must be a scalar of length 1 or a list of length ",m, " consisting of vectors of lengths ", paste(misalign.n, collapse=" "), sep=""))
  }
  
  storage.mode(beta.starting) <- "double"
  storage.mode(phi.starting) <- "double"
  storage.mode(A.starting) <- "double"
  storage.mode(nu.starting) <- "double"
  storage.mode(w.starting) <- "double"

  ####################################################
  ##Priors
  ####################################################
  beta.Norm <- 0
  beta.prior <- "flat"
  K.IW <- 0
  nu.Unif <- 0
  phi.Unif <- 0

  if(missing(priors)) {stop("error: prior list for the parameters must be specified")}
    
  names(priors) <- tolower(names(priors))

  if("beta.norm" %in% names(priors)){
    beta.Norm <- priors[["beta.normal"]]
    if(!is.list(beta.Norm) || length(beta.Norm) != 2){stop("error: beta.Norm must be a list of length 2")}
    if(length(beta.Norm[[1]]) != p ){stop(paste("error: beta.Norm[[1]] must be a vector of length, ",p, " with elements corresponding to betas' mean",sep=""))}
    if(length(beta.Norm[[2]]) != p ){stop(paste("error: beta.Norm[[2]] must be a vector of length, ",p, " with elements corresponding to betas' sd",sep=""))}
    beta.prior <- "normal"
  }
  
  if(!"k.iw" %in% names(priors)){stop("error: K.IW must be specified")}
  K.IW <- priors[["k.iw"]]
  if(!is.list(K.IW) || length(K.IW) != 2){stop("error: K.IW must be a list of length 2")}
  if(length(K.IW[[1]]) != 1 ){stop("error: K.IW[[1]] must be of length 1 (i.e., the IW df hyperparameter)")}   
  if(length(K.IW[[2]]) != m^2 ){stop(paste("error: K.IW[[2]] must be a vector or matrix of length, ",m^2, ", (i.e., the IW scale matrix hyperparameter)",sep=""))}
  
  if(!"phi.unif" %in% names(priors)){stop("error: phi.Unif must be specified")}
  phi.Unif <- priors[["phi.unif"]]
  if(!is.list(phi.Unif) || length(phi.Unif) != 2){stop("error: phi.Unif must be a list of length 2")}
  if(length(phi.Unif[[1]]) != m){stop(paste("error: phi.Unif[[1]] must be a vector of length, ",m, "",sep=""))}
  if(length(phi.Unif[[2]]) != m){stop(paste("error: phi.Unif[[2]] must be a vector of length, ",m, "",sep=""))}
  if(any(phi.Unif[[2]]-phi.Unif[[1]] <= 0)){stop("error: phi.Unif has zero support")}
  phi.Unif <- as.vector(t(cbind(phi.Unif[[1]],phi.Unif[[2]])))
  
  if(cov.model == "matern"){
    
    if(!"nu.unif" %in% names(priors)){stop("error: nu.Unif must be specified")}
    nu.Unif <- priors[["nu.unif"]]
    if(!is.list(nu.Unif) || length(nu.Unif) != 2){stop("error: nu.Unif must be a list of length 2")}
    if(length(nu.Unif[[1]]) != m){stop(paste("error: nu.Unif[[1]] must be a vector of length, ",m, "",sep=""))}
    if(length(nu.Unif[[2]]) != m){stop(paste("error: nu.Unif[[2]] must be a vector of length, ",m, "",sep=""))}
    if(any(nu.Unif[[2]]-nu.Unif[[1]] <= 0)){stop("error: nu.Unif has zero support")}
    nu.Unif <- as.vector(t(cbind(nu.Unif[[1]],nu.Unif[[2]])))
  }
    
  storage.mode(K.IW[[1]]) <- "double"; storage.mode(K.IW[[2]]) <- "double"
  storage.mode(nu.Unif) <- "double"
  storage.mode(phi.Unif) <- "double"

  ####################################################
  ##Tuning values
  ####################################################
  beta.tuning <- 0
  phi.tuning <- 0
  A.tuning <- 0
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
      
      if(!is.amcmc){
        beta.tuning <- diag(beta.tuning)
      }
      
    }else{
      stop("error: beta tuning is misspecified")
    }
    
    if(!"a" %in% names(tuning)){stop("error: A must be specified in tuning value list")}
    A.tuning <- as.vector(tuning[["a"]])
    if(length(A.tuning) != n.ltr){stop(paste("error: A must be of length ",n.ltr," in tuning value list",sep=""))}
    
    if(!"phi" %in% names(tuning)){stop("error: phi must be specified in tuning value list")}
    phi.tuning <- tuning[["phi"]]
    if(length(phi.tuning) != m){stop(paste("error: phi must be of length ",m," in tuning value list",sep=""))}
    
    if(cov.model == "matern"){
      if(!"nu" %in% names(tuning)){stop("error: nu must be specified in tuning value list")}
      nu.tuning <- tuning[["nu"]]
      if(length(nu.tuning) != m){stop(paste("error: nu must be of length ",m," in tuning value list",sep=""))}
    }    
    
    if(!"w" %in% names(tuning)){stop("error: w must be specified in tuning value list")}
    w.tuning <- tuning[["w"]]

    if(is.numeric(w.tuning) && length(w.tuning) == 1){
      w.tuning <- rep(w.tuning, n)
    }else if(is.list(w.tuning) && length(w.tuning) == m && all(sapply(w.tuning, length) == misalign.n)){
      w.tuning <- do.call(c, w.tuning)
    }else{
      stop(paste("error: w in the tuning value list must be a scalar of length 1 or a list of length ",m, " consisting of vectors of lengths ", paste(misalign.n, collapse=" "), sep=""))
    }
 
  }else{##no tuning provided
    
    if(!is.amcmc){
      stop("error: tuning value list must be specified")
    }
    
    beta.tuning <- rep(0.01,p)
    phi.tuning <- rep(0.01,m)
    A.tuning <- rep(0.01,m*(m-1)/2+m)
    nu.tuning <- rep(0.01,m)
    w.tuning <- rep(0.01,n)
  }

  storage.mode(beta.tuning) <- "double"
  storage.mode(phi.tuning) <- "double"
  storage.mode(A.tuning) <- "double"
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

  out <- .Call("spGLMMisalign_AMCMC", Y, X, misalign.p, misalign.n, m, coords.D, family, weights,
               beta.prior, beta.Norm,
               K.IW, nu.Unif, phi.Unif,
               phi.starting, A.starting, nu.starting, beta.starting, w.starting,
               phi.tuning, A.tuning, nu.tuning, beta.tuning, w.tuning,
               cov.model, n.batch, batch.length, accept.rate, verbose, n.report)
 
  run.time <- proc.time() - ptm

  ##parameter names
  out$p.beta.theta.samples <-  mcmc(t(out$p.beta.theta.samples))

  col.names <- rep("null",ncol(out$p.beta.theta.samples))
  
  if(cov.model != "matern"){
    col.names <- c(x.names, rep("K",n.ltr), paste("phi[",1:m,"]",sep=""))
  }else{
    col.names <- c(x.names, rep("K",n.ltr), paste("phi[",1:m,"]",sep=""), paste("nu[",1:m,"]",sep=""))
  }
    
  colnames(out$p.beta.theta.samples) <- col.names
  
  AtA <- function(x, m){
    A <- matrix(0, m, m)
    A[lower.tri(A, diag=TRUE)] <- x
    (A%*%t(A))[lower.tri(A, diag=TRUE)]
  }
  
  K.names <- paste("K[",matrix(apply(cbind(expand.grid(1:m,1:m)), 1, function(x) paste(x, collapse=",")),m,m)[lower.tri(matrix(0,m,m), diag=TRUE)],"]",sep="")
  
  colnames(out$p.beta.theta.samples)[colnames(out$p.beta.theta.samples)%in%"K"] <- K.names
  out$p.beta.theta.samples[,K.names] <- t(apply(out$p.beta.theta.samples[,K.names,drop=FALSE], 1, AtA, m))
    
  out$weights <- weights
  out$family <- family
  out$Y <- Y
  out$X <- X
  out$m <- m
  out$misalign.p <- misalign.p
  out$misalign.n <- misalign.n
  out$coords <- coords
  out$x.names <- x.names
  out$cov.model <- cov.model
  out$run.time <- run.time

  class(out) <- "spMisalignGLM"
  out  

}

