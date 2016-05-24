spDynLM <- function(formula, data = parent.frame(), coords, knots,
                    starting, tuning, priors, cov.model, get.fitted=FALSE,
                    n.samples,  verbose=TRUE, n.report=100, ...){
  
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
  ##Formula
  ####################################################
  if(missing(formula)){stop("error: formula must be specified")}
  
  if(class(formula) == "list"){
    
    holder <- mkspDynMats(formula, data)
    Y <- holder[[1]]
    X <- as.matrix(holder[[2]])
    x.names <- holder[[3]]

  }else{
    stop("error: formula is misspecified")
  }

  p <- ncol(X)
  n <- nrow(Y)
  N.t <- length(formula)
  miss <- as.vector(is.na(Y))

  ##make sure storage mode is correct
  storage.mode(Y) <- "double"
  storage.mode(X) <- "double"
  storage.mode(p) <- "integer"
  storage.mode(n) <- "integer"
  storage.mode(N.t) <- "integer"
  storage.mode(miss) <- "integer"

  ####################################################
  ##Distance matrices
  ####################################################
  
  ####################
  ##Coords
  #################### 
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
  coords.knots.D <- 0
  
  if(is.pp){
    knots.D <- iDist(knot.coords)
    m <- nrow(knots.D)
    coords.knots.D <- iDist(coords, knot.coords)

    if(min(coords.knots.D) == 0){

      stop("error: knots and observation coordinates cannot coincide. At least one knot location coincides with an observed coordinate.")
    }
    
  }

  storage.mode(m) <- "integer"
  storage.mode(knots.D) <- "double"
  storage.mode(coords.knots.D) <- "double"
  
  ####################################################
  ##Covariance model
  ####################################################
  if(missing(cov.model)){stop("error: cov.model must be specified")}
  if(!cov.model%in%c("gaussian","exponential","matern","spherical"))
    {stop("error: specified cov.model '",cov.model,"' is not a valid option; choose, from gaussian, exponential, matern, spherical.")}
  
  ####################################################
  ##Priors
  ####################################################
  if(missing(priors)) {stop("error: prior list for the parameters must be specified")}
    
  names(priors) <- tolower(names(priors))

  if(!"beta.0.norm" %in% names(priors)){stop("error: beta.0.norm must be specified")}
  beta.0.Norm <- priors[["beta.0.norm"]]
  if(!is.list(beta.0.Norm) || length(beta.0.Norm) != 2){stop("error: beta.0.Norm must be a list of length 2")}
  if(length(beta.0.Norm[[1]]) != p ){stop(paste("error: beta.0.Norm[[1]] must be a vector of length, ",p, "",sep=""))}
  if(length(beta.0.Norm[[2]]) != p^2 ){stop(paste("error: beta.0.Norm[[2]] must be a ",p,"x",p," covariance matrix",sep=""))}

  if(!"sigma.sq.ig" %in% names(priors)){stop("error: sigma.sq.IG must be specified")}
  sigma.sq.IG <- priors[["sigma.sq.ig"]]
  if(!is.list(sigma.sq.IG) || length(sigma.sq.IG) != 2){stop("error: sigma.sq.IG must be a list of length 2")}
  if(length(sigma.sq.IG[[1]]) != N.t){stop(paste("error: sigma.sq.IG[[1]] must be a vector of length, ",N.t, "",sep=""))}
  if(length(sigma.sq.IG[[2]]) != N.t){stop(paste("error: sigma.sq.IG[[2]] must be a vector of length, ",N.t, "",sep=""))}
  sigma.sq.IG <- as.vector(t(cbind(sigma.sq.IG[[1]],sigma.sq.IG[[2]])))

  if(!"tau.sq.ig" %in% names(priors)){stop("error: tau.sq.IG must be specified")}
  tau.sq.IG <- priors[["tau.sq.ig"]]
  if(!is.list(tau.sq.IG) || length(tau.sq.IG) != 2){stop("error: tau.sq.IG must be a list of length 2")}
  if(length(tau.sq.IG[[1]]) != N.t){stop(paste("error: tau.sq.IG[[1]] must be a vector of length, ",N.t, "",sep=""))}
  if(length(tau.sq.IG[[2]]) != N.t){stop(paste("error: tau.sq.IG[[2]] must be a vector of length, ",N.t, "",sep=""))}
  tau.sq.IG <- as.vector(t(cbind(tau.sq.IG[[1]],tau.sq.IG[[2]])))
  
  if(!"phi.unif" %in% names(priors)){stop("error: phi.Unif must be specified")}
  phi.Unif <- priors[["phi.unif"]]
  if(!is.list(phi.Unif) || length(phi.Unif) != 2){stop("error: phi.Unif must be a list of length 2")}
  if(length(phi.Unif[[1]]) != N.t){stop(paste("error: phi.Unif[[1]] must be a vector of length, ",N.t, "",sep=""))}
  if(length(phi.Unif[[2]]) != N.t){stop(paste("error: phi.Unif[[2]] must be a vector of length, ",N.t, "",sep=""))}
  if(any(phi.Unif[[2]]-phi.Unif[[1]] <= 0)){stop("error: phi.Unif has zero support")}
  phi.Unif <- as.vector(t(cbind(phi.Unif[[1]],phi.Unif[[2]])))

  nu.Unif <- 0
  if(cov.model == "matern"){
    if(!"nu.unif" %in% names(priors)){stop("error: nu.Unif must be specified")}
    nu.Unif <- priors[["nu.unif"]]
    if(!is.list(nu.Unif) || length(nu.Unif) != 2){stop("error: nu.Unif must be a list of length 2")}
    if(length(nu.Unif[[1]]) != N.t){stop(paste("error: nu.Unif[[1]] must be a vector of length, ",N.t, "",sep=""))}
    if(length(nu.Unif[[2]]) != N.t){stop(paste("error: nu.Unif[[2]] must be a vector of length, ",N.t, "",sep=""))}
    if(any(nu.Unif[[2]]-nu.Unif[[1]] <= 0)){stop("error: nu.Unif has zero support")}
    nu.Unif <- as.vector(t(cbind(nu.Unif[[1]],nu.Unif[[2]])))
  }
   
  if(!"sigma.eta.iw" %in% names(priors)){stop("error: Sigma.eta.IW must be specified")}
  sigma.eta.IW <- priors[["sigma.eta.iw"]]
  if(!is.list(sigma.eta.IW) || length(sigma.eta.IW) != 2){stop("error: Sigma.eta.IW must be a list of length 2")}
  if(length(sigma.eta.IW[[1]]) != 1){stop("error: Sigma.eta.IW[[1]] must be of length 1 (i.e., the IW df hyperparameter)")}   
  if(length(sigma.eta.IW[[2]]) != p^2){stop(paste("error: Sigma.eta.IW[[2]] must be a vector or matrix of length, ",p^2, ", (i.e., the IW scale matrix hyperparameter)",sep=""))}

  storage.mode(sigma.sq.IG) <- "double"
  storage.mode(tau.sq.IG) <- "double"
  storage.mode(phi.Unif) <- "double"
  storage.mode(nu.Unif) <- "double"
  storage.mode(sigma.eta.IW[[1]]) <- "double"; storage.mode(sigma.eta.IW[[2]]) <- "double"
  
  ####################################################
  ##Starting values
  ####################################################
  if(missing(starting)){stop("error: starting value list for the parameters must be specified")}
  
  names(starting) <- tolower(names(starting))   

  if(!"beta" %in% names(starting)){stop("error: beta must be specified in starting value list")}
  beta.starting <- starting[["beta"]]
  if(length(beta.starting) != N.t*p){stop(paste("error: beta starting must be of length ",N.t,"*",p,sep=""))}

  if(!"sigma.sq" %in% names(starting)){stop("error: sigma.sq must be specified in starting value list")}
  sigma.sq.starting <- starting[["sigma.sq"]]
  if(length(sigma.sq.starting) != N.t){stop(paste("error: sigma.sq starting must be a vector of length, ",N.t, "",sep=""))}

  if(!"tau.sq" %in% names(starting)){stop("error: a prior was spcified for tau.sq therefore tau.sq must be specified in starting value list")}
  tau.sq.starting <- starting[["tau.sq"]]
  if(length(tau.sq.starting) != N.t){stop(paste("error: tau.sq starting must be a vector of length, ",N.t, "",sep=""))}

  if(!"phi" %in% names(starting)){stop("error: phi must be specified in starting value list")}
  phi.starting <- starting[["phi"]]
  if(length(phi.starting) != N.t){stop(paste("error: phi starting must be a vector of length, ",N.t, "",sep=""))}
  
  nu.starting <- 0
  if(cov.model == "matern"){
    if(!"nu" %in% names(starting)){stop("error: nu must be specified in starting value list")}
      nu.starting <- starting[["nu"]]
    if(length(nu.starting) != N.t){stop(paste("error: nu starting must be a vector of length, ",N.t, "",sep=""))}
  }

  if(!"sigma.eta" %in% names(starting)){stop("error: Sigma.eta must be specified in starting value list")}
  sigma.eta.starting <- as.vector(starting[["sigma.eta"]])
  if(length(sigma.eta.starting) != p^2){stop(paste("error: Sigma.eta must be a positive definite matrix of length, ",p^2, sep=""))}
   
  storage.mode(beta.starting) <- "double"
  storage.mode(phi.starting) <- "double"
  storage.mode(sigma.sq.starting) <- "double"
  storage.mode(tau.sq.starting) <- "double"
  storage.mode(nu.starting) <- "double"
  storage.mode(sigma.eta.starting) <- "double"

  ####################################################
  ##Tuning values
  #################################################### 
  if(missing(tuning)){stop("error: tuning value vector for the spatial parameters must be specified")}
  
  names(tuning) <- tolower(names(tuning))
   
  if(!"phi" %in% names(tuning)){stop("error: phi must be specified in tuning value list")}
  phi.tuning <- tuning[["phi"]]
  if(length(phi.tuning) != N.t){stop(paste("error: phi tuning must be a vector of length, ",N.t, "",sep=""))}
  
  nu.tuning <- 0
  if(cov.model == "matern"){
    if(!"nu" %in% names(tuning)){stop("error: nu must be specified in tuning value list")}
    nu.tuning <- tuning[["nu"]]
    if(length(nu.tuning) != N.t){stop(paste("error: nu tuning must be a vector of length, ",N.t, "",sep=""))}
  }
    
  storage.mode(phi.tuning) <- "double"
  storage.mode(nu.tuning) <- "double"

  ####################################################
  ##Other stuff
  ####################################################
  if(missing(n.samples)){stop("error: n.samples needs to be specified")}
  storage.mode(n.samples) <- "integer"

  storage.mode(get.fitted) <- "integer"
  storage.mode(n.report) <- "integer"
  storage.mode(verbose) <- "integer"

  ####################################################
  ##Pack it up and off it goes
  ####################################################
  ptm <- proc.time()

  if(is.pp){
    out <- .Call("spPPDynLM", Y, t(X), p, n, m, N.t, coords.D, knots.D, coords.knots.D,
                 beta.0.Norm, sigma.sq.IG, tau.sq.IG, nu.Unif, phi.Unif, sigma.eta.IW,
                 beta.starting, phi.starting, sigma.sq.starting, tau.sq.starting, nu.starting, sigma.eta.starting,
                 phi.tuning, nu.tuning,
                 cov.model, n.samples, miss, get.fitted, verbose, n.report)   
  }else{
    out <- .Call("spDynLM", Y, t(X), p, n, N.t, coords.D,
                 beta.0.Norm, sigma.sq.IG, tau.sq.IG, nu.Unif, phi.Unif, sigma.eta.IW,
                 beta.starting, phi.starting, sigma.sq.starting, tau.sq.starting, nu.starting, sigma.eta.starting,
                 phi.tuning, nu.tuning,
                 cov.model, n.samples, miss, get.fitted, verbose, n.report)   
  }
  
  run.time <- proc.time() - ptm

  out$p.beta.0.samples <- mcmc(t(out$p.beta.0.samples))
  colnames(out$p.beta.0.samples) <- x.names
  
  out$p.beta.samples <- mcmc(t(out$p.beta.samples))
  colnames(out$p.beta.samples) <- as.vector(t(sapply(paste(x.names,".t",sep=""),paste,1:N.t,sep="")))
  
  out$p.theta.samples <- mcmc(t(out$p.theta.samples))
  
  if(cov.model != "matern"){
    colnames(out$p.theta.samples) <- as.vector(t(sapply(paste(c("sigma.sq", "tau.sq", "phi"),".t",sep=""),paste,1:N.t,sep="")))
  }else{
    colnames(out$p.theta.samples) <- as.vector(t(sapply(paste(c("sigma.sq", "tau.sq", "phi", "nu"),".t",sep=""),paste,1:N.t,sep="")))
  }

  out$p.sigma.eta.samples <- mcmc(t(out$p.sigma.eta.samples))
  colnames(out$p.sigma.eta.samples) <- paste("eta[",matrix(apply(cbind(expand.grid(x.names,x.names)), 1, function(x) paste(x, collapse=",")),p,p),"]",sep="")
    
  out$Y <- Y
  out$X <- X
  out$coords <- coords
  out$cov.model <- cov.model
  out$x.names <- x.names
  out$run.time <- run.time
  out$missing.indx <- miss

  if(is.pp){
    out$knot.coords <- knot.coords
  }
 
  class(out) <- "spDynLM"
  out  
}

