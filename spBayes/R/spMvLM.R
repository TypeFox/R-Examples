spMvLM <- function(formula, data = parent.frame(), coords, knots,
                   starting, tuning, priors, cov.model,
                   modified.pp=TRUE, amcmc, n.samples,
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

  if(class(formula) == "formula"){
    
    holder <- parseFormula(formula, data)
    Y <- holder[[1]]
    X <- holder[[2]]
    x.names <- holder[[3]]
    m <- 1
    
  }else if(class(formula) == "list"){
    
    mv.mats <- mkMats(formula, data)
    Y <- mv.mats[[1]]
    X <- mv.mats[[2]]
    x.names <- mv.mats[[3]]
    m <- length(formula)
    
  }else{
    stop("error: formula is misspecified")
  }
 
  p <- ncol(X)
  n <- nrow(X)/m 
  n.ltr <- m*(m+1)/2

  ##make sure storage mode is correct
  storage.mode(Y) <- "double"
  storage.mode(X) <- "double"
  storage.mode(m) <- "integer"
  storage.mode(p) <- "integer"
  storage.mode(n) <- "integer"

  ####################################################
  ##sampling method
  ####################################################
  n.batch <- 0
  batch.length <- 0
  accept.rate <- 0
  is.amcmc <- TRUE
  
  if(missing(amcmc)){
   
    if(missing(n.samples)){stop("error: n.samples needs to be specified")}
    
    n.batch <- n.samples
    batch.length <- 1
    is.amcmc <- FALSE
        
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
  
  g <- 0
  knots.D <- 0
  knots.coords.D <- 0
  
  if(is.pp){
    knots.D <- iDist(knot.coords)
    g <- nrow(knots.D)
    knots.coords.D <- iDist(knot.coords, coords)
  }

  storage.mode(modified.pp) <- "integer"
  storage.mode(g) <- "integer"
  storage.mode(knots.D) <- "double"
  storage.mode(knots.coords.D) <- "double"
  
  ####################################################
  ##Covariance model
  ####################################################
  if(missing(cov.model)){stop("error: cov.model must be specified")}
  if(!cov.model%in%c("gaussian","exponential","matern","spherical"))
    {stop("error: specified cov.model '",cov.model,"' is not a valid option; choose, from gaussian, exponential, matern, spherical.")}

  ####################################################
  ##Priors
  ####################################################
  beta.Norm <- 0
  beta.prior <- "flat"
  K.prior <- 0
  K.prior.name <- 0
  Psi.prior <- 0
  Psi.prior.name <- 0
  Psi.diag <- FALSE
  nu.Unif <- 0
  phi.Unif <- 0
  nugget <- FALSE
  
  if(missing(priors)) {stop("error: prior list for the parameters must be specified")}
    
  names(priors) <- tolower(names(priors))

  if("beta.norm" %in% names(priors)){
    beta.Norm <- priors[["beta.norm"]]
    if(!is.list(beta.Norm) || length(beta.Norm) != 2){stop("error: beta.Norm must be a list of length 2")}
    if(length(beta.Norm[[1]]) != p){stop(paste("error: beta.Norm[[1]] must be a vector of length, ",p, "",sep=""))}
    if(length(beta.Norm[[2]]) != p^2){stop(paste("error: beta.Norm[[2]] must be a ",p,"x",p," covariance matrix",sep=""))}
    beta.prior <- "normal"
  }
  
  if("k.iw" %in% names(priors)){
    K.prior <- priors[["k.iw"]]
    if(!is.list(K.prior) || length(K.prior) != 2){stop("error: K.IW must be a list of length 2")}
    if(length(K.prior[[1]]) != 1){stop("error: K.IW[[1]] must be of length 1 (i.e., the IW df hyperparameter)")}   
    if(length(K.prior[[2]]) != m^2){stop(paste("error: K.IW[[2]] must be a vector or matrix of length, ",m^2, ", (i.e., the IW scale matrix hyperparameter)",sep=""))}
    K.prior.name <- "IW"
  }else if("a.norm" %in% names(priors)){
    K.prior <- priors[["a.norm"]]
    if(!is.list(K.prior) || length(K.prior) != 2){stop("error: A.Norm must be a list of length 2")}
    if(length(K.prior[[1]]) != n.ltr){stop(paste("error: A.Norm[[1]] must be a vector of length, ",n.ltr, "",sep=""))}
    if(length(K.prior[[2]]) != n.ltr){stop(paste("error: A.Norm[[2]] must be a vector of length, ",n.ltr, "",sep=""))}
    K.prior.name <- "normal"
  }else{
    stop("error: a valid prior for K or A must be specified")
  }
  
  if("psi.iw" %in% names(priors)){
    Psi.prior <- priors[["psi.iw"]]
    if(!is.list(Psi.prior) || length(Psi.prior) != 2){stop("error: Psi.IW must be a list of length 2")}
    if(length(Psi.prior[[1]]) != 1){stop("error: Psi.IW[[1]] must be of length 1 (i.e., the IW df hyperparameter)")}   
    if(length(Psi.prior[[2]]) != m^2){stop(paste("error: Psi.IW[[2]] must be a vector or matrix of length, ",m^2, ", (i.e., the IW scale matrix hyperparameter)",sep=""))}
    Psi.prior.name <- "IW"
    nugget <- TRUE
  }
 
  if("l.norm" %in% names(priors)){
    Psi.prior <- priors[["l.norm"]]
    if(!is.list(Psi.prior) || length(Psi.prior) != 2){stop("error: L.Norm must be a list of length 2")}
    if(length(Psi.prior[[1]]) != n.ltr){stop(paste("error: L.Norm[[1]] must be a vector of length, ",n.ltr, "",sep=""))}
    if(length(Psi.prior[[2]]) != n.ltr){stop(paste("error: L.Norm[[2]] must be a vector of length, ",n.ltr, "",sep=""))}
    Psi.prior.name <- "normal"
    nugget <- TRUE
  }

  if("psi.ig" %in% names(priors)){
    Psi.prior <- priors[["psi.ig"]]
    if(!is.list(Psi.prior) || length(Psi.prior) != 2){stop("error: Psi.IG must be a list of length 2")}
    if(length(Psi.prior[[1]]) != m){stop(paste("error: Psi.IG[[1]] must be a vector of length, ",m, "",sep=""))}
    if(length(Psi.prior[[2]]) != m){stop(paste("error: Psi.IG[[2]] must be a vector of length, ",m, "",sep=""))}
    Psi.prior.name <- "IG"
    Psi.diag <- TRUE
    nugget <- TRUE
  }
   
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
   
  storage.mode(K.prior[[1]]) <- "double"; storage.mode(K.prior[[2]]) <- "double"
  if(nugget){
    storage.mode(Psi.prior[[1]]) <- "double"; storage.mode(Psi.prior[[2]]) <- "double"
  }
  storage.mode(nu.Unif) <- "double"
  storage.mode(phi.Unif) <- "double"
  storage.mode(nugget) <- "integer"
  storage.mode(Psi.diag) <- "integer"
 
  ####################################################
  ##Starting values
  ####################################################
  beta.starting <- 0
  A.starting <- 0
  L.starting <- 0
  phi.starting <- 0
  nu.starting <- 0
  
  if(missing(starting)){stop("error: starting value list for the parameters must be specified")}
  
  names(starting) <- tolower(names(starting))   

  if(is.pp){
    if(!"beta" %in% names(starting)){
      beta.starting <- as.vector(coefficients(lm(Y~X-1)))
    }else{
      beta.starting <- starting[["beta"]]
    }
  }
  
  if(!"a" %in% names(starting)){stop("error: A must be specified in starting")}
  A.starting <- starting[["a"]]
  if(length(A.starting) != n.ltr){stop(paste("error: A must be of length ",n.ltr," in starting value list",sep=""))}

  if(nugget){
    if(Psi.diag){
      if(!"psi" %in% names(starting)){stop("error: Psi is specified as diagonal so Psi must be specified in starting value list")}
      L.starting <- as.vector(starting[["psi"]])
      if(length(L.starting) != m){stop(paste("error: Psi is specified as diagonal so Psi must be of length ",m," in starting value list",sep=""))}
    }else{
      if(!"l" %in% names(starting)){stop("error: L must be specified in starting value list")}
      L.starting <- as.vector(starting[["l"]])
      if(length(L.starting) != n.ltr){stop(paste("error: L must be of length ",n.ltr," in starting value list",sep=""))}
    }
  }
  
  if(!"phi" %in% names(starting)){stop("error: phi must be specified in starting")}
  phi.starting <- starting[["phi"]]
  if(length(phi.starting) != m){stop(paste("error: phi must be of length ",m," in starting value list",sep=""))}
  
  if(cov.model == "matern"){
    if(!"nu" %in% names(starting)){stop("error: nu must be specified in starting")}
    nu.starting <- starting[["nu"]]
    if(length(nu.starting) != m){stop(paste("error: nu must be of length ",m," in starting value list",sep=""))}
  }

  storage.mode(beta.starting) <- "double"
  storage.mode(phi.starting) <- "double"
  storage.mode(A.starting) <- "double"
  storage.mode(L.starting) <- "double"
  storage.mode(nu.starting) <- "double"
      
  ####################################################
  ##Tuning values
  ####################################################
  phi.tuning <- 0
  A.tuning <- 0
  L.tuning <- 0
  nu.tuning <- 0
  
  if(missing(tuning)){stop("error: tuning value vector for the spatial parameters must be specified")}
  
  names(tuning) <- tolower(names(tuning))
  
  if(!"a" %in% names(tuning)){stop("error: A must be specified in tuning value list")}
  A.tuning <- as.vector(tuning[["a"]])
  if(length(A.tuning) != n.ltr){stop(paste("error: A must be of length ",n.ltr," in tuning value list",sep=""))}
  
  if(nugget){
    if(Psi.diag){
      if(!"psi" %in% names(tuning)){stop("error: Psi is specified as diagonal so Psi must be specified in tuning value list")}
      L.tuning <- as.vector(tuning[["psi"]])
      if(length(L.tuning) != m){stop(paste("error: Psi is specified as diagonal so Psi must be of length ",m," in tuning value list",sep=""))}
    }else{
      if(!"l" %in% names(tuning)){stop("error: L must be specified in tuning value list")}
      L.tuning <- as.vector(tuning[["l"]])
      if(length(L.tuning) != n.ltr){stop(paste("error: L must be of length ",n.ltr," in tuning value list",sep=""))}
    }
  }
  
  if(!"phi" %in% names(tuning)){stop("error: phi must be specified in tuning value list")}
  phi.tuning <- tuning[["phi"]]
  if(length(phi.tuning) != m){stop(paste("error: phi must be of length ",m," in tuning value list",sep=""))}
    
  if(cov.model == "matern"){
    if(!"nu" %in% names(tuning)){stop("error: nu must be specified in tuning value list")}
    nu.tuning <- tuning[["nu"]]
    if(length(nu.tuning) != m){stop(paste("error: nu must be of length ",m," in tuning value list",sep=""))}
  }    
    
  storage.mode(phi.tuning) <- "double"
  storage.mode(A.tuning) <- "double"
  storage.mode(L.tuning) <- "double"
  storage.mode(nu.tuning) <- "double"

  ####################################################
  ##Other stuff
  ####################################################
  storage.mode(n.report) <- "integer"
  storage.mode(verbose) <- "integer"
  
  ####################################################
  ##Fix for no nugget pp model
  ####################################################
  if(is.pp && !nugget){
    
    if(!modified.pp){
      modified.pp <- TRUE
      warning("Because the specified predictive process model has no nugget, a modified predictive process was used.")
    }

    if(any(knots.coords.D == 0)){

      jit <- 1000*.Machine$double.eps
      warning(paste("Because the specified predictive process model has no nugget, a shift of distance ",jit," was added to those ", sum(knots.coords.D == 0), " knot locations that coincide with observed locations.",sep=""))
      
      while(any(knots.coords.D == 0)){
        ind <- which(knots.coords.D == 0, arr.ind = TRUE)
        knot.coords[ind[1,]] <- knot.coords[ind[1,]]+jit       
        knots.coords.D <- iDist(knot.coords, coords)
        knots.D <- iDist(knot.coords)
      }
    }
     
  }

  ####################################################
  ##Pack it up and off it goes
  ####################################################
  ptm <- proc.time()

  if(is.pp){

    ##transpose X to simplify recovery of beta
    XX <- t(X)
    storage.mode(XX) <- "double"
    
    out <- .Call("spPPMvLM", Y, XX, p, n, m, g, coords.D, knots.D, knots.coords.D, modified.pp, 
                 beta.prior, beta.Norm,
                 K.prior, K.prior.name,
                 Psi.prior, Psi.prior.name, Psi.diag,
                 nu.Unif, phi.Unif,
                 beta.starting, phi.starting, A.starting, L.starting, nu.starting, 
                 phi.tuning, A.tuning, L.tuning, nu.tuning,
                 nugget, cov.model, is.amcmc, n.batch, batch.length, accept.rate, verbose, n.report)
  }else{
    out <- .Call("spMvLM", Y, X, p, n, m, coords.D,
                 beta.prior, beta.Norm,
                 K.prior, K.prior.name,
                 Psi.prior, Psi.prior.name, Psi.diag,
                 nu.Unif, phi.Unif,
                 phi.starting, A.starting, L.starting, nu.starting, 
                 phi.tuning, A.tuning, L.tuning, nu.tuning,
                 nugget, cov.model, is.amcmc, n.batch, batch.length, accept.rate, verbose, n.report)
  }
  
  run.time <- proc.time() - ptm

  ##parameter names
  if(is.pp){
    out$p.beta.samples <- mcmc(t(out$p.beta.samples))
    colnames(out$p.beta.samples) <- x.names
  }
  
  out$p.theta.samples <- mcmc(t(out$p.theta.samples))

  col.names <- rep("null",ncol(out$p.theta.samples))
  
  if(!nugget && cov.model != "matern"){
    col.names <- c(rep("K",n.ltr), paste("phi[",1:m,"]",sep=""))
  }else if(nugget && cov.model != "matern"){
    if(Psi.diag){
      col.names <- c(rep("K",n.ltr), rep("Psi",m), paste("phi[",1:m,"]",sep=""))
    }else{
      col.names<- c(rep("K",n.ltr), rep("Psi",n.ltr), paste("phi[",1:m,"]",sep=""))
    }
  }else if(!nugget && cov.model == "matern"){
    col.names <- c(rep("K",n.ltr), paste("phi[",1:m,"]",sep=""), paste("nu[",1:m,"]",sep=""))
  }else{
    if(Psi.diag){
      col.names<- c(rep("K",n.ltr), rep("Psi",m), paste("phi[",1:m,"]",sep=""), paste("nu[",1:m,"]",sep=""))
    }else{
      col.names <- c(rep("K",n.ltr), rep("Psi",n.ltr), paste("phi[",1:m,"]",sep=""), paste("nu[",1:m,"]",sep=""))
    }
  }
    
  colnames(out$p.theta.samples) <- col.names
  
  AtA <- function(x, m){
    A <- matrix(0, m, m)
    A[lower.tri(A, diag=TRUE)] <- x
    (A%*%t(A))[lower.tri(A, diag=TRUE)]
  }
  
  K.names <- paste("K[",matrix(apply(cbind(expand.grid(1:m,1:m)), 1, function(x) paste(x, collapse=",")),m,m)[lower.tri(matrix(0,m,m), diag=TRUE)],"]",sep="")
  
  colnames(out$p.theta.samples)[colnames(out$p.theta.samples)%in%"K"] <- K.names
  out$p.theta.samples[,K.names] <- t(apply(out$p.theta.samples[,K.names,drop=FALSE], 1, AtA, m))
  
  if(nugget){
    if(Psi.diag){
      Psi.names <- paste("Psi[",1:m,",",1:m,"]",sep="")
      colnames(out$p.theta.samples)[colnames(out$p.theta.samples)%in%"Psi"] <- Psi.names
    }else{
      Psi.names <- paste("Psi[",matrix(apply(cbind(expand.grid(1:m,1:m)), 1, function(x) paste(x, collapse=",")),m,m)[lower.tri(matrix(0,m,m), diag=TRUE)],"]",sep="")
      colnames(out$p.theta.samples)[colnames(out$p.theta.samples)%in%"Psi"] <- Psi.names
      out$p.theta.samples[,Psi.names] <- t(apply(out$p.theta.samples[,Psi.names,drop=FALSE], 1, AtA, m))
    }
  }

  out$Y <- Y
  out$X <- X
  out$m <- m
  out$coords <- coords
  out$is.pp <- is.pp
  if(is.pp){
    out$knot.coords <- knot.coords
    out$modified.pp <- modified.pp
  }  
  out$cov.model <- cov.model
  out$nugget <- nugget
  out$Psi.diag <- Psi.diag
  out$beta.prior <- beta.prior
  out$beta.Norm <- beta.Norm
  out$x.names <- x.names
  out$run.time <- run.time
  
  class(out) <- "spMvLM"
  out  

}

