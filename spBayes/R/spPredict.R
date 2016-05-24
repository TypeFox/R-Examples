spPredict <- function(sp.obj, pred.coords, pred.covars, start=1, end, thin=1, verbose=TRUE, n.report=100, ...){
  
  ####################################################
  ##Check for unused args
  ####################################################
  formal.args <- names(formals(sys.function(sys.parent())))
  elip.args <- names(list(...))
  for(i in elip.args){
    if(! i %in% formal.args)
      warning("'",i, "' is not an argument")
  }
  
  if(missing(sp.obj)){stop("error: spPredict expects sp.obj\n")}
  if(!class(sp.obj) %in% c("spLM", "spMvLM", "spGLM", "spMvGLM","bayesLMRef","bayesGeostatExact","nonSpGLM","nonSpMvGLM","spMisalignLM","spMisalignGLM")){
    stop("error: requires an output object of class spLM, spMvLM, spGLM, spMvGLM, bayesGeostatExact, bayesLMRef, nonSpGLM, nonSpMvGLM, spMisalignLM, or spMisalignGLM\n")
  }

  obj.class <- class(sp.obj)

  ##
  ##non spatial model prediction
  ##
  if(obj.class %in% c("nonSpGLM", "nonSpMvGLM")){

    X <- sp.obj$X
    Y <- sp.obj$Y
    family <- sp.obj$family
    p.beta.samples <- sp.obj$p.beta.samples
    n.samples <- nrow(p.beta.samples)
    
    ##subsamples
    if(missing(end)){end <- n.samples}
    
    if(!is.numeric(start) || start >= n.samples)
      stop("error: invalid start")
    if(!is.numeric(end) || end > n.samples) 
      stop("error: invalid end")
    if(!is.numeric(thin) || thin >= n.samples) 
      stop("error: invalid thin")
       
    s.indx <- seq(start, end, by=as.integer(thin))
    
    p.beta.samples <- p.beta.samples[s.indx,,drop=FALSE]
    n.samples <- nrow(p.beta.samples)

    out <- list()
    
    if(family == "binomial"){
      out$p.predictive.samples <- apply(p.beta.samples, 1,  function(s){1/(1+exp(-pred.covars%*%s))})          
    }else{##poisson
       out$p.predictive.samples <- apply(p.beta.samples, 1,  function(s){exp(pred.covars%*%s)}) 
     }

    return(out)
  }

  if(obj.class == "bayesLMRef"){

    X <- sp.obj$X
    Y <- sp.obj$Y
    p <- ncol(X)
    n <- nrow(X)
    p.beta.tauSq.samples <- sp.obj$p.beta.tauSq.samples
    n.samples <- nrow(p.beta.tauSq.samples)

    ##subsamples
    if(missing(end)){end <- n.samples}
    
    if(!is.numeric(start) || start >= n.samples)
      stop("error: invalid start")
    if(!is.numeric(end) || end > n.samples) 
      stop("error: invalid end")
    if(!is.numeric(thin) || thin >= n.samples) 
      stop("error: invalid thin")
       
    s.indx <- seq(start, end, by=as.integer(thin))
    
    p.beta.tauSq.samples <- p.beta.tauSq.samples[s.indx,]
    n.samples <- nrow(p.beta.tauSq.samples)
    
    ##check covars
    if(missing(pred.covars)){stop("error: pred.covars must be specified\n")}
    if(!any(is.data.frame(pred.covars), is.matrix(pred.covars))){stop("error: pred.covars must be a data.frame or matrix\n")}
    if(ncol(pred.covars) != ncol(X)){ stop(paste("error: pred.covars must have ",p," columns\n"))}

    out <- list()
    out$p.predictive.samples <- apply(p.beta.tauSq.samples, 1, function(s){rnorm(n, pred.covars%*%s[1:p], sqrt(s[p+1]))})
    
    return(out)
  }

  ##
  ##bayesGeostatExact
  ##
  if(obj.class == "bayesGeostatExact"){
    
    X <- sp.obj$X
    n <- sp.obj$n
    p <- sp.obj$p
    Y <- sp.obj$Y
    coords <- sp.obj$coords
    cov.model <- sp.obj$cov.model
    samples <- sp.obj$p.samples
    phi <- sp.obj$phi
    n.samples <- sp.obj$n.samples
    
    if(cov.model == "matern")
      nu <- sp.obj$nu
    
    ##check covars
    if(missing(pred.covars)){stop("error: pred.covars must be specified\n")}
    
    if(!any(is.data.frame(pred.covars), is.matrix(pred.covars))){stop("error: pred.covars must be a data.frame or matrix\n")}
    
    if(ncol(pred.covars) != ncol(X))
      stop(paste("error: pred.covars must have ",p," columns\n"))
    
    ##thin
    if(missing(end))
      end <- n.samples
    
    if(!is.numeric(start) || start >= n.samples)
      stop("error: invalid start")
    if(!is.numeric(end) || end > n.samples) 
      stop("error: invalid end")
    if(!is.numeric(thin) || thin >= n.samples) 
      stop("error: invalid thin")
    
    samples <- samples[seq(start, end, by=as.integer(thin)),]
    n.samples <- nrow(samples)
    
    ##get samples
    beta <- as.matrix(samples[,1:p])
    tau.sq <- samples[,"tau.sq"]
    sigma.sq <- samples[,"sigma.sq"]    
    
    ##make R
    D <- as.matrix(dist(coords))
    
    if(cov.model == "exponential"){
      R <- exp(-phi*D)
    }else if(cov.model == "matern"){
      R <- (D*phi)^nu/(2^(nu-1)*gamma(nu))*besselK(x=D*phi, nu=nu)
      diag(R) <- 1
    }else if(cov.model == "gaussian"){
      R <- exp(-1*((phi*D)^2))
    }else if(cov.model == "spherical"){
      R <- D
      R[TRUE] <- 1
      R[D > 0 & D < 1/phi] <- 1-1.5*phi*D[D > 0 & D <= 1/phi]+0.5*((phi*D[D > 0 & D <= 1/phi])^3)
      R[D >= 1/phi] <- 0   
    }else{
      stop("error: in spPredict, specified cov.model '",cov.model,"' is not a valid option; choose, from gaussian, exponential, matern, spherical.")
    }
    
    n.pred <- nrow(pred.coords)
    
    y.pred <- matrix(0, n.pred, n.samples)
    
    R.eigen <- eigen(R)
    R.vals <- R.eigen$values
    R.vecs <- R.eigen$vectors
    R.vects.t <- t(R.vecs)
    
    if(verbose)
      cat("Predicting ...\n")
    
    report <- 1
    
    ##for each pred point by each sample
    for(i in 1:n.pred){
      
      D.pred <- sqrt((pred.coords[i,1]-coords[,1])^2 + (pred.coords[i,2]-coords[,2])^2)
      
      if(cov.model == "exponential"){
        gamma <- exp(-phi*D.pred)
      }else if(cov.model == "matern"){
        gamma <- (D.pred*phi)^nu/(2^(nu-1)*gamma(nu))*besselK(x=D.pred*phi, nu=nu)
      }else if(cov.model == "gaussian"){
        gamma <- exp(-1*((phi*D.pred)^2))
      }else if(cov.model == "spherical"){
        gamma <- D.pred
        gamma[TRUE] <- 1
        gamma[D.pred > 0 & D.pred < 1/phi] <- 1-1.5*phi*D.pred[D.pred > 0 & D.pred <= 1/phi]+0.5*((phi*D.pred[D.pred > 0 & D.pred <= 1/phi])^3)
        gamma[D.pred >= 1/phi] <- 0   
      }else{
        stop("error: in spPredict, specified cov.model '",cov.model,"' is not a valid option; choose, from gaussian, exponential, matern, spherical.")
      }
      
      gamma <- as.matrix(gamma)
      
      for(s in 1:n.samples){
        
        R.inv <- R.vecs%*%diag(1/(R.vals+tau.sq[s]/sigma.sq[s]))%*%t(R.vecs)
        
        mu <- pred.covars[i,]%*%beta[s,]+t(gamma)%*%R.inv%*%(Y-X%*%beta[s,])
        S <- sigma.sq[s]*(1-t(gamma)%*%R.inv%*%gamma)+tau.sq[s]
        
        y.pred[i,s] <- rnorm(1, mu, sqrt(S))
        
      }

      
      if(verbose){
        if(report == 10){
          cat(paste("Percent complete: ",100*i/n.pred,"\n",sep=""))
          report <- 0
        }
        report <- report+1
      }
    }
    
    out <- list()
    out$p.predictive.samples <- y.pred
    
    return(out)
  }
 
  
  ##
  ##spatial model prediction
  ##
  if(obj.class %in% c("spGLM", "spMvGLM")){

    if(missing(pred.coords)){stop("error: pred.coords must be specified\n")}
    if(!any(is.data.frame(pred.coords), is.matrix(pred.coords))){stop("error: pred.coords must be a data.frame or matrix\n")}
    if(!ncol(pred.coords) == 2){stop("error: pred.coords must have two columns (assumed to be X, Y)\n")}
    
    family <- sp.obj$family
    X <- sp.obj$X
    Y <- sp.obj$Y
    p <- ncol(X)
    m <- 1 ##for spGLM
    if(obj.class == "spMvGLM"){m <- sp.obj$m}
    obs.coords <- sp.obj$coords
    n <- nrow(obs.coords)
    cov.model <- sp.obj$cov.model
    p.beta.theta.samples <- sp.obj$p.beta.theta.samples
    n.samples <- nrow(p.beta.theta.samples)
    is.pp <- sp.obj$is.pp
    q <- nrow(pred.coords)
    
    ##check covars
    if(missing(pred.covars)){stop("error: pred.covars must be specified\n")}
    if(!any(is.data.frame(pred.covars), is.matrix(pred.covars))){stop("error: pred.covars must be a data.frame or matrix\n")}
    if(ncol(pred.covars) != ncol(X)){ stop(paste("error: pred.covars must have ",p," columns\n"))}
    
    ##subsamples
    if(missing(end)){end <- n.samples}
    
    if(!is.numeric(start) || start >= n.samples)
      stop("error: invalid start")
    if(!is.numeric(end) || end > n.samples) 
      stop("error: invalid end")
    if(!is.numeric(thin) || thin >= n.samples) 
      stop("error: invalid thin")
       
    s.indx <- seq(start, end, by=as.integer(thin))

    p.beta.theta.samples <- t(p.beta.theta.samples[s.indx,,drop=FALSE])
    n.samples <- ncol(p.beta.theta.samples)

    knots.obs.D  <- NULL
    knots.D <- NULL
    knots.pred.D <- NULL
    obs.pred.D <- NULL
    obs.D <- NULL
    pred.D <- NULL
    p.w.samples <- NULL
    
    if(is.pp){
      knot.coords <- sp.obj$knot.coords
      g <- nrow(knot.coords)
      p.w.samples <- sp.obj$p.w.knots.samples[,s.indx,drop=FALSE]
      knots.D <- iDist(knot.coords)
      pred.knots.D <- iDist(pred.coords, knot.coords)
    }else{
      p.w.samples <- sp.obj$p.w.samples[,s.indx,drop=FALSE]
      obs.pred.D <- iDist(obs.coords, pred.coords)
      obs.D <- iDist(obs.coords)
    }
    
    storage.mode(X) <- "double"
    storage.mode(Y) <- "double"
    storage.mode(n) <- "integer"
    storage.mode(p) <- "integer"
    storage.mode(m) <- "integer"
    storage.mode(pred.covars) <- "double"
    storage.mode(q) <- "integer"
    storage.mode(p.beta.theta.samples) <- "double"
    storage.mode(n.samples) <- "integer"
    storage.mode(p.w.samples) <- "double"   
    storage.mode(verbose) <- "integer"
    storage.mode(n.report) <- "integer"
    
    if(is.pp){
      
      storage.mode(g) <- "integer"
      storage.mode(knots.D) <- "double"
      storage.mode(pred.knots.D) <- "double"

      out <- .Call("spPPMvGLMPredict", family, X, Y, n, m, g, p, pred.covars, q, knots.D, pred.knots.D, 
                   p.beta.theta.samples, p.w.samples, n.samples,
                   cov.model, verbose, n.report)
    }else{
      
      storage.mode(obs.pred.D) <- "double"
      storage.mode(obs.D) <- "double"
       
      out <- .Call("spMvGLMPredict", family, X, Y, n, m, p, pred.covars, q, obs.D, obs.pred.D, 
                   p.beta.theta.samples, p.w.samples, n.samples, cov.model,
                   verbose, n.report)
    }

    out
        
  }else if(obj.class == "spLM"){

    if(missing(pred.coords)){stop("error: pred.coords must be specified\n")}
    if(!any(is.data.frame(pred.coords), is.matrix(pred.coords))){stop("error: pred.coords must be a data.frame or matrix\n")}
    if(!ncol(pred.coords) == 2){stop("error: pred.coords must have two columns (assumed to be X, Y)\n")}
        
    X <- sp.obj$X
    Y <- sp.obj$Y
    p <- ncol(X)
    n <- nrow(X)
    obs.coords <- sp.obj$coords
    cov.model <- sp.obj$cov.model
    p.theta.samples <- sp.obj$p.theta.samples
    n.samples <- nrow(p.theta.samples)
    is.pp <- sp.obj$is.pp
    nugget <- sp.obj$nugget
    beta.prior <- sp.obj$beta.prior
    beta.Norm <- sp.obj$beta.Norm
    ##r.indx <- sp.obj$r.indx
    get.beta <- sp.obj$get.beta
    
    ##check covars
    if(missing(pred.covars)){stop("error: pred.covars must be specified\n")}
    if(!any(is.data.frame(pred.covars), is.matrix(pred.covars))){stop("error: pred.covars must be a data.frame or matrix\n")}
    if(ncol(pred.covars) != ncol(X)){ stop(paste("error: pred.covars must have ",p," columns\n"))}
    
    ##subsamples
    if(missing(end)){end <- n.samples}
    
    if(!is.numeric(start) || start >= n.samples)
      stop("error: invalid start")
    if(!is.numeric(end) || end > n.samples) 
      stop("error: invalid end")
    if(!is.numeric(thin) || thin >= n.samples) 
      stop("error: invalid thin")
    
    s.indx <- seq(start, end, by=as.integer(thin))
    
    beta <- NULL;
    
    if(is.pp){
      beta <- sp.obj$p.beta.samples[s.indx,,drop=FALSE]
      knot.coords <- sp.obj$knot.coords
      m <- nrow(knot.coords)
      modified.pp <- sp.obj$modified.pp
    }   
   
    p.theta.samples <- p.theta.samples[s.indx,,drop=FALSE]##note, for flat we could use the p.theta.recover.samples
    
    n.samples <- nrow(p.theta.samples)
    
    ##recover beta if needed
    if(!is.pp && beta.prior == "flat"){
      ## if(all(s.indx %in% r.indx)){
      ##   beta <- sp.obj$p.beta.samples[s.indx,,drop=FALSE]
      ## }else{
        beta <- spRecover(sp.obj, get.beta=TRUE, get.w=FALSE, start=start, end=end, thin=thin)$p.beta.recover.samples
      ## }
    }
    
    sigma.sq.indx <- 0; tau.sq.indx <- 0; phi.indx <- 0; nu.indx <- 0
    
    if(!nugget && cov.model != "matern"){
      sigma.sq.indx <- 0; phi.indx <- 1
    }else if(nugget && cov.model != "matern"){
      sigma.sq.indx <- 0; tau.sq.indx <- 1; phi.indx <- 2
    }else if(!nugget && cov.model == "matern"){
      sigma.sq.indx <- 0; phi.indx <- 1; nu.indx <- 2
    }else{
      sigma.sq.indx <- 0; tau.sq.indx <- 1; phi.indx <- 2; nu.indx <- 3
    }
    
    knots.obs.D  <- NULL
    knots.D <- NULL
    knots.pred.D <- NULL
    obs.pred.D <- NULL
    obs.D <- NULL
    
    if(is.pp){
      knots.obs.D <- iDist(knot.coords, obs.coords)
      knots.D <- iDist(knot.coords)
      knots.pred.D <- iDist(knot.coords, pred.coords)
    }else{
      obs.pred.D <- iDist(obs.coords, pred.coords)
      obs.D <- iDist(obs.coords)
    }
    
    n.pred <- nrow(pred.coords)
    
    storage.mode(X) <- "double"
    storage.mode(Y) <- "double"
    storage.mode(n) <- "integer"
    storage.mode(p) <- "integer"
    storage.mode(pred.covars) <- "double"
    storage.mode(n.pred) <- "integer"
    storage.mode(p.theta.samples) <- "double"
    storage.mode(n.samples) <- "integer"
    storage.mode(beta) <- "double"
    storage.mode(sigma.sq.indx) <- "integer"
    storage.mode(tau.sq.indx) <- "integer"
    storage.mode(phi.indx) <- "integer"
    storage.mode(nu.indx) <- "integer"
    storage.mode(verbose) <- "integer"
    storage.mode(n.report) <- "integer"
    
    if(is.pp){
      
      storage.mode(m) <- "integer"
      storage.mode(knots.D) <- "double"
      storage.mode(knots.obs.D) <- "double"
      storage.mode(knots.pred.D) <- "double"
      storage.mode(modified.pp) <- "integer"
      
      out <- .Call("spPPLMPredict", X, Y, n, p, m, pred.covars, n.pred,
                   p.theta.samples, n.samples,
                   beta, sigma.sq.indx, tau.sq.indx, phi.indx, nu.indx,
                   nugget, knots.D, knots.obs.D, knots.pred.D, cov.model, modified.pp, verbose, n.report)
    }else{
      
     storage.mode(obs.pred.D) <- "double"
     storage.mode(obs.D) <- "double"
     
     out <- .Call("spLMPredict", X, Y, n, p, pred.covars, n.pred,
                  p.theta.samples, n.samples,
                  beta.prior, beta.Norm, beta, sigma.sq.indx, tau.sq.indx, phi.indx, nu.indx,
                  obs.D, obs.pred.D, cov.model, nugget, verbose, n.report)
   }
    
    out
  
  }else if(obj.class == "spMvLM"){

    if(missing(pred.coords)){stop("error: pred.coords must be specified\n")}
    if(!any(is.data.frame(pred.coords), is.matrix(pred.coords))){stop("error: pred.coords must be a data.frame or matrix\n")}
    if(!ncol(pred.coords) == 2){stop("error: pred.coords must have two columns (assumed to be X, Y)\n")}
        
    X <- sp.obj$X
    Y <- sp.obj$Y
    p <- ncol(X)
    m <- sp.obj$m
    obs.coords <- sp.obj$coords
    n <- nrow(obs.coords)
    cov.model <- sp.obj$cov.model
    p.theta.samples <- sp.obj$p.theta.samples
    n.samples <- nrow(p.theta.samples)
    is.pp <- sp.obj$is.pp
    nugget <- sp.obj$nugget
    Psi.diag <- sp.obj$Psi.diag
    beta.prior <- sp.obj$beta.prior
    beta.Norm <- sp.obj$beta.Norm
    ##r.indx <- sp.obj$r.indx
    get.beta <- sp.obj$get.beta
    
    q <- nrow(pred.coords)
    
    ##check covars
    if(missing(pred.covars)){stop("error: pred.covars must be specified\n")}
    if(!any(is.data.frame(pred.covars), is.matrix(pred.covars))){stop("error: pred.covars must be a data.frame or matrix\n")}
    if(ncol(pred.covars) != ncol(X)){ stop(paste("error: pred.covars must have ",p," columns\n"))}
    
    ##subsamples
    if(missing(end)){end <- n.samples}
    
    if(!is.numeric(start) || start >= n.samples)
      stop("error: invalid start")
    if(!is.numeric(end) || end > n.samples) 
      stop("error: invalid end")
    if(!is.numeric(thin) || thin >= n.samples) 
      stop("error: invalid thin")
    
    s.indx <- seq(start, end, by=as.integer(thin))
    
    beta <- NULL
    
    if(is.pp){
      beta <- sp.obj$p.beta.samples[s.indx,,drop=FALSE]
      knot.coords <- sp.obj$knot.coords
      g <- nrow(knot.coords)
      modified.pp <- sp.obj$modified.pp
    }  
       
    p.theta.samples <- t(p.theta.samples[s.indx,,drop=FALSE])##note, for flat we could use the p.theta.recover.samples
    n.samples <- ncol(p.theta.samples)
    
    ##recover beta if needed (note, beta samples not needed for beta normal)
    if(!is.pp && beta.prior == "flat"){
      ## if(all(s.indx %in% r.indx)){
      ##   beta <- sp.obj$p.beta.samples[s.indx,,drop=FALSE]
      ## }else{
        beta <- spRecover(sp.obj, get.beta=TRUE, get.w=FALSE, start=start, end=end, thin=thin)$p.beta.recover.samples
      ## }
    }
        
    knots.obs.D  <- NULL
    knots.D <- NULL
    knots.pred.D <- NULL
    obs.pred.D <- NULL
    obs.D <- NULL
    
    if(is.pp){
      knots.obs.D <- iDist(knot.coords, obs.coords)
      knots.D <- iDist(knot.coords)
      knots.pred.D <- iDist(knot.coords, pred.coords)
    }else{
      obs.pred.D <- iDist(obs.coords, pred.coords)
      obs.D <- iDist(obs.coords)
    }
    
    storage.mode(X) <- "double"
    storage.mode(Y) <- "double"
    storage.mode(n) <- "integer"
    storage.mode(p) <- "integer"
    storage.mode(m) <- "integer"
    storage.mode(pred.covars) <- "double"
    storage.mode(q) <- "integer"
    storage.mode(p.theta.samples) <- "double"
    storage.mode(n.samples) <- "integer"
    storage.mode(beta) <- "double"
    storage.mode(nugget) <- "integer"
    storage.mode(Psi.diag) <- "integer"
    storage.mode(verbose) <- "integer"
    storage.mode(n.report) <- "integer"
    
    if(is.pp){
      
      storage.mode(g) <- "integer"
      storage.mode(knots.D) <- "double"
      storage.mode(knots.obs.D) <- "double"
      storage.mode(knots.pred.D) <- "double"
      storage.mode(modified.pp) <- "integer"

      out <- .Call("spPPMvLMPredict", X, Y, n, m, g, p, pred.covars, q,
                   knots.D, knots.obs.D, knots.pred.D, 
                   p.theta.samples, beta, n.samples,
                   nugget, Psi.diag, cov.model,
                   modified.pp, verbose, n.report)
    }else{
      
     storage.mode(obs.pred.D) <- "double"
     storage.mode(obs.D) <- "double"


     Z <- t(pred.covars)
     storage.mode(Z) <- "double"
     
     out <- .Call("spMvLMPredict", X, Y, n, m, p, Z, q, obs.D, obs.pred.D, 
                  p.theta.samples, beta, n.samples,
                  beta.prior, beta.Norm,
                  nugget, Psi.diag, cov.model,
                  verbose, n.report)
   }
    
    out
    
  }else if(obj.class == "spMisalignLM"){
    
    X <- sp.obj$X
    Y <- sp.obj$Y
    m <- sp.obj$m
    obs.coords <- sp.obj$coords
    misalign.p <- sp.obj$misalign.p
    misalign.n <- sp.obj$misalign.n
    cov.model <- sp.obj$cov.model
    p.theta.samples <- sp.obj$p.theta.samples
    n.samples <- nrow(p.theta.samples)
    nugget <- sp.obj$nugget
    beta.prior <- sp.obj$beta.prior
    beta.Norm <- sp.obj$beta.Norm
     
    ##get prediction covariates and coordinates
    if(missing(pred.covars)){stop("error: pred.covars must be specified\n")}
    if(any(!is.list(pred.covars), length(pred.covars) != m)){stop(paste("error: pred.covars must be a list of length have ",m,"\n"))}
    if(!all(unlist(lapply(pred.covars, is.matrix)))){stop("error: pred.covars must be a list of matrices\n")}

    misalign.n.pred <- unlist(lapply(pred.covars, nrow))
    misalign.p.pred <- unlist(lapply(pred.covars, ncol))
    X.pred <- do.call(adiag, pred.covars)

    if(!identical(misalign.p, misalign.p.pred)){stop("error: the number of columns in the model matrices do not match those listed in pred.covars\n")}

    if(missing(pred.coords)){stop("error: pred.coords must be specified\n")}
    if(any(!is.list(pred.coords), length(pred.coords) != m)){stop(paste("error: pred.coords must be a list of length have ",m,"\n"))}
    if(!all(unlist(lapply(pred.coords, is.matrix)))){stop("error: pred.coords must be a list of matrices\n")}

    if(!identical(unlist(lapply(pred.coords, nrow)), misalign.n.pred)){stop("error: number of rows in the pred.coords matrices do not match those listed in pred.covars\n")}
    if(any(unlist(lapply(pred.coords, ncol)) != 2)){stop("error: all matrices in pred.coords must have two columns\n")}
    
    pred.coords <- as.matrix(do.call(rbind, pred.coords))
 
    ##subsamples
    if(missing(end)){end <- n.samples}
    
    if(!is.numeric(start) || start >= n.samples)
      stop("error: invalid start")
    if(!is.numeric(end) || end > n.samples) 
      stop("error: invalid end")
    if(!is.numeric(thin) || thin >= n.samples) 
      stop("error: invalid thin")
    
    s.indx <- seq(start, end, by=as.integer(thin))
           
    p.theta.samples <- t(p.theta.samples[s.indx,,drop=FALSE])##note, for flat we could use the p.theta.recover.samples
    n.samples <- ncol(p.theta.samples)
    
    ##recover beta if needed (note, beta samples not needed for beta normal)
    beta <- NULL
    if(beta.prior == "flat"){
      beta <- t(spRecover(sp.obj, get.beta=TRUE, get.w=FALSE, start=start, end=end, thin=thin)$p.beta.recover.samples)
    }
        
    pred.obs.D <- iDist(pred.coords, obs.coords)
    obs.D <- iDist(obs.coords)
        
    storage.mode(X) <- "double"
    storage.mode(Y) <- "double"
    storage.mode(m) <- "integer"
    storage.mode(misalign.n) <- "integer"
    storage.mode(misalign.p) <- "integer"
    storage.mode(p.theta.samples) <- "double"
    storage.mode(n.samples) <- "integer"
    storage.mode(beta) <- "double"
    storage.mode(nugget) <- "integer"
    storage.mode(verbose) <- "integer"
    storage.mode(n.report) <- "integer"
    storage.mode(misalign.n.pred) <- "integer"
    storage.mode(misalign.p.pred) <- "integer"
    storage.mode(pred.obs.D) <- "double"
    storage.mode(obs.D) <- "double"
    
    Z <- X.pred
    storage.mode(Z) <- "double"
    
     out <- .Call("spMisalignPredict", Y, X, m, misalign.n, misalign.p,
                  Z, misalign.n.pred, misalign.p.pred,
                  obs.D, pred.obs.D, 
                  p.theta.samples, beta, n.samples,
                  beta.prior, beta.Norm,
                  nugget, cov.model,
                  verbose, n.report)


    out$p.beta.samples.recover <- mcmc(t(beta))
    out
    
  }else if(obj.class == "spMisalignGLM"){

    family <- sp.obj$family
    X <- sp.obj$X
    Y <- sp.obj$Y
    m <- sp.obj$m
    obs.coords <- sp.obj$coords
    misalign.p <- sp.obj$misalign.p
    misalign.n <- sp.obj$misalign.n
    cov.model <- sp.obj$cov.model
    p.w.samples <- sp.obj$p.w.samples
    p.beta.theta.samples <- sp.obj$p.beta.theta.samples
    n.samples <- nrow(p.beta.theta.samples)
     
    ##get prediction covariates and coordinates
    if(missing(pred.covars)){stop("error: pred.covars must be specified\n")}
    if(any(!is.list(pred.covars), length(pred.covars) != m)){stop(paste("error: pred.covars must be a list of length have ",m,"\n"))}
    if(!all(unlist(lapply(pred.covars, is.matrix)))){stop("error: pred.covars must be a list of matrices\n")}

    misalign.n.pred <- unlist(lapply(pred.covars, nrow))
    misalign.p.pred <- unlist(lapply(pred.covars, ncol))
    X.pred <- do.call(adiag, pred.covars)

    if(!identical(misalign.p, misalign.p.pred)){stop("error: the number of columns in the model matrices do not match those listed in pred.covars\n")}

    if(missing(pred.coords)){stop("error: pred.coords must be specified\n")}
    if(any(!is.list(pred.coords), length(pred.coords) != m)){stop(paste("error: pred.coords must be a list of length have ",m,"\n"))}
    if(!all(unlist(lapply(pred.coords, is.matrix)))){stop("error: pred.coords must be a list of matrices\n")}

    if(!identical(unlist(lapply(pred.coords, nrow)), misalign.n.pred)){stop("error: number of rows in the pred.coords matrices do not match those listed in pred.covars\n")}
    if(any(unlist(lapply(pred.coords, ncol)) != 2)){stop("error: all matrices in pred.coords must have two columns\n")}
    
    pred.coords <- as.matrix(do.call(rbind, pred.coords))
 
    ##subsamples
    if(missing(end)){end <- n.samples}
    
    if(!is.numeric(start) || start >= n.samples)
      stop("error: invalid start")
    if(!is.numeric(end) || end > n.samples) 
      stop("error: invalid end")
    if(!is.numeric(thin) || thin >= n.samples) 
      stop("error: invalid thin")
    
    s.indx <- seq(start, end, by=as.integer(thin))

    p.w.samples <- p.w.samples[,s.indx,drop=FALSE]
 
    p.samples <- t(p.beta.theta.samples[s.indx,,drop=FALSE])
    n.samples <- ncol(p.samples)

    pred.obs.D <- iDist(pred.coords, obs.coords)
    obs.D <- iDist(obs.coords)
        
    storage.mode(X) <- "double"
    storage.mode(Y) <- "double"
    storage.mode(m) <- "integer"
    storage.mode(misalign.n) <- "integer"
    storage.mode(misalign.p) <- "integer"
    storage.mode(p.samples) <- "double"
    storage.mode(p.w.samples) <- "double"
    storage.mode(n.samples) <- "integer"
    storage.mode(verbose) <- "integer"
    storage.mode(n.report) <- "integer"
    storage.mode(misalign.n.pred) <- "integer"
    storage.mode(misalign.p.pred) <- "integer"
    storage.mode(pred.obs.D) <- "double"
    storage.mode(obs.D) <- "double"
    
    Z <- X.pred
    storage.mode(Z) <- "double"

    out <- .Call("spMisalignGLMPredict", family, Y, X, m, misalign.n, misalign.p,
                 Z, misalign.n.pred, misalign.p.pred,
                 obs.D, pred.obs.D, 
                 p.samples, p.w.samples, n.samples,
                 cov.model,
                 verbose, n.report)
    out
    
  }else{
    stop("error: wrong class\n")
  }
  
}
