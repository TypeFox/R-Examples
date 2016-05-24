spRecover <- function(sp.obj, get.beta=TRUE, get.w=TRUE, start=1, end, thin=1, verbose=TRUE, n.report=100, ...){
  
  ####################################################
  ##Check for unused args
  ####################################################
  formal.args <- names(formals(sys.function(sys.parent())))
  elip.args <- names(list(...))
  for(i in elip.args){
    if(! i %in% formal.args)
      warning("'",i, "' is not an argument")
  }
  
  if(missing(sp.obj)){stop("error: spRecover expects sp.obj\n")}
  if(!class(sp.obj) %in% c("spLM", "spMvLM", "spMisalignLM")){
    stop("error: requires an output object of class spLM, spMvLM, or spMisalignLM\n")
  }

  obj.class <- class(sp.obj)

  ##need beta to get.w
  if(get.w){
    get.beta <- TRUE
  }
  
  if(obj.class == "spLM"){
    
    Y <-sp.obj$Y
    X <- sp.obj$X
    p <- ncol(X)
    n <- nrow(X)
    coords <- sp.obj$coords
    p.theta.samples <- sp.obj$p.theta.samples
    n.samples <- nrow(p.theta.samples)
    n.params <- ncol(p.theta.samples)
    cov.model <- sp.obj$cov.model
    nugget <- sp.obj$nugget
    beta.prior <- sp.obj$beta.prior
    beta.Norm <- sp.obj$beta.Norm
    x.names <- sp.obj$x.names
    is.pp <- sp.obj$is.pp
    modified.pp <- sp.obj$modified.pp
  
    coords.D <- iDist(coords)
    
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
    
    if(missing(end))
      end <- n.samples
    
    if(!is.numeric(start) || start >= n.samples)
      stop("error: invalid start")
    if(!is.numeric(end) || end > n.samples) 
      stop("error: invalid end")
    if(!is.numeric(thin) || thin >= n.samples) 
      stop("error: invalid thin")

    p.theta.samples <- p.theta.samples[seq(start, end, by=as.integer(thin)),,drop=FALSE]
    n.samples <- nrow(p.theta.samples)
    
    ##if pp
    p.beta.samples <- NULL
    knots.D <- NULL
    knots.obs.D <- NULL
    m <- NULL
    if(is.pp){
      p.beta.samples <- sp.obj$p.beta.samples[seq(start, end, by=as.integer(thin)),,drop=FALSE]
      knots.D <- iDist(sp.obj$knot.coords)
      knots.obs.D <- iDist(sp.obj$knot.coords, coords)
      m <- nrow(knots.D)
    }
    
    storage.mode(Y) <- "double"
    storage.mode(X) <- "double"
    storage.mode(p) <- "integer"
    storage.mode(n) <- "integer"
    storage.mode(coords.D) <- "double"
    storage.mode(p.theta.samples) <- "double"
    storage.mode(n.samples) <- "integer"
    storage.mode(n.params) <- "integer"
    storage.mode(sigma.sq.indx) <- "integer"
    storage.mode(tau.sq.indx) <- "integer"
    storage.mode(phi.indx) <- "integer"
    storage.mode(nu.indx) <- "integer"
    storage.mode(nugget) <- "integer"
    storage.mode(get.beta) <- "integer"
    storage.mode(get.w) <- "integer"
    storage.mode(verbose) <- "integer"
    storage.mode(n.report) <- "integer"
    
    storage.mode(m) <- "integer"
    storage.mode(is.pp) <- "integer"
    storage.mode(modified.pp) <- "integer"
    storage.mode(knots.D) <- "double"
    storage.mode(knots.obs.D) <- "double"
    storage.mode(p.beta.samples) <- "double"
    
    if(is.pp){
      out <- .Call("spPPLMRecover", X, Y, n, p, m, 
                   p.theta.samples, n.samples, 
                   p.beta.samples, sigma.sq.indx, tau.sq.indx, phi.indx, nu.indx,
                   nugget, knots.D, knots.obs.D, cov.model, modified.pp, verbose, n.report)

      sp.obj$p.beta.recover.samples <- mcmc(p.beta.samples)
      sp.obj$p.theta.recover.samples <- mcmc(p.theta.samples)
      
      sp.obj$p.w.recover.samples <- out$p.w.samples
      sp.obj$p.wStr.recover.samples <- out$p.wStr.samples
      
      
    }else{
      out <- .Call("spLMRecover", Y, X, p, n, coords.D,
                   p.theta.samples, n.samples, 
                   sigma.sq.indx, tau.sq.indx, phi.indx, nu.indx,
                   beta.prior, beta.Norm, 	   
                   nugget, cov.model,
                   get.beta, get.w,
                   verbose, n.report)
      
      rownames(out$p.beta.samples) <- x.names
      sp.obj$p.beta.recover.samples <- mcmc(t(out$p.beta.samples))
      sp.obj$p.theta.recover.samples <- mcmc(p.theta.samples)
      
      if(get.w){
        sp.obj$p.w.recover.samples <- out$p.w.samples
      }
           
    }
        
    class(sp.obj) <- "spLM"

    sp.obj
    
  }else if(obj.class == "spMvLM"){
    
    Y <-sp.obj$Y
    X <- sp.obj$X
    p <- ncol(X)
    m <- sp.obj$m
    coords <- sp.obj$coords
    n <- nrow(coords)
    p.theta.samples <- sp.obj$p.theta.samples
    n.samples <- nrow(p.theta.samples)
    cov.model <- sp.obj$cov.model
    nugget <- sp.obj$nugget
    Psi.diag <- sp.obj$Psi.diag
    beta.prior <- sp.obj$beta.prior
    beta.Norm <- sp.obj$beta.Norm
    x.names <- sp.obj$x.names
    is.pp <- sp.obj$is.pp
    modified.pp <- sp.obj$modified.pp

    coords.D <- iDist(coords)
    
    if(missing(end))
      end <- n.samples
    
    if(!is.numeric(start) || start >= n.samples)
      stop("error: invalid start")
    if(!is.numeric(end) || end > n.samples) 
      stop("error: invalid end")
    if(!is.numeric(thin) || thin >= n.samples) 
      stop("error: invalid thin")
    
    p.theta.samples <- t(p.theta.samples[seq(start, end, by=as.integer(thin)),,drop=FALSE])
    n.samples <- ncol(p.theta.samples)

    ##if pp
    p.beta.samples <- NULL
    knots.D <- NULL
    knots.obs.D <- NULL
    if(is.pp){
      p.beta.samples <- sp.obj$p.beta.samples[seq(start, end, by=as.integer(thin)),,drop=FALSE]
      knots.D <- iDist(sp.obj$knot.coords)
      knots.obs.D <- iDist(sp.obj$knot.coords, coords)
      g <- nrow(knots.D)
    }
    
    storage.mode(Y) <- "double"
    storage.mode(X) <- "double"
    storage.mode(p) <- "integer"
    storage.mode(n) <- "integer"
    storage.mode(m) <- "integer"
    storage.mode(coords.D) <- "double"
    storage.mode(p.theta.samples) <- "double"
    storage.mode(n.samples) <- "integer"
    storage.mode(nugget) <- "integer"
    storage.mode(Psi.diag) <- "integer"
    storage.mode(get.beta) <- "integer"
    storage.mode(get.w) <- "integer"
    storage.mode(verbose) <- "integer"
    storage.mode(n.report) <- "integer"
    
    if(is.pp){
      storage.mode(g) <- "integer"
      storage.mode(knots.D) <- "double"
      storage.mode(knots.obs.D) <- "double"
      storage.mode(modified.pp) <- "integer"
      
      out <- .Call("spPPMvLMRecover", X, Y, n, m, g, p,
                   knots.D, knots.obs.D,
                   p.theta.samples, p.beta.samples, n.samples,   
                   nugget, Psi.diag, cov.model,
                   modified.pp, verbose, n.report)

      sp.obj$p.beta.recover.samples <- mcmc(p.beta.samples)
      sp.obj$p.theta.recover.samples <- mcmc(t(p.theta.samples))
      
      sp.obj$p.w.recover.samples <- out$p.w.samples
      sp.obj$p.wStr.recover.samples <- out$p.wStr.samples
      
    }else{
      
      out <- .Call("spMvLMRecover", Y, X, p, n, m, coords.D,
                   p.theta.samples, n.samples,
                   beta.prior, beta.Norm, 	   
                   nugget, Psi.diag, cov.model,
                   get.beta, get.w,
                   verbose, n.report)
      
      rownames(out$p.beta.samples) <- x.names
      sp.obj$p.beta.recover.samples <- mcmc(t(out$p.beta.samples))
      sp.obj$p.theta.recover.samples <- mcmc(t(p.theta.samples))
      
      if(get.w){
        sp.obj$p.w.recover.samples <- out$p.w.samples
      }
      
    }
    
    class(sp.obj) <- "spMvLM"

    sp.obj
    
  }else if(obj.class == "spMisalignLM"){
    
    Y <-sp.obj$Y
    X <- sp.obj$X
    m <- sp.obj$m
    coords <- sp.obj$coords
    misalign.p <- sp.obj$misalign.p
    misalign.n <- sp.obj$misalign.n
    p.theta.samples <- sp.obj$p.theta.samples
    n.samples <- nrow(p.theta.samples)
    cov.model <- sp.obj$cov.model
    nugget <- sp.obj$nugget
    beta.prior <- sp.obj$beta.prior
    beta.Norm <- sp.obj$beta.Norm
    x.names <- sp.obj$x.names
    
    coords.D <- iDist(coords)
    
    if(missing(end))
      end <- n.samples
    
    if(!is.numeric(start) || start >= n.samples)
      stop("error: invalid start")
    if(!is.numeric(end) || end > n.samples) 
      stop("error: invalid end")
    if(!is.numeric(thin) || thin >= n.samples) 
      stop("error: invalid thin")
    
    p.theta.samples <- t(p.theta.samples[seq(start, end, by=as.integer(thin)),,drop=FALSE])
    n.samples <- ncol(p.theta.samples)

    storage.mode(Y) <- "double"
    storage.mode(X) <- "double"
    storage.mode(misalign.p) <- "integer"
    storage.mode(misalign.n) <- "integer"
    storage.mode(m) <- "integer"
    storage.mode(coords.D) <- "double"
    storage.mode(p.theta.samples) <- "double"
    storage.mode(n.samples) <- "integer"
    storage.mode(nugget) <- "integer"
    storage.mode(get.beta) <- "integer"
    storage.mode(get.w) <- "integer"
    storage.mode(verbose) <- "integer"
    storage.mode(n.report) <- "integer"
    
    out <- .Call("spMisalignRecover", Y, X, misalign.p, misalign.n, m, coords.D,
                 p.theta.samples, n.samples,
                 beta.prior, beta.Norm, 	   
                 nugget, cov.model,
                 get.beta, get.w,
                 verbose, n.report)
    
    rownames(out$p.beta.samples) <- x.names
    sp.obj$p.beta.recover.samples <- mcmc(t(out$p.beta.samples))
    sp.obj$p.theta.recover.samples <- mcmc(t(p.theta.samples))
    
    if(get.w){
      sp.obj$p.w.recover.samples <- out$p.w.samples
    }
    
    class(sp.obj) <- "spMisalignLM"
    
    sp.obj
    
  }else{
    stop("error: wrong class\n")
  }
  
}
