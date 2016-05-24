spDiag <- function(sp.obj, start=1, end, thin=1, verbose=TRUE, n.report=100, ...){
  
  ####################################################
  ##Check for unused args, thin, etc.
  ####################################################
  formal.args <- names(formals(sys.function(sys.parent())))
  elip.args <- names(list(...))
  for(i in elip.args){
    if(! i %in% formal.args)
      warning("'",i, "' is not an argument")
  }

  if(missing(sp.obj)){stop("error: spDiag expects sp.obj\n")}
  if(!class(sp.obj) %in% c("spLM","spMvLM", "spGLM", "spMvGLM","bayesLMRef","nonSpGLM","nonSpMvGLM")){
    stop("error: spDiag requires an output object of class spLM, spMvLM, spGLM, spMvGLM, bayesLMRef, nonSpGLM, or nonSpMvGLM\n")}
  if(!is.logical(verbose)){stop("error: verbose must be of type logical\n")}
  
  rmvn <- function(n, mu=0, V = matrix(1)){
    p <- length(mu)
    if(any(is.na(match(dim(V),p))))
      stop("Dimension problem!")
    D <- chol(V)
    t(matrix(rnorm(n*p), ncol=p) %*% D + rep(mu,rep(n,p)))
  }
    
  lower2full <- function(x, m){
    A <- matrix(0, m, m)
    A[lower.tri(A, diag=TRUE)] <- x
    A[upper.tri(A, diag=FALSE)] <- t(A)[upper.tri(A, diag=FALSE)]
    A
  }

  GP <- function(Y.rep, y){
    mu.rep <- apply(Y.rep, 1, mean)
    var.rep <- apply(Y.rep, 1, var)
    G <- sum((y-mu.rep)^2)
    P <- sum(var.rep)
    D <- G+P

    GPD <- matrix(0,3,1)
    rownames(GPD) <- c("G","P","D")
    colnames(GPD) <- "value"
    GPD[1,1] <- G
    GPD[2,1] <- P
    GPD[3,1] <- D
    GPD
  }
  
  GR <- function(Y.rep, y){
    mu.rep <- apply(Y.rep, 1, mean)
    var.rep <- apply(Y.rep, 1, var)
    GRS <- -sum(((y-mu.rep)/sqrt(var.rep))^2) - sum(log(var.rep))
    GRS
  }
    
  if(class(sp.obj) == "bayesLMRef"){

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
    
    beta <- p.beta.tauSq.samples[s.indx,1:p,drop=FALSE]
    tau.sq <- p.beta.tauSq.samples[s.indx,p+1,drop=FALSE]
    
    beta.mu <- apply(beta, 2, mean)
    tau.sq.mu <- mean(tau.sq)

    n.samples <- nrow(beta)

    status <- 0
    
    d <- rep(0, n.samples)
    DIC <- matrix(0,4,1)
    
    Y.rep <- matrix(0, n, n.samples)
    
    if(verbose){
      cat("-------------------------------------------------\n\t\tCalculating scores\n-------------------------------------------------\n")
    }
       
    for(s in 1:n.samples){
      
      Q <- Y-X%*%beta[s,]
      
      d[s] <- n*log(tau.sq[s])+(t(Q)%*%Q)/tau.sq[s]
      
      status <- 0
      
      ##GP
      mu <- X%*%beta[s,]
      Y.rep[,s] <- rnorm(n, mu, sqrt(tau.sq[s]))
      
      if(verbose){
        if(status == n.report){
          cat(paste("Sampled: ",s," of ",n.samples,", ",round(100*s/n.samples,2),"%\n", sep=""))
          status <- 0
        }
        status <- status+1
      }
      
    }
    
    d.bar <- mean(d)
    
    ##Get d.bar.omega
    Q <- Y-X%*%beta.mu
    
    d.bar.omega <- n*log(tau.sq.mu)+(t(Q)%*%Q)/tau.sq.mu
    
    pd <- d.bar - d.bar.omega
    dic <- d.bar + pd
    
    rownames(DIC) <- c("bar.D", "D.bar.Omega", "pD", "DIC")
    colnames(DIC) <- c("value")
    DIC[1,1] <- d.bar
    DIC[2,1] <- d.bar.omega
    DIC[3,1] <- pd
    DIC[4,1] <- dic
    
    out <- list()
    out$DIC <- DIC
    out$GP <- GP(Y.rep,Y)
    out$GRS <- GR(Y.rep,Y)
    
    return(out)
  }
    
  if(class(sp.obj) %in% c("nonSpGLM", "nonSpMvGLM")){

    X <- sp.obj$X
    Y <- sp.obj$Y
    family <- sp.obj$family
    weights <- sp.obj$weights

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
    
    beta <- p.beta.samples[s.indx,,drop=FALSE]
   
    beta.mu <- apply(beta, 2, mean)

    n.samples <- nrow(beta)

    status <- 0
    
    d <- rep(0, n.samples)
    DIC <- matrix(0,4,1)

    if(verbose){
      cat("-------------------------------------------------\n\t\tCalculating scores\n-------------------------------------------------\n")
    }
    
    for(s in 1:n.samples){
      
      if(family == "poisson"){
        d[s] <- -2*sum(Y*(X%*%beta[s,])-exp(X%*%beta[s,])+log(weights))
      }else if(family == "binomial"){
        pp <- 1/(1+exp(-X%*%beta[s,]))
        d[s] <- -2*sum(Y*log(pp)+(weights-Y)*log(1-pp))
      }else{
        stop("error: family is misspecified")
      }

    }
    
    d.bar <- mean(d)
    
    if(family == "poisson"){
      d.bar.omega <- -2*sum(Y*(X%*%beta.mu)-exp(X%*%beta.mu)+log(weights))
    }else if(family == "binomial"){
      pp <- 1/(1+exp(-X%*%beta.mu))
      d.bar.omega <- -2*sum(Y*log(pp)+(weights-Y)*log(1-pp))
    }else{
      stop("error: family is misspecified")
    }
    
    pd <- d.bar - d.bar.omega
    dic <- d.bar + pd
    
    DIC <- matrix(0,4,1)
    
    rownames(DIC) <- c("bar.D", "D.bar.Omega", "pD", "DIC")
    colnames(DIC) <- c("value")
    DIC[1,1] <- d.bar
    DIC[2,1] <- d.bar.omega
    DIC[3,1] <- pd
    DIC[4,1] <- dic
    
    out <- list()    
    out$DIC <- DIC
    
    return(out)
  }

  
  if(class(sp.obj) %in% c("spLM", "spMvLM")){
    
    ##recover beta and/or w (burn-in and thinning before spDiag)
    is.pp <- sp.obj$is.pp
    
    if(!all(c("p.beta.recover.samples", "p.theta.recover.samples", "p.w.recover.samples") %in% names(sp.obj))){
      if(!is.pp){
        stop("error: run spRecover(sp.obj, get.beta=T, get.w=T, ...) before calling spDiag")
      }else{
        stop("error: run spRecover(sp.obj, get.w=T, ...) before calling spDiag")
      }
    }

    theta <- sp.obj$p.theta.recover.samples
    n.samples <- nrow(theta)
    
    beta <- sp.obj$p.beta.recover.samples
    
    w <- sp.obj$p.w.recover.samples
    
    if(class(sp.obj) == "spMvLM"){
      
      if(!sp.obj$nugget){
        stop("DIC cannot be computed for a no nugget model.")
      }
      
      Y <- sp.obj$Y
      X <- sp.obj$X
      m <- sp.obj$m ##number of outcomes
      p <- sp.obj$p
      n <- nrow(sp.obj$coords)
      nltr <- m*(m+1)/2
      Psi.diag <- sp.obj$Psi.diag
      modified.pp <- sp.obj$modified.pp
      cov.model <- sp.obj$cov.model
      
      if(Psi.diag){
        Psi.names <- paste("Psi[",1:m,",",1:m,"]",sep="")
      }else{
        Psi.names <- paste("Psi[",matrix(apply(cbind(expand.grid(1:m,1:m)), 1, function(x) paste(x, collapse=",")),m,m)[lower.tri(matrix(0,m,m), diag=TRUE)],"]",sep="")
      }
      Psi <- as.matrix(theta[,Psi.names])
      Psi.mu <- apply(Psi, 2, mean)
      
      K.names <- paste("K[",matrix(apply(cbind(expand.grid(1:m,1:m)), 1, function(x) paste(x, collapse=",")),m,m)[lower.tri(matrix(0,m,m), diag=TRUE)],"]",sep="")  
      K <- as.matrix(theta[,K.names])
      K.mu <- apply(K, 2, mean)
      
      phi.names <- paste("phi[",1:m,"]",sep="")
      phi <- as.matrix(theta[,phi.names])
      phi.mu <- apply(phi, 2, mean);
      
      nu.names <- paste("nu[",1:m,"]",sep="")
      nu <- nu.mu <- NULL;
      if(cov.model=="matern"){
        nu <- as.matrix(theta[,nu.names])
        nu.mu <- apply(nu, 2, mean);
      }
      
      beta.mu <- apply(beta, 2, mean)
      
      w.mu <- as.matrix(apply(w,1,mean))
      
      status <- 0
      
      d <- rep(0, n.samples)
      DIC <- matrix(0,4,1)
      
      Y.rep <- matrix(0, n*m, n.samples)
      
      if(verbose)
        cat("-------------------------------------------------\n\t\tCalculating scores\n-------------------------------------------------\n")
      
      ##if non-pp or non-modified pp
      if(is.pp && modified.pp){
        
        knots.D <- iDist(sp.obj$knot.coords)
        q <- nrow(knots.D) ##number of knots
        obs.D <- iDist(sp.obj$coords)
        obs.knots.D <- iDist(sp.obj$coords, sp.obj$knot.coords)
        cov.model <- sp.obj$cov.model
        
        ##note, C_e_r is the mnxmn mxm block diag covariance matrix marginalized over \tile{\eps} use to compute GP
        C.eps.r <- matrix(0, m, n*m)
        
        for(s in 1:n.samples){
          
          KK <- lower2full(K[s,], m)
          
          if(Psi.diag){
            PP <- diag(Psi[s,], m)
          }else{
            PP <- lower2full(Psi[s,], m)
          }
          
          Q <- Y-X%*%beta[s,]-w[,s]
          
          theta.1 <- phi[s,]
          theta.2 <- NULL
          if(cov.model=="matern"){
            theta.2 <- nu[s,]
          }
          
          storage.mode(Q) <- "double"
          storage.mode(KK) <- "double"
          storage.mode(PP) <- "double"
          storage.mode(theta.1) <- "double"
          storage.mode(theta.2) <- "double"
          storage.mode(knots.D) <- "double"
          storage.mode(obs.knots.D) <- "double"
          storage.mode(q) <- "integer"
          storage.mode(n) <- "integer"
          storage.mode(m) <- "integer"
          
          d[s] <- .Call("spMPPMvDIC", Q, knots.D, obs.knots.D, n, m, q, PP, KK, theta.1, theta.2, cov.model, C.eps.r)
          
          ##GP
          mu <- X%*%beta[s,]+w[,s]
          for(i in 1:n){
            Y.rep[((i-1)*m+1):((i-1)*m+m),s] <- rmvn(1, mu[((i-1)*m+1):((i-1)*m+m)], C.eps.r[,((i-1)*m+1):((i-1)*m+m)])
          }
          
          if(verbose){
            if(status == n.report){
              cat(paste("Sampled: ",s," of ",n.samples,", ",round(100*s/n.samples,2),"%\n", sep=""))
              status <- 0
            }
            status <- status+1
          }
          
        }
        
        d.bar <- mean(d)
        
        ##Get d.bar.omega  
        KK <- lower2full(K.mu, m)
        
        if(Psi.diag){
          PP <- diag(Psi.mu, m)
        }else{
          PP <- lower2full(Psi.mu, m)
        }
        
        Q <- Y-X%*%beta.mu-w.mu
        
        storage.mode(KK) <- "double"
        storage.mode(PP) <- "double"
        storage.mode(phi.mu) <- "double"
        storage.mode(nu.mu) <- "double"
        storage.mode(Q) <- "double"
        
        d.bar.omega <- .Call("spMPPMvDIC", Q, knots.D, obs.knots.D, n, m, q, PP, KK, phi.mu, nu.mu, cov.model, C.eps.r)
        
        pd <- d.bar - d.bar.omega
        dic <- d.bar + pd
        
        rownames(DIC) <- c("bar.D", "D.bar.Omega", "pD", "DIC")
        colnames(DIC) <- c("value")
        DIC[1,1] <- d.bar
        DIC[2,1] <- d.bar.omega
        DIC[3,1] <- pd
        DIC[4,1] <- dic
        
        out <- list()
        out$DIC <- DIC
        out$GP <- GP(Y.rep,Y)
        out$GRS <- GR(Y.rep,Y)
        
        out
        
      }else{
        
        for(s in 1:n.samples){
          
          if(Psi.diag){
            PP <- diag(Psi[s,], m)
          }else{
            PP <- lower2full(Psi[s,], m)
          }
          
          Q <- Y-X%*%beta[s,]-w[,s]
          
          d[s] <- n*determinant(PP)$modulus+t(Q)%*%(diag(n)%x%chol2inv(chol(PP)))%*%Q
          
          ##GP
          mu <- X%*%beta[s,]+w[,s]
          for(i in 1:n){
            Y.rep[((i-1)*m+1):((i-1)*m+m),s] <- rmvn(1, mu[((i-1)*m+1):((i-1)*m+m)], PP)
          }
          
          if(verbose){
            if(status == n.report){
              cat(paste("Sampled: ",s," of ",n.samples,", ",round(100*s/n.samples,2),"%\n", sep=""))
              status <- 0
            }
            status <- status+1
          }
          
        }
        
        d.bar <- mean(d)
        
        ##Get d.bar.omega
        if(Psi.diag){
          PP <- diag(Psi.mu, m)
        }else{
          PP <- lower2full(Psi.mu, m)
        }
        
        Q <- Y-X%*%beta.mu-w.mu
        
        d.bar.omega <- n*determinant(PP)$modulus+t(Q)%*%(diag(n)%x%chol2inv(chol(PP)))%*%Q
        
        pd <- d.bar - d.bar.omega
        dic <- d.bar + pd
        
        rownames(DIC) <- c("bar.D", "D.bar.Omega", "pD", "DIC")
        colnames(DIC) <- c("value")
        DIC[1,1] <- d.bar
        DIC[2,1] <- d.bar.omega
        DIC[3,1] <- pd
        DIC[4,1] <- dic
        
        out <- list()
        out$DIC <- DIC
        out$GP <- GP(Y.rep,Y)
        out$GRS <- GR(Y.rep,Y)
        
        out      
      }            
      
    }else if(class(sp.obj) == "spLM"){
      
      if(!sp.obj$nugget){
        stop("DIC cannot be computed for a no nugget model.")
      }
      
      is.pp <- sp.obj$is.pp
      modified.pp <- sp.obj$modified.pp
      
      Y <- sp.obj$Y
      X <- sp.obj$X
      n <- nrow(X)
      
      ##get samples
      tau.sq <- theta[,"tau.sq"]
      
      ##get sample means
      beta.mu <- apply(beta, 2, mean)
      tau.sq.mu <- mean(tau.sq)
      w.mu <- as.matrix(apply(w,1,mean))
      
      status <- 0
      
      d <- rep(0, n.samples)
      DIC <- matrix(0,4,1)
      
      Y.rep <- matrix(0, n, n.samples)

      if(verbose){
        cat("-------------------------------------------------\n\t\tCalculating scores\n-------------------------------------------------\n")
      }
      
      ##if non-pp or non-modified pp
      if(is.pp && modified.pp){
        
        knots.D <- iDist(sp.obj$knot.coords)
        obs.knots.D <- iDist(sp.obj$coords, sp.obj$knot.coords)
        cov.model <- sp.obj$cov.model
        
        sigma.sq <- as.matrix(theta[,"sigma.sq"])
        phi <- as.matrix(theta[,"phi"])
        
        sigma.sq.mu <- mean(sigma.sq) 
        phi.mu <- mean(phi);
        
        if(cov.model=="matern"){
          nu <- as.matrix(theta[,"nu"])
          nu.mu <- mean(nu);
        }
        
        for(s in 1:n.samples){
          
          Q <- Y-X%*%beta[s,]-w[,s]
          
          if(cov.model != "matern"){
            ct <- sigma.sq[s]*spCor(obs.knots.D, cov.model, phi[s])
            C.str <- sigma.sq[s]*spCor(knots.D, cov.model, phi[s])
          }else{
            ct <- sigma.sq[s]*spCor(obs.knots.D, cov.model, phi[s], nu[s])
            C.str <- sigma.sq[s]*spCor(knots.D, cov.model, phi[s], nu[s])
          }
          
          
          ##ct C^{*-1} c
          C <- ct%*%chol2inv(chol(C.str))%*%t(ct)
          
          e <- tau.sq[s]+sigma.sq[s]-diag(C)
          
          d[s] <- sum(log(e))+t(Q)%*%(Q/e)
          
          ##GP
          mu <- X%*%beta[s,]+w[,s]
          Y.rep[,s] <- rnorm(n, mu, sqrt(e))
          
          if(verbose){
            if(status == n.report){
              cat(paste("Sampled: ",s," of ",n.samples,", ",round(100*s/n.samples,2),"%\n", sep=""))
              status <- 0
            }
            status <- status+1
          }
          
        }
        
        d.bar <- mean(d)
        
        ##Get d.bar.omega 
        Q <- Y-X%*%beta.mu-w.mu
        
        if(cov.model != "matern"){
          ct <- sigma.sq.mu*spCor(obs.knots.D, cov.model, phi.mu)
          C.str <- sigma.sq.mu*spCor(knots.D, cov.model, phi.mu)
        }else{
          ct <- sigma.sq.mu*spCor(obs.knots.D, cov.model, phi.mu, nu.mu)
          C.str <- sigma.sq.mu*spCor(knots.D, cov.model, phi.mu, nu.mu)
        }
        
        ##ct C^{*-1} c
        C <- ct%*%chol2inv(chol(C.str))%*%t(ct)
        
        e <- tau.sq.mu+sigma.sq.mu-diag(C)
        
        d.bar.omega <-  sum(log(e))+t(Q)%*%(Q/e)
        
        pd <- d.bar - d.bar.omega
        dic <- d.bar + pd
        
        rownames(DIC) <- c("bar.D", "D.bar.Omega", "pD", "DIC")
        colnames(DIC) <- c("value")
        DIC[1,1] <- d.bar
        DIC[2,1] <- d.bar.omega
        DIC[3,1] <- pd
        DIC[4,1] <- dic
        
        out <- list()
        out$DIC <- DIC
        out$GP <- GP(Y.rep,Y)
        out$GRS <- GR(Y.rep,Y)
        
        out
        
      }else{
        
        for(s in 1:n.samples){
          
          Q <- Y-X%*%beta[s,]-w[,s]
          
          d[s] <- n*log(tau.sq[s])+(t(Q)%*%Q)/tau.sq[s]
          
          status <- 0
          
          ##GP
          mu <- X%*%beta[s,]+w[,s]
          Y.rep[,s] <- rnorm(n, mu, sqrt(tau.sq[s]))
          
          if(verbose){
            if(status == n.report){
              cat(paste("Sampled: ",s," of ",n.samples,", ",round(100*s/n.samples,2),"%\n", sep=""))
              status <- 0
            }
            status <- status+1
          }
          
        }
        
        d.bar <- mean(d)
        
        ##Get d.bar.omega
        Q <- Y-X%*%beta.mu-w.mu
        
        d.bar.omega <- n*log(tau.sq.mu)+(t(Q)%*%Q)/tau.sq.mu
        
        pd <- d.bar - d.bar.omega
        dic <- d.bar + pd
        
        rownames(DIC) <- c("bar.D", "D.bar.Omega", "pD", "DIC")
        colnames(DIC) <- c("value")
        DIC[1,1] <- d.bar
        DIC[2,1] <- d.bar.omega
        DIC[3,1] <- pd
        DIC[4,1] <- dic
        
        out <- list()
        out$DIC <- DIC
        out$GP <- GP(Y.rep,Y)
        out$GRS <- GR(Y.rep,Y)
        
        out
      }
    }
  } else if(class(sp.obj) %in% c("spGLM","spMvGLM")){ 
    
    family <- sp.obj$family
    X <- sp.obj$X
    Y <- sp.obj$Y
    p <- ncol(X)
    m <- 1 ##for spGLM
    if(class(sp.obj) == "spMvGLM"){m <- sp.obj$m}
    obs.coords <- sp.obj$coords
    n <- nrow(obs.coords)
    cov.model <- sp.obj$cov.model
    p.beta.theta.samples <- sp.obj$p.beta.theta.samples
    n.samples <- nrow(p.beta.theta.samples)
    weights <- as.vector(t(sp.obj$weights))
        
    ##subsamples
    if(missing(end)){end <- n.samples}
    
    if(!is.numeric(start) || start >= n.samples)
      stop("error: invalid start")
    if(!is.numeric(end) || end > n.samples) 
      stop("error: invalid end")
    if(!is.numeric(thin) || thin >= n.samples) 
      stop("error: invalid thin")
    
    s.indx <- seq(start, end, by=as.integer(thin))

    beta <- sp.obj$p.beta.theta.samples[s.indx,1:p,drop=FALSE]
    w <- sp.obj$p.w.samples[,s.indx,drop=FALSE]

    beta.mu <- apply(beta, 2, mean)
    w.mu <- as.matrix(apply(w,1,mean))
    
    n.samples <- nrow(beta)

    status <- 0
    
    d <- rep(0, n.samples)
    DIC <- matrix(0,4,1)

    if(verbose){
      cat("-------------------------------------------------\n\t\tCalculating scores\n-------------------------------------------------\n")
    }
    
    for(s in 1:n.samples){
      
      if(family == "poisson"){
        d[s] <- -2*sum(Y*(X%*%beta[s,]+w[,s])-exp(X%*%beta[s,]+w[,s])+log(weights))
      }else if(family == "binomial"){
        pp <- 1/(1+exp(-X%*%beta[s,]-w[,s]))
        d[s] <- -2*sum(Y*log(pp)+(weights-Y)*log(1-pp))
      }else{
        stop("error: family is misspecified")
      }

    }
    
    d.bar <- mean(d)
    
    if(family == "poisson"){
      d.bar.omega <- -2*sum(Y*(X%*%beta.mu+w.mu)-exp(X%*%beta.mu+w.mu)+log(weights))
    }else if(family == "binomial"){
      pp <- 1/(1+exp(-X%*%beta.mu-w.mu))
      d.bar.omega <- -2*sum(Y*log(pp)+(weights-Y)*log(1-pp))
    }else{
      stop("error: family is misspecified")
    }
    
    pd <- d.bar - d.bar.omega
    dic <- d.bar + pd
    
    DIC <- matrix(0,4,1)
    
    rownames(DIC) <- c("bar.D", "D.bar.Omega", "pD", "DIC")
    colnames(DIC) <- c("value")
    DIC[1,1] <- d.bar
    DIC[2,1] <- d.bar.omega
    DIC[3,1] <- pd
    DIC[4,1] <- dic
    
    out <- list()    
    out$DIC <- DIC
    
    out
    
  }else{
    stop("error: wrong class of input object.\n")
  }
}
