## __________________________________________________________
##
## InferEdges
##
## INPUT
##	X	: n*p empirical data matrix
##	penalty : penalty mask to use (numeric matrix or scalar)
##
## OPTIONAL INPUT
##      ctrl.inference : 
##
## OUTPUT
##     Theta         : inferred network
##     loglik        : loglikelihood for the corresponding model
##     loglik.l1     : loglikelihhod + |penalty * Theta|_1
##     BIC           : loglikelihood - 1/2 df * log(n-1)
##     AIC           : loglikelihood - df
##     method        : the method used to perform the inference
##
## Estimate the inverse covariance matrix from the empirical
## covariance matrix solving p independent l1-penalized regressions
## (Meinshausen and Bulhman) or the a l1-penalized likelihood of
## Gaussian multivariate sample (Banerjee et al, Friedman et al) or
## the l1-penalized likelihood of a first order, vector autoregressive
## (VAR1) process (Charbonnier et al.).
## __________________________________________________________
##
InferEdges <- function(X, penalty, ctrl.inf = OptInference()) {
  
  ## Initializing
  n <- nrow(X)
  p <- ncol(X)
  
  if (is.matrix(penalty)) {
    Rho <- penalty
  } else if (length(penalty)==1) {
    Rho <- matrix(penalty, p, p)
  }
  
  method        <- ctrl.inf$method
  initial.guess <- ctrl.inf$initial.guess
  sym.rule      <- ctrl.inf$sym.rule
  tasks         <- ctrl.inf$tasks
  coupling      <- ctrl.inf$coupling
  solver        <- ctrl.inf$solver
  eps           <- ctrl.inf$eps
  maxIt         <- ctrl.inf$maxIt
  alpha         <- ctrl.inf$alpha
    
  ## PENALIZED IID GAUSSIAN PSEUDO-LIKELIHOOD
  if (method == "neighborhood.selection") {
    out <- neighborhood.selection(X,Rho,initial.guess,solver,eps)
  } else
  ## PENALIZED IID GAUSSIAN LIKELIHOOD WITH GLASSO
  if (method == "graphical.lasso") {
    out <- graphical.lasso(X,Rho,initial.guess,solver,eps,maxIt)
  } else 
  ## PENALIZED VAR1-LIKELIHOOD
  if (method == "var1.inference") {
    out <- var1.inference(X,Rho,initial.guess,solver,eps)
  }
  ## MULTIPLE TASKS INFERENCE
  if (method == "multi.gaussian") {
    out <- multiple.static(X,Rho,tasks=tasks,coupling=coupling,alpha=alpha,
                           initial.guess=initial.guess,eps=eps)
  }
  if (method == "multi.var1") {
    out <- multiple.dynamic(X,Rho,tasks=tasks,coupling=coupling,alpha=alpha,
                            initial.guess=initial.guess,eps=eps)
  }
  Theta <- switch(sym.rule, "AND" = out$Theta.and, "OR" = out$Theta.or, out$Theta)

  ## Get the number of edges inferred
  if (sum(is.na(Theta))>0) {
    n.edges <- NA
  } else {
    n.edges <- switch(method,
                      "neighborhood.selection" = nb.edges(Theta,FALSE),
                      "graphical.lasso"        = nb.edges(Theta,FALSE),
                      "var1.inference"         = nb.edges(Theta,TRUE),
                      "multi.gaussian"         = sapply(Theta,nb.edges,FALSE),
                      "multi.var1"             = sapply(Theta,nb.edges,TRUE))
  }
  
  return(list(Theta=Theta, Beta=out$Beta, n.edges=n.edges, loglik=out$loglik,
              loglik.l1=out$loglik.l1, BIC=out$BIC, AIC=out$AIC))
}

OptInference <- function(method = "neighborhood.selection", initial.guess =NULL,
                         sym.rule = NULL, tasks =NULL, coupling="intertwined",
                         solver = NULL, eps = 1e-8, maxIt = 50, alpha = 0.5) {
  
  ## Enforce multiple inference if factors are specified
  if (nlevels(tasks) > 1 & !(method %in% c("multi.gaussian","multi.var1"))) {
    method   <- "multi.gaussian"
    coupling <- "intertwined"
  }
  
  ## Enforce no symmetrization rule for the var1 inference model
  if (method %in% c("var1.inference","multi.var1")) {
    sym.rule <- "NO"
  } else if (is.null(sym.rule)) {
    sym.rule <- "AND"
  }
  
  ## Setting default Lasso solver according to the inference method
  if (is.null(solver)) {
    if (method == "graphical.lasso") {
      solver <- "shooting"
    } else {
      solver <- "active.set"
    }
  }
  
  return(list(method = method, initial.guess = initial.guess,
              sym.rule = sym.rule, tasks = tasks, coupling = coupling,
              solver = solver, eps = eps, maxIt = maxIt, alpha = alpha))
  
}

nb.edges <- function(x,directed) {

  edges <- sum(abs(x)>0)    
  if (!directed) {
    edges <- edges - sum(abs(diag(x))>0)
    edges <- edges/2
  }
 
  return(edges)
}

neighborhood.selection <- function(X, Rho, initial.guess = NULL,
                                   solver = "active.set", eps = 1e-8) {

  ## INITIALIZING
  n <- nrow(X)
  p <- ncol(X)
  S <- var(X,na.rm=TRUE)

  ## WARM RESTART IF POSSIBLE
  if (is.null(initial.guess)) {
    Beta  <- matrix(0,p,p)
    Theta <- matrix(0,p,p)
    diag(Theta) <- diag(1/S)
  } else {
    Beta <- initial.guess
    Theta <- matrix(0,p,p)
    ## Reconstruction of the full matrix Theta  
    for(k in 1:p) {
      Theta[k,k]  <- 1/(S[k,k])
      Theta[-k,k] <- - Beta[-k,k]*Theta[k,k]
    }
  }
  dimnames(Theta) <- list(colnames(X), colnames(X))

  ## Handling case where Rho = 0 or when initial guess already has
  ## more non zero entries than possibly activated
  if ((max(Rho[upper.tri(Rho)]) == 0) |
      (sum(abs(Theta) > 0)-sum(diag(abs(Theta))>0)) >= n*(p-1) ) {
    
    if (max(Rho[upper.tri(Rho)]) == 0) {
      Sm1 <- try(solve(S),silent=TRUE)
      if (is.matrix(Sm1)) {
        Theta <- Sm1
      }
    }
    
    Theta.and <- sign(Theta) * pmin(abs(Theta),t(abs(Theta)))
    Theta.or <- pmax(Theta,t(Theta)) - pmax(-Theta,-t(Theta)) 
    diag(Theta.or) <- diag(Theta.or)/2
    D  <- diag(Theta)
    Theta.tilde <- Theta %*% diag(D^(-1/2))
    loglik <- (n/2)*(log(prod(D))-sum(t(Theta.tilde)%*%var(X) %*% Theta.tilde))
    loglik.l1 <- loglik - (n/2) * sum(abs(Rho * Theta))
    df        <- sum(abs(Theta)>0) - sum(abs(diag(Theta))>0) 
    BIC       <- loglik - .5* df * log(n)
    AIC       <- loglik - df            
    
    return(list(Theta=Theta,Theta.or=Theta.or,Theta.and=Theta.and, Beta=Beta,
                loglik=loglik, loglik.l1=loglik.l1, BIC=BIC, AIC=AIC))
    
  }
  
  if ("shooting" %in% solver) {
    ## SHOOTING ALGORITHM
    Rho[is.infinite(Rho)] <- -1 
    out <- .C("Lasso",
              as.integer (p),
              as.double  (S),
              as.double  (Rho),
              as.double  (eps),
              Beta =  as.double(Beta),
              PACKAGE="simone")
    Beta <- matrix(out$Beta, ncol=p)    
  } else {
    ## ACTIVE SET ALGORITHM
    for(k in 1:p) {
      ## Solving the Lasso problem on the kth column
      out <- LassoConstraint(S[-k,-k],S[-k,k],Rho[-k,k],
                             beta=Beta[-k,k],eps=eps,maxSize=min(n,p-1))
      if (out$converged) {
        Beta[-k,k] <- out$beta
      } else {
        cat("out of convergence... stopping here")
        return(list(Theta=NA,Theta.or=NA,Theta.and=NA,
                    loglik=NA, loglik.l1=NA, BIC=NA, AIC=NA))
      }      
    }
  }

  ## Reconstruction of the full matrix Theta  
  for(k in 1:p) {
    Theta[k,k]  <- 1/S[k,k]
    Theta[-k,k] <- -Beta[-k,k]*Theta[k,k]
  }
  
  ## POST-SYMETRIZATION WITH AND RULE
  Theta.and <- sign(Theta) * pmin(abs(Theta),t(abs(Theta)))
  
  ## POST-SYMETRIZATION WITH OR RULE
  Theta.or <- pmax(Theta,t(Theta)) - pmax(-Theta,-t(Theta)) 
  diag(Theta.or) <- diag(Theta.or)/2
  
  ## Compute various criterion (Pseudo-likelihood of a Gaussian vector)
  D  <- diag(Theta)
  Theta.tilde <- Theta %*% diag(D^(-1/2))
  loglik <- (n/2) * (log(prod(D)) - sum(t(Theta.tilde) %*% S %*% Theta.tilde))
  loglik.l1 <- loglik - (n/2) * sum(abs(Rho * Theta))
  df        <- sum(abs(Theta)>0) - sum(abs(diag(Theta))>0)
  BIC       <- loglik - .5* df * log(n)
  AIC       <- loglik - df    

  return(list(Theta=Theta,Theta.or=Theta.or,Theta.and=Theta.and,Beta=Beta,
              loglik=loglik, loglik.l1=loglik.l1, BIC=BIC, AIC=AIC))
}

graphical.lasso <- function(X, Rho, initial.guess = NULL,
                            solver = "shooting", eps = 1e-6, maxIt = 25) {

  ## INITIALIZING
  n <- nrow(X)
  p <- ncol(X)
  S <- var(X,na.rm=TRUE)

  Sigma <- matrix(S + diag(diag(Rho)),nrow=p, ncol=p)
  ## Test the positive definiteness of the starting matrix
  if (sum(abs(eigen(Sigma)$values) < eps) > 0) {
    diag(Rho) <- diag( rowSums(abs(S)) + eps)
    Sigma     <- matrix(S + diag(diag(Rho)),nrow=p, ncol=p)
    cat("diag(Rho) too small to start with positive definite matrix - Enforcing brutaly Sigma0 to be diagonal dominant\n")  
  }   

  ## WARM RESTART IF POSSIBLE
  if (is.null(initial.guess)) {
    Beta  <- matrix(0,p,p)
    Theta <- matrix(0,p,p) 
    diag(Theta) <- diag(1/S)
    initial.guess <- Theta
  } else {
    Beta <- initial.guess
    Theta <- matrix(0,p,p) 
    ## Computing Theta.hat by blockwise invertion of Sigma ...
    for(k in 1:p) { # let's work on the kth column of Sigma
      ## Reconstruction of the full matrix Theta
      Theta[k,k]  <- 1/(Sigma[k,k] - Sigma[k,-k] %*% Beta[-k,k]/2)
      Theta[-k,k] <- - Theta[k,k] * Beta[-k,k]/2
    }
  }
  dimnames(Theta) <- list(colnames(X), colnames(X))
  
  if ("shooting" %in% solver) {
    ## SHOOTING ALGORITHM
    Rho[is.infinite(Rho)] <- -1 
    out <- .C("GLasso",
              as.integer (p),
              as.double  (S),
              as.double  (Rho),
              as.double  (eps),
              as.integer (maxIt),
              Sigma   = as.double(Sigma),
              Beta    = as.double(Beta),
              finalIt = as.integer(1),
              PACKAGE ="simone")
    
    Sigma  <- matrix(out$Sigma, ncol=p)
    Beta   <- matrix(out$Beta, nrow=p, ncol=p)
    n.iter <- out$finalIt

    ## Computing Theta.hat by blockwise invertion of Sigma ...
    for(k in 1:p) { # let's work on the kth column of Sigma
      ## Reconstruction of the full matrix Theta
      Theta[k,k]  <- 1/(Sigma[k,k] - Sigma[k,-k] %*% Beta[-k,k]/2)
      Theta[-k,k] <- - Theta[k,k] * Beta[-k,k]/2
    }
    ## The dual gap (Should  be zero)
    PD.gap <- abs(sum(diag(Theta %*% S)) + sum(abs(Rho * Theta)) - p)
    
  } else {
    ## ACTIVE SET ALGORITHM
    n.iter <- 0

    ## Initial value of Theta (Sigma is positive definite)
    threshold <- eps * sum(abs(Sigma)) / p
    ##PD.gap <- abs(sum(diag(Theta %*% S)) + sum(abs(Rho * Theta)) - p)
    PD.gap <- Inf    
    while((PD.gap > threshold) & (n.iter < maxIt)) {

      n.iter <- n.iter + 1

      for(k in 1:p) { # let's work on the kth column of Sigma
        ## Solving the Lasso problem on the kth column
        out <- LassoConstraint(Sigma[-k,-k]/2,-S[-k,k],Rho[-k,k],
                               beta=Beta[-k,k],eps=eps,maxSize=min(n,p-1))
        if (out$converged) {
          Beta[-k,k] <- out$beta
        } else {
          cat("out of convergence... stopping here")
          return(list(Sigma=NA,Theta=NA,Theta.and=NA,Theta.or=NA,
                      Beta=NA,PD.gap=NA,n.iter=NA,
                      loglik=NA, loglik.l1=NA, BIC=NA, AIC=NA))
        }
        
        Sigma[-k,k] <- Sigma[-k,-k] %*% Beta[-k,k]/2
        Sigma[k,-k] <- Sigma[-k,k]
        
        ## Reconstruction of the full matrix Theta
        Theta[k,k]  <- 1/(Sigma[k,k] - Sigma[k,-k] %*% Beta[-k,k]/2)
        Theta[-k,k] <- -Theta[k,k] * Beta[-k,k]/2
      }
      PD.gap <- abs(sum(diag(Theta %*% S)) + sum(abs(Rho * Theta)) - p)      
    }
  }
     
  ## POST-SYMETRIZATION WITH AND RULE
  Theta.and <- sign(Theta) * pmin(abs(Theta),t(abs(Theta)))
  
  ## POST-SYMETRIZATION WITH OR RULE
  Theta.or <- pmax(Theta,t(Theta)) - pmax(-Theta,-t(Theta)) 
  diag(Theta.or) <- diag(Theta.or)/2
  
  ## Compute selection criterion
  loglik <- (n/2)*(log(det(Theta)) - sum(diag(Theta %*% S)))
  loglik.l1   <- loglik - (n/2) * sum(abs(Rho * Theta))
  df     <- sum(abs(Theta)>0) - sum(abs(diag(Theta))>0)
  BIC    <- loglik - .5* df * log(n)
  AIC    <- loglik - df    

  return(list(Sigma=Sigma,Theta=Theta,Theta.and=Theta.and,Theta.or=Theta.or,
              Beta=Beta,PD.gap=PD.gap,n.iter=n.iter,
              loglik=loglik, loglik.l1=loglik.l1, BIC=BIC, AIC=AIC))
}

var1.inference <- function(X, Rho, initial.guess = NULL,
                           solver="active.set", eps=1e-8) {
  
  ## Initializing
  n <- nrow(X)
  p <- ncol(X)

  if (is.null(initial.guess)) {
    Theta <- matrix(0,p,p)
    Beta  <- matrix(0,p,p)
  } else {
    Theta <- initial.guess
    Beta  <- initial.guess
  }
  dimnames(Theta) <- list(colnames(X), colnames(X))

  ## Empirical covariances
  S <- 1/(n-1) * t(X[-n,]) %*% X[-n,]
  V <- 1/(n-1) * t(X[-n,]) %*% X[-1,]
 
  ## Handling case where Rho = 0
  if (max(Rho) == 0) {
    Sm1 <- try(solve(S),silent=TRUE)
    if (is.matrix(Sm1)) {
      Theta <- Sm1 %*% V
    }
    loglik <- (n/2)*(2*sum(diag(t(V)%*% Theta))-sum(diag(t(Theta) %*% S %*% Theta)))
    loglik.l1   <- loglik - (n/2) * sum(abs(Rho * Theta))
    df     <- sum(abs(Theta)>0)
    BIC    <- loglik - .5* df * log(n-1)
    AIC    <- loglik - df        
    return(list(Theta=Theta, loglik=loglik, loglik.l1=loglik.l1, BIC=BIC, AIC=AIC))
  }
    
  ## Solving the Lasso problem
  ## "active.set" et "shooting" methods are available
  if ("shooting" %in% solver) {
    Rho[is.infinite(Rho)] <- -1 
    out <- .C("ARLasso",
              as.integer (p),
              as.double  (S), 
              as.double  (V),
              as.double  (Rho),
              as.double  (eps),
              Theta = as.double(Theta),
              PACKAGE = "simone")
    Theta <- matrix(out$Theta, ncol=p)
  } else {
    for(k in 1:p) { # let's work on the kth column of Wxs
      ## Solving the Lasso problem on the kth column
      out <- LassoConstraint(S,-V[,k],Rho[,k],beta=Theta[,k],eps=eps,
                             maxSize=min(n,p))
      if (out$converged) {
        Theta[,k] <- out$beta
      } else {
        cat("out of convergence... stopping here")
        return(list(Theta=NA, Beta=NA, loglik=NA, loglik.l1=NA, BIC=NA, AIC=NA))
      }
    }
  }

  ## Compute selection criterion
  loglik <- (n/2)*(2*sum(diag(t(V)%*% Theta)) - sum(diag(t(Theta) %*% S %*% Theta)))
  loglik.l1   <- loglik - (n/2) * sum(abs(Rho * Theta))
  df     <- sum(abs(Theta)>0)
  BIC    <- loglik - .5* df * log(n-1)
  AIC    <- loglik - df    
  
  return(list(Theta=Theta, Beta=Beta, loglik=loglik, loglik.l1=loglik.l1,
              BIC=BIC, AIC=AIC))
}

multiple.static <- function(X, Rho, tasks=factor(rep(1,nrow(X))),
                            coupling="intertwined", alpha=0.5,
                            initial.guess = NULL, eps=1e-8) {
 
  ## Basic parameters
  T  <- nlevels(tasks)
  p  <- ncol(X)
  nt <- table(tasks)
  n  <- sum(nt)
  
  ## __________________________________________________
  ## 0. BUILDING-UP THE MATRICES

  ## S.t stocked the variance of each task as a list
  S.t <- by(X,tasks,var,na.rm=TRUE) 

  if (coupling %in% c("intertwined")) {
    S <- matrix(0,p,p)
    for (t in 1:T) {
      S <- S + S.t[[t]] *  nt[t] / n
    }
    ## add S to each S.t to get Stilde.t
    S.t <- lapply(S.t,function(x) alpha * x + (1-alpha) * S )
  }
  
  
  ## __________________________________________________
  ## 1. SOLVING THE UNDERLYING MATRICIAL LASSO PROBLEMS
  if (is.null(initial.guess)) {
    Beta <- matrix(0,(p-1)*T, p)
  } else {
    Beta <- initial.guess
  }
  
  coupling.method <- switch(coupling,
                            "intertwined" = LassoConstraint,
                            "groupLasso" = GroupLassoConstraint,
                            "coopLasso"  = CoopLassoConstraint)

  pen.norm <- switch(coupling, "intertwined" = 1, sqrt(p))
  
  for (k in 1:p) {

    ## The penalty parameters
    Lambda <- cbind(rep(Rho[-k,k],T)) * pen.norm
    ## Build the current C.11, C.12 (according to the T kth columns)
    C.11 <- matrix(0,(p-1)*T,(p-1)*T)
    C.12 <- NULL
    for (t in 1:T) {
      current.block <- ((t-1)*(p-1)+1):(t*(p-1))
      C.11[current.block,current.block] <- S.t[[t]][-k,-k]
      C.12 <- cbind(c(C.12,S.t[[t]][-k,k]))      
    }
    
    out <- coupling.method(C.11, C.12, Lambda, T, beta=Beta[,k], eps=eps,
                           maxSize=min(n,p*T))
    if (out$converged) {
      Beta[,k] <- out$beta
    } else {
      cat("out of convergence... stopping here")
      return(list(Theta=NA,Theta.and=NA,Theta.or=NA, Beta=NA,
                  loglik=NA, loglik.l1=NA, BIC=NA, AIC=NA))
    }
    
  }

  ## __________________________________________________
  ## 2. FULL RECONSTRUCTION OF THE K(t)s MATRICES
  Theta     <- list()
  Theta.and <- list()
  Theta.or  <- list()
  loglik    <- 0
  loglik.l1 <- 0
  BIC       <- 0
  AIC       <- 0
  for (t in 1:T) {
    Theta.c <- matrix(0,p,p)
    dimnames(Theta.c) <- list(colnames(X), colnames(X))

    indices <- ((t-1)*(p-1)+1):(t*(p-1))
    for (k in 1:p) {
      Theta.c[k,k]  <- 1/(S.t[[t]][k,k])
      Theta.c[-k,k] <- Beta[indices,k] * Theta.c[k,k] 
    }

    ## Compute various criterion (Pseudo-likelihood of a Gaussian vector)
    D  <- diag(Theta.c)
    Theta.tilde <- Theta.c %*% diag(D^(-1/2))
    loglik.c    <- (nt[t]/2)*(log(prod(D)) -
                              sum(t(Theta.tilde) %*% var(X)%*%Theta.tilde))
    loglik.l1.c <- loglik.c - (nt[t]/2) * sum(abs(Rho * Theta.c))
    df.c        <- sum(abs(Theta.c)>0)- sum(abs(D)>0)
    BIC.c       <- loglik.c - .5* df.c * log(nt[t])
    AIC.c       <- loglik.c - df.c
    
    ## Compute selection criterion for the whole problem
    loglik    <- loglik + loglik.c
    loglik.l1 <- loglik.l1 +  loglik.l1.c
    BIC       <- BIC + BIC.c
    AIC       <- AIC + AIC.c

    ## Matrix Theta and various criteria
    Theta[[t]] <- Theta.c
    
    ## POST-SYMETRIZATION WITH AND RULE
    Theta.and[[t]] <- sign(Theta.c) * pmin(abs(Theta.c),t(abs(Theta.c)))
    
    ## POST-SYMETRIZATION WITH OR RULE
    Theta.or[[t]] <- pmax(Theta.c,t(Theta.c)) - pmax(-Theta.c,-t(Theta.c))
    diag(Theta.or[[t]]) <- diag(Theta.or[[t]])/2
    
  }

  return(list(Theta=Theta,Theta.and=Theta.and,Theta.or=Theta.or, Beta=Beta,
              loglik=loglik, loglik.l1=loglik.l1, BIC=BIC, AIC=AIC))

}



multiple.dynamic <- function(X, Rho, tasks=factor(rep(1,nrow(X))),
                             coupling="intertwined",  alpha=0.5,
                             initial.guess = NULL,eps=1e-8) {
    
  ## Basic parameters
  T  <- nlevels(tasks)
  p  <- ncol(X)
  nt <- table(tasks)
  n  <- sum(nt)
  
  ## __________________________________________________
  ## 0. BUILDING-UP THE MATRICES

  ## S.t stocked the variance of each task as a list
  S.t <- list()
  V.t <- list()
  for (t in 1:T) {
    x <- as.matrix(X[tasks == levels(tasks)[t],])
    S.t[[t]] <- 1/(nt[t]-1) * t(x[-nt[t],]) %*% x[-nt[t],]
    V.t[[t]] <- 1/(nt[t]-1) * t(x[-nt[t],]) %*% x[-1,]
  }
  
  if (coupling %in% c("intertwined")) {
    S <- matrix(0,p,p)
    V <- matrix(0,p,p)
    for (t in 1:T) {
      S <- S + S.t[[t]] *  nt[t] / n
      V <- V + V.t[[t]] *  nt[t] / n
    }
    ## add S to each S.t to get Stilde.t
    S.t <- lapply(S.t,function(x) alpha * x + (1-alpha) * S )
    V.t <- lapply(V.t,function(x) alpha * x + (1-alpha) * V )
  }
  
  
  ## __________________________________________________
  ## 1. SOLVING THE UNDERLYING MATRICIAL LASSO PROBLEMS
  if (is.null(initial.guess)) {
    Beta <- matrix(0,p*T, p)
  } else {
    Beta <- initial.guess
  }
  
  coupling.method <- switch(coupling,
                            "intertwined" = LassoConstraint,
                            "groupLasso" = GroupLassoConstraint,
                            "coopLasso"  = CoopLassoConstraint)

  pen.norm <- switch(coupling, "intertwined" = 1, sqrt(T))
  
  for (k in 1:p) {

    ## The penalty parameters
    Lambda <- cbind(rep(Rho[,k],T)) * pen.norm
    ## Build the current C.11, C.12 (according to the T kth columns)
    C <- matrix(0,p*T,p*T)
    D <- NULL
    for (t in 1:T) {
      current.block <- ((t-1)*p+1):(t*p)
      C[current.block,current.block] <- S.t[[t]]
      D <- cbind(c(D,V.t[[t]][,k]))
    }
    
    out <- coupling.method(C, -D, Lambda, T, beta=Beta[,k], eps=eps,
                           maxSize=min(n,p*T))
    if (out$converged) {
      Beta[,k] <- out$beta
    } else {
      cat("out of convergence... stopping here")
      return(list(Theta=NA, Beta=NA, loglik=NA, loglik.l1=NA, BIC=NA, AIC=NA))
    }    
  }

  ## __________________________________________________
  ## 2. FULL RECONSTRUCTION OF THE K(t)s MATRICES
  Theta     <- list()
  loglik    <- 0
  loglik.l1 <- 0
  BIC       <- 0
  AIC       <- 0
  for (t in 1:T) {
    Theta.c <- matrix(0,p,p)
    dimnames(Theta.c) <- list(colnames(X), colnames(X))

    indices <- ((t-1)*p+1):(t*p)
    for (k in 1:p) {
      Theta.c[,k] <- Beta[indices,k]
    }
    
    ## Compute selection criterion for the current task
    loglik.c    <- (nt[t]/2)*(2*sum(diag(t(V.t[[t]])%*% Theta.c)) -
                              sum(diag(t(Theta.c) %*% S.t[[t]] %*% Theta.c)))
    loglik.l1.c <- loglik.c - (nt[t]/2) * sum(abs(Rho * Theta.c))
    df.c        <- sum(abs(Theta.c)>0)
    BIC.c       <- loglik.c - .5* df.c * log(nt[t]-1)
    AIC.c       <- loglik.c - df.c
    
    ## Compute selection criterion for the whole problem
    loglik    <- loglik + loglik.c
    loglik.l1 <-  loglik.l1 +  loglik.l1.c
    BIC       <- BIC + BIC.c
    AIC       <- AIC + AIC.c
    
    ## Matrix Theta and various criteria
    Theta[[t]] <- Theta.c    
  }

  return(list(Theta=Theta,Beta=Beta,loglik=loglik,loglik.l1=loglik.l1,
              BIC=BIC,AIC=AIC))

}

