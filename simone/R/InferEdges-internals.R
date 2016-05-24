LassoConstraint <- function(S,s,rho,T=NULL,beta=rep(0,ncol(S)),eps=1e-6,
                             maxSize=ncol(S),maxIt=10*maxSize) {
  
  ## Initialisation
  p <- ncol(S)
  sigma <- abs(beta) > 0
  theta <- sign(beta)
  it <- 0
  
  abs.dL.min <- abs(dL1(beta,S,s,rho))
  ## Check if optimality is reached at starting point
  if (all(abs.dL.min < eps) | sum(sigma) >= maxSize){
    return(list(beta=beta,converged=TRUE))
  } else {
    ## Look for the most violating entry
    l <- which.max(abs.dL.min)
    nabla.f <- S %*% beta + s
    sigma[l] <- TRUE
    theta[l] <- -sign(nabla.f[l])
  }
  
  while (1) {
    it <- it+1

    sign.feasible <- FALSE
    
    while (!sign.feasible) {
    
      ## (1) OPTIMIZATION OVER THE ACTIVE SET
      ## ________________________________________________________________
      h   <- rep(0,p)
      x <- try(solve(S[sigma,sigma],cbind(rho*theta+s)[sigma]),silent=TRUE)
      if (is.vector(x)) {
        h[sigma] <- - beta[sigma] - x 
      } else {
        out <- optim(beta[sigma],method="BFGS",fn=L1,gr=dL1,
                     S=S[sigma,sigma],s=s[sigma],Lambda=rho[sigma])
        h[sigma] <- out$par-beta[sigma]
        ##if (out$convergence != 0) {cat('\nconvergence pb in optim:',out$convergence)}
      }
      
      ## (2) TEST LAGRAGIAN VIOLATION AND GROUP DELETION IF APPLICABLE
      ## ________________________________________________________________
      ## Check for sign feasibility
      sign.feasible <- all(sign(beta+h)[sigma]==theta[sigma])
      if (!sign.feasible) {
        ind <- which(sign(beta+h)[sigma]!=theta[sigma])
        gamma <- -beta[sigma][ind]/h[sigma][ind]
        gamma[is.nan(gamma)] <- Inf
        k <- which(gamma == min(gamma))
        if (gamma[k[1]] == Inf) {
          return(list(beta=beta,converged=FALSE))
        }
        
        ## Remove the kth active constraint, which make beta sign feasible
        beta <- beta + gamma[k[1]]*h
        sigma[sigma][ind][k] <- FALSE
      } else {
        beta <- beta+h
      }
      theta[sigma] <- sign(beta)[sigma]
    }

    ## (3) OPTIMALITY TESTING AND GROUP ACTIVATION IF APPLICABLE
    ## ________________________________________________________________
    abs.dL.min <- abs(dL1(beta,S,s,rho))
    if (all(abs.dL.min[!sigma] < eps) | it > maxIt | sum(sigma) >= maxSize){
      if (it > maxIt & any(abs.dL.min[!sigma] > eps)) {
        return(list(beta=beta,converged=FALSE))
      } else {
        return(list(beta=beta,converged=TRUE))
      }
    } else {
      ## Look for the most violating entry
      l <- which.max(abs.dL.min)
      nabla.f <- S %*% beta + s
      sigma[l] <- TRUE
      theta[l] <- -sign(nabla.f[l])
    }
  }
}

L1 <- function(beta,S,s,Lambda) {
  L <- 0.5*t(beta) %*% S %*% beta + t(beta) %*% s + sum(Lambda*abs(beta))
  return(L)
}

dL1 <- function(beta,S,s,Lambda) {
  dL <- S %*% beta + s
  zero <- beta==0
  ## terms whose beta is different from zero
  dL[!zero] <- dL[!zero] + Lambda[!zero] * sign(beta[!zero])
  ## terms whose beta is zero
  can.be.null <- abs(dL[zero]) < Lambda[zero]
  ## which ones can be null ?
  dL[zero][can.be.null] <- 0
  dL[zero][!can.be.null] <- dL[zero][!can.be.null] -
    Lambda[zero][!can.be.null]*sign(dL[zero][!can.be.null])  
  return(dL)
}

GroupLassoConstraint <- function(C11,C12,Lambda,T,beta=rep(0,nrow(C11)),
                                 eps=1e-6, maxSize=ncol(C11),
                                 maxIt=10*maxSize) {
  
  ## INITIALIZATION
  p <- nrow(C11)/T
  active.set <- rep(rowSums(matrix(abs(beta)>0,p,T)) > 0,T)
  it <- 0

  ## Check if optimality is reached at starting point
  dL.min <- dL2(beta,C11,C12,Lambda,T,eps)
  abs.dL.min <- abs(dL.min)
  n2.dL <- sqrt(rowSums(matrix(dL.min,ncol=T)^2))
  if (all(abs.dL.min < eps) | sum(active.set) >= maxSize) {
    return(list(beta=beta,converged=TRUE))
  } else {
    ## UPDATE THE ACTIVE SET / ACTIVE GROUP
    i <- which.max(n2.dL)
    new <- seq(from=i,to=i+(T-1)*p,by=p) 
    active.set[new] <- TRUE
  }
  
  while (1) {
    
    it <- it+1
    
    ## (1) OPTIMIZATION OVER A
    ## _____________________________________________________________
    n2.beta <- rep(sqrt(rowSums(matrix(beta,ncol=T)^2)),T)
    res <- optim(beta[active.set], method="L-BFGS-B", fn=L2, gr=dL2,
                 C=C11[active.set,active.set],
                 c=C12[active.set],
                 rho=Lambda[active.set],T=T, eps,
                 control=list(maxit=200))
    beta[active.set] <- res$par
##    if (res$convergence != 0) {
      #cat('\nconvergence pb in optim:',res$convergence,"",res$message)
##      return(beta)
##    }

    ## (2) GROUP DELETION IF APPLICABLE
    ## _____________________________________________________________
    ## Check if some betas vanished during optimization
    dL.min <- dL2(beta,C11,C12,Lambda,T,eps)
    n2.beta <- sqrt(rowSums(matrix(beta,ncol=T)^2))
    n2.dL <- sqrt(rowSums(matrix(dL.min,ncol=T)^2))
    for (i in intersect(which(n2.beta < eps),which(active.set[1:p]))  ){
      if (n2.dL[i] < eps) {
        ind <- seq(from=i,to=i+(T-1)*p,by=p)
        active.set[ind] <- FALSE
      }
    }
    
    ## (3) OPTIMALITY TESTING AND GROUP ACTIVATION IF APPLICABLE
    ## _____________________________________________________________
    ## Get the entry for which the constraint is the most violated    
    abs.dL.min <- abs(dL.min)
    if (all(abs.dL.min[!active.set] < eps) | it > maxIt | sum(active.set) >= maxSize){
      if (it > maxIt & any(abs.dL.min[!active.set] > eps)) {
        return(list(beta=beta,converged=FALSE))
      } else {
        return(list(beta=beta,converged=TRUE))
      }
    } else {
      ## UPDATE THE ACTIVE SET / ACTIVE GROUP
      i <- which.max(n2.dL)
      new <- seq(from=i,to=i+(T-1)*p,by=p)
      active.set[new] <- TRUE
    }
  }
}

L2 <- function(beta,C,c,rho,T,eps) {
  n2 <- sqrt(rowSums(matrix(beta,ncol=T)^2))
  L <- .5*t(beta) %*% C %*% beta + t(c) %*% beta + sum(rho[1:length(n2)]*n2)
  return(L)
}

dL2 <- function(beta,C,c,rho,T,eps) {
  n2 <- rep(sqrt(rowSums(matrix(beta,ncol=T)^2)),T)
  zero <- n2 < eps

  ## terms whose norm2 is different from zero
  dL <- rep(0,length(beta))
  nabla.f <- C %*% beta + c
  dL[!zero] <- nabla.f[!zero] + rho[!zero] * beta[!zero] / n2[!zero]

  ## terms whose norm2 is zero
  n2.nabla.f <- rep(sqrt(rowSums(matrix(nabla.f,ncol=T)^2)),T)
  to.zero <- n2.nabla.f[zero] < rho[zero] + eps

  ## which ones can be null ?
  dL[zero][to.zero] <- 0
  dL[zero][!to.zero] <- nabla.f[zero][!to.zero] * (1 - rho[zero][!to.zero]/n2.nabla.f[zero][!to.zero])
  
  return(dL)
}  

CoopLassoConstraint <- function(C11,C12,Lambda,T,beta=rep(0,nrow(C11)),
                                 eps=1e-6, maxSize=ncol(C11),
                                 maxIt=10*maxSize){
  
  ## INITIALIZATION
  p <- nrow(C11)/T
  active.set <- rep(rowSums(matrix(abs(beta)>0,p,T)) > 0,T)
  h <- rep(0,T*p)
  it <- 0

  ## Check if optimality is reached at the starting point
  dL.min <- dL3(beta,C11,C12,Lambda,T,eps)
  if (all(abs(dL.min) < eps) | sum(active.set) >= maxSize) {
    return(list(beta=beta,converged=TRUE))
  } else {
    ## UPDATE THE ACTIVE SET / ACTIVE GROUP
    n2.dL.pos <- sqrt(rowSums(matrix(pmax(0, dL.min),ncol=T)^2))
    n2.dL.neg <- sqrt(rowSums(matrix(pmax(0,-dL.min),ncol=T)^2))
    i <- which.max(pmax(n2.dL.pos,n2.dL.neg))
    new <- seq(from=i,to=i+(T-1)*p,by=p)
    active.set[new] <- TRUE
  }
  
  while (1) {
    
    it <- it+1
    
    ## (1) OPTIMIZATION OVER A
    ## _____________________________________________________________
    a.pos <- rep(sqrt(rowSums(matrix(pmax(0, beta),ncol=T)^2)),T)
    a.neg <- rep(sqrt(rowSums(matrix(pmax(0,-beta),ncol=T)^2)),T)
    upper.bound <-rep(Inf,T*p)
    lower.bound <-rep(-Inf,T*p)
    lower.bound[a.neg==0] <- 0
    upper.bound[a.pos==0] <- 0
    if (-sign(n2.dL.pos[new][1]-n2.dL.neg[new][1]) == -1) {
      lower.bound[new] <- -Inf
    } else {
      upper.bound[new] <- Inf
    }
    res <- optim(beta[active.set], fn=L3, gr=dL3, method="L-BFGS-B",
                 C=C11[active.set,active.set], c=C12[active.set],
                 rho=Lambda[active.set], T=T, eps,
                 lower=lower.bound[active.set],upper=upper.bound[active.set],
                 control=list(maxit=200))
    beta[active.set] <- res$par
    beta <- beta+h
    
    ## (2) TEST LAGRAGIAN VIOLATION AND GROUP DELETION IF APPLICABLE
    ## _____________________________________________________________
    
    ## DESACTIVATION ON active.set = not(A+ cap A-) 
    dL.min <- dL3(beta,C11,C12,Lambda,T,eps)
    n2.beta <- sqrt(rowSums(matrix(beta,ncol=T)^2))
    n2.dL   <- sqrt(rowSums(matrix(dL.min,ncol=T)^2))
    for (i in intersect(which(n2.beta < eps),which(active.set[1:p]))  ){
      if (n2.dL[i] < eps) {
        ind <- seq(from=i,to=i+(T-1)*p,by=p)
        active.set[ind] <- FALSE
      }
    }
    
    ## (3) OPTIMALITY TESTING AND GROUP ACTIVATION IF APPLICABLE
    ## _____________________________________________________________
    ## Get the entry for which the constraint is the most violated
    if (all(abs(dL.min[!active.set]) < eps) | it > maxIt | sum(active.set) >= maxSize){
      if (it > maxIt & any(abs(dL.min[!active.set]) > eps)) {
        return(list(beta=beta,converged=FALSE))
      } else {
        return(list(beta=beta,converged=TRUE))
      }
    } else {
      ## UPDATE THE ACTIVE SET / ACTIVE GROUP
      n2.dL.pos <- sqrt(rowSums(matrix(pmax(0, dL.min),ncol=T)^2))
      n2.dL.neg <- sqrt(rowSums(matrix(pmax(0,-dL.min),ncol=T)^2))
      i <- which.max(pmax(n2.dL.pos,n2.dL.neg))
      new <- seq(from=i,to=i+(T-1)*p,by=p)
      active.set[new] <- TRUE
    }
  }
}

L3 <- function(beta,C,c,rho,T,eps) {
  a.pos <- sqrt(rowSums(matrix(pmax(0, beta),ncol=T)^2))
  a.neg <- sqrt(rowSums(matrix(pmax(0,-beta),ncol=T)^2))
  
  L <- 0.5*t(beta) %*% C %*% beta + t(c) %*% beta + sum(rho[1:length(a.pos)]*(a.pos+a.neg))
  return(L)
}

dL3 <- function(beta,C,c,rho,T,eps) {
  a.pos <- rep(sqrt(rowSums(matrix(pmax(0, beta),ncol=T)^2)),T)
  a.neg <- rep(sqrt(rowSums(matrix(pmax(0,-beta),ncol=T)^2)),T)
  
  ## manage cases where alpha+ > 0 or alpha- < 0
  dL <- C %*% beta + c

  neg <- beta < -eps
  pos <- beta >  eps
  dL[pos] <- dL[pos] + rho[pos]*pmax(rep(0,sum(pos)), beta[pos]) / a.pos[pos]
  dL[neg] <- dL[neg] - rho[neg]*pmax(rep(0,sum(neg)),-beta[neg]) / a.neg[neg]

  ## manage cases where beta = 0 and (alpha- = 0 and/or alpha+ = 0)
  n2.theta.pos <- rep(sqrt(rowSums(matrix(pmax(0, dL),ncol=T)^2)),T)
  n2.theta.neg <- rep(sqrt(rowSums(matrix(pmax(0,-dL),ncol=T)^2)),T)

  zero.pos <- abs(beta) < eps & dL >  eps & a.neg < eps
  zero.neg <- abs(beta) < eps & dL < -eps & a.pos < eps

  pos.to.zero <- n2.theta.pos[zero.pos] < rho[zero.pos] + eps
  neg.to.zero <- n2.theta.neg[zero.neg] < rho[zero.neg] + eps

  dL[zero.pos][!pos.to.zero] <- dL[zero.pos][!pos.to.zero] * (1-rho[zero.pos][!pos.to.zero]/n2.theta.pos[zero.pos][!pos.to.zero])
  dL[zero.neg][!neg.to.zero] <- dL[zero.neg][!neg.to.zero] * (1-rho[zero.neg][!neg.to.zero]/n2.theta.neg[zero.neg][!neg.to.zero])
  
  dL[zero.pos][pos.to.zero] <- 0
  dL[zero.neg][neg.to.zero] <- 0
  
  return(dL)
}
