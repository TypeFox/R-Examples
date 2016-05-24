#------------------
# Benjamin Risk
# 25 October 2012
# functions for conducting infomax ICA
# Does not constrain W to be orthogonal 

infomaxICA <- function(X, n.comp, W.list = NULL, whiten = FALSE, maxit = 500, eps = 1e-08, alpha.eps = 1e-08, verbose = FALSE, restarts=0) {

  ##Internal functions:
  lInfomax <- function(W,Z) {
    n = ncol(Z)
    S = W%*%Z
    n*log(abs(det(W)))+sum(-S-2*log(1+exp(-S)))
  }
  
  gradInfomax <- function(W,Z) {
    n = ncol(Z)
    S = W%*%Z
    y = 1 / (1+exp(-S))
    solve(t(W)) + ((1 - 2*y) %*% t(Z)) / n
  }
  
  myMixmat <- function (p = 2) {
    a <- matrix(rnorm(p * p), p, p)
    sa <- svd(a)
    d <- sort(runif(p,min=1,max=10))
    mat <- sa$u %*% (sa$v * d)
    attr(mat, "condition") <- d[p]/d[1]
    mat
  }
  
  if(nrow(X) < ncol(X)) stop("Use X = S M parameterization")
  
  p = ncol(X)
  if(is.null(n.comp)) d = ncol(X) else d = n.comp
  
  if(p!=d && whiten==FALSE) stop('Use whitening if p!=n.comp')
  
  runs = restarts + 1
  dsq = d^2

  if(whiten) {
    zData <- whitener(X,n.comp=d)
    Z = zData$Z
    whitener = zData$whitener
  } else {
    Z = scale(X)
    whitener = diag(p)
  }
  Z = t(Z)
  
  
  if(is.null(W.list)) {
    if(whiten)  {
      theta.list = lapply(rep(choose(d,2),runs),runif,min=0,max=2*pi)
      W.list = lapply(theta.list,theta2W) 
    } else {
      W.list = lapply(rep(p,runs),myMixmat)    
    }   
  }
  loglik.v=numeric(runs)
  out.list = NULL
  for(k in 1:runs) {
    w.init = W.list[[k]]
    alpha = 1
    curF = lInfomax(W = w.init, Z = Z)
    deltaW <- gradInfomax(W = w.init, Z = Z)
    normGrad = sqrt(sum(deltaW^2)/(dsq))
    iter = 1
    Table = NULL
    oldW <- w.init
    while (iter < maxit) {
      if ((normGrad < eps) || (alpha<alpha.eps)) 
        break
      alpha <- 2 * alpha
      newW <- oldW + alpha * deltaW
      newF <- lInfomax(newW, Z)
      while (newF <= curF) { # Searches for alpha that reduces f.
        alpha <- alpha/2
        if(alpha < alpha.eps) { 
          warning("alpha is less than alpha.eps -- if norm of gradient is still large, then try a different w.init")
          break
        }
        newW <- oldW + alpha * deltaW
        newF <- lInfomax(newW, Z)
      }  
      deltaW = gradInfomax(W = newW, Z = Z)
      normGrad <- sqrt(sum(deltaW^2)/dsq)
      curF <- newF
      rowTable <- c(iter, curF, log10(normGrad), alpha)
      if(verbose) message('iter: ',rowTable[1],'; newF: ',round(rowTable[2],6),'; log10||grad||: ',round(rowTable[3],6),'; alpha: ',alpha)
      Table <- rbind(Table, rowTable)
      oldW = newW
      iter = iter + 1
    }
    colnames(Table) = c("Iter","f","||Grad||","alpha")
    convergence <- 1*(normGrad < eps)
    if(alpha < alpha.eps && convergence == 0) convergence = 2
    if (convergence==0) 
      warning("convergence not obtained in ", maxit, 
              " iterations used.")
    #else if (convergence==2)
     # warning("check convergence: alpha is less than alpha.eps, so the norm of the gradient is greater than eps but probably sufficiently small")
    colnames(Table)=c('Iter','lInfomax','log10||Grad||','alpha') 
    S = oldW%*%Z
    ##STANDARDIZE:
    dMat <- diag(1/apply(S,1,sd))
    S <- dMat%*%S
    oldW <- dMat%*%oldW
    out.list[[k]] = list(S = t(S), W = t(oldW), M = solve(t(oldW))%*%ginv(whitener),f = curF, Table = Table, convergence = convergence)
    loglik.v[k] = mean(dlogis(S,location=0,scale=sqrt(3)/pi,log=TRUE))
  }
  out.list[[which.max(loglik.v)]]
}
