
kVB <- 0x001000
kVBMacKay <- 0x001001
kPXVB <- 0x001002

kDefaultMethod <- kVB

kDefaultMaxIt <- 100000
kDefaultPruning <- 1e+8
kDefaultReltol <- sqrt(.Machine$double.eps)

## faster than 1 / (1 + exp(-x))
sigmoid <- function(x) return(0.5 * tanh(0.5 * x) + 0.5)

## Sparse logistic regression with automatic relevance determination
SlrArd <- function(T, X, bias=TRUE, method=c("VB", "VBMacKay", "PX-VB"),
                   control=list(), check.lb=TRUE, a0=0, b0=0, mu0=0, xi0=1){
  ## Check method
  if(any("VB" == method))
    nmethod <- kVB
  else if(any("VBMacKay" == method))
    nmethod <- kVBMacKay
  else if(any("PX-VB" == method))
    nmethod <- kPXVB
  else
    nmethod <- kDefaultMethod

  ## Check input
  if(is.matrix(T))
    T <- c(T)

  N <- length(T)
  if(N != nrow(X))
    stop("length(T) and nrow(X) must be same.")

  if(bias){
    X <- cbind("Bias"=1, X)
  }
  M <- ncol(X)

  cnames <- colnames(X)

  ## Check control options
  valid <- match(names(control), c("maxit", "reltol", "trace", "REPORT", "pruning"))
  if(any(is.na(valid)))
    stop(cat("unknown control options:",
             sprintf("%s", names(control)[is.na(valid)])))

  if(is.null(control$maxit))
    control$maxit <- kDefaultMaxIt
  if(is.null(control$reltol))
    control$reltol <- kDefaultReltol
  if(is.null(control$trace))
    control$trace <- TRUE
  if(is.null(control$REPORT))
    control$REPORT <- floor(control$maxit / 20)
  if(is.null(control$pruning))
    control$pruning <- kDefaultPruning

  ## Initialize parameters
  mu <- rep(mu0, M)
  S <- matrix(nrow=M, ncol=M, 0)
  diag(S) <- 1
  xi <- rep(xi0, N)
  lambda <- 0.5 * (sigmoid(xi) - 0.5) / xi

  a <- a0 + 0.5
  if(b0 == 0)
    eta <- rep(1, M)
  else
    eta <- rep(a0 / b0, M)

  ## Indices of active parameters
  is <- 1:M

  ## Pre-processing
  tt <- (T - 0.5) * X
  tmp <- colSums(tt)
  XX <- t(apply(X=X, MARGIN=1, FUN=function(x) cbind(x) %*% x))

  if(check.lb){
    lb.div <- FALSE
    lb.diff <- rep(0, control$maxit)
  }

  converged <- FALSE
  LB <- -Inf
  for(i in 1:control$maxit){
    ## Show trace
    if(control$trace){
      if(i %% control$REPORT == 0){
        print(sprintf("Iter#%d: %d parameters remain, lower bound = %f + const.",
                      i, length(is), LB))
      }
    }

    ## Update parameters
    lXX <- colSums(lambda*XX)
    S <- solve(diag(eta) + 2 * lXX)
    mu <- S %*% tmp
    if(nmethod == kVB){
      ## VB
      b <- b0 + 0.5 * (diag(S) + c(mu^2))
      eta <- a / b
    }
    else if(nmethod == kVBMacKay){
      ## VB with effective number of well-determined parameters (MacKay, 1992)
      eta <- (1 - eta*(diag(S)) + 2*a0) / (c((mu)^2) + 2*b0)
      b <- a / eta
    }
    else if(nmethod == kPXVB){
      ## PX-VB
      b <- b0 + 0.5 * (diag(S) + c(mu^2))
      eta <- a / b

      pl <- polyroot(c(-2*b0*sum(eta), 0, 2*M*a0, -c(mu) %*% tmp,
                       2*sum((S + mu%*%t(mu))*lXX)))
      pl <- Re(pl[Re(pl) > 0 & abs(Im(pl)) < 1e-8])
      if(length(pl) == 0){
        stop("No real roots found. Please use VB instead.")
      }
      else if(length(pl) > 1){
        warning("Multiple roots found. First one is used.")
      }
      c2 <- pl[1]

      b <- c2 * b
      eta <- eta / c2

      S <- c2 * S
      mu <- sqrt(c2) * mu
    }
    else{
      stop("never reach")
    }

    ## Prune ineffective parameters
    eff <- eta < control$pruning

    if(!all(eff)){
      if(control$trace){
        print(sprintf("Iter#%d: pruned %s", i, do.call(paste, as.list(is[!eff]))))
      }

      eta <- eta[eff]
      tt <- tt[,eff]
      tmp <- colSums(tt)
      X <- X[,eff]
      mu <- mu[eff]
      b <- b[eff]
      is <- is[eff]

      XX <- XX[,as.logical(c(cbind(eff) %*% eff))]
      S <- matrix(nrow=length(is), S[as.logical(c(cbind(eff) %*% eff))])

##      LB <- -Inf
    }

    xi <- sqrt(rowSums((X %*% (S + mu %*% t(mu))) * X))
    sxi <- sigmoid(xi)
    lambda <- 0.5 * (sxi - 0.5) / xi

    ## Check lower-bound
    nLB <- sum(log(sxi) + tt %*% mu - 0.5 * xi) - a * sum(log(b)) - 0.5*log(det(S))

    if(check.lb){
      if(is.infinite(nLB))
        lb.div <- TRUE
      lb.diff[i] <- nLB - LB
    }

    if(!is.infinite(nLB) &&
       abs(nLB - LB) < control$reltol * (abs(nLB) + control$reltol)){
      converged <- TRUE
      LB <- nLB
      break
    }

    LB <- nLB
  }

  ## Check error
  if(check.lb){
    if(lb.div)
      warning(sprintf("Sometimes the lower bound diverges. Check lb.diff."))
    if(any(lb.diff < -1e-12)){
      warning(sprintf("Sometimes the lower bound decreases. Check lb.diff."))
    }
  }

  ## Return results
  res <- list(coefficients=rep(0, M), irrelevance=rep(Inf, M),
              iterations=i, converged=converged, lower.bound=LB, method=method,
              fitted.values=c(sigmoid(X %*% mu)))
  res$coefficients[is] <- mu
  res$irrelevance[is] <- eta
  res$residuals <- T - res$fitted.values
  names(res$coefficients) <- names(res$irrelevance) <- cnames

  if(check.lb)
    res$lb.diff <- lb.diff

  return(res)
}
