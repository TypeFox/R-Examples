checkArgs <- function(prior, start, niter, pMoves, thin, prop, y){

  # default prior distributions
  ind <- TRUE
  pr <- getPrior(d = (-3+ind))
  if (!missing(prior)) {
    pr[names(prior)] <- prior
  }
  prior <- pr

  pp <- getProp()
  if(!missing(prop)){
    pp[names(prop)] <- prop
  }
  prop <- pp
  
  # find start values
  if(is.null(start)){
    start <- list()
    nJ <- max(round(prior$lambda),3)
    start$jl <- seq(0.001,0.999, length=nJ)
    start$jv <- rep(nJ, nJ)
    start$jh <- rep(1/nJ, nJ)
  }
  
  # default move probabilities
  if(is.null(pMoves)){
    pMoves <- c(UPDT = 1/3, ADD = 1/3, RMV = 1/3)
  }
  if(niter != round(niter)){
    warning("niter needs to be an integer.")
    niter <- round(niter)
  }
  if(length(pMoves) != 3){
    stop("pMoves needs to have length 3.")
  }
  if(sum(pMoves) != 1){
    stop("pMoves needs to sum to 1.")
  }
  if(thin != round(thin)){
    warning("thin needs to be an integer.")
    thin <- round(thin)
  }
  list(prior = prior, niter = niter, start = start,
       prop = prop, pMoves = pMoves, thin = thin,
       niter = niter)
}

getPrior <- function(V = NULL, m = NULL , a = 0, 
                     d, alpha = 1, lambda = 1, vL = 1,
                     vU = 70, la = 1, lb = 1){
  list(V=V, m=m, a=a, d=d, alpha=alpha, lambda=lambda,
      vL=vL, vU=vU, la=la, lb=lb) 
}

getProp <- function(hVar=100, locPr=10, varPr=5){
  list(hVar = hVar, locPr = locPr, varPr = varPr)
}

bnpmr <- function(y, x, prior = NULL, start = NULL, niter = 10000,
                 pMoves = NULL, thin = 1, burnIn = 0, prop = NULL, 
                 seed=1, size = 50){
  x <- (x-min(x))/(max(x) - min(x))
  cargs <- checkArgs(prior, start, niter, pMoves, thin, prop, y)
  prior <- cargs$prior
  niter <- cargs$niter
  pMoves <- cargs$pMoves
  thin <- cargs$thin
  prop <- cargs$prop
  start <- cargs$start
  quadInd <- 0
  
  #size <-ceiling(prior$lambda + 6*sqrt(prior$lambda)) # choose size
  jlOut <- double(length = size*floor((niter-burnIn)/thin))
  jvOut <- double(length = size*floor((niter-burnIn)/thin))
  jhOut <- double(length = size*floor((niter-burnIn)/thin))
  betaOut <- double(length = (2+quadInd)*floor((niter-burnIn)/thin))
  s2Out <- double(length = floor((niter-burnIn)/thin))
  jl <- double(length = size)
  jv <- double(length = size)
  jh <- double(length = size)
  jl[1:length(start$jl)] <- start$jl
  jv[1:length(start$jv)] <- start$jv
  jh[1:length(start$jh)] <- start$jh
  dimCount <- integer(length = floor((niter-burnIn)/thin))

  if(is.null(prior$V)){
    V <- Vinv <- rep(-1, (2+quadInd)^2)
    m <- rep(0,2+quadInd)
  } else {
    Vinv <- c(solve(prior$V))
    m <- c(prior$m)
  }
  
  #double jl[], double jv[], double jh[], int *sizetrans, int *dimtrans, 
  #       int dimCount[], double x[], double ytrans[], int *Ntrans,
  #       double pMoves[], int *oR, double prior[], double Vtrans[], double Vinvtrans[], double mtrans[],
  #       double prop[], int *niter,  int *thin, int *burnIn, int *seed, double jlOut[], double jvOut[], double jhOut[],
  #      double betaOut[], double *s2Out
         
  priorvec <- c(prior$vL, prior$vU, prior$la, prior$lb, prior$alpha, prior$lambda, prior$a, prior$d)
  propvec <- c(prop$hVar, prop$varPr, prop$locPr)

  res <- .C("MH", as.double(jl), as.double(jv), as.double(jh), as.integer(size), as.integer(length(start$jl)),
            dimcount=as.integer(dimCount), as.double(x), as.double(y), as.integer(length(x)), as.double(pMoves),
            as.integer(quadInd),
            as.double(priorvec), as.double(V), as.double(Vinv), as.double(m), as.double(propvec),            
            as.integer(niter), as.integer(thin), as.integer(burnIn), as.integer(seed), 
            jl=as.double(jlOut), jv=as.double(jvOut), jh=as.double(jhOut), beta=as.double(betaOut), s2=as.double(s2Out),
            PACKAGE="bnpmr")
  res
}

pred.bnpmr <- function(x, res){

  jl <- res$jl[res$jl!=0]
  jv <- res$jv[res$jv!=0]
  jh <- res$jh[res$jh!=0]
  s2 <- res$s2
  beta <- res$beta
  beta <- matrix(beta, ncol=2, byrow=TRUE)
  dimcount <- res$dimcount
  N <- length(dimcount)
  cs <- cumsum(c(0,dimcount))
  out <- matrix(ncol = length(x), nrow = N)
  for(n in 1:N){
    val <- rep(0, length(x))
    for(i in (cs[n]+1):cs[n+1]){
      val <- val + jh[i]*ptsp(x,jl[i],jv[i])
    }
    out[n,] <- beta[n,1] + beta[n,2]*val
  }
  out
}

ptsp <- function(x, m, n){
  res <- numeric(length(x))
  ind <- x <= m
  res[ind] <- m*(x[ind]/m)^n
  res[!ind] <- 1-(1-m)*((1-x[!ind])/(1-m))^n
  res
}

