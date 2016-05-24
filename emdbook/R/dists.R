rchibarsq <- function(n,df=1,mix=0.5) {
  ifelse(runif(n)>mix,
         rchisq(n,df),
         if (df==1) 0 else rchisq(n,df-1))
}

dchibarsq <- function(x,df=1,mix=0.5,log=FALSE) {
  df <- rep(df,length.out=length(x))
  mix <- rep(mix,length.out=length(x))
  c1 <- ifelse(df==1,0,dchisq(x,df-1))
  c2 <- dchisq(x,df)
  r <- mix*c1+(1-mix)*c2
  zeros <- (x==0 & df==1)
  if (any(zeros)) {
    r[zeros] <- Inf
  }
  if (log) log(r) else r
}

pchibarsq <- function(p,df=1,mix=0.5,lower.tail=TRUE,log.p=FALSE) {
  df <- rep(df,length.out=length(p))
  mix <- rep(mix,length.out=length(p))
  c1 <- ifelse(df==1,if (lower.tail) 1 else 0,
               pchisq(p,df-1,lower.tail=lower.tail))
  c2 <- pchisq(p,df,lower.tail=lower.tail)
  r <- mix*c1+(1-mix)*c2
  if (log.p) log(r) else r
}

qchibarsq <- function(q,df=1,mix=0.5) {
  n <- max(length(q),length(df),length(mix))
  df <- rep(df,length.out=n)
  mix <- rep(mix,length.out=n)
  q <- rep(q,length.out=n)
  tmpf2 <- function(q,df,mix) {
    if (df>1) {
      tmpf <- function(x) {
        pchibarsq(x,df,mix)-q
      }
      uniroot(tmpf,lower=qchisq(q,df-1),upper=qchisq(q,df))$root
    } else {
      newq <- (q-mix)/(1-mix)
      ifelse(newq<0,0,qchisq(newq,df=1))
    }
  }
  mapply(tmpf2,q,df,mix)
  ##  if (any(df>1)) stop("df>1 not implemented yet")
  ## ?? need uniroot() solution to find quantiles of
  ##   mixtures?
}


dmvnorm <- function (x, mu, Sigma, log = FALSE, tol = 1e-06) {
    if (is.vector(x)) 
        x = t(as.matrix(x))
    n = length(mu)
    if (is.vector(mu)) {
        p <- length(mu)
        if (is.matrix(x)) {
            mu <- matrix(rep(mu, nrow(x)), ncol = p, byrow = TRUE)
        }
    }
    else {
        p <- ncol(mu)
    }
    if (!all(dim(Sigma) == c(p, p)) || nrow(x) != nrow(mu)) 
        stop("incompatible arguments")
    eS <- eigen(Sigma, symmetric = TRUE) 
    ev <- eS$values
    if (!all(ev >= -tol * abs(ev[1]))) 
        stop("Sigma is not positive definite")
    z = t(x - mu)
    logdetS = try(determinant(Sigma, logarithm = TRUE)$modulus)
    attributes(logdetS) <- NULL
    iS = try(solve(Sigma))
    if (class(iS) == "try-error" || class(logdetS) == "try-error") {
        warning("difficulty inverting/taking determinant of Var-Cov matrix")
        return(NA)
    }
    ssq = diag(t(z) %*% iS %*% z)
    loglik = -(n * (log(2*pi)) +  logdetS + ssq)/2
    if (log) loglik else exp(loglik)
  }


nonint <- function(x) (abs((x) - floor((x)+0.5)) > 1e-7)

dbetabinom <- function(x,prob,size,theta,shape1,shape2,log=FALSE) {
  if (missing(prob) && !missing(shape1) && !missing(shape2)) {
    prob <- shape1/(shape1+shape2)
    theta <- shape1+shape2
  }
  v <- lfactorial(size)-lfactorial(x)-lfactorial(size-x)-
      lbeta(theta*(1-prob),theta*prob)+lbeta(size-x+theta*(1-prob),x+theta*prob)
  if (any(n <- nonint(x))) {
      warning("non-integer x detected; returning zero probability")
      v[n] <- -Inf
  }
  if (log) v else exp(v)
}

rbetabinom <- function(n,prob,size,theta,shape1,shape2) {
  if (!missing(prob) && !missing(size) && missing(shape1) && missing(shape2)) {
    shape1 <- theta*prob
    shape2 <- theta*(1-prob)
  }
  rbinom(n,size=size,prob=rbeta(n,shape1,shape2))
}

## could implement distribution function as:
## D(x) = 1 - (n B(b+n-x-1,a+x+1) Gamma(n) 3_F_2(1,a+x+1,-n+x+1;
##   x+2, -b-n+x+2;1))/(B(a,b) B(n-x,x+2) Gamma(n+2))
## (from Mathworld) -- would need to get 3F2 somewhere?
## not in GSL!
## http://tolstoy.newcastle.edu.au/R/e2/help/07/02/10988.html
##  suggests Davies package could be adapted
##  or ?? SuppDists ??

## zero-inflated negative binomial

dzinbinom <- function(x,mu,size,zprob,log=FALSE) {
  logv = log(1-zprob) + dnbinom(x,mu=mu,size=size,log=TRUE)
  logv = ifelse(x==0,log(zprob+exp(logv)),logv)
  if (log) logv else exp(logv)
}

rzinbinom <- function(n,mu,size,zprob) {
  ifelse(runif(n)<zprob,0,rnbinom(n,mu=mu,size=size))
}

## add q, p functions for zinbinom?  where else is zinbinom
## implemented?

