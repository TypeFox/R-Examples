##%%%%%%%%%%  Binomial Family %%%%%%%%%%

## Make pseudo data for logistic regression
mkdata.binomial <- function(y,eta,wt,offset)
{
    if (is.vector(y)) y <- as.matrix(y)
    if (is.null(wt)) wt <- rep(1,dim(y)[1])
    if (is.null(offset)) offset <- rep(0,dim(y)[1])
    if (dim(y)[2]==1) {
        if ((max(y)>1)|(min(y)<0))
            stop("gss error: binomial responses should be between 0 and 1")
    }
    else {
        if (min(y)<0)
            stop("gss error: paired binomial response should be nonnegative")
        wt <- wt * (y[,1]+y[,2])
        y <- y[,1]/(y[,1]+y[,2])
    }
    odds <- exp(eta)
    p <- odds/(1+odds)
    u <- p - y
    w <- p/(1+odds)
    ywk <- eta-u/w-offset
    wt <- w*wt
    list(ywk=ywk,wt=wt,u=u*wt)
}

## Calculate deviance residuals for logistic regression
dev.resid.binomial <- function(y,eta,wt)
{
    if (is.vector(y)) y <- as.matrix(y)
    if (is.null(wt)) wt <- rep(1,dim(y)[1])
    if (dim(y)[2]>1) {
        wt <- wt * (y[,1]+y[,2])
        y <- y[,1]/(y[,1]+y[,2])
    }
    odds <- exp(eta)
    as.vector(2*wt*(y*log(ifelse(y==0,1,y*(1+odds)/odds))
                    +(1-y)*log(ifelse(y==1,1,(1-y)*(1+odds)))))
}

## Calculate null deviance for logistic regression
dev.null.binomial <- function(y,wt,offset)
{
    if (is.vector(y)) y <- as.matrix(y)
    if (is.null(wt)) wt <- rep(1,dim(y)[1])
    if (dim(y)[2]>1) {
        wt <- wt * (y[,1]+y[,2])
        y <- y[,1]/(y[,1]+y[,2])
    }
    p <- sum(wt*y)/sum(wt)
    odds <- p/(1-p)
    if (!is.null(offset)) {
        eta <- log(odds) - mean(offset)
        repeat {
            odds <- exp(eta+offset)
            p <- odds/(1+odds)
            u <- p - y
            w <- p/(1+odds)
            eta.new <- eta-sum(wt*u)/sum(wt*w)
            if (abs(eta-eta.new)/(1+abs(eta))<1e-7) break
            eta <- eta.new    
        }
    }
    sum(2*wt*(y*log(ifelse(y==0,1,y*(1+odds)/odds))
              +(1-y)*log(ifelse(y==1,1,(1-y)*(1+odds)))))
}


##%%%%%%%%%%  Poisson Family %%%%%%%%%%

## Make pseudo data for Poisson regression
mkdata.poisson <- function(y,eta,wt,offset)
{
    if (is.null(wt)) wt <- rep(1,length(y))
    if (is.null(offset)) offset <- rep(0,length(y))
    if (min(y)<0)
        stop("gss error: Poisson response should be nonnegative")
    lambda <- exp(eta)
    u <- lambda - y
    w <- lambda
    ywk <- eta-u/w-offset
    wt <- w*wt
    list(ywk=ywk,wt=wt,u=u*wt)
}

## Calculate deviance residuals for Poisson regression
dev.resid.poisson <- function(y,eta,wt)
{
    if (is.null(wt)) wt <- rep(1,length(y))
    lambda <- exp(eta)
    as.vector(2*wt*(y*log(ifelse(y==0,1,y/lambda))-(y-lambda)))
}

## Calculate null deviance for Poisson regression
dev.null.poisson <- function(y,wt,offset)
{
    if (is.null(wt)) wt <- rep(1,length(y))
    lambda <- sum(wt*y)/sum(wt)
    if (!is.null(offset)) {
        eta <- log(lambda) - mean(offset)
        repeat {
            lambda <- exp(eta+offset)
            u <- lambda - y
            w <- lambda
            eta.new <- eta-sum(wt*u)/sum(wt*w)
            if (abs(eta-eta.new)/(1+abs(eta))<1e-7) break
            eta <- eta.new    
        }
    }
    sum(2*wt*(y*log(ifelse(y==0,1,y/lambda))-(y-lambda)))
}


##%%%%%%%%%%  Gamma Family %%%%%%%%%%

## Make pseudo data for Gamma regression
mkdata.Gamma <- function(y,eta,wt,offset)
{
    if (is.null(wt)) wt <- rep(1,length(y))
    if (is.null(offset)) offset <- rep(0,length(y))
    if (min(y)<=0)
        stop("gss error: gamma responses should be positive")
    mu <- exp(eta)
    u <- 1-y/mu
    ywk <- eta-u-offset
    list(ywk=ywk,wt=wt,u=u*wt)
}

## Calculate deviance residuals for Gamma regression
dev.resid.Gamma <- function(y,eta,wt)
{
    if (is.null(wt)) wt <- rep(1,length(y))
    mu <- exp(eta)
    as.vector(2*wt*(-log(y/mu)+(y-mu)/mu))
}

## Calculate null deviance for Gamma regression
dev.null.Gamma <-
function(y,wt,offset) {
  if (is.null(wt)) wt <- rep(1,length(y))
  mu <- sum(wt*y)/sum(wt)
  if (!is.null(offset)) {
    eta <- log(mu)-mean(offset)
    repeat {
      mu <- exp(eta+offset)
      u <- 1-y/mu
      eta.new <- eta-sum(wt*u)/sum(wt)
      if (abs(eta-eta.new)/(1+abs(eta))<1e-7) break
      eta <- eta.new    
    }
  }
  sum(2*wt*(-log(y/mu)+(y-mu)/mu))
}


##%%%%%%%%%%  Inverse Gaussian Family %%%%%%%%%%

## Make pseudo data for IG regression
mkdata.inverse.gaussian <- function(y,eta,wt,offset)
{
    if (is.null(wt)) wt <- rep(1,length(y))
    if (is.null(offset)) offset <- rep(0,length(y))
    if (min(y)<=0)
        stop("gss error: inverse gaussian responses should be positive")
    mu <- exp(eta)
    u <- (1-y/mu)/mu
    w <- 1/mu
    ywk <- eta-u/w-offset
    wt <- w*wt
    list(ywk=ywk,wt=wt,u=u*wt)
}

## Calculate deviance residuals for IG regression
dev.resid.inverse.gaussian <- function(y,eta,wt)
{
    if (is.null(wt)) wt <- rep(1,length(y))
    mu <- exp(eta)
    as.vector(wt*((y-mu)^2/(y*mu^2)))
}

## Calculate null deviance for IG regression
dev.null.inverse.gaussian <- function(y,wt,offset)
{
    if (is.null(wt)) wt <- rep(1,length(y))
    mu <- sum(wt*y)/sum(wt)
    if (!is.null(offset)) {
        eta <- log(mu)-mean(offset)
        repeat {
            mu <- exp(eta+offset)
            u <- (1-y/mu)/mu
            w <- 1/mu
            eta.new <- eta-sum(wt*u)/sum(wt*w)
            if (abs(eta-eta.new)/(1+abs(eta))<1e-7) break
            eta <- eta.new    
        }
    }
    sum(wt*((y-mu)^2/(y*mu^2)))
}


##%%%%%%%%%%  Negative Binomial Family %%%%%%%%%%

## Make pseudo data for NB regression
mkdata.nbinomial <- function(y,eta,wt,offset,nu)
{
    if (is.vector(y)) y <- as.matrix(y)
    if (is.null(wt)) wt <- rep(1,dim(y)[1])
    if (is.null(offset)) offset <- rep(0,dim(y)[1])
    if (dim(y)[2]==2) {
        if (min(y[,1])<0)
            stop("gss error: negative binomial response should be nonnegative")
        if (min(y[,2])<=0)
            stop("gss error: negative binomial size should be positive")
        odds <- exp(eta)
        p <- odds/(1+odds)
        q <- 1/(1+odds)
        u <- y[,1]*p-y[,2]*q
        w <- y[,2]*q
        ywk <- eta-u/w-offset
        wt <- w*wt
        list(ywk=ywk,wt=wt,u=u*wt)
    }
    else {
        if (min(y)<0)
            stop("gss error: negative binomial response should be nonnegative")
        odds <- exp(eta)
        p <- odds/(1+odds)
        q <- 1/(1+odds)
        if (is.null(nu)) log.nu <- log(mean(y*odds))
        else log.nu <- log(nu)
        repeat {
            nu <- exp(log.nu)
            ua <- sum(digamma(y+nu)-digamma(nu)+log(p))*nu
            wa <- sum(trigamma(y+nu)-trigamma(nu))*nu*nu+ua
            log.nu.new <- log.nu - ua/wa
            if (abs(log.nu-log.nu.new)/(1+abs(log.nu))<1e-7) break
            log.nu <- log.nu.new
        }
        u <- y*p-nu*q
        w <- nu*q
        ywk <- eta-u/w-offset
        wt <- w*wt
        list(ywk=ywk,wt=wt,nu=nu,u=u*wt)
    }
}

## Calculate deviance residuals for NB regression
dev.resid.nbinomial <- function(y,eta,wt)
{
    if (is.null(wt)) wt <- rep(1,dim(y)[1])
    odds <- exp(eta)
    p <- odds/(1+odds)
    q <- 1/(1+odds)
    as.vector(2*wt*(y[,1]*log(ifelse(y[,1]==0,1,y[,1]/(y[,1]+y[,2])/q))
                    +y[,2]*log(y[,2]/(y[,1]+y[,2])/p)))
}

## Calculate null deviance for NB regression
dev.null.nbinomial <- function(y,wt,offset)
{
    if (is.null(wt)) wt <- rep(1,dim(y)[1])
    p <- sum(wt*y[,2])/sum(wt*y)
    if (!is.null(offset)) {
        eta <- log(p/(1-p)) - mean(offset)
        repeat {
            odds <- exp(eta+offset)
            p <- odds/(1+odds)
            q <- 1/(1+odds)
            u <- y[,1]*p-y[,2]*q
            w <- y[,2]*q
            eta.new <- eta-sum(wt*u)/sum(wt*w)
            if (abs(eta-eta.new)/(1+abs(eta))<1e-7) break
            eta <- eta.new    
        }
    }
    sum(2*wt*(y[,1]*log(ifelse(y[,1]==0,1,y[,1]/(y[,1]+y[,2])/q))
              +y[,2]*log(y[,2]/(y[,1]+y[,2])/p)))
}
