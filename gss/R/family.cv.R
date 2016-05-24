##%%%%%%%%%%  Binomial Family %%%%%%%%%%

## Calculate CV score for binomial regression
cv.binomial <- function(y,eta,wt,hat,alpha)
{
    if (is.vector(y)) y <- as.matrix(y)
    if (is.null(wt)) wt <- rep(1,dim(y)[1])
    if (dim(y)[2]==1) {
        if ((max(y)>1)|(min(y)<0))
            stop("gss error: binomial responses should be between 0 and 1")
        m <- rep(1,dim(y)[1])
    }
    else {
        if (min(y)<0)
            stop("gss error: paired binomial response should be nonnegative")
        m <- y[,1]+y[,2]
        y <- y[,1]/m
    }
    wtt <- wt * m
    odds <- exp(eta)
    p <- odds/(1+odds)
    w <- p/(1+odds)
    lkhd <- -sum(wtt*(y*eta-log(1+odds)))/sum(wtt)
    aux1 <- sum(hat/w)/(sum(wtt)-sum(hat))
    aux2 <- sum(wtt*y/(1+odds))/sum(wtt)
    list(score=lkhd+abs(alpha)*aux1*aux2,varht=1,w=as.vector(wtt*w))
}


##%%%%%%%%%%  Poisson Family %%%%%%%%%%

## Calculate CV score for Poisson regression
cv.poisson <- function(y,eta,wt,hat,alpha,sr,q)
{
    if (is.null(wt)) wt <- rep(1,length(y))
    if (min(y)<0)
        stop("gss error: Poisson response should be nonnegative")
    nxi <- ncol(q)
    nn <- ncol(sr)
    nnull <- nn-nxi
    lambda <- exp(eta)
    w <- as.vector(lambda)
    lkhd <- -sum(wt*(y*eta-lambda))/sum(wt*y)
    ## matrix H
    mu <- apply(wt*w*sr,2,sum)/sum(wt*w)
    v <- t(sr)%*%(wt*w*sr)/sum(wt*w)-outer(mu,mu)
    v[(nnull+1):nn,(nnull+1):nn] <- v[(nnull+1):nn,(nnull+1):nn]+q/sum(wt*y)
    ## Cholesky decomposition of H
    z <- chol(v,pivot=TRUE)
    v <- z
    rkv <- attr(z,"rank")
    while (v[rkv,rkv]<v[1,1]*sqrt(.Machine$double.eps)) rkv <- rkv-1
    if (rkv<nn) v[(rkv+1):nn,(rkv+1):nn] <- diag(v[1,1],nn-rkv)
    ## trace
    mu <- apply(wt*y*sr,2,sum)/sum(wt*y)
    sr <- sqrt(wt*y)*t(t(sr)-mu)
    sr <- backsolve(v,t(sr[,attr(z,"pivot")]),transpose=TRUE)
    aux1 <- sum(sr^2)
    aux2 <- 1/sum(wt*y)/(sum(wt*y)-1)
    list(score=lkhd+abs(alpha)*aux1*aux2,varht=1,w=as.vector(wt*w))
}


##%%%%%%%%%%  Gamma Family %%%%%%%%%%

## Calculate CV score for Gamma regression
cv.Gamma <- function(y,eta,wt,hat,rss,alpha)
{
    if (is.null(wt)) wt <- rep(1,length(y))
    if (min(y)<=0)
        stop("gss error: gamma responses should be positive")
    mu <- exp(eta)
    u <- 1-y/mu
    lkhd <- sum(wt*(y/mu+eta))/sum(wt)
    aux1 <- sum(hat)/(sum(wt)-sum(hat))
    aux2 <- -sum(wt*y*u/mu)/sum(wt)
    list(score=lkhd+alpha*aux1*aux2,varht=rss/(1-mean(hat)),w=as.vector(wt))
}


##%%%%%%%%%%  Inverse Gaussian Family %%%%%%%%%%

## Calculate CV score for inverse gaussian regression
cv.inverse.gaussian <- function(y,eta,wt,hat,rss,alpha)
{
    if (is.null(wt)) wt <- rep(1,length(y))
    if (min(y)<=0)
        stop("gss error: inverse gaussian responses should be positive")
    mu <- exp(eta)
    u <- (1-y/mu)/mu
    w <- 1/mu
    eta1 <- eta+hat/(1-hat)*u/w
    mu1 <- exp(eta1)
    lkhd <- sum(wt*((y/mu/2-1)/mu))/sum(wt)
    aux1 <- sum(hat/w)/(sum(wt)-sum(hat))
    aux2 <- -sum(wt*y*u/mu/mu)/sum(wt)
#    aux3 <- -sum(wt*y*(1/2/mu/mu-1/2/mu1/mu1))/sum(wt)
    list(score=lkhd+alpha*aux1*aux2,varht=rss/(1-mean(hat)),w=as.vector(wt*w))
}


##%%%%%%%%%%  Negative Binomial Family %%%%%%%%%%

## Calculate CV score for NB regression
cv.nbinomial <- function(y,eta,wt,hat,alpha)
{
    if (is.null(wt)) wt <- rep(1,dim(y)[1])
    if (is.null(offset)) offset <- rep(0,dim(y)[1])
    if (min(y[,1])<0)
        stop("gss error: negative binomial response should be nonnegative")
    if (min(y[,2])<=0)
        stop("gss error: negative binomial size should be positive")
    odds <- exp(eta)
    p <- odds/(1+odds)
    q <- 1/(1+odds)
    u <- y[,1]*p-y[,2]*q
    w <- y[,2]*q
    lkhd <- sum(wt*(-(y[,1]+y[,2])*log(q)-y[,2]*eta))/sum(wt)
    lkhd <- lkhd+sum(wt*(lgamma(y[,2])-lgamma(y[,1]+y[,2])))/sum(wt)
    aux1 <- sum(hat/w)/(sum(wt)-sum(hat))
    aux2 <- sum(wt*y[,1]*p*u)/sum(wt)
    list(score=lkhd+alpha*aux1*aux2,varht=1,w=as.vector(wt*w))
}


##%%%%%%%%%%  Weibull Family %%%%%%%%%%

## Calculate CV score for Weibull regression
cv.weibull <- function(y,eta,wt,hat,nu,alpha)
{
    if (is.vector(y)) stop("gss error: missing censoring indicator")
    if (is.null(wt)) wt <- rep(1,dim(y)[1])
    if (is.null(offset)) offset <- rep(0,dim(y)[1])
    xx <- y[,1]
    delta <- as.logical(y[,2])
    if (dim(y)[2]>=3) zz <- y[,3]
    else zz <- rep(0,length(xx))
    if (any(zz<0)|any(zz>=xx))
        stop("gss error: inconsistent life time data")
    u <- nu*(delta-(xx^nu-zz^nu)*exp(-nu*eta))
    w <- nu^2*(xx^nu-zz^nu)*exp(-nu*eta)
    lkhd <- sum(wt*((xx^nu-zz^nu)*exp(-nu*eta)-delta*(nu*(log(xx)-eta)+log(nu))))/sum(wt)
    aux1 <- sum(hat/w)/(sum(wt)-sum(hat))
    aux2 <- sum(wt*nu*delta*abs(u))/sum(wt)
    list(score=lkhd+alpha*aux1*aux2,varht=1,w=as.vector(wt*w))
}


##%%%%%%%%%%  Log Normal Family %%%%%%%%%%

## Calculate CV score for log normal regression
cv.lognorm <- function(y,eta,wt,hat,nu,alpha)
{
    if (is.vector(y)) stop("gss error: missing censoring indicator")
    if (is.null(wt)) wt <- rep(1,dim(y)[1])
    if (is.null(offset)) offset <- rep(0,dim(y)[1])
    xx <- y[,1]
    delta <- as.logical(y[,2])
    if (dim(y)[2]>=3) zz <- y[,3]
    else zz <- rep(0,length(xx))
    if (any(zz<0)|any(zz>=xx))
        stop("gss error: inconsistent life time data")
    xx <- nu*(log(xx)-eta)
    zz <- nu*(log(zz)-eta)
    s.xx <- ifelse(xx<7,dnorm(xx)/(1-pnorm(xx)),xx+1/xx)
    s.zz <- ifelse(zz<7,dnorm(zz)/(1-pnorm(zz)),zz+1/zz)
    s.xx <- pmax(s.xx,s.zz)
    u <- nu*(delta*(s.xx-xx)-(s.xx-s.zz))
    w <- (s.xx^2/2-xx*s.xx+xx^2/2+log(s.xx)+log(2*pi)/2)
    w <- nu^2*(w-ifelse(s.zz==0,0,(s.zz^2/2-zz*s.zz+zz^2/2+log(s.zz)+log(2*pi)/2)))
    w <- ifelse(w<1e-6,1e-6,w)
    aux1 <- sum(hat/w)/(sum(wt)-sum(hat))
    aux2 <- sum(wt*nu*delta*abs((s.xx-xx)*u))/sum(wt)
    s.xx <- ifelse(xx<7,log(1-pnorm(xx)),-xx^2/2-log(xx+1/xx)-log(2*pi)/2)
    s.zz <- ifelse(zz<7,log(1-pnorm(zz)),-zz^2/2-log(zz+1/zz)-log(2*pi)/2)
    s.xx <- pmin(s.xx,s.zz)
    lkhd <- sum(wt*(delta*(xx^2/2+s.xx-log(nu))+s.zz-s.xx))/sum(wt)
    list(score=lkhd+alpha*aux1*aux2,varht=1,w=as.vector(wt*w))
}


##%%%%%%%%%%  Log Logistic Family %%%%%%%%%%

## Calculate CV score for log logistic regression
cv.loglogis <- function(y,eta,wt,hat,nu,alpha)
{
    if (is.vector(y)) stop("gss error: missing censoring indicator")
    if (is.null(wt)) wt <- rep(1,dim(y)[1])
    if (is.null(offset)) offset <- rep(0,dim(y)[1])
    xx <- y[,1]
    delta <- as.logical(y[,2])
    if (dim(y)[2]>=3) zz <- y[,3]
    else zz <- rep(0,length(xx))
    if (any(zz<0)|any(zz>=xx))
        stop("gss error: inconsistent life time data")
    xx <- 1/(1+exp(nu*(log(xx)-eta)))
    zz <- 1/(1+exp(nu*(log(zz)-eta)))
    u <- nu*(delta*xx-(zz-xx))
    w <- nu^2/2*(zz^2-xx^2)
    lkhd <- sum(wt*(delta*(-log(1-xx)-log(nu))+log(zz)-log(xx)))/sum(wt)
    aux1 <- sum(hat/w)/(sum(wt)-sum(hat))
    aux2 <- sum(wt*nu*delta*xx*abs(u))/sum(wt)
    list(score=lkhd+alpha*aux1*aux2,varht=1,w=as.vector(wt*w))
}
