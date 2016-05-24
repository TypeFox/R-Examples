##%%%%%%%%%%  Weibull Family %%%%%%%%%%

## Make pseudo data for Weibull regression
mkdata.weibull <- function(y,eta,wt,offset,nu)
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
    if (nu[[2]]) {
        lkhd <- function(log.nu) {
            nu <- exp(log.nu)
            -sum(wt*(delta*(nu*(log(xx)-eta)+log.nu)
                 -(xx^nu-zz^nu)*exp(-nu*eta)))
        }
        if (is.null(nu[[1]])) nu[[1]] <- 1
        nu[[1]] <- exp(nlm(lkhd,log(nu[[1]]),stepmax=.5)$est)
    }
    u <- nu[[1]]*(delta-(xx^nu[[1]]-zz^nu[[1]])*exp(-nu[[1]]*eta))
    w <- nu[[1]]^2*(xx^nu[[1]]-zz^nu[[1]])*exp(-nu[[1]]*eta)
    ywk <- eta-u/w-offset
    wt <- w*wt
    list(ywk=ywk,wt=wt,nu=nu,u=u*wt)
}

## Calculate deviance residuals for Weibull regression
dev.resid.weibull <- function(y,eta,wt,nu)
{
    if (is.null(wt)) wt <- rep(1,dim(y)[1])
    xx <- y[,1]
    delta <- as.logical(y[,2])
    if (dim(y)[2]>=3) zz <- y[,3]
    else zz <- rep(0,length(xx))
    z <- -delta*nu*(log(xx)-eta)+(xx^nu-zz^nu)*exp(-nu*eta)
    as.numeric(2*wt*(z+delta*(log(xx^nu)-log(xx^nu-zz^nu)-1)))
}

## Calculate null deviance for Weibull regression
dev.null.weibull <- function(y,wt,offset,nu)
{
    if (is.null(wt)) wt <- rep(1,dim(y)[1])
    xx <- y[,1]
    delta <- as.logical(y[,2])
    if (dim(y)[2]>=3) zz <- y[,3]
    else zz <- rep(0,length(xx))
    if (is.null(offset)) offset <- rep(0,length(xx))
    eta <- log(sum(wt*(xx^nu-zz^nu)*exp(-nu*offset))/sum(wt*delta))/nu
    eta <- eta + offset    
    z <- -delta*nu*(log(xx)-eta)+(xx^nu-zz^nu)*exp(-nu*eta)
    sum(2*wt*(z+delta*(log(xx^nu)-log(xx^nu-zz^nu)-1)))
}


##%%%%%%%%%%  Log Normal Family %%%%%%%%%%

## Make pseudo data for log normal regression
mkdata.lognorm <- function(y,eta,wt,offset,nu)
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
    if (nu[[2]]) {
        lkhd <- function(log.nu) {
            nu <- exp(log.nu)
            xx.wk <- nu*(log(xx)-eta)
            zz.wk <- nu*(log(zz)-eta)
            -sum(wt*(delta*(-xx.wk^2/2-log(1-pnorm(xx.wk))+log.nu)
                 +log((1-pnorm(xx.wk))/(1-pnorm(zz.wk)))))
        }
        if (is.null(nu[[1]])) nu[[1]] <- 1
        nu[[1]] <- exp(nlm(lkhd,log(nu[[1]]),stepmax=.5)$est)
    }
    xx <- nu[[1]]*(log(xx)-eta)
    zz <- nu[[1]]*(log(zz)-eta)
    s.xx <- ifelse(xx<7,dnorm(xx)/(1-pnorm(xx)),xx+1/xx)
    s.zz <- ifelse(zz<7,dnorm(zz)/(1-pnorm(zz)),zz+1/zz)
    s.xx <- pmax(s.xx,s.zz)
    u <- nu[[1]]*(delta*(s.xx-xx)-(s.xx-s.zz))
    w <- (s.xx^2/2-xx*s.xx+xx^2/2+log(s.xx)+log(2*pi)/2)
    w <- nu[[1]]^2*(w-ifelse(s.zz==0,0,(s.zz^2/2-zz*s.zz+zz^2/2+log(s.zz)+log(2*pi)/2)))
    w <- ifelse(w<1e-6,1e-6,w)
    ywk <- eta-u/w-offset
    wt <- w*wt
    list(ywk=ywk,wt=wt,nu=nu,u=u*wt)
}

## Calculate deviance residuals for log normal regression
dev.resid.lognorm <- function(y,eta,wt,nu)
{
    if (is.null(wt)) wt <- rep(1,dim(y)[1])
    xx <- y[,1]
    delta <- as.logical(y[,2])
    if (dim(y)[2]>=3) zz <- y[,3]
    else zz <- rep(0,length(xx))
    dev0 <- NULL
    for (i in 1:length(xx)) {
        if (!delta[i]|!zz[i]) dev0 <- c(dev0,0)
        else {
            fun.wk <- function(eta) {
                (nu*(log(xx[i])-eta))^2/2+log(1-pnorm(nu*(log(zz[i])-eta)))
            }
            dev0 <- c(dev0,nlm(fun.wk,log(xx[i]),stepmax=1)$min)
        }
    }
    xx <- nu*(log(xx)-eta)
    zz <- nu*(log(zz)-eta)
    s.xx <- log(1-pnorm(xx))
    s.zz <- log(1-pnorm(zz))
    z <- -delta*(-xx^2/2-s.xx)-s.xx+s.zz
    as.numeric(2*wt*(z-dev0))
}

dev0.resid.lognorm <- function(y,eta,wt,nu)
{
    if (is.null(wt)) wt <- rep(1,dim(y)[1])
    xx <- y[,1]
    delta <- as.logical(y[,2])
    if (dim(y)[2]>=3) zz <- y[,3]
    else zz <- rep(0,length(xx))
    xx <- nu*(log(xx)-eta)
    zz <- nu*(log(zz)-eta)
    s.xx <- ifelse(xx<7,log(1-pnorm(xx)),-xx^2/2-log(xx+.15)-log(2*pi)/2)
    s.zz <- ifelse(zz<7,log(1-pnorm(zz)),-zz^2/2-log(zz+.15)-log(2*pi)/2)
    s.xx <- pmin(s.xx,s.zz)
    z <- -delta*(-xx^2/2-s.xx)-s.xx+s.zz
    as.numeric(2*wt*z)
}

## Calculate null deviance for log normal regression
dev.null.lognorm <- function(y,wt,offset,nu)
{
    if (is.null(wt)) wt <- rep(1,dim(y)[1])
    xx <- y[,1]
    delta <- as.logical(y[,2])
    if (dim(y)[2]>=3) zz <- y[,3]
    else zz <- rep(0,length(xx))
    dev0 <- NULL
    for (i in 1:length(xx)) {
        if (!delta[i]|!zz[i]) dev0 <- c(dev0,0)
        else {
            fun.wk <- function(eta) {
                (nu*(log(xx[i])-eta))^2/2+log(1-pnorm(nu*(log(zz[i])-eta)))
            }
            dev0 <- c(dev0,nlm(fun.wk,log(xx[i]),stepmax=1)$min)
        }
    }
    if (is.null(offset)) offset <- rep(0,length(xx))
    lkhd <- function(eta) {
        eta <- eta + offset
        xx.wk <- nu*(log(xx)-eta)
        zz.wk <- nu*(log(zz)-eta)
        -sum(wt*(delta*(-xx.wk^2/2-log(1-pnorm(xx.wk)))
                 +log((1-pnorm(xx.wk))/(1-pnorm(zz.wk)))))
    }
    eta <- nlm(lkhd,mean(log(xx)-offset),stepmax=1)$est + offset
    xx <- nu*(log(xx)-eta)
    zz <- nu*(log(zz)-eta)
    z <- -delta*(-xx^2/2-log(1-pnorm(xx)))-log((1-pnorm(xx))/(1-pnorm(zz)))
    sum(2*wt*(z-dev0))
}


##%%%%%%%%%%  Log Logistic Family %%%%%%%%%%

## Make pseudo data for log logistic regression
mkdata.loglogis <- function(y,eta,wt,offset,nu)
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
    if (nu[[2]]) {
        lkhd <- function(log.nu) {
            nu <- exp(log.nu)
            xx.wk <- nu*(log(xx)-eta)
            zz.wk <- nu*(log(zz)-eta)
            -sum(wt*(delta*(xx.wk-log(1+exp(xx.wk))+log.nu)
                 -log((1+exp(xx.wk))/(1+exp(zz.wk)))))
        }
        if (is.null(nu[[1]])) nu[[1]] <- 1
        nu[[1]] <- exp(nlm(lkhd,log(nu[[1]]),stepmax=.5)$est)
    }
    xx <- 1/(1+exp(nu[[1]]*(log(xx)-eta)))
    zz <- 1/(1+exp(nu[[1]]*(log(zz)-eta)))
    u <- nu[[1]]*(delta*xx-(zz-xx))
    w <- nu[[1]]^2/2*(zz^2-xx^2)
    w <- pmax(w,1e-6)
    ywk <- eta-u/w-offset
    wt <- w*wt
    list(ywk=ywk,wt=wt,nu=nu,u=u*wt)
}

## Calculate deviance residuals for log logistic regression
dev.resid.loglogis <- function(y,eta,wt,nu)
{
    if (is.null(wt)) wt <- rep(1,dim(y)[1])
    xx <- y[,1]
    delta <- as.logical(y[,2])
    if (dim(y)[2]>=3) zz <- y[,3]
    else zz <- rep(0,length(xx))
    dev0 <- NULL
    for (i in 1:length(xx)) {
        if (!delta[i]) dev0 <- c(dev0,0)
        else {
            if (!zz[i]) dev0 <- c(dev0,2*log(2))
            else {
                if ((xx[i]/zz[i])^nu<=2) dev0 <- c(dev0,nu*log(xx[i]/zz[i]))
                else dev0 <- c(dev0,2*log(2)-log(xx[i]^nu/(xx[i]^nu-zz[i]^nu)))
            }
        }
    }
    xx <- nu*(log(xx)-eta)
    zz <- nu*(log(zz)-eta)
    z <- -delta*(xx-log(1+exp(xx)))+log((1+exp(xx))/(1+exp(zz)))
    as.numeric(2*wt*(z-dev0))
}

dev0.resid.loglogis <- function(y,eta,wt,nu)
{
    if (is.null(wt)) wt <- rep(1,dim(y)[1])
    xx <- y[,1]
    delta <- as.logical(y[,2])
    if (dim(y)[2]>=3) zz <- y[,3]
    else zz <- rep(0,length(xx))
    xx <- nu*(log(xx)-eta)
    zz <- nu*(log(zz)-eta)
    z <- -delta*(xx-log(1+exp(xx)))+log((1+exp(xx))/(1+exp(zz)))
    as.numeric(2*wt*z)
}

## Calculate null deviance for log logistic regression
dev.null.loglogis <- function(y,wt,offset,nu)
{
    if (is.null(wt)) wt <- rep(1,dim(y)[1])
    xx <- y[,1]
    delta <- as.logical(y[,2])
    if (dim(y)[2]>=3) zz <- y[,3]
    else zz <- rep(0,length(xx))
    dev0 <- NULL
    for (i in 1:length(xx)) {
        if (!delta[i]) dev0 <- c(dev0,0)
        else {
            if (!zz[i]) dev0 <- c(dev0,2*log(2))
            else {
                if ((xx[i]/zz[i])^nu<=2) dev0 <- c(dev0,nu*log(xx[i]/zz[i]))
                else dev0 <- c(dev0,2*log(2)-log(xx[i]^nu/(xx[i]^nu-zz[i]^nu)))
            }
        }
    }
    if (is.null(offset)) offset <- rep(0,length(xx))
    lkhd <- function(eta) {
        eta <- eta + offset
        xx.wk <- nu*(log(xx)-eta)
        zz.wk <- nu*(log(zz)-eta)
        -sum(wt*(delta*(xx.wk-log(1+exp(xx.wk)))
                 -log((1+exp(xx.wk))/(1+exp(zz.wk)))))
    }
    eta <- nlm(lkhd,mean(log(xx)-offset))$est + offset
    xx <- nu*(log(xx)-eta)
    zz <- nu*(log(zz)-eta)
    z <- -delta*(xx-log(1+exp(xx)))+log((1+exp(xx))/(1+exp(zz)))
    sum(2*wt*(z-dev0))
}
