#*************************************************************************
# Initializes random number generator seed. Argument should be a
# positive integer
rngseed <- function(seed)
{
    seed <- as.integer(seed)
    if(seed<=0)stop("'seed' must be a positive integer")
    .Fortran("rngs", seed, PACKAGE = "mix")
    invisible()
}
#*************************************************************************
# Changes NA's to single precision missing value code
.na.to.snglcode <- function(x,mvcode)
{
    x[is.na(x)] <- mvcode
    x
}
#*************************************************************************
# Changes NA's to integer missing value code
.na.to.icode <- function(x,mvcode)
{
    x[is.na(x)] <- as.integer(mvcode)
    x
}
#*************************************************************************
# Changes missing value code to NA
.code.to.na <- function(x,mvcode)
{
    x[x==mvcode] <- NA
    x
}
#*************************************************************************
#  Perform preliminary manipulations on matrix of mixed data.
#  The categorical variables must be in columns 1:p of x,
#  and must be coded as positive integers. The other columns are
#  assumed to be continuous.
prelim.mix <- function(x,p)
{
    x <- as.matrix(x)
# get dimensions of x, separate into w and z
    n <- nrow(x); p <- as.integer(p); q <- as.integer(ncol(x)-p)
    w <- x[,(1:p)]; storage.mode(w) <- "integer"
    z <- x[,(p+1):(p+q)]; storage.mode(z) <- "double"
# store names
   if(!is.null(dimnames(x)[[2]])){
       wnames <- dimnames(x)[[2]][1:p]
       znames <- dimnames(x)[[2]][(p+1):(p+q)]
   } else{
       wnames <- NULL
       znames <- NULL
   }
    rnames <- dimnames(x)[[1]]
    w <- matrix(w,n,p); z <- matrix(z,n,q)
# get dimensions of contingency table
    d <- integer(p)
    for(j in 1:p) d[j] <- as.integer(max(w[!is.na(w[,j]),j]))
# missingness indicators for z
    rz <- 1*is.na(z)
    nmisz <- as.integer(apply(rz,2,sum))
    mdpz <- as.integer((rz%*%(2^((1:q)-1)))+1)
    rz <- 1-rz; storage.mode(rz) <- "integer"
# get missing data patterns for w
    rw <- 1*is.na(w)
    nmisw <- as.integer(apply(rw,2,sum))
    nmis <- c(nmisw,nmisz)
    names(nmis) <- c(wnames,znames)
    mdpw <- as.integer((rw%*%(2^((1:p)-1)))+1)
    rw <- 1-rw; storage.mode(rw) <- "integer"
# calculate the known part of the cell number for rows of w
    cumd <- as.integer(round(exp(cumsum(log(d)))))
    mobs <- 1+((.na.to.icode(w,1)-1)*rw)%*%(cumd/d)
# do row sort
    ro <- order(mdpz,mdpw,mobs)
    w <- matrix(w[ro,],n,p); mdpw <- mdpw[ro]; mobs <- mobs[ro]; rw <- matrix(rw[ro,],n,p)
    z <- matrix(z[ro,],n,q); mdpz <- mdpz[ro]; rz <- matrix(rz[ro,],n,q)
    ro <- order(ro)
# compress missing data patterns
    mdpzst <- as.integer(seq(along=mdpz)[!duplicated(mdpz)])
    mdpz <- unique(mdpz); npattz <- length(mdpz)
    mdpzfin <- as.integer(c(mdpzst[2:npattz]-1,n))
    if(npattz==1)mdpzfin <- n
    mdpzgrp <- numeric(); mdpwst <- numeric(); mdpwtmp <- numeric()
    for(i in 1:npattz){
        tmp <- mdpw[mdpzst[i]:mdpzfin[i]]
        mdpwst <- c(mdpwst,
                    seq(along=tmp)[!duplicated(tmp)]+mdpzst[i]-1)
        tmp <- unique(tmp)
        mdpzgrp <- c(mdpzgrp,length(tmp))
        mdpwtmp <- c(mdpwtmp,tmp)}
    mdpw <- as.integer(mdpwtmp); npattw <- length(mdpw)
    storage.mode(mdpzgrp) <- "integer"
    storage.mode(mdpwst) <- "integer"
    mdpwfin <- as.integer(c(mdpwst[2:npattw]-1,n))
    mdpwgrp <- numeric(); mobsst <- numeric(); mobstmp <- numeric()
    for(i in 1:npattw){
        tmp <- mobs[mdpwst[i]:mdpwfin[i]]
        mobsst <- c(mobsst,
                    seq(along=tmp)[!duplicated(tmp)]+mdpwst[i]-1)
        tmp <- unique(tmp)
        mdpwgrp <- c(mdpwgrp,length(tmp))
        mobstmp <- c(mobstmp,tmp)}
    mobs <- as.integer(mobstmp); ngrp <- length(mobs)
    storage.mode(mdpwgrp) <- "integer"
    storage.mode(mobsst) <- "integer"
# create r-matrix for display purposes
    r <- cbind(rw,rz)[mdpwst,]
    tmp <- as.character(c(mdpwst[2:npattw],n+1)-mdpwst)
    dimnames(r) <- list(tmp,c(wnames,znames))
    ncells <- cumd[p]; jmp <- as.integer(cumd/d)
    rz <- matrix(rz[mdpzst,],npattz,ncol(rz)); rw <- rw[mdpwst,]
# form matrix of packed storage indices
    npsi <- as.integer(q*(q+1)/2)
    psi <- .Fortran("mkpsi",as.integer(q-1),matrix(as.integer(0),q,q),
                    PACKAGE="mix")[[2]]
# center and scale the columns of z
    mvcode <- max(z[!is.na(z)])+1000
    z <- .na.to.snglcode(z,mvcode)
    tmp <- .Fortran("ctrsc",z,n,q,numeric(q),numeric(q),mvcode, PACKAGE="mix")
    z <- tmp[[1]]; xbar <- tmp[[4]]; sdv <- tmp[[5]]
    z <- .code.to.na(z,mvcode)
# return list
    nmobs <- as.integer(c(mobsst[-1],n+1)-mobsst)
    list(w=w,n=n,p=p,d=d,jmp=jmp,z=z,q=q,
         r=r,rz=rz,rw=rw,nmis=nmis,ro=ro,
         mdpzgrp=mdpzgrp,mdpwgrp=mdpwgrp,
         mobs=mobs,mobsst=mobsst,nmobs=nmobs,ncells=ncells,ngrp=ngrp,
         npattz=npattz,npattw=npattw,rnames=rnames,npsi=npsi,psi=psi,
         xbar=xbar,sdv=sdv,wnames=wnames,znames=znames,rnames=rnames)
}
#*************************************************************************
em.mix <- function(s,start,prior=1,maxits=1000,showits=TRUE,eps=.0001)
{
    if(length(prior)==1) prior <- rep(prior,s$ncells)
    prior <- as.double(prior)
    w <- !is.na(prior)
    prior[!w] <- -999.0
    s$z <- .na.to.snglcode(s$z,999)
    tp <- integer(s$p); tq <- integer(s$q)
    sigma <- numeric(s$npsi); pii <- numeric(s$ncells)
    mu <- matrix(0,s$q,s$ncells)
    tmp <- .Fortran("tobsm",s$q,s$psi,s$npsi,sigma,s$ncells,mu,pii,
                    s$npattz,s$rz,s$mdpzgrp,s$npattw,s$p,s$rw,s$mdpwgrp,
                    s$ngrp,s$mobs,s$mobsst,s$nmobs,s$n,s$z,tp,tq,
                    PACKAGE="mix")
    kn1 <- tmp[[4]];kn2 <- tmp[[6]];kn3 <- tmp[[7]]
    if(missing(start)){
        tmp <- .Fortran("stvlm",s$q,s$psi,s$npsi,sigma,s$ncells,mu,
                        PACKAGE="mix")
        sigma <- tmp[[4]];mu <- tmp[[6]]
        pii <- rep(1,s$ncells)
        pii[!w] <- 0}
    if(!missing(start)){
        sigma <- start$sigma; mu <- start$mu; pii <- start$pi
        if(any(pii[w]==0)){
            warning("Starting value on the boundary")}
        if(any(!w)){
            if(any(pii[!w]!=0)){
                stop("starting value has nonzero elements for structural zeros")}}}
    converged <- FALSE; it <- 0
    if(showits) cat(paste("Steps of EM:","\n"))
    t1 <- sigma; t2 <- mu; t3 <- pii
    while((!converged)&(it<maxits)){
        it <- it+1
        if(showits) cat(paste(format(it),"...",sep=""))
        if(it>1){sigma <- t1; mu <- t2; pii <- t3}
        tmp <- .Fortran("estepm",s$q,s$psi,s$npsi,s$ncells,sigma,mu,pii,kn1,
                        kn2,kn3,t1,t2,t3,s$npattz,s$rz,tq,tq,s$mdpzgrp,s$npattw,
                        s$p,s$rw,tp,s$mdpwgrp,s$ngrp,
                        s$mobs,s$mobsst,s$nmobs,s$n,s$z,s$d,s$jmp,tp,kn3,
                        PACKAGE="mix")
        t1 <- tmp[[11]];t2 <- tmp[[12]];t3 <- tmp[[13]]
        tmp <- .Fortran("mstepm",s$q,s$psi,s$npsi,s$ncells,t1,t2,t3,s$n,prior,
                        PACKAGE="mix")
        t1 <- tmp[[5]]; t2 <- tmp[[6]]; t3 <- tmp[[7]]
        if(any(t3<0))
            stop("Estimate outside the parameter space. Check prior.")
        t3[w][t3[w]<(.0000001/sum(w))] <- 0
        c3 <- all(abs(t3-pii)<=(eps*abs(pii)))
        c2 <- all(abs(t2-mu)<=(eps*abs(mu)))
        c1 <- all(abs(t1-sigma)<=(eps*abs(sigma)))
        converged <- c3&c2&c1}
    if(showits) cat("\n")
    list(sigma=t1,mu=t2,pi=t3)
}
#*************************************************************************
loglik.mix <- function(s,theta)
{
    s$z <- .na.to.snglcode(s$z,999)
    tp <- integer(s$p); tq <- integer(s$q)
    ll <- .Fortran("lobsm",s$q,s$psi,s$npsi,s$ncells,theta$sigma,
                   theta$mu,theta$pi,s$npattz,s$rz,tq,tq,s$mdpzgrp,s$npattw,s$p,
                   s$rw,tp,s$mdpwgrp,s$ngrp,s$mobs,s$mobsst,s$nmobs,
                   s$n,s$z,s$d,s$jmp,integer(s$p),ll=numeric(1),
                   PACKAGE="mix")$ll
    ll
}
#*************************************************************************
# Retrieves list of parameters from theta. If corr=F, returns a list
# containing an array of cell probabilities, a matrix of cell means, and
# a covariance matrix. If corr=T, returns a list containing an array of
# cell probabilities, a matrix of cell means, a vector of standard
# deviations, and a correlation matrix.
getparam.mix <- function(s,theta,corr=FALSE)
{
    pii <- array(theta$pi,s$d)
    if(!is.null(s$wnames)){
        pinames <- as.list(1:s$p)
        for(i in 1:s$p){
            pinames[[i]] <- paste(s$wnames[i],"=",format(1:s$d[i]),sep="")}
        dimnames(pii) <- pinames}
    mu <- theta$mu*s$sdv + s$xbar
    dimnames(mu) <- list(s$znames,NULL)
    sigma <- theta$sigma[s$psi]
    sigma <- matrix(sigma,s$q,s$q)
    tmp <- matrix(s$sdv,s$q,s$q)
    sigma <- sigma*tmp*t(tmp)
    dimnames(sigma) <- list(s$znames,s$znames)
    mu[,pii==0] <- NA
    if(corr){
        sdv <- sqrt(diag(sigma))
        names(sdv) <- s$znames
        tmp <- matrix(sdv,s$q,s$q)
        r <- sigma/(tmp*t(tmp)); dimnames(r) <- list(s$znames,s$znames)
        result <- list(pi=pii,mu=mu,sdv=sdv,r=r)}
    else result <- list(pi=pii,mu=mu,sigma=sigma)
    result
}
#*************************************************************************
da.mix <- function(s,start,steps=1,prior=.5,showits=FALSE)
{
    if(length(prior)==1) prior <- rep(prior,s$ncells)
    prior <- as.double(prior)
    w <- !is.na(prior)
    prior[!w] <- -999.0
    s$z <- .na.to.snglcode(s$z,999)
    s$w <- .na.to.icode(s$w,999)
    tp <- integer(s$p); tq <- integer(s$q)
    tmp <- .Fortran("tobsm",s$q,s$psi,s$npsi,numeric(s$npsi),s$ncells,
                    matrix(0,s$q,s$ncells),numeric(s$ncells),
                    s$npattz,s$rz,s$mdpzgrp,s$npattw,s$p,s$rw,s$mdpwgrp,
                    s$ngrp,s$mobs,s$mobsst,s$nmobs,s$n,s$z,tp,tq,
                    PACKAGE="mix")
    kn1 <- tmp[[4]];kn2 <- tmp[[6]];kn3 <- tmp[[7]]
    sigma <- start$sigma; pii <- start$pi; mu <- start$mu
    t1 <- sigma; t2 <- mu; t3 <- pii
    if(showits) cat(paste("Steps of Data Augmentation:","\n"))
    for(i in 1:steps){
        if(showits) cat(paste(format(i),"...",sep=""))
        tmp <- .Fortran("istepm",s$q,s$psi,s$npsi,s$ncells,sigma,mu,pii,kn1,
                        kn2,kn3,t1=t1,t2=t2,t3=t3,s$npattz,s$rz,tq,tq,s$mdpzgrp,
                        s$npattw,s$p,s$rw,tp,s$mdpwgrp,s$ngrp,
                        s$mobs,s$mobsst,s$nmobs,s$n,s$z,s$d,s$jmp,tp,kn3,kn1,s$w,
                        numeric(s$q), PACKAGE="mix")
        t1 <- tmp$t1; t2 <- tmp$t2; t3 <- tmp$t3
        tmp <- .Fortran("pstepm",s$q,s$psi,s$npsi,s$ncells,sigma=t1,mu=t2,pi=t3,
                        s$n,s$p,as.double(prior),numeric(s$npsi),matrix(0,s$q,s$q),
                        numeric(s$q),tq,err=numeric(1), PACKAGE="mix")
    {if(tmp$err==1)
         stop("Improper posterior--empty cells")
    else{
        sigma <- tmp$sigma; mu <- tmp$mu; pii <- tmp$pi}}}
    if(showits) cat("\n")
    mu[,!w] <- NA
    list(sigma=sigma,mu=mu,pi=pii)
}
#*************************************************************************
imp.mix <- function(s,theta,x)
{
    sigma <- theta$sigma
    mu <- theta$mu
    pii <- theta$pi
    z <- .na.to.snglcode(s$z,999)
    w <- .na.to.icode(s$w,999)
    tp <- integer(s$p); tq <- integer(s$q)
    tmp <- .Fortran("istepm",s$q,s$psi,s$npsi,s$ncells,sigma,mu,pii,
                    sigma,mu,pii,sigma,mu,pii,s$npattz,s$rz,tq,tq,s$mdpzgrp,
                    s$npattw,s$p,s$rw,tp,s$mdpwgrp,s$ngrp,
                    s$mobs,s$mobsst,s$nmobs,s$n,z=z,s$d,s$jmp,tp,pii,sigma,w=w,
                    numeric(s$q), PACKAGE="mix")
    w <- tmp$w[s$ro,]
    z <- tmp$z*matrix(s$sdv,s$n,s$q,TRUE)+matrix(s$xbar,s$n,s$q,TRUE)
    z <- z[s$ro,]
    if(!missing(x)){
        zorig <- x[,(s$p+1):(s$p+s$q)]
        z[!is.na(zorig)] <- zorig[!is.na(zorig)]}
    ximp <- cbind(w,z)
    dimnames(ximp) <- list(s$rnames,c(s$wnames,s$znames))
    ximp
}
#*************************************************************************
ecm.mix <- function(s,margins,design,start,prior=1,maxits=1000,
                    showits=TRUE, eps=.0001)
{
    if(length(prior)==1) prior <- rep(prior,s$ncells)
    prior <- as.double(prior)
    w <- !is.na(prior)
    prior[!w] <- -999.0
    storage.mode(design) <- "double"; storage.mode(margins) <- "integer"
    s$z <- .na.to.snglcode(s$z,999)
    tp <- integer(s$p); tq <- integer(s$q)
    sigma <- numeric(s$npsi); pii <- numeric(s$ncells)
    mu <- matrix(0,s$q,s$ncells)
    tmp <- .Fortran("tobsm",s$q,s$psi,s$npsi,sigma,s$ncells,mu,pii,
                    s$npattz,s$rz,s$mdpzgrp,s$npattw,s$p,s$rw,s$mdpwgrp,
                    s$ngrp,s$mobs,s$mobsst,s$nmobs,s$n,s$z,tp,tq,
                    PACKAGE="mix")
    kn1 <- tmp[[4]];kn2 <- tmp[[6]];kn3 <- tmp[[7]]
    if(missing(start)){
        tmp <- .Fortran("stvlm",s$q,s$psi,s$npsi,sigma,s$ncells,mu,
                        PACKAGE="mix")
        sigma <- tmp[[4]];mu <- tmp[[6]]
        pii <- rep(1,s$ncells)
        pii[!w] <- 0}
    if(!missing(start)){
        sigma <- start$sigma; mu <- start$mu; pii <- start$pi
        if(any(pii[w]==0)){
            warning("Starting value on the boundary")}
        if(any(!w)){
            if(any(pii[!w]!=0)){
                stop("Starting value has nonzero elements for structural zeros")}}}
    eps1 <- .0000001*s$n/sum(w)
    r <- ncol(design); npsir <- as.integer(r*(r+1)/2); wk <- numeric(npsir)
    tr <- integer(r); wkr <- numeric(r); wkd <- numeric(s$ncells)
    beta <- matrix(0,r,s$q)
    psir <- .Fortran("mkpsi",as.integer(r-1),matrix(as.integer(0),r,r),
                     PACKAGE="mix")[[2]]
    converged <- FALSE; it <- 0
    if(showits) cat(paste("Steps of ECM:","\n"))
    t1 <- sigma; t2 <- mu; t3 <- pii
    while((!converged) && (it < maxits) )
    {
        it <- it+1
        if(showits) cat(paste(format(it),"...",sep=""))
        if(it>1){sigma <- t1; mu <- t2; pii <- t3}
        tmp <- .Fortran("estepm",s$q,s$psi,s$npsi,s$ncells,sigma,mu,pii,kn1,
                        kn2,kn3,t1,t2,t3,s$npattz,s$rz,tq,tq,s$mdpzgrp,s$npattw,
                        s$p,s$rw,tp,s$mdpwgrp,s$ngrp,
                        s$mobs,s$mobsst,s$nmobs,s$n,s$z,s$d,s$jmp,tp,kn3,
                        PACKAGE="mix")
        t1 <- tmp[[11]];t2 <- tmp[[12]];t3 <- tmp[[13]]
        tmp <- .Fortran("mstepcm",s$q,s$psi,s$npsi,s$ncells,t1,t2,t3,sigma=t1,
                        mu=t2,s$n,r,design,wk,tr,psir,npsir,wkr,wkd,beta,
                        PACKAGE="mix")
        t1 <- tmp$sigma; t2 <- tmp$mu
        t3[w] <- t3[w]+prior[w]-1
        if(any(t3[w]<0))
            stop("Estimate outside the parameter space. Check prior.")
        t3 <- .Fortran("ipf",s$ncells,t3,pii,length(margins),margins,
                       s$p,s$d,s$jmp,tp,tp,tp,eps1, PACKAGE="mix")[[3]]
        t3[w] <- t3[w]/sum(t3[w])
        c3 <- all(abs(t3-pii)<=(eps*abs(pii)))
        c2 <- all(abs(t2-mu)<=(eps*abs(mu)))
        c1 <- all(abs(t1-sigma)<=(eps*abs(sigma)))
        converged <- c3&c2&c1}
    if(showits) cat("\n")
    list(sigma=t1,mu=t2,pi=t3)
}
#*************************************************************************
dabipf.mix <- function(s,margins,design,start,steps=1,prior=.5,showits=FALSE)
{
    if(length(prior)==1) prior <- rep(prior,s$ncells)
    prior <- as.double(prior)
    w <- !is.na(prior)
    prior[!w] <- -999.0
    storage.mode(design) <- "double"; storage.mode(margins) <- "integer"
    s$z <- .na.to.snglcode(s$z,999)
    s$w <- .na.to.icode(s$w,999)
    tp <- integer(s$p); tq <- integer(s$q)
    r <- ncol(design); npsir <- as.integer(r*(r+1)/2); wk <- numeric(npsir)
    tr <- integer(r); wkr <- numeric(r); wkd <- numeric(s$ncells)
    beta <- matrix(0,r,s$q)
    psir <- .Fortran("mkpsi",as.integer(r-1),matrix(as.integer(0),r,r),
                     PACKAGE="mix")[[2]]
    tmp <- .Fortran("tobsm",s$q,s$psi,s$npsi,numeric(s$npsi),s$ncells,
                    matrix(0,s$q,s$ncells),numeric(s$ncells),
                    s$npattz,s$rz,s$mdpzgrp,s$npattw,s$p,s$rw,s$mdpwgrp,
                    s$ngrp,s$mobs,s$mobsst,s$nmobs,s$n,s$z,tp,tq,
                    PACKAGE="mix")
    kn1 <- tmp[[4]];kn2 <- tmp[[6]];kn3 <- tmp[[7]]
    sigma <- start$sigma; pii <- start$pi; mu <- start$mu
    t1 <- sigma; t2 <- mu; t3 <- pii
    if(showits) cat(paste("Steps of Data Augmentation-Bayesian IPF:","\n"))
    for(i in 1:steps){
        if(showits) cat(paste(format(i),"...",sep=""))
        tmp <- .Fortran("istepm",s$q,s$psi,s$npsi,s$ncells,sigma,mu,pii,kn1,
                        kn2,kn3,t1=t1,t2=t2,t3=t3,s$npattz,s$rz,tq,tq,s$mdpzgrp,
                        s$npattw,s$p,s$rw,tp,s$mdpwgrp,s$ngrp,
                        s$mobs,s$mobsst,s$nmobs,s$n,s$z,s$d,s$jmp,tp,kn3,kn1,s$w,
                        numeric(s$q), PACKAGE="mix")
        t1 <- tmp$t1; t2 <- tmp$t2; t3 <- tmp$t3
        tmp <- .Fortran("pstepcm",s$q,s$psi,s$npsi,s$ncells,t1,t2,t3,
                        sigma=t1,mu=t2,s$n,r,design,wk,tr,psir,npsir,
                        wkr,wkd,tq,numeric(s$npsi),beta,matrix(0,s$q,s$q),
                        PACKAGE="mix")
        sigma <- tmp$sigma; mu <- tmp$mu
        tmp <- .Fortran("bipf",s$ncells,t3,pii=pii,prior,
                        length(margins),margins,s$p,s$d,s$jmp,tp,tp,tp,
                        err=as.integer(0), PACKAGE="mix")
        if(tmp$err==1) stop("Improper posterior - Check prior.")
        pii <- tmp$pii; pii[w] <- pii[w]/sum(pii[w])}
    if(showits) cat("\n")
    list(sigma=sigma,mu=mu,pi=pii)
}
#*************************************************************************
# multiple imputation inference
mi.inference <- function(est,std.err,confidence=.95)
{
    qstar <- est[[1]]
    for(i in 2:length(est)){qstar <- cbind(qstar,est[[i]])}
    qbar <- apply(qstar,1,mean)
    u <- std.err[[1]]
    for(i in 2:length(std.err)){u <- cbind(u,std.err[[i]])}
    dimnames(u)[[1]] <- dimnames(qstar)[[1]]
    u <- u^2
    ubar <- apply(u,1,mean)
    bm <- apply(qstar,1,var)
    m <- dim(qstar)[2]
    tm <- ubar+((1+(1/m))*bm)
    rem <- (1+(1/m))*bm/ubar
    nu <- (m-1)*(1+(1/rem))**2
    alpha <- 1-(1-confidence)/2
    low <- qbar-qt(alpha,nu)*sqrt(tm)
    up <- qbar+qt(alpha,nu)*sqrt(tm)
    pval <- 2*(1-pt(abs(qbar/sqrt(tm)),nu))
    fminf <- (rem+2/(nu+3))/(rem+1)
    result <- list(est=qbar,std.err=sqrt(tm),df=nu,signif=pval,lower=low,
                   upper=up,r=rem,fminf=fminf)
    result
}
#***********************************************************************
