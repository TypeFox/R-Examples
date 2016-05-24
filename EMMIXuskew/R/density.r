#
#  EM algorithm for Mixture of Unrestricted Multivariate Skew t-distributioins
#  Package: EMMIX-uskew
#  Version: 0.11-5
#
#  Code by S. X. Lee
#  Updated on 14 Oct 2013
#
# Lee S. and Mclachlan, G.J. (2010) On the fitting of finite mixtures of
#   multivariate skew t-distributions via the EM algorithm
#

################################################################################
#  SECTION 2
#                 Multivariate Skew-t Mixture Densities
#
################################################################################
#

dmsn <- function(dat, mu=NULL, sigma=NULL, delta=NULL, known=NULL) {
    if(missing(dat)) stop("dat must be provided")
    n <- nrow(dat);  p <- ncol(dat);                
    if(is.null(p)) {p<-length(dat); n<-1; dat <- matrix(as.numeric(dat),1,p)}
    if(!is.matrix(dat)) dat <- as.matrix(dat)
    if(!is.null(known)) {
        mu <- known$mu
        sigma <- known$sigma
        delta <- known$delta
    } 
    if(is.null(mu)) mu <- rep(0,p)
    if(is.null(delta)) delta <- rep(0,p)
    if(is.null(sigma)) sigma<-diag(p)
    mu <- as.numeric(mu)
    delta <- as.numeric(delta)
    Delta <- diag(p); diag(Delta) <- delta            
    Ydiff <- t(dat) - mu%*% matrix(1,1,n)             
    Omega <- sigma + Delta %*% t(Delta)               
    InvOmega <- solve(Omega)                          
    Lambda <- diag(p) - Delta %*% InvOmega %*% Delta  
    Q <- Delta %*% InvOmega %*% Ydiff                 
    N1 <- dmn(dat, mu=as.numeric(mu), Sigma=Omega, log=FALSE);  N2 <- N1     
    for (i in 1:n) {
        N2[i] <- sadmvn(lower=rep(-Inf,p), upper=Q[,i], mean=matrix(rep(0,p),p), varcov=Lambda) 
    }
    SN <- 2^p * N1 * N2
    return(SN)
}


dmst <- function(dat, mu=NULL, sigma=NULL, delta=NULL, dof=1, known=NULL, tmethod=1) {
    if(dof==Inf) return(dmsn(dat, mu, sigma, delta, known))
    if(missing(dat)) stop("dat must be provided")
    n <- nrow(dat);  p <- ncol(dat);                   
    if(is.null(p)) {p<-length(dat); n<-1; dat <- matrix(as.numeric(dat),1,p)}
    if(!is.matrix(dat)) dat <- as.matrix(dat)
    if(!is.null(known)) {
        mu <- known$mu
        sigma <- known$sigma
        delta <- known$delta
        dof <- known$dof
    } 
    if(is.null(mu)) mu <- rep(0,p)
    if(is.null(delta)) delta <- rep(0,p)
    if(is.null(sigma)) sigma<-diag(p)
    mu <- as.numeric(mu)
    delta <- as.numeric(delta)
    Delta <- diag(p); diag(Delta) <- delta            
    Ydiff <- t(dat) - mu%*% matrix(1,1,n)             
    Omega <- sigma + Delta %*% t(Delta)               
    InvOmega <- solve(Omega)                          
    Lambda <- diag(p) - Delta %*% InvOmega %*% Delta  
    Q <- Delta %*% InvOmega %*% Ydiff                 
    eta <- mahalanobis(dat, as.numeric(mu), InvOmega, T)  
    N <- sqrt((dof+p)/abs(dof+eta))                   
    T1 <- dmt(dat, as.numeric(mu), Omega, dof);  T2 <- T1     
    for (i in 1:n) T2[i] <- pmt(N[i]*Q[,i], matrix(rep(0,p),p), Lambda, round(dof+p), tmethod)
    ST <- 2^p * T1 * T2
    return(ST)
}

dfmmst <- function(dat, mu=NULL, sigma=NULL, delta=NULL, dof=NULL, pro=NULL, known=NULL, tmethod=1) {
    n <- nrow(dat);  p <- ncol(dat)
    if(is.null(p)) {p<-length(dat); n<-1; dat <- matrix(as.numeric(dat),1,p)}
    if(!is.matrix(dat)) dat <- as.matrix(dat)  
    if(!is.null(known)) {
        mu <- known$mu
        sigma <- known$sigma
        delta <- known$delta
        dof <- known$dof
        pro <- known$pro
    }
    if(missing(dof)) stop("dof is missing.")
    g <- length(dof)
    if(g==1) return(dmst(dat, known=as.fmmst(1, mu, sigma, delta, dof), tmethod=tmethod)) 
    if(is.null(dof)) dof <- rep(1,g)
    if(is.null(pro)) pro <- rep(1/g,g)
    ST <- rep(0, length=n)
    for (i in 1:g) ST <- ST + pro[i]*dmst(dat, as.numeric(mu[[i]]), sigma[[i]], as.numeric(delta[[i]]), round(dof[i]), tmethod=tmethod)
    return(ST)
}


dmt <- function(dat, mu, sigma, dof = Inf, log = FALSE) {
    if (dof == Inf)  return(dmn(dat, mu, sigma, log = log))
    d  <- if(is.matrix(sigma)) ncol(sigma) else 1
    n  <- if(d==1) length(dat) else nrow(dat)
    X <- t(matrix(dat, nrow = n, ncol = d)) - mu
    Q <- apply((solve(sigma) %*% X) * X, 2, sum)
    logDet <- sum(logb(abs(diag(qr(sigma)$qr))))
    logPDF <- (lgamma((dof + d)/2) - 0.5 * (d * logb(pi * dof) + logDet)
               - lgamma(dof/2) - 0.5 * (dof + d) * logb(1 + Q/dof))
    if(log) logPDF else exp(logPDF)
}

pmt <- function(dat, mu=rep(0,length(dat)), sigma=diag(length(dat)), dof=Inf, method=1, ...){
  if (dof == Inf)   return(sadmvn(lower=rep(-Inf,length(dat)), upper=dat, mean=mu, varcov=sigma))      
  if(dof%%1 != 0) {
      if(method==2) return(mtcdf2(dat, mu, sigma, dof, ...))
      if(method==3) return(mtcdf(dat, mu, sigma, dof, ...))
  } 
  dof <- round(dof)
  if(length(dat) == 2) 
    pmt.biv(dof, lower=rep(-Inf, 2), upper=dat, mu, sigma) 
  else  
    pmt.mtv(dof, lower=rep(-Inf, length(dat)), upper=dat, mu, sigma)  
}


dfmmt <- function(dat, mu=NULL, sigma=NULL, dof=NULL, pro=NULL, known=NULL) {
    n <- nrow(dat);  p <- ncol(dat)
    if(is.null(p)) {p<-length(dat); n<-1; dat <- matrix(as.numeric(dat),1,p)}
    if(!is.matrix(dat)) dat <- as.matrix(dat)  
    if(!is.null(known)) {
        mu <- known$mu
        sigma <- known$sigma
        dof <- known$dof
        pro <- known$pro
    }
    if(missing(dof)) stop("dof is missing.")
    g <- length(dof)
    if(g==1) return(dmt(dat, mu, sigma, dof)) 
    if(is.null(dof)) dof <- rep(1,g)
    if(is.null(pro)) pro <- rep(1/g,g)
    MT <- rep(0, length=n)
    for (i in 1:g) MT <- MT + pro[i]*dmt(dat, as.numeric(mu[[i]]), sigma[[i]], round(dof[i]))
    return(MT)
}    

