#
#  EM algorithm for Mixture of Unrestricted Multivariate Skew t-distributioins
#  Package: EMMIX-uskew
#  Version: 0.11-5
#
#  Code by S. X. Lee
#  Updated on 31 Jan, 2012
#
# Lee S. and Mclachlan, G.J. (2010) On the fitting of finite mixtures of
#   multivariate skew t-distributions via the EM algorithm
#
# The following code are adapted from other packages

################################################################################
#  SECTION 6
#                           Miscellaneous Tools
#
################################################################################

#skewness
skewness <- function (x, na.rm = FALSE){
    if (is.matrix(x))
        apply(x, 2, skewness, na.rm = na.rm)
    else if (is.vector(x)) {
    if (na.rm) x <- x[!is.na(x)]
	       n <- length(x)
        (sum((x-mean(x))^3)/n)/(sum((x-mean(x))^2)/n)^(3/2)
    }
    else if (is.data.frame(x))
        sapply(x, skewness, na.rm = na.rm)
    else skewness(as.vector(x), na.rm = na.rm)
}

#mahalanobis
mahalanobis. <- function(x, center, cov, inverted=FALSE, ...) {
    x <- if(is.vector(x)) matrix(x, ncol=length(x)) else as.matrix(x)
    x <- t(sweep(x, 2, center))# = (x - center)
    retval <- colSums(x * if(inverted) cov%*%x else solve(cov, x, ...))
    names(retval) <- rownames(x)
    retval
}

mahalanobis <- function(x, center, cov, inverted=FALSE, ...){
    x <- if(is.vector(x)) matrix(x, ncol=length(x)) else as.matrix(x)
    x <- sweep(x, 2, center)# = (x - center)
    if(!inverted)
	cov <- solve(cov, ...)
    retval <- rowSums((x%*%cov) * x)
    names(retval) <- rownames(x)
    retval
}

#
isPosDef <- function(M) { 
    if ( all(M == t(M) ) ) { 
        if (  all(eigen(M)$values > 0) ) {TRUE}
        else {FALSE} 
    }else {FALSE}  
} 

#
is.whole <- function(a) { 
   (is.numeric(a) && floor(a)==a) ||
   (is.complex(a) && floor(Re(a)) == Re(a) && floor(Im(a)) == Im(a))
}

#bivariate t-distribution
pmt.biv <- function(dof, lower, upper, mu, Sigma){
    if(any(dim(Sigma) != c(2,2))) stop("dimensions mismatch")
    if(length(mu) != 2) stop("dimensions mismatch") 
    if(round(dof) != dof) warning("non integer dof is rounded to integer") 
    nu <- if(dof<Inf) as.integer(round(dof)) else 0
    if(dof==Inf) nu <- 0
    sd <- sqrt(diag(Sigma))
    rho <- cov2cor(Sigma)[1,2]
    lower <- as.double((lower-mu)/sd)
    upper <- as.double((upper-mu)/sd)
    if(any(lower > upper)) stop("lower>upper integration limits")
    if(any(lower == upper)) return(0)
    infin <- c(2,2)
    infin <- replace(infin, (upper == Inf) & (lower > -Inf), 1)
    infin <- replace(infin, (upper < Inf) & (lower == -Inf), 0)
    infin <- replace(infin, (upper == Inf) & (lower == -Inf), -1)
    infin <- as.integer(infin)
    lower <- replace(lower, lower == -Inf, 0)
    upper <- replace(upper, upper == Inf, 0)
    rho   <- as.double(rho)
    prob  <- as.double(0)
    result <- .Fortran("smvbvt", prob, nu, lower, upper, infin, rho)
    return(result[[1]])
}   
   
#multivariate t-distribtuion (for p> 2)
pmt.mtv <- function(dof, lower, upper, mu, Sigma, maxpts=2000*d, abseps=1e-6, releps=0) {
    if(dof == Inf) return(sadmvn(lower, upper, mu, Sigma, maxpts, abseps, releps))
    if(any(lower > upper)) stop("lower>upper integration limits")
    if(any(lower == upper)) return(0)
    if(round(dof) != dof) warning("non integer dof is rounded to integer") 
    dof <- as.integer(round(dof))
    d  <- as.integer(if(is.matrix(Sigma)) ncol(Sigma) else 1)
    Sigma  <- matrix(Sigma, d, d)
    sd  <- sqrt(diag(Sigma))
    rho <- cov2cor(Sigma)
    lower <- as.double((lower-mu)/sd)
    upper <- as.double((upper-mu)/sd)
    if(d == 1) return(pt(upper, dof) - pt(lower, dof))
    infin <- rep(2,d)
    infin <- replace(infin, (upper == Inf) & (lower > -Inf), 1)
    infin <- replace(infin, (upper < Inf) & (lower == -Inf), 0)
    infin <- replace(infin, (upper == Inf) & (lower == -Inf), -1)
    infin <- as.integer(infin)
    lower <- replace(lower, lower == -Inf, 0)
    upper <- replace(upper, upper == Inf, 0)
    correl <- rho[upper.tri(rho, diag=FALSE)]
    maxpts <- as.integer(maxpts)
    abseps <- as.double(abseps)
    releps <- as.double(releps)
    error  <- as.double(0)
    value  <- as.double(0)
    inform <- as.integer(0)
    result <- .Fortran("sadmvt", d, dof, lower, upper, infin, correl, maxpts,
                     abseps, releps, error, value, inform)
    prob <- result[[11]]
    attr(prob,"error")  <- result[[10]]
    attr(prob,"status") <- switch(1+result[[12]], 
                  "normal completion", "accuracy non achieved", "oversize")
    return(prob)
}

#multivariate normal distribtuion
sadmvn <- function(lower, upper, mean, varcov, maxpts=2000*d, abseps=1e-6, releps=0){
      if(any(lower > upper)) stop("lower>upper integration limits")
      if(any(lower == upper)) return(0)
      d <- as.integer(if(is.matrix(varcov)) ncol(varcov) else 1)
      varcov <- matrix(varcov, d, d)
      sd  <- sqrt(diag(varcov))
      rho <- cov2cor(varcov)
      lower <- as.double((lower-mean)/sd)
      upper <- as.double((upper-mean)/sd)
      if(d == 1) return(pnorm(upper)-pnorm(lower))
      infin <- rep(2,d)
      infin <- replace(infin, (upper == Inf) & (lower > -Inf), 1)
      infin <- replace(infin, (upper < Inf) & (lower == -Inf), 0)
      infin <- replace(infin, (upper == Inf) & (lower == -Inf), -1)
      infin <- as.integer(infin)
      lower <- replace(lower, lower == -Inf, 0)
      upper <- replace(upper, upper == Inf, 0)
      correl <- as.double(rho[upper.tri(rho, diag=FALSE)])
      maxpts <- as.integer(maxpts)
      abseps <- as.double(abseps)
      releps <- as.double(releps)
      error  <- as.double(0)
      value  <- as.double(0)
      inform <- as.integer(0)
      result <- .Fortran("sadmvn", d, lower, upper, infin, correl, maxpts,
                   abseps, releps, error, value, inform)
      prob <- result[[10]]
      attr(prob,"error")  <- result[[9]]
      attr(prob,"status") <- switch(1+result[[11]], 
                    "normal completion", "accuracy non achieved", "oversize")
      return(prob)
}

#multivariate normal density
dmn <- function(x, mu=rep(0,d), Sigma, log=FALSE) {
      d  <- if(is.matrix(Sigma)) ncol(Sigma) else 1
      if(d>1 & is.vector(x)) x <- matrix(x, 1, d)
      n  <- if(d==1)  length(x) else nrow(x) 
      X  <- t(matrix(x, nrow=n, ncol=d)) - mu
      Q  <- apply((solve(Sigma)%*% X)* X, 2, sum) 
      logDet <- sum(logb(abs(diag(qr(Sigma)[[1]]))))
      logPDF <- as.vector(Q + d*logb(2*pi)+logDet)/(-2)
      if(log) logPDF else exp(logPDF)
}

#error rate
error.rate <- function(clust1,clust2){
    if((n=length(clust1))!=length(clust2))  stop("error: length not equal")
    if( (g=length(table(clust1)))!=length(table(clust2))) stop("the number of clusters are not equal")
    permute<-function(a){
        n<-length(a)
        if(n==1) f<-a
        else{
            nm<-gamma(n)
            f<-array(0,c(n,n*nm))
            j<-1
        
            for(i in a){
                f[1, (j-1)*nm+1:nm]<-i
                f[-1,(j-1)*nm+1:nm]<-permute(setdiff(a,i))
                j<-j+1
            }
        }    
        f
    }
       
    id<-1:n
    cmb<-permute(1:g)
    nperm<-ncol(cmb)
    rate<-rep(0,nperm)
    
    for(i in 1:nperm){
        tmp<-rep(0,g)
        tc<-rep(0,n)
        for(j in 1:g)
            tc[clust2==j]=cmb[j,i]
    
        for(j in 1:g){  
            tmp1<-0 
            for(k in (1:g)[-j])
                tmp1<-tmp1+length(intersect(id[clust1==j],id[tc==k]))
            tmp[j]<-tmp1
        }
        rate[i]<-sum(tmp)/n
    }
    min(rate)
}

