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

################################################################################
#  SECTION 3
#                   Random Multivariate Skew-t Variables
#
################################################################################

rmst <- function(n=1, mu = NULL, sigma=NULL, delta=NULL, dof=1, known=NULL) {
    if(n==0) stop("n must be at least 1.")
    if(!is.null(known)) {
        mu <- known$mu
        sigma <- known$sigma
        delta <- known$delta
        dof <- known$dof
    }      
    if(is.null(mu)) stop("mu is missing") else p <- length(as.numeric(mu))
    if(is.null(delta)) delta <- rep(0,p) else if(length(as.numeric(delta))!=p) stop(paste("delta should be a",p,"by 1 matrix"))
    if(is.null(sigma)) sigma<-diag(p) else if(length(c(sigma))!=(p*p)) stop(paste("sigma should be a",p,"by",p,"matrix"))        
    Delta <- diag(p); diag(Delta)<- delta              
    w <- rgamma(n, dof/2, dof/2)                       
    w <- (sqrt(w))^(-1)                                
    U0 <- mvrnorm(n, rep(0,p), sigma)                  
    U1 <- abs(mvrnorm(n, rep(0,p), diag(p)))           
    Z <- Delta %*% t(U1) + t(U0)                       
    Y <- matrix(mu,p,1)%*%rep(1,n) + Z/w               
#    Y <- rbind(Y, rep(1,n))                           
    return(t(Y))                                    
}

rfmmst <- function(g, n, mu, sigma, delta, dof=rep(10,g), pro=rep(1/g,g), known=NULL) {
    if(g==1) return(cbind(rmst(n, mu, sigma, delta, dof, known), rep(1,n)))
    if(!is.null(known)) {
        mu <- known$mu
        sigma <- known$sigma
        delta <- known$delta
        dof <- known$dof
        pro <- known$pro
    }
    if(is.null(mu[[1]])) stop("mu is missing") else p <- length(as.numeric(mu[[1]]))
    for(i in 1:g) {
        if(length(as.numeric(mu[[i]]))!=p) stop("dimension mismatch in component means.") 
        if(length(as.numeric(delta[[i]]))!=p) stop("dimension mismatch in component skewness parameters.")
        if(length(as.numeric(sigma[[i]]))!=(p*p)) stop("dimension mismatch in component scale matrices.")
    }    
    if(length(as.numeric(n))!=1) {CompSize <- cs <- n}
    else {CompSize <- cs <- table(sample(1:g, n, replace = TRUE, prob = pro))}
    if(length(CompSize) < g) {
        CompSize <- rep(0,g)
        for(i in as.numeric(names(cs))) CompSize[i] <- cs[paste(i)]
    }
    names(CompSize) <- NULL
    Y <- array(0, c(10,p))      
    for (i in 1:g) {
        if(CompSize[i]>0) Y <- rbind(Y, rmst(CompSize[i], as.numeric(mu[[i]]), sigma[[i]], as.numeric(delta[[i]]), dof[i]))
    }
    Y <- Y[-(1:10),] 
    Y <- cbind(Y,rep(1:g,CompSize))                            
    return(Y)
}


