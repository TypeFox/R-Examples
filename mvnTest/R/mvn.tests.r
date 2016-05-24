#########################################
## the Anderson-Darling test                
#########################################

AD.test<-function(data, qqplot=FALSE)
 {  ## This function calculates the value of the Anderson-Darling test
    ## and approximate p-value
    ## data - as n times p matrix or data frame, where p is dimension
    if (!is.data.frame(data) && !is.matrix(data)) stop('data supplied must be either of class \"data frame\" or \"matrix\"')
    if (dim(data)[2] < 2 || is.null(dim(data))) {stop("data dimesion has to be more than 1")}
    if (dim(data)[1] < 3) {stop("not enough data for assessing mvn")}
    data.name <- deparse(substitute(data))
    xp <- as.matrix(data)
    p <- dim(xp)[2]      
    n <- dim(xp)[1]
    ## getting MLEs...
    s.mean <- colMeans(xp)
    s.cov <- (n-1)/n*cov(xp)
    s.cov.inv <- solve(s.cov) # inverse matrix of S (matrix of sample covariances)
    D <- rep(NA,n) # vector of (Xi-mu)'S^-1(Xi-mu)...
    for (j in 1:n) 
       D[j] <- t(xp[j,]-s.mean)%*%(s.cov.inv%*%(xp[j,]-s.mean))
    D.or <- sort(D) ## get ordered statistics
    Gp <- pchisq(D.or,df=p)
    ## getting the value of A-D test... 
    ind <- c(1:n)
    an <- (2*ind-1)*(log(Gp[ind])+log(1 - Gp[n+1-ind]))
    AD <- -n - sum(an) / n

    if (qqplot) {   ## not checked yet....
            xp <- scale(xp, scale = FALSE)
            D <- xp %*%s.cov.inv%*% t(xp)
            d <- diag(D)    
            r <- rank(d)  
            chi2q <- qchisq((r-0.5)/n,p)
            plot(d, chi2q, pch = 19, main = "Chi-Square Q-Q Plot", xlab = "Squared Mahalanobis Distance", ylab = "Chi-Square Quantile")
            abline(0, 1,lwd = 2, col = "black")
        }

    ## getting the p-value...
    N <- 1e4
    U <- rep(0,N) ## initializing values of the AD test
    for (i in 1:N) { ## loop through N reps
        dat <- rmvnorm(n=n,mean=s.mean,sigma=s.cov) # generating n p-variate normal rv   
        mean1 <- colMeans(dat)
        cov1 <- (n-1)/n*cov(dat)
        cov.inv <- solve(cov1) # inverse matrix of S (matrix of sample covariances)
        D <- rep(NA,n) # vector of (Xi-mu)'S^-1(Xi-mu)...
        for (j in 1:n) 
           D[j] <- t(dat[j,]-mean1)%*%(cov.inv%*%(dat[j,]-mean1))
        Gp <- pchisq(sort(D),df=p)
        ## getting the value of A-D test... 
        an <- (2*ind-1)*(log(Gp[ind])+log(1 - Gp[n+1-ind]))
        U[i] <- -n - sum(an) / n
     } 
    
    p.value <- (sum(U >= AD)+1)/(N+1)
    result <- new("ad",AD=AD, p.value=p.value, data.name=data.name)
    result
} ## AD.test



#################################################
##      Cramer-von Mises test                  ##
#################################################

CM.test<-function(data, qqplot=FALSE)
 {  ##  This function calculates the value of the Cramer-von Mises test and
    ## approximate  p-value
    ## data - as n times p matrix or data frame, where p is dimension
    if (!is.data.frame(data) && !is.matrix(data)) stop('data supplied must be either of class \"data frame\" or \"matrix\"')
    if (dim(data)[2] < 2 || is.null(dim(data))) {stop("data dimesion has to be more than 1")}
    if (dim(data)[1] < 3) {stop("not enough data for assessing mvn")}
    data.name <- deparse(substitute(data))
    xp <- as.matrix(data)
    p <- dim(xp)[2]      
    n <- dim(xp)[1]
    ## getting MLEs...
    s.mean <- colMeans(xp)
    s.cov <- (n-1)/n*cov(xp)
    s.cov.inv <- solve(s.cov) # inverse matrix of S (matrix of sample covariances)
    D <- rep(NA,n) # vector of (Xi-mu)'S^-1(Xi-mu)...
    for (j in 1:n) 
       D[j] <- t(xp[j,]-s.mean)%*%(s.cov.inv%*%(xp[j,]-s.mean))
    D.or <- sort(D) ## get ordered statistics
    Gp <- pchisq(D.or,df=p)
    ## getting the value of CM test...
    ind <- c(1:n) 
    w <- (Gp - (2 * ind - 1) / (2 * n)) ^ 2
    CM <- sum(w) + 1 / (12 * n)

    if (qqplot) {   ## not checked yet....
            xp <- scale(xp, scale = FALSE)
            D <- xp %*%s.cov.inv%*% t(xp)
            d <- diag(D)    
            r <- rank(d)  
            chi2q <- qchisq((r-0.5)/n,p)
            plot(d, chi2q, pch = 19, main = "Chi-Square Q-Q Plot", xlab = "Squared Mahalanobis Distance", ylab = "Chi-Square Quantile")
            abline(0, 1,lwd = 2, col = "black")
        }

    ## getting the p-value...
       N <- 1e4
       U <- rep(0,N) ## initializing values of the AD test
       for (i in 1:N) { ## loop through N reps
           dat <- rmvnorm(n=n,mean=s.mean,sigma=s.cov) # generating n p-variate normal rv   
           mean1 <- colMeans(dat)
           cov1 <- (n-1)/n*cov(dat)
           cov.inv <- solve(cov1)
           D <- rep(NA,n) # vector of (Xi-mu)'S^-1(Xi-mu)...
           for (j in 1:n)
               D[j] <- t(dat[j,]-mean1)%*%(cov.inv%*%(dat[j,]-mean1))
           Gp <- pchisq(sort(D),df=p)
           w <- (Gp - (2 * ind - 1) / (2 * n)) ^ 2
           U[i] <-  sum(w) + 1 / (12 * n)        
       }       
       p.value <- (sum(U >= CM)+1)/(N+1)
    
    result <- new("cm",CM=CM, p.value=p.value, data.name=data.name)
    result
 } ## CM.test



####################################################
## function to calculate the Doornik-Hansen (1994) test (Doornik and Hansen 1994)  
####################################################

DH.test<-function(data,qqplot=FALSE){   
    ##  This function calculates the Doornik-Hansen (1994) test statistic
    ## returns the test statistic and the approximate p-value
    ## data - as n times p matrix or data frame, where p is dimension
    if (!is.data.frame(data) && !is.matrix(data)) stop('data supplied must be either of class \"data frame\" or \"matrix\"')
    if (dim(data)[2] < 2 || is.null(dim(data))) {stop("data dimesion has to be more than 1")}
    if (dim(data)[1] < 3) {stop("not enough data for assessing mvn")}
    data.name <- deparse(substitute(data))
    data <- as.matrix(data)

    transform<-function(data)
    {
    x<-scale(data, scale=F)
    x<-t(x)
    s<-var(data)
    v<-diag(s)^(-1/2)
    v <- diag(v)
    c.data<-v %*% s %*% v
    values<- (tmp <- eigen(c.data, symmetric=T))$values
    vectors<-tmp$vectors
    lambda<-diag(values^(-.5))
    trans<- vectors %*% lambda %*% t(vectors) %*% v %*% x
    return(t(trans))
    }


find.z1<-function(n,trans)

    {  
      get.skew<-function(data)
    {  
            get.m2<-function(data)
        {
         d<-dim(data)
         n<-d[1]
         p<-d[2]
         mean.data<-apply(data,2,mean)
         m2<-rep(0,p)
        for(i in 1:p)
              m2[i]<-(sum((data[,i]-mean.data[i])^2))/n
         return(m2)
         }

        get.m3<-function(data)
         {
            d<-dim(data)
            n<-d[1]
            p<-d[2]
            mean.data<-apply(data,2,mean)
            m3<-rep(0,p)
            for(i in 1:p)
               m3[i]<- sum((data[,i]-mean.data[i])^3)/n
            return(m3)
        }

        d<- dim(data)
        n<-d[1]
        p<-d[2]
        m2<-get.m2(data)
        m3<-get.m3(data)
        b1<-rep(0,p)
        for(i in 1:p)
             b1[i]<-m3[i]/m2[i]^1.5
        return(b1)
    }
        skew <-get.skew(trans)
        B<- 3*(n^2+27*n-70)*(n+1)*(n+3)/((n-2)*(n+5)*(n+7)*(n+9))
        w.2<- -1+ (2*(B-1))^.5 
        delta <- 1/log(sqrt(w.2))^.5
        y <- skew * ((w.2-1)*(n+1)*(n+3)/(12*(n-2)))^.5
        z1<-delta*log(y+(y^2+1)^.5)
        return(z1)
}

find.z2<-function(n,trans)
    {    
    get.kurt<-function(data)
    {
           get.m2<-function(data)
        {
         d<-dim(data)
         n<-d[1]
         p<-d[2]
         mean.data<-apply(data,2,mean)
         m2<-rep(0,p)
        for(i in 1:p)
              m2[i]<-sum((data[,i]-mean.data[i])^2)/n            
        return(m2)
         }         
   get.m4<-function(data)
    {
        d<-dim(data)
        n<-d[1]
        p<-d[2]
        mean.data<-apply(data,2,mean)
        m4<-rep(0,p)
        for(i in 1:p)
              m4[i]<-sum((data[,i]-mean.data[i])^4)/n
        return(m4)
    }
        d<- dim(data)
        n<-d[1]
        p<-d[2]
        m2<-get.m2(data)
        m4<-get.m4(data)
        b2<-rep(0,p)
        for(i in 1:p)
                 b2[i]<-m4[i]/m2[i]^2
        return(b2)
}
    get.skew<-function(data)
    {  
      get.m2<-function(data)
        {
         d<-dim(data)
         n<-d[1]
         p<-d[2]
         mean.data<-apply(data,2,mean)
         m2<-rep(0,p)
        for(i in 1:p)        
            m2[i]<-sum((data[,i]-mean.data[i])^2)/n
         return(m2)
         }

        get.m3<-function(data)
         {
            d<-dim(data)
            n<-d[1]
            p<-d[2]
            mean.data<-apply(data,2,mean)
            m3<-rep(0,p)
            for(i in 1:p)
                      m3[i]<-sum((data[,i]-mean.data[i])^3)/n
            return(m3)
        }

        d<- dim(data)
        n<-d[1]
        p<-d[2]
        m2<-get.m2(data)
        m3<-get.m3(data)

        b1<-rep(0,p)
        for(i in 1:p)
               b1[i]<-m3[i]/m2[i]^1.5
        return(b1)
    }

        skew<-get.skew(trans)
        kurt<-get.kurt(trans)
        delta<-(n-3)*(n+1)*((n^2)+(15*n)-4)
        a<- (n-2)*(n+5)*(n+7)*(n^2+27*n-70)/(6*delta)
        c<- (n-7)*(n+5)*(n+7)*(n^2+2*n-5)/(6*delta)
        k<-(n+5)*(n+7)*(n^3+37*n^2+11*n-313)/(12*delta)
        alpha<-a + skew^2*c
        chi<-(kurt-1-skew^2)*2*k
        z2<-((chi/(2*alpha))^(1/3)-1+ 1/(9*alpha))*3*alpha^0.5
	return(z2)
    }


    trans<-transform(data)
    d<-dim(trans)
    n<-d[1]
    p<-d[2]
    z1<-find.z1(n,trans)
    z2<-find.z2(n,trans)
    DH<-sum(t(z1)*z1)+sum(t(z2)*z2)
    p.value<- 1-pchisq(DH,df=2*p)

    if (qqplot) {   ## NOT checked yet....
            xp <- scale(data, scale = FALSE)
            Sa <- cov(xp)
            D <- xp %*%solve(Sa)%*% t(xp)
            d <- diag(D)    
            r <- rank(d)  
            chi2q <- qchisq((r-0.5)/n,p)
            plot(d, chi2q, pch = 19, main = "Chi-Square Q-Q Plot", xlab = "Squared Mahalanobis Distance", ylab = "Chi-Square Quantile")
            abline(0, 1,lwd = 2, col = "black")
        }
    result <- new("dh",DH = DH, p.value=p.value, data.name=data.name)
    result
} ## DH.test



####################################################
## function to calculate the Royston test (Royston 1992)
####################################################

R.test<-function(data, qqplot=FALSE)
{   ## This function finds the equivalent degrees of freedom and calculates 
    ## the Royston (1992) test and approximate p-value
    ## data - as n times p matrix or data frame, where p is dimension
    if (!is.data.frame(data) && !is.matrix(data)) stop('data supplied must be either of class \"data frame\" or \"matrix\"')
    if (dim(data)[2] < 2 || is.null(dim(data))) {stop("data dimesion has to be more than 1")}
    if (dim(data)[1] < 3) {stop("not enough data for assessing mvn")}
    data.name <- deparse(substitute(data))
    data <- as.matrix(data)

    find.coeff.new<-function(n)
    {
    ## This function finds the coefficients used to obtain the values after applying Royston's normalizing transformation.
        ## n is the number of observations
        ## p is the number of variables
        if(n <= 3) {
            stop(" n is too small!")
        }
        else if(n >= 4 && n <=11) {
            x = n
        }
        else if (n >= 12 && n <=2000){
            x = log(n)
        }
        else{
            stop("n is too large!")
        }


        if(n >= 4 && n <=11) {
            gamma = -2.273 + 0.459 * x
        }
        else if(n >= 12 && n <= 2000) {
            gamma = 0

        }
        else {
            stop("n is not in the proper size range")
        }
        ret <- list(gamma = gamma)
        if(n >= 4 && n <= 11) {
            mu = 0.5440 - 0.39978 * x + 0.025054 * x^2 - 0.0006714 * x^3
        }
        else if(n >= 12 && n <= 2000) {
            mu = -1.5861 - 0.31082 * x - 0.083751 * x^2 + 0.0038915 * x^3
        }
        else {
            stop(" n is not in the proper size range")
        }
        ret$mu <- mu
        if(n >= 4 && n <= 11) {
            sigma = exp(1.3822 - 0.77857 *x + 0.062767 * x^2 - 0.0020322 *x^3)
        }
        else if(n >= 12 && n <= 2000) {
            sigma = exp(-0.4803 -0.082676 *x + 0.0030302 * x^2)
        }
        else {
            stop("n is not in the proper size range")
        }
        ret$sigma <- sigma
        return(ret)
    }

    our.test.shapiro<-function(data, n, p)
    {
    ## This function finds the values of Zj after the Royston transformation is applied
        ## data is the data file being analyzed, it is in matrix form
        ## n is the number of observations
        ## is the number of variables

        W <- rep(0, p)

        for(i in 1:p)   W[i] <- shapiro.test(data[, i])$statistic
        return(W)
    }

    r.trans.new<-function( n, p, a, b)
    {
    ## This function normalized the Shapiro-Wilks statistic under Royston's Transformation
        ## data is the data set used, in matrix form
        ## n is the number of observations
        ## p is the number of variables

        Z <- rep(0, p)
        if (n >= 4 && n <= 11)
        {
           for(i in 1:p) Z[i] <- (-log(a$gamma-log(1-b[i]))-a$mu)/a$sigma
        }
        else if (n >= 12 && n <= 2000)
        {
           for (i in 1:p) Z[i] <- (log(1-b[i]) + a$gamma - a$mu)/ a$sigma
        }
        else {
            stop ("WARNING : n is not in the proper range")
        }
        return(Z)
    }

    Royston<-function( n, p, z)
    {
    ## This function finds the statistics to be compared for multivariate normality
        ## data is the data set to be analysed
        ## n is the number of observations
        ## p is the number of variables

        G <- rep(0, p)
        for(i in 1:p) G[i] <- qnorm((pnorm( - z[i]))/2)^2
        return(G)
    }

    d <- dim(data)
    n <- d[1]
    p <- d[2]
    a <- find.coeff.new(n)
    b <- our.test.shapiro(data, n, p)
    z <- r.trans.new( n, p, a, b)
    r <- Royston( n, p, z)
    u <- 0.715
    v <- 0.21364+0.015124*log(n)^2-0.0018034*log(n)^3
    la <- 5
    corr.data <- cor(data)
    new.corr <- corr.data^la * (1 - u  * (1 - corr.data)^u/v)
    total <- sum(new.corr) - p
    avg.corr <- total/(p^2 - p)
    est.e <- p/(1 + (p - 1) * avg.corr)
    ## equivalent degrees of freedom
    res <- list()
    R <- est.e *sum(r)/p
    p.value <- 1 - (pchisq(R, df = est.e))

    if (qqplot) {   ## not checked yet....
            xp <- scale(data, scale = FALSE)
            Sa <- cov(xp)
            D <- xp %*%solve(Sa)%*% t(xp)
            d <- diag(D)    
            r <- rank(d)  
            chi2q <- qchisq((r-0.5)/n,p)
            plot(d, chi2q, pch = 19, main = "Chi-Square Q-Q Plot", xlab = "Squared Mahalanobis Distance", ylab = "Chi-Square Quantile")
            abline(0, 1,lwd = 2, col = "black")
        }
    result <- new("r",R=R,p.value=p.value, data.name=data.name)
    result
} ## R.test


####################################################
## function to calculate the Henze-Zirkler test (Henze and Zirkler 1990)
####################################################

HZ.test<-function(data, qqplot=FALSE) 
{   ## This function calculates the Henze-Zirkler test 
    ## statistic for assessing  multivariate normality of a data set to be analyzed
    ## returns the value of the test and approximate p-value
    ## r is the desired alpha level of significance
    ## data - as n times p matrix or data frame, where p is dimension
    if (!is.data.frame(data) && !is.matrix(data)) stop('data supplied must be either of class \"data frame\" or \"matrix\"')
    if (dim(data)[2] < 2 || is.null(dim(data))) {stop("data dimesion has to be more than 1")}
    if (dim(data)[1] < 3) {stop("not enough data for assessing mvn")}
    data.name <- deparse(substitute(data))
    data <- as.matrix(data)

    loop4<- function(data)
    {
    n <- dim(data)[1]
    s.inv <- solve(var(data))
    X <- X.2 <- data %*% s.inv %*% t(data)
    diag(X) <- 0
    for(i in 1:(n-1)) {
        for(j in (i+1):n) 
            X[j,i] <- X[i, j] <- X.2[i,i] - 2 * X.2[i,j] + X.2[j,j]
    }
    return(X)
    }
    
    Square<-function(data)
    {
        S <- var(data)
        S.inv <- solve(S)
        x.diff <- scale(data, scale = F)
        square <- diag(x.diff %*% S.inv %*% t(x.diff))
        return(square)
    }


    d <- dim(data)
    n <- d[1]
    p <- d[2]

    Beta <- 1/sqrt(2) * ((2 * p + 1)/4)^(1/(p + 4)) * (n^(1/(p + 4)))

    a <- var(data)
    m <- dim(a)[1]
    ranka <- qr(a)$rank

    if( ranka == m){

    s <- loop4(data)
    square <- Square(data)

    T.beta <- n * (1/n^2 * sum(exp( -Beta^2/2 * s)) - 2 * ((1 + (
        Beta^2))^( - p/2)) * (1/n) * (sum(exp( - (Beta^2/(2 * (1 +
        Beta^2))) * square))) + ((1 + (2 * Beta^2))^( - p/2)))
    }
    else
    { 
        T.beta <- n * 4
    }
    ans <- list(T = T.beta)
    mu <- 1 - ((1 + 2 * Beta^2)^( - p/2)) * (1 + (p * Beta^2)/(1 + 2 *
        Beta^2) + (p * (p + 2) *Beta^4)/(2 * (1 + 2 * Beta^2)^2))
   
    w.beta <- (1 + Beta^2) * (1 + 3 * Beta^2)
    sigma.sq <- 2 * (1 + 4 * Beta^2)^( - p/2) + 2 * (1 + 2 * Beta^2)^(- p) *
        (1 + (2 * p * Beta^4)/(1 + (2 * Beta^2))^2 + (
        3 * p * (p + 2) * Beta^8)/(4 * (1 + 2 * Beta^2)^4)) - 4 *
        (w.beta^( - p/2)) * (1 + (3 * p * Beta^4)/(2 * w.beta) + (
        p * (p + 2) * Beta^8)/(2 * w.beta^2))
   

    p.mu<-log(sqrt(mu^4/(sigma.sq + mu^2)))

    p.sigma<-sqrt(log((sigma.sq + mu^2)/mu^2))

    p.value<- 1 - (plnorm(T.beta,p.mu,p.sigma))

    if (qqplot) {   ## not checked yet....
            xp <- scale(data, scale = FALSE)
            Sa <- cov(xp)
            D <- xp %*%solve(Sa)%*% t(xp)
            d <- diag(D)    
            r <- rank(d)
            chi2q <- qchisq((r-0.5)/n,p)
            plot(d, chi2q, pch = 19, main = "Chi-Square Q-Q Plot", xlab = "Squared Mahalanobis Distance", ylab = "Chi-Square Quantile")
            abline(0, 1,lwd = 2, col = "black")
        }
    result <- new("hz", HZ=T.beta, p.value=p.value, data.name=data.name)
    result    
} ## HZ.test




