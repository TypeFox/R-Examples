iind<-function(X,Y,method,n.simu,n, p.X, p.Y)
    {
    # W.star <- sum(diag(solve(cov(X)) %*%  %*% solve(cov(Y)) %*% cov(Y,X)))
    covXY <- cov(X,Y)
    ch.X <- chol(cov(X))
    ch.Y <- chol(cov(Y))
    p1 <- backsolve(ch.X, forwardsolve(ch.X, covXY))
    p2 <- backsolve(ch.Y, forwardsolve(ch.Y, t(covXY)))
    W.star <- sum(diag(p1 %*% p2))
    
    STATISTIC <- n*W.star
    
    names(STATISTIC) <- "Q.2"
    dfs <- p.X * p.Y
    
    METHOD<-"Test of independence using Pillai's trace"
    
    
    res <- switch(method,
     "approximation" = {PVAL <- 1-pchisq(STATISTIC,df=dfs)
                        PARAMETER <- dfs
                        names(PARAMETER)<-"df"
                        list(statistic=STATISTIC,p.value=as.numeric(PVAL),method=METHOD,parameter=PARAMETER)}
     ,
     "permutation" = {
                    W.star.simu<-function(X,Y,index, covXi = syminv(cov(X)), covYi= syminv(cov(Y))) 
                        {
                        # covXi <- syminv(cov(X))
                        # covYi <- syminv(cov(Y))
                        covXY <- cov(X,Y[index,])
                        sum( diag( covXi %*% covXY %*% covYi %*% t(covXY)))
                        }
                    
                    covXi = syminv(cov(X))
                    covYi = syminv(cov(Y))
                    statistics <- replicate (n.simu, W.star.simu(X,Y,sample(1:n), covXi=covXi, covYi=covYi))
                    PVAL<-mean(statistics>W.star)
                    PARAMETER <- n.simu
                    names(PARAMETER)<-"replications"
                    list(statistic=STATISTIC,p.value=as.numeric(PVAL),method=METHOD,parameter=PARAMETER)
                    })
    res
    }



ssind.inner<-function(X, Y, method, n.simu, n, p.X, p.Y, p)
    {
    U.x <- spatial.sign(X)
    U.y <- spatial.sign(Y)
    U.xy <- crossprod(U.x,U.y)
    Q.stat <- sum((U.xy)^2)
    dfs <- p.X * p.Y
    STATISTIC <- dfs /n * Q.stat 
    names(STATISTIC) <- "Q.2"
    
    
    METHOD<-"Spatial sign test of independence using inner standardization"
    
    res <- switch(method,
     "approximation" = {PVAL <- 1-pchisq(STATISTIC, df = dfs)
                        PARAMETER <- dfs
                        names(PARAMETER)<-"df"
                        list(statistic=STATISTIC,p.value=as.numeric(PVAL),method=METHOD,parameter=PARAMETER)}
    ,
     "permutation" = {
                    Q2simu<-function(U.x,U.y,index) 
                        {
                        U.xy <- crossprod(U.x,U.y[index,])
                        sum((U.xy)^2)
                        }
                    
                    statistics <- replicate (n.simu, Q2simu(U.x,U.y,sample(1:n)))
                    PVAL<-mean(statistics>Q.stat)
                    PARAMETER <- n.simu
                    names(PARAMETER)<-"replications"
                    list(statistic=STATISTIC,p.value=as.numeric(PVAL),method=METHOD,parameter=PARAMETER)
                    })
    res
    }


srind.inner <- function(X, Y, method, n.simu, n, p.X, p.Y, p)
    {
    R.x <- spatial.rank(X)
    R.y <- spatial.rank(Y)
    R.xy <- crossprod(R.x,R.y)
    Q.stat <- sum((R.xy)^2)
    dfs <- p.X * p.Y
    STATISTIC <- dfs * n * Q.stat / (sum(R.x^2) * sum(R.y^2))
    names(STATISTIC) <- "Q.2"

    METHOD<-"Spatial rank test of independence using inner standardization"
    
    res <- switch(method,
     "approximation" = {PVAL <- 1-pchisq(STATISTIC, df = dfs)
                        PARAMETER <- dfs
                        names(PARAMETER)<-"df"
                        list(statistic=STATISTIC,p.value=as.numeric(PVAL),method=METHOD,parameter=PARAMETER)}
    ,
     "permutation" = {
                    Q2simu<-function(R.x,R.y,index) 
                        {
                        R.xy <- crossprod(R.x,R.y[index,])
                        sum((R.xy)^2)
                        }
                    
                    statistics <- replicate (n.simu, Q2simu(R.x,R.y,sample(1:n)))
                    PVAL<-mean(statistics>Q.stat)
                    PARAMETER <- n.simu
                    names(PARAMETER)<-"replications"
                    list(statistic=STATISTIC,p.value=as.numeric(PVAL),method=METHOD,parameter=PARAMETER)
                    })
    res
    }


symmsind.inner <- function(X, Y, method, n.simu, n, p.X, p.Y, p)
    {
    Rcov.x <- rank.shape(X)
    Rcov.y <- rank.shape(Y)
    sq.Rcov.x <- mat.sqrt(solve(Rcov.x))
    sq.Rcov.y <- mat.sqrt(solve(Rcov.y))
    X.inner <- X %*% sq.Rcov.x
    Y.inner <- Y %*% sq.Rcov.y
    R.x <- spatial.rank(X.inner, shape = FALSE)
    R.y <- spatial.rank(Y.inner, shape = FALSE)
    
    US.X <- spatial.sign(pair.diff(X.inner), center = FALSE, shape = FALSE)
    US.Y <- spatial.sign(pair.diff(Y.inner), center = FALSE, shape = FALSE)
    
    dfs <- p.X * p.Y
    
    US.xy <- crossprod(US.X,US.Y)
    
    Q.stat <- sum(US.xy^2) 
    # Note the factor is omitted because here not all pairwise differences but only choose(n,2)
    # differences are used.
    STATISTIC <- (dfs*n) / ((n-1)^2) * Q.stat / (sum(R.x^2)*(sum(R.y^2)))
    names(STATISTIC) <- "Q.2"

    METHOD<-"Symmetrized spatial sign test of independence using inner standardization"
    
    res <- switch(method,
     "approximation" = {PVAL <- 1-pchisq(STATISTIC,df=dfs)
                        PARAMETER <- dfs
                        names(PARAMETER)<-"df"
                        list(statistic=STATISTIC,p.value=as.numeric(PVAL),method=METHOD,parameter=PARAMETER)}
    ,
     "permutation" = {
                    stop("'method = 'permutation'' is not implemented for symmetrized sign scores")
                    })
    res
    }
