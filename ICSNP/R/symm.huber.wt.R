symm.huber.wt <- function(X, wt =  rep(1, nrow(X)), qg = 0.9, init = NULL, eps = 1e-06, maxiter = 100, na.action=na.fail)
    {
    if (length(wt) != nrow(X)) 
            stop("length of 'wt' must equal the number of rows in 'x'")
    x <- data.frame(wt=wt)
    x$X <- as.matrix(X)
    x <- na.action(x)
    if (!all(sapply(x, is.numeric))) 
        stop("'X' and 'wt' must be numeric")
    X <- x$X
    wt<- x$wt

    p <- ncol(X)
    
    c.square <- 2 * qchisq(qg, p)
    sigma.square <- 2 * pchisq(c.square / 2, p + 2) + (c.square / p) * (1 - qg)
    
        
    if (any(wt < 0) || (sum(wt)) == 0) 
            stop("weights must be non-negative and not all zero")
    
    if (p<2) stop("'X' must be at least bivariate")  
    
    data2 <- pair.diff(X)
    pwt <- as.vector(pair.prod(matrix(wt,ncol=1)))
    sum.pwt <- sum(pwt)
    w.data2 <- sqrt(pwt) * data2
    
    if (is.numeric(init)) V.0 <- init else V.0<- crossprod(w.data2) / sum.pwt
    
    iter<-0
    differ<-Inf
    while(differ>eps)
        {
        if (iter==maxiter)
            {
             stop("maxiter reached without convergence")
            }
        iter<-iter+1
        d2<-(mahalanobis(data2,0,V.0))
        w2<-ifelse(d2<=c.square, 1/sigma.square, (c.square/d2)/sigma.square)
        w.data3 <- sqrt(w2) * w.data2
        V.new<- (1/sum.pwt) * crossprod(w.data3)
        differ <- frobenius.norm(V.new-V.0)    
        V.0<-V.new
        }
    colnames(V.new) <- colnames(X)
    rownames(V.new) <- colnames(X)    
    V.new
    }
