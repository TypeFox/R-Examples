duembgen.shape.wt <- function(X, wt =  rep(1, nrow(X)), init = NULL, eps = 1e-6, maxiter = 100, na.action = na.fail)
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
    
        
        if (any(wt < 0) || (sum(wt)) == 0) 
            stop("weights must be non-negative and not all zero")
    
    if (p<2) stop("'X' must be at least bivariate")  
    
    data2 <- pair.diff(X)
    center.ind<-apply(data2,1,setequal,y=rep(0,p))
    n.del<-sum(center.ind)
    
    if (n.del!= 0)
        {
        data2<-data2[center.ind==F,]
        pwt <- as.vector(pair.prod(matrix(wt,ncol=1)))[center.ind==F]
        if (n.del>1)
            {warning(paste(n.del ,"one pairwise difference equal to the origin was removed"))
        }else
            {warning("One pairwise difference equal to the origin was removed")}
        } else
        {pwt <- as.vector(pair.prod(matrix(wt,ncol=1)))[center.ind==F]}
    
    
    
    sum.pwt <- sum(pwt)
    w.data2 <- pwt * data2
    
    iter <-0
    if (is.numeric(init)) V.0<-solve(init)
    else V.0<-solve(t(w.data2)%*%w.data2)
    differ<-Inf
    while(TRUE)
        {
        if (differ<eps) break
        if (iter>=maxiter)
            {
             stop("maxiter reached without convergence")
            }
        # print(iter)
            V.new<-.wt.duembgen.step(V.0,data2,pwt,p,sum.pwt)
            differ<-frobenius.norm(V.new-V.0)
            V.0<-V.new
            iter=iter+1 
            }
    
    
    V.shape<-solve(V.new)
    V<-V.shape/det(V.shape)^(1/p)
    
    colnames(V) <- colnames(X)
    rownames(V) <- colnames(X)
    return(V)
    
    }

#
# internal function of the weighted duembgen shape matrix
#

.wt.duembgen.step<-function(V.old,datas,pwt,p,sum.pwt)
        {
        sqrt.V.old <- mat.sqrt(V.old)
        r <- sqrt(rowSums((datas %*% sqrt.V.old)^2))
        datas2 <- sqrt(pwt)* (1/r) * datas
        datas3 <- datas2%*%sqrt.V.old
        M.V.old<-p/sum.pwt*crossprod(datas3)
        M.V.old.inv <- solve(M.V.old)
        V.new<-sum(diag(V.old %*% M.V.old.inv))^(-1)*(sqrt.V.old %*% M.V.old.inv %*% sqrt.V.old)
        return(V.new)
        }
