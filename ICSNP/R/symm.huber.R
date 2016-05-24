symm.huber <- function(X, qg=0.9, init=NULL, eps = 1e-06, maxiter = 100, na.action=na.fail)
    {
    X <- na.action(X)
    X <- as.matrix(X)
    CNAMES <- colnames(X)
    p <- ncol(X)
    
    c.square <- 2 * qchisq(qg, p)
    sigma.square <- 2 * pchisq(c.square / 2, p + 2) + (c.square / p) * (1 - qg)
    
    if (p<2) stop("'X' must be at least bivariate")  

    n<-dim(X)[1]

    if (is.numeric(init)) V.0 <- init else V.0<- cov(X)
 
    iter<-0
    differ<-Inf
    
    while(differ>eps)
        {
        if (iter==maxiter)
            {
             stop("maxiter reached without convergence")
            }
        iter<-iter+1
        
        V.new<- SSCov.hub(X,solve(V.0),c.square,sigma.square)
        differ <- frobenius.norm(V.new-V.0)    
        V.0<-V.new
        }
    colnames(V.new) <- CNAMES
    rownames(V.new) <- CNAMES
    V.new
    }

# Internal function for symm.huber

SSCov.hub<-function(X,V,c.s,sig.s)
{
X<-as.matrix(X)
d<-dim(X)
matrix(.C("symm_huber", as.double(X), as.double(V), as.integer(d), as.double(c.s), as.double(sig.s),res=double(d[2]^2))$res, ncol=d[2],byrow=T)/(d[1]*(d[1]-1)/2)
}
