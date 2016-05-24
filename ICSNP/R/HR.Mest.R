`HR.Mest` <-
function(X,maxiter=100,eps.scale=1e-6,eps.center=1e-6,na.action=na.fail)
    {
    X<-na.action(X)
    if(!all(sapply(na.action(X), is.numeric))) stop("'X' must be numeric") 
    X<-as.matrix(X)

    n <- dim(X)[1]
    p <- dim(X)[2]
               
    if (p<2) stop("'X' must be at least bivariate")


    # starting value is marginal median vector
    theta0 <- apply(X,2,median) 
    V0 <- tyler.shape(X,location=theta0,eps=eps.scale,maxiter=maxiter,print.it=FALSE)
    A0 <- solve(mat.sqrt(V0))                                                 
    differ <- Inf 
    iter <- 0
     while(differ>eps.scale)
        { 
        if (iter==maxiter)
            {
             stop("maxiter reached without convergence")
            } 
        iter <- iter + 1
        theta.k1 <- center.step(X,A0,maxiter,eps.center) 
        V.k1 <- tyler.shape(X,location=theta.k1,eps=eps.scale,maxiter=maxiter,print.it=FALSE)
        A.k1 <- solve(mat.sqrt(V.k1))                                                
        theta.k2 <- center.step(X,A.k1,maxiter,eps.center)   
        V.k2 <- tyler.shape(X,location=theta.k2,eps=eps.scale,maxiter=maxiter,print.it=FALSE)                                           
        A.k2 <- solve(mat.sqrt(V.k2))
        A0 <- A.k2 
        differ <- sqrt(sum((theta.k1-theta.k2)^2))
        }
    scatter <- solve(A.k2 %*% t(A.k2))
    theta.k2 <- as.vector(theta.k2)
    names(theta.k2) <- colnames(X)
    res <- list(center = theta.k2, scatter = V.k2)
    return(res)
    }
