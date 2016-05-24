ae.hl.estimate<- function(X, init=NULL, shape=TRUE, steps=Inf, maxiter=500,  eps=1e-06, na.action=na.fail)
{ 
    X <- na.action(X)
    X <- as.matrix(X)
    n <- dim(X)[1]
    p <- dim(X)[2]

    if(p==1) return(median(c(pairsum(X)/2,X)))
    if(is.finite(steps)) maxiter<-Inf
    
    if (is.matrix(shape)) {
        X <- X %*% solve(mat.sqrt(shape))
        res <- as.vector(mat.sqrt(shape)%*%hl.location(X, init=init, steps=steps, eps=eps))
        attr(res, "shape") <- shape
        return(res)
    }else if (is.logical(shape)) {
        if (!shape) {
            res <- hl.location(X, init=init, steps=steps, eps=eps)
            attr(res, "shape") <- diag(p)
            return(res)
        }
    }
    if(shape)
    {
    if(is.function(init)) init<-init(X)
    if(is.null(init)) {
        init <- spat.median(X)
    }else if (any(!is.vector(init), !is.numeric(init))) { 
        stop("'init' must be a numeric vector or NULL")
    }else if (length(init) != p) 
        stop("'init' is of wrong dimension")

    mu <- init
    V <- signrank.shape(X,fixed.loc=TRUE,location=mu)
    differ <- Inf
    iter <- 0
    while(TRUE)
    {
     if(iter>=steps) return(res)
     if(iter>=maxiter) stop("maxiter reached")
      iter<-iter+1
      sqrtV<-mat.sqrt(V)
      A <- solve(sqrtV)
      Z<- sweep(X,2,mu)%*%t(A)  
      mu.new <- mu+sqrtV%*%hl.center.step(Z,rep(0,p)) 
      V.new<-sqrtV%*%SRCov(Z,rep(0,p))%*%sqrtV 
     V.new<-to.shape(V.new)
     res<-as.vector(mu.new)
     attr(res,"shape")<-V.new
     if(all(is.infinite(steps),mat.norm(V.new-V)<eps,sqrt(sum((mu.new -mu)^2))<eps))  return(res)
     mu<-mu.new
     V<-V.new
    }
    }
}

