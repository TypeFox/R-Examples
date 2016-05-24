`spatial.median` <-
function(X,init=NULL,maxiter=500,eps=1e-6,print.it=FALSE,na.action=na.fail)
    {
        X<-na.action(X)
        if(!all(sapply(X, is.numeric))) stop("'X' must be numeric") 
        X<-as.matrix(X)
        
        if (dim(X)[2]==1) return(median(X))

        if (is.null(init)) y <- apply(X,2,median)
        else y <- init
        iter <-0 
        differ <- Inf 
        while(differ>eps)
            {
            if (iter==maxiter) stop("maxiter reached without convergence")
            y.new <- spatial.median.step(y,X)
            differ <- sqrt(sum((y-y.new)^2))
            iter <- iter+1
            y <- y.new
            }

        if (print.it==T) print(paste("convergence was reached after",iter, "iterations"))
   
        return(y)
    }
