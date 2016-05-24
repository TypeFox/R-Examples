pBayes <- function(x, method="m.ind", const=NULL){
# computes pseudo-Bayes estimates of cell counts in x
    n <- sum(x)
    
    if(is.null(dim(x))){
        d <- length(x)
        p <- 1  
    }  
    else{ 
        d <- prod(dim(x))
        p <- length(dim(x))  
    }
    
    if(method=="jeffreys" | method=="Jeffreys"){
        const <- 1/2
        if(const*d > 0.20*n) warning("The total of the added constant is greater \n 
                                      than the 10% of the total observed cells counts")
        xx <- x + const
        K <- d*const
        if(p>1) ee <- array(1/d, dim(x))
        else ee <- rep(1/d, d)
    }
    if(method=="minimax"){
        const <- sqrt(n)/d
        if(const*d > 0.20*n) warning("The total of the added constant is greater \n 
                                      than the 10% of the total observed cells counts")
        xx <- x + const
        K <- d*const
        if(p>1) ee <- array(1/d, dim(x))
        else ee <- rep(1/d, d)
    }
    if(method=="invcat"){
        const <- 1/d
        if(const*d > 0.20*n) warning("The total of the added constant is greater \n 
                                      than the 10% of the total observed cells counts")
        xx <- x + const
        K <- d*const
        if(p>1) ee <- array(1/d, dim(x))
        else ee <- rep(1/d, d)
    }
    if(method=="user"){
        if(const*d > 0.20*n) warning("The total of the added constant is greater \n 
                                      than the 10% of the total observed cells counts")
        xx <- x + const
        K <- d*const
        if(p>1) ee <- array(1/d, dim(x))
        else ee <- rep(1/d, d)
    }
    if(method=="m.ind"){
        if(p==1) stop("m.ind can be used with at least 2-ways tables")
        pp <- prop.table(x)
        ee <- loglin(pp, as.list(1:p), fit=T)$fit
        K <- (1^2 - sum(pp^2))/sum((pp-ee)^2)
        tt1 <- n/(n+K)*pp + K/(n+K)*ee
        xx <- tt1*n
    }
    if(method=="h.assoc"){
        if(p==1) stop("h.assoc can be used with at least 2-ways tables")
        pp <- prop.table(x)
        l1 <- as.list(1:p)
        l2 <- as.list(data.frame(combn(1:p, 2)))
        ll <- c(l1, l2)
        names(ll) <- NULL       
        ee <- loglin(pp, ll, fit=T)$fit
        K <- (1^2 - sum(pp^2))/sum((pp-ee)^2)
        tt1 <- n/(n+K)*pp + K/(n+K)*ee
        xx <- tt1*n
    }
    vv <- c(n=n, no.cells=d, av.cfr=n/d, no.0s=sum(x==0), const=const, K=K, rel.K=K/(n+K))
    list(info=vv, prior=ee*n, pseudoB=xx)
}