mvhuberM <- function(X, qg=0.9, fixed.loc=FALSE, location=NULL, init=NULL, steps=Inf, eps=1e-06,  maxiter=100, na.action=na.fail)
{ 
    X <- na.action(X)
    X <- as.matrix(X)
    n <- dim(X)[1]
    p <- dim(X)[2]
    c1 <- sqrt(qchisq(qg, p))
    c2 <- pchisq(c1^2, p + 2) + ((c1^2 * (1-qg))/p)
    
    if(is.null(init)){
       V <- cov(X)
    }else {
     if(is.function(init)) init <- init(X)
     V <- init
     }

    if(is.null(location)){
       mu <- colMeans(X)
    }else {
     if(is.function(location)) location <- location(X)
     mu <- location
     }

    if(!fixed.loc){
 
     iter <- 0
     while(TRUE) {
        if(iter==steps) break
        if(iter>=maxiter) stop("maxiter reached")
        iter <- iter + 1
        r <- sqrt(mahalanobis(X, center = mu, cov = V))
        w1 <- (r <= c1) + (c1 * (r > c1))/r
        w2 <- w1^2/c2
        mu.new <- apply(w1 * X, 2, sum)/sum(w1)   
        y <- sqrt(w2) * scale(X, center = mu.new, scale = F)
        V.new <- ((n - 1) * var(y))/n
        
        if(all(is.infinite(steps),mat.norm(V.new-V)<eps,sqrt(sum((mu.new -mu)^2))<eps)){
        V <- V.new
        mu <- mu.new
        break
        }
        V <- V.new
        mu <- mu.new
     }
    }else{
     iter <- 0
     while(TRUE) {
        if(iter==steps) break
        if(iter>=maxiter) stop("maxiter reached")
        iter <- iter + 1
        r <- sqrt(mahalanobis(X, center = mu, cov = V))
        w1 <- (r <= c1) + (c1 * (r > c1))/r
        w2 <- w1^2/c2
        y <- sqrt(w2) * scale(X, center = mu, scale = F)
        V.new <- ((n - 1) * var(y))/n
        
        if(all(is.infinite(steps),mat.norm(V.new-V)<eps)) break
        V <- V.new
     }

    }
    list(location = mu, scatter = V)
}





