Dynamic.Ridge <-
function(time,death,X,R,Ft,lambda=0,fun=c("dynamic","weighted","simple"),eps=1e-06, iter.max=200, theta, mon=FALSE, lambdaFixed=FALSE) {
    call <- match.call()
    if (missing(fun))
        stop("fun is missing with no default!") 
    lambda <- lambda
    ord <- order(time)
    Ftime <- Ft[ord, ]
    n <- nrow(X); d <- death; p <- ncol(X); q <- ncol(Ftime)
    if(missing(R)){v <- 0} else{v <- ncol(R)}
    dl <- rep(0,v); dm <- p*q+v
    if(missing(theta)){  theta <- matrix(rep(0, p * q), ncol = q)}
    X <- apply(X, 2, function(x) { x - mean(x)})
    means <- matrix(apply(matrix(apply(X, 2, function(x) x *Ft), ncol = p * q), 2, mean), ncol = q, byrow = TRUE)
    eventtimes <- (1:n)[d == 1]
    k <- length(eventtimes)
    hi <- rep(1,k)
    JJ <- Ftime
    Ftime <- Ftime[eventtimes, ]
    it <- 0; itl <- 0; dlt <- 5
    while(dlt>eps & itl <iter.max){     
        delta <- 5;  iter <- 0; 
        while (delta > eps) {
            th <- theta
            it <- it + 1
            if(v>0){            
                tmp <- sapply(eventtimes, sumeventsf, X, R, Ftime, theta,dl,n,p,v,eventtimes)}
            else{tmp <- sapply(eventtimes, sumevents, X, Ftime, theta, n, p, eventtimes)}
            hi <- unlist(tmp[1, ])
            
            if(fun=="weighted"){
                W <- diag(hi)} else {W <- diag(k)}
            
            if(fun=="simple"){
                scores <- matrix(unlist(tmp[2, ]), ncol = k) %*% rep.int(1,k) + matrix(t(X[eventtimes, ]) %*% Ftime, ncol = 1) - lambda * as.vector(theta)}
            else {scores <- matrix(unlist(tmp[2, ]), ncol = k) %*% rep.int(1,k) + matrix(t(X[eventtimes, ]) %*% Ftime, ncol = 1) - lambda * as.vector(theta %*%t(Ftime)%*%W %*% Ftime)}
            
            if(v>0){scores <- rbind(scores, matrix(unlist(tmp[3,]),ncol=k)%*%rep.int(1,k)+t(R[eventtimes,]) %*%rep.int(1,k))}
            
            if(fun=="simple"){	          	     
                pen <- lambda*diag(dm)
                if(v>0){pen[(p*q+1):dm,(p*q+1):dm] <- 0}}
            else{pen<- lambda * kronecker( t(Ftime)%*%W %*% Ftime, diag(rep(1,p)))
                 
                 if(v>0){
                     pen <-  0*diag(dm)
                     p2 <-  lambda * kronecker( t(Ftime)%*%W %*% Ftime, diag(p))
                     dmp <- dim(p2)[2]
                     pen[1:dmp,1:dmp] <- p2}}
            
            if(v>0){  h2 <--(matrix(matrix(unlist(tmp[4,]),ncol=k) %*% rep(1,k), ncol=dm, byrow=FALSE)) }
            else{h2 <--(matrix(matrix(unlist(tmp[3,]),ncol=k) %*% rep(1,k), ncol=dm, byrow=FALSE)) }
            
            hess <- h2-pen
            shess <- solve(hess)
            coef <- c(matrix(theta,ncol=1),dl)-(shess%*%matrix(t(scores),ncol=1))
            theta <- matrix(coef[1:(p*q)],ncol=q)
            if(v>0){dl <- coef[(p*q+1):dm]}
            delta <- max(abs(theta - th))
            if(v>0){ Hat <- (matrix(matrix(unlist(tmp[4,]),ncol=k) %*% rep(1,k), ncol=dm, byrow=FALSE))%*% shess}
            else{Hat <- (matrix(matrix(unlist(tmp[3,]),ncol=k) %*% rep(1,k), ncol=dm, byrow=FALSE))%*% shess}
        }
        
        if(fun=="simple")
            lambdaA <- as.numeric(-sum(diag(Hat[1:(p*q),1:(p*q)]))/sum( theta%*% (t(theta))))
        else    lambdaA <- as.numeric(-sum(diag(Hat[1:(p*q),1:(p*q)]))/sum( theta%*%t(Ftime)%*%W%*%Ftime%*% (t(theta))))
        
        dlt <- sum(abs(lambda-lambdaA))
        if(lambdaFixed==TRUE){dlt= eps}
        if(lambdaFixed==FALSE){lambda <- lambdaA}
        itl <- itl+1
        
        if(mon==TRUE){cat("iteration",itl,". Penalty weight=",lambda,"\n")}
        
    }
    if(fun=="simple"){
        lik <- sum(log(hi)) + sum(X[eventtimes, ] * t(theta %*%t(Ftime)))-.5*lambda*sum(diag( theta%*%t(theta) ))}
    else{lik <- sum(log(hi)) + sum(X[eventtimes, ] * t(theta %*%t(Ftime)))-.5*lambda*sum(diag(theta %*%t(Ftime)%*%W %*% Ftime%*%t(theta)) ) }
    
    if(v>0){
        fit <- list(call = call, theta=theta, fixed.coef=dl, hess=hess, loglik=lik, time = time,death=death, X=X, R=R,Ft=Ft ,iter =it,inter.it =itl, lambda=lambda,h2=h2,hat=Hat)}
    else {fit <- list(call = call, theta=theta, hess=hess, loglik=lik, time = time,death=death, X=X, Ft=Ft ,iter =it,inter.it =itl, lambda=lambda,h2=h2,hat=Hat)}
    
    class(fit) <- "cox.dynamic.ridge"
    fit}
