robust <- function(x, w, eta, etanew, grad){
    dx <- c(0, diff(x))
    loglik <- Lhat_eta(x, w, eta)$ll
    liknew <- Lhat_eta(x, w, etanew)$ll
    #dirder <- as.numeric(t(grad) %*% (etanew - eta))
    dirder <- crossprod(grad, etanew - eta)
    iter <- 0
    
    while ((liknew < loglik) && (iter < 20)){
        iter <- iter + 1
        etanew <- (eta + etanew)/2
        liknew <- Lhat_eta(x, w, etanew)$ll
        dirder <- dirder/2
    }
    
    t0 <- (2 - 2 * (liknew - loglik)/dirder)^(-1)
    
    if (t0 < 1){eta <- (1 - t0) * eta + t0 * etanew} else {eta <- etanew}
    
    return(eta)
}
