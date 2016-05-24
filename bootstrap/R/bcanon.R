"bcanon" <- function(x,nboot,theta,...,alpha =
                     c(.025,.05,.1,.16,.84,.9,.95,.975)) { 
    if (!all(alpha < 1) || !all(alpha > 0))
      stop("All elements of alpha must be in (0,1)")

    alpha_sorted <- sort(alpha)
    if (nboot <= 1/min(alpha_sorted[1],1-alpha_sorted[length(alpha_sorted)]))
      warning("nboot is not large enough to estimate your chosen alpha.")

    call <- match.call()
    n <- length(x)
    thetahat <- theta(x,...)
    bootsam<- matrix(sample(x,size=n*nboot,replace=TRUE),nrow=nboot)

    thetastar <- apply(bootsam,1,theta,...)
    z0 <- qnorm(sum(thetastar<thetahat)/nboot)
    
    u <- rep(0,n)
    for(i in 1:n){
        u[i] <- theta(x[-i],...)
    }
    uu <- mean(u)-u
    acc <- sum(uu*uu*uu)/(6*(sum(uu*uu))^1.5)
    
    zalpha <- qnorm(alpha)
    
    tt <- pnorm(z0+ (z0+zalpha)/(1-acc*(z0+zalpha)))
    
    confpoints <- quantile(x=thetastar,probs=tt,type=1)
    names(confpoints) <- NULL
    confpoints <- cbind(alpha,confpoints)
    dimnames(confpoints)[[2]] <- c("alpha","bca point")
    return(list(confpoints=confpoints, 
                z0=z0, 
                acc=acc, 
                u=u, 
                call=call))
}
