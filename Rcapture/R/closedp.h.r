"closedp.h" <- function(X, dfreq=FALSE, m="Mh", h="Poisson", a=2)
{

        X<-as.matrix(X)
        t <- ifelse(dfreq,dim(X)[2]-1,dim(X)[2])

    #####################################################################################################################################
    # Validation des arguments fournis en entrée
    
    # Argument dfreq
    if(!is.logical(dfreq)||!isTRUE(all.equal(length(dfreq),1))) stop("'dfreq' must be a logical object of length 1")

    # Argument X
    X <- as.matrix(X)
    if (dfreq)
    {
        if (any(X[,1:t]!=1&X[,1:t]!=0)) stop("every columns of 'X' but the last one must contain only zeros and ones")
        if (any((X[,t+1]%%1)!=0)) stop("the last column of 'X' must contain capture history frequencies, therefore integers")
    } else {
        if(any(X!=1&X!=0)) stop("'X' must contain only zeros and ones")
    }
    
    # Argument m
    if(!isTRUE(all.equal(length(m),1))) stop("'m' must be of length 1")
    if(!identical(m,"Mh")&&!identical(m,"Mth"))
              stop("'m' can only take the values \"Mh\" and \"Mth\"")
    
    # Argument h
    if(!isTRUE(all.equal(length(h),1))) stop("'h' must be of length 1")
    if(!is.function(h)&&!identical(h,"Poisson"))
    stop("'h' must be a function or a character string taking the value \"Poisson\"")
    
    # Argument a
    if(!isTRUE(all.equal(length(a),1))) stop("'a' must be of length 1")
    if (!is.numeric(a)) stop("'a' must be a numeric value")
    
    #####################################################################################################################################

        histpos <- histpos.t(t)
        Y <- histfreq.t(X,dfreq=dfreq)
        nbcap <- apply(histpos, 1, sum)

        mX2 <- if (identical(h,"Poisson")) a^nbcap - 1 else if (is.function(h)) h(nbcap)
        mX <- if (identical(m,"Mh")) cbind(nbcap,mX2) else if(identical(m,"Mth")) cbind(histpos,mX2)
        anaM <- glm(Y~mX,family=poisson)
        NM <- sum(Y)+exp(anaM$coef[1]) # calcul de la taille de la population N
        varcovM <- summary(anaM)$cov.unscaled
        erreurtypeM <- sqrt(exp(anaM$coef[1])+(exp(2*anaM$coef[1]))*varcovM[1,1])
        M <- matrix(c(NM,erreurtypeM,anaM$dev,anaM$df.residual,anaM$aic),nrow=1)


        # Préparation des sorties
        closedp.call<-match.call()
        modelname <- if (is.function(h)) paste(m,closedp.call$h) else if(identical(h,"Poisson")) paste(m,paste(h,a,sep=""))
        dimnames(M) <- list(modelname,c("abundance","stderr","deviance","df","AIC"))
        ans <- list(n=sum(Y),results=M,glm=anaM)
        class(ans) <- "closedp.custom"
        ans

}
