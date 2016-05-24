"closedp.mX" <- function(X,dfreq=FALSE,mX,mname="Customized model")
{

        X<-as.matrix(X)
        t <- ifelse(dfreq,dim(X)[2]-1,dim(X)[2])

    #####################################################################################################################################
    # Validation des arguments fournis en entrée

    # Argument dfreq
    if(!is.logical(dfreq)||!isTRUE(all.equal(length(dfreq),1))) stop("'dfreq' must be a logical object of length 1")
    
    # Argument X
    if (dfreq)
    {
        if (any(X[,1:t]!=1&X[,1:t]!=0)) stop("every columns of 'X' but the last one must contain only zeros and ones")
        if (any((X[,t+1]%%1)!=0)) stop("the last column of 'X' must contain capture history frequencies, therefore integers")
    } else {
        if(any(X!=1&X!=0)) stop("'X' must contain only zeros and ones")
    }
    
    # Argument mX
    mX<-as.matrix(mX)
    if (!isTRUE(all.equal(2^t-1,dim(mX)[1]))) stop("'mX' must have 2^t-1 rows")
    
    # Argument mname
    if(!is.character(mname)) stop("'mname' must be a character string specifying the model's name")

    #####################################################################################################################################

        Y <- histfreq.t(X,dfreq=dfreq)

        anaM <- glm(Y~mX,family=poisson)
        NM <- sum(Y)+exp(anaM$coef[1])
        varcovM <- summary(anaM)$cov.unscaled
        erreurtypeM <- sqrt(exp(anaM$coef[1])+(exp(2*anaM$coef[1]))*varcovM[1,1])
        M <- matrix(c(NM,erreurtypeM,anaM$dev,anaM$df.residual,anaM$aic),nrow=1)

         
        # Préparation des sorties
        dimnames(M) <- list(mname,c("abundance","stderr","deviance","df","AIC"))
        ans <- list(n=sum(Y),results=M,glm=anaM)
        class(ans) <- "closedp.custom"
        ans

}


print.closedp.custom <- function(x, ...) {
        cat("\nNumber of captured units:",x$n,"\n\n")
        cat("Abundance estimation and model fit:\n")
        tableau <- x$results
        tableau[,c(1,2)] <- round(tableau[,c(1,2)],1)
        tableau[,4] <- round(tableau[,4],0)
        tableau[,c(3,5)] <- round(tableau[,c(3,5)],3)       
        print.default(tableau, print.gap = 2, quote = FALSE, right=TRUE)
        cat("\n")
        invisible(x)
}

boxplot.closedp.custom <- function(x,...) {
        boxplot((x$glm$y-fitted(x$glm))/sqrt(fitted(x$glm)),main="Boxplot of Pearson Residuals for the customized model")     
}
