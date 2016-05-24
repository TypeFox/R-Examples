

####
## PRINT METHOD
####
##

print.mvnBvs <- function(x, digits=3, ...)
{
	nChain = length(x)-1
    
    nChain = length(x)-1
    p <- dim(x$chain1$B.p)[1]
    q <- dim(x$chain1$B.p)[2]
    
    nS <- dim(x$chain1$B.p)[3]
    value <- list(model=class(x)[2])
    cov.names <- x$chain1$covNames
    out.names <- colnames(x$chain1$B.p[,,1], do.NULL=FALSE, prefix = "Outcome.")
    
    cat("\nMultivariate Bayesian Variable Selection \n")
    nChain = x$setup$nChain
    
    if(class(x)[2] == "factor-analytic")
    {
        ##
        cat("\nCovariance Structure: Factor-analytic \n")
    }
    if(class(x)[2] == "unstructured")
    {
        ##
        cat("\nCovariance Structure: Unstructured \n")
    }
    ##
    cat("\nNumber of chains:    ", nChain,"\n")
    ##
    cat("Number of scans:     ", x$setup$numReps,"\n")
    ##
    cat("Thinning:            ", x$setup$thin,"\n")
    ##
    cat("Percentage of burnin: ", x$setup$burninPerc*100, "%\n", sep = "")
    
    cat("\n#####")
    
    #B and Gamma
    B <- array(NA, c(p,q, nS*nChain))
    Gamma <- array(NA, c(p,q, nS*nChain))
    B[,,1:nS] <- x$chain1$B.p
    Gamma[,,1:nS] <- x$chain1$gamma.p
    
    if(nChain > 1){
        for(i in 2:nChain){
            nam <- paste("chain", i, sep="")
            B[,,(nS*(i-1)+1):(nS*i)] <- x[[nam]]$B.p
            Gamma[,,(nS*(i-1)+1):(nS*i)] <- x[[nam]]$gamma.p
        }
    }
    
    Gamma.Mean	<- apply(Gamma, c(1,2), mean)
    B.Mean <- matrix(NA, p, q)
    B.Sd <-	matrix(NA, p, q)
    for(k in 1:p)
    {
        for(j in 1:q)
        {
            if(any(B[k,j,]!=0))
            {
                B.Mean[k, j] <- mean(B[k,j,B[k,j,]!=0])
                B.Sd[k, j] <- sd(B[k,j,B[k,j,]!=0])
            }else
            {
                B.Mean[k, j] <- 0
                B.Sd[k, j] <- 0
            }
        }
    }
    rownames(Gamma.Mean) <- cov.names
    rownames(B.Mean) <- cov.names
    rownames(B.Sd) <- cov.names
    
    colnames(Gamma.Mean) <- out.names
    colnames(B.Mean) <- out.names
    colnames(B.Sd) <- out.names
    
    output.BG <- list()
    output.BG[[cov.names[1]]] <- cbind(B.Mean[1,], B.Sd[1,], Gamma.Mean[1,])
    colnames(output.BG[[cov.names[1]]]) <- c("beta|gamma=1", "SD", "gamma")
    
    if(p > 1){
        for(i in 2:p){
            nam <- cov.names[i]
            output.BG[[nam]] <- cbind(B.Mean[i,], B.Sd[i,], Gamma.Mean[i,])
            colnames(output.BG[[nam]]) <- c("beta|gamma=1", "SD", "gamma")
        }
    }
    
    ##
    cat("\nRegression Coefficients and Inclusion Probabilities \n")
    print(output.BG, digits=digits)
}


####
## SUMMARY METHOD
####
##

summary.mvnBvs <- function(object, digits=3, ...)
{
    x <- object
	nChain = length(x)-1
    p <- dim(x$chain1$B.p)[1]
    q <- dim(x$chain1$B.p)[2]
    
    nS <- dim(x$chain1$B.p)[3]
    value <- list(model=class(x)[2])
    cov.names <- x$chain1$covNames
    out.names <- colnames(x$chain1$B.p[,,1], do.NULL=FALSE, prefix = "Outcome.")

    # estimates
    ##
    
    #B and Gamma
    B <- array(NA, c(p,q, nS*nChain))
    Gamma <- array(NA, c(p,q, nS*nChain))
    B[,,1:nS] <- x$chain1$B.p
    Gamma[,,1:nS] <- x$chain1$gamma.p
    
    if(nChain > 1){
        for(i in 2:nChain){
            nam <- paste("chain", i, sep="")
            B[,,(nS*(i-1)+1):(nS*i)] <- x[[nam]]$B.p
            Gamma[,,(nS*(i-1)+1):(nS*i)] <- x[[nam]]$gamma.p
        }
    }
    
    Gamma.Mean	<- apply(Gamma, c(1,2), mean)
    B.Mean <- matrix(NA, p, q)
    B.Sd <-	matrix(NA, p, q)
    for(k in 1:p)
    {
        for(j in 1:q)
        {
            if(any(B[k,j,]!=0))
            {
                B.Mean[k, j] <- mean(B[k,j,B[k,j,]!=0])
                B.Sd[k, j] <- sd(B[k,j,B[k,j,]!=0])
            }else
            {
                B.Mean[k, j] <- 0
                B.Sd[k, j] <- 0
            }
        }
    }
    rownames(Gamma.Mean) <- cov.names
    rownames(B.Mean) <- cov.names
    rownames(B.Sd) <- cov.names
    
    colnames(Gamma.Mean) <- out.names
    colnames(B.Mean) <- out.names
    colnames(B.Sd) <- out.names
    
    output.BG <- list()
    output.BG[[cov.names[1]]] <- cbind(B.Mean[1,], B.Sd[1,], Gamma.Mean[1,])
    colnames(output.BG[[cov.names[1]]]) <- c("beta|gamma=1", "SD", "gamma")
    
    if(p > 1){
        for(i in 2:p){
            nam <- cov.names[i]
            output.BG[[nam]] <- cbind(B.Mean[i,], B.Sd[i,], Gamma.Mean[i,])
            colnames(output.BG[[nam]]) <- c("beta|gamma=1", "SD", "gamma")
        }
    }
    
    #beta0
    beta0 <- x$chain1$beta0.p
    if(nChain > 1){
        for(i in 2:nChain){
            nam <- paste("chain", i, sep="")
            beta0 <- rbind(beta0, x[[nam]]$beta0.p)
        }
    }

    beta0.Med <- apply(beta0, 2, median)
    beta0.Mean <- apply(beta0, 2, mean)
    beta0.Sd <- apply(beta0, 2, sd)
    beta0.Ub <- apply(beta0, 2, quantile, prob = 0.975)
    beta0.Lb <- apply(beta0, 2, quantile, prob = 0.025)
    
    output.beta0 <- cbind(beta0.Mean, beta0.Lb, beta0.Ub)
    dimnames(output.beta0) <- list(out.names, c("beta0", "LL", "UL"))
    
    value$BetaGamma <- output.BG
    value$beta0 <- output.beta0
    
    
    if(class(x)[2] == "unstructured")
    {
        
        #Sigma
        Sigma <- array(NA, c(q,q, nS*nChain))
        Sigma[,,1:nS] <- x$chain1$Sigma.p
        
        if(nChain > 1){
            for(i in 2:nChain){
                nam <- paste("chain", i, sep="")
                Sigma[,,(nS*(i-1)+1):(nS*i)] <- x[[nam]]$Sigma.p
            }
        }
        
        Sigma.Med <- apply(Sigma, c(1,2), median)
        Sigma.Sd <- apply(Sigma, c(1,2), sd)
        Sigma.Ub <- apply(Sigma, c(1,2), quantile, prob = 0.975)
        Sigma.Lb <- apply(Sigma, c(1,2), quantile, prob = 0.025)
        
        dimnames(Sigma.Med) <- list(c(rep("", q)), c("Sigma-PM", rep("", q-1)))
        dimnames(Sigma.Sd) <- list(c(rep("", q)), c("Sigma-SD", rep("", q-1)))
        dimnames(Sigma.Lb) <- list(c(rep("", q)), c("Sigma-LL", rep("", q-1)))
        dimnames(Sigma.Ub) <- list(c(rep("", q)), c("Sigma-UL", rep("", q-1)))
        
        value$Sigma.PM <- Sigma.Med
        value$Sigma.SD <- Sigma.Sd
        value$Sigma.LL <- Sigma.Lb
        value$Sigma.UL <- Sigma.Ub

    }
    
    if(class(x)[2] == "factor-analytic")
    {
        
        #lambda
        lambda <- x$chain1$lambda.p
        if(nChain > 1){
            for(i in 2:nChain){
                nam <- paste("chain", i, sep="")
                lambda <- rbind(lambda, x[[nam]]$lambda.p)
            }
        }
        
        lambda.Med <- apply(lambda, 2, median)
        lambda.Mean <- apply(lambda, 2, mean)
        lambda.Sd <- apply(lambda, 2, sd)
        lambda.Ub <- apply(lambda, 2, quantile, prob = 0.975)
        lambda.Lb <- apply(lambda, 2, quantile, prob = 0.025)
        
        output.lambda <- cbind(lambda.Mean, lambda.Lb, lambda.Ub)
        dimnames(output.lambda) <- list(out.names, c("lambda", "LL", "UL"))
        
        value$lambda <- output.lambda
        
        sigSq.p <- x$chain1$sigSq.p
        
        if(nChain > 1){
            for(i in 2:nChain){
                nam <- paste("chain", i, sep="")
                sigSq.p <- rbind(sigSq.p, x[[nam]]$sigSq.p)
            }
        }
        
        sigSq.pMed <- apply(sigSq.p, 2, median)
        sigSq.pUb <- apply(sigSq.p, 2, quantile, prob = 0.975)
        sigSq.pLb <- apply(sigSq.p, 2, quantile, prob = 0.025)
        
        output.sigSq <- cbind(sigSq.pMed, sigSq.pLb, sigSq.pUb)
        dimnames(output.sigSq) <- list("", c( "sigmaSq", "LL", "UL"))
        
        value$sigmaSq <- output.sigSq
    }
    
    
    value$setup <- x$setup
    class(value) <- "summ.mvnBvs"
    return(value)
}



####
## PRINT.SUMMARY METHOD
####
##


print.summ.mvnBvs <- function(x, digits=3, ...)
{
    
    cat("\nMultivariate Bayesian Variable Selection \n")
    nChain = x$setup$nChain
    
    if(x$model == "factor-analytic")
    {
        ##
        cat("\nCovariance Structure: Factor-analytic \n")
    }
    if(x$model == "unstructured")
    {
        ##
        cat("\nCovariance Structure: Unstructured \n")
    }
    cat("\n#####")
    
    ##
    cat("\nRegression Coefficients and Inclusion Probabilities \n")
    print(x$BetaGamma, digits=digits)
    
    cat("\n#####")
    
    ##
    cat("\nIntercepts \n")
    print(x$beta0, digits=digits)
    
    if(x$model == "factor-analytic")
    {
        cat("\n#####")
        ##
        cat("\nFactor Loadings \n")
        print(x$lambda, digits=digits)
        
        cat("\n#####")
        ##
        cat("\nResidual Variance \n")
        print(x$sigmaSq, digits=digits)
    }
    
}














