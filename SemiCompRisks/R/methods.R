

####
## PRINT METHOD
####
##
print.Freq <- function(x, digits=3, ...)
{   
    obj <- x
  ##
  logEst <- obj$estimate
  logSE  <- sqrt(diag(obj$Finv))
  value  <- cbind(logEst, logSE, logEst - 1.96*logSE, logEst + 1.96*logSE)
  ## 
  dimnames(value) <- list(obj$myLabels, c( "Estimate", "SE", "LL", "UL"))
  
  ##
  if(class(obj)[2] == "Surv")
  {
    ##
    cat("\nAnalysis of independent univariate time-to-event data \n")
    ## 
    #cat("\nBaseline hazard function components:\n")
    #print(round(value[c(1:2),], digits=digits))
    ##
    cat("\nRegression coefficients:\n")
    print(round(value[-c(1:2),], digits=digits))
  }

  ##
  if(class(obj)[2] == "ID")
  {
    ##
    cat("\nAnalysis of independent semi-competing risks data \n")
    cat(class(obj)[5], "assumption for h3\n")
    ##
    #cat("\nBaseline hazard function components:\n")
    #print(round(value[c(1:6),], digits=digits))
    ##
    value_theta <- matrix(exp(value[7,]), ncol = 4)
    dimnames(value_theta) <- list("", c( "Estimate", "SE", "LL", "UL"))
    value_theta[1,2] <- value[7,2] * exp(value[7,1])
    
    cat("\nVariance of frailties, theta:\n")
    if(obj$frailty == TRUE) print(round(value_theta, digits=digits))
    if(obj$frailty == FALSE) print("NA")    
    ##
    cat("\nRegression coefficients:\n")
    if(obj$frailty == TRUE) print(round(value[-c(1:7),], digits=digits))
    if(obj$frailty == FALSE) print(round(value[-c(1:6),], digits=digits))  
  }
  
  ##
  invisible()
}

print.Bayes <- function(x, digits=3, ...)
{
	nChain = x$setup$nChain    
    
    if(class(x)[2] == "ID")
    {
        if(class(x)[3] == "Cor")
        {
            ##
            cat("\nAnalysis of cluster-correlated semi-competing risks data \n")
        }
        if(class(x)[3] == "Ind")
        {
            ##
            cat("\nAnalysis of independent semi-competing risks data \n")
        }
        ##
        cat(x$setup$model, "assumption for h3\n")
    }
    if(class(x)[2] == "Surv")
    {
        if(class(x)[3] == "Cor")
        {
            ##
            cat("\nAnalysis of cluster-correlated univariate time-to-event data \n")
        }
        if(class(x)[3] == "Ind")
        {
            ##
            cat("\nAnalysis of independent univariate time-to-event data \n")
        }
    }
    
    ##
    cat("\nNumber of chains:    ", nChain,"\n")
    ##
    cat("Number of scans:     ", x$setup$numReps,"\n")
    ##
    cat("Thinning:            ", x$setup$thin,"\n")
    ##
    cat("Percentage of burnin: ", x$setup$burninPerc*100, "%\n", sep = "")
    
  
    # convergence diagnostics
    
    if(nChain > 1){
        
        cat("\n######\n")
        cat("Potential Scale Reduction Factor\n")
        
        if(class(x)[2] == "ID")
        {

            theta <- x$chain1$theta.p
            for(i in 2:nChain){
                nam <- paste("chain", i, sep = "")
                theta <- cbind(theta, x[[nam]]$theta.p)
            }
            psrftheta <- matrix(calcPSR(theta), 1, 1)
            dimnames(psrftheta) <- list("", "")
            
            
            cat("\nVariance of frailties, theta:")
            print(round(psrftheta, digits=digits))
            

            beta.names <- unique(c(x$chain1$covNames1, x$chain1$covNames2, x$chain1$covNames3))
            nP         <- length(beta.names)
            output <- matrix(NA, nrow=nP, ncol=3)
            dimnames(output) <- list(beta.names, c("beta1", "beta2", "beta3"))
            
            if(length(x$chain1$beta1.p) != 0){
                
                #beta1
                
                p1	= dim(x$chain1$beta1.p)[2]
                
                psrfBeta1 <- rep(NA, p1)
                for(j in 1:p1){
                    
                    #namPara = paste("beta_", j, sep = "")
                    
                    beta1 <- x$chain1$beta1[,j]
                    for(i in 2:nChain){
                        nam <- paste("chain", i, sep = "")
                        beta1 <- cbind(beta1, x[[nam]]$beta1[,j])
                    }
                    psrfBeta1[j] <- calcPSR(beta1)
                }
                
                for(i in 1:nP)
                {
                    for(k in 1:p1) if(x$chain1$covNames1[k] == beta.names[i]) output[i,1] <- psrfBeta1[k]
                }
                
            }
            
            if(length(x$chain1$beta2.p) != 0){
                
                #beta2
                
                p2	= dim(x$chain1$beta2.p)[2]
                
                psrfBeta2 <- rep(NA, p2)
                for(j in 1:p2){
                    
                    #namPara = paste("beta_", j, sep = "")
                    
                    beta2 <- x$chain1$beta2[,j]
                    for(i in 2:nChain){
                        nam <- paste("chain", i, sep = "")
                        beta2 <- cbind(beta2, x[[nam]]$beta2[,j])
                    }
                    psrfBeta2[j] <- calcPSR(beta2)
                }
                for(i in 1:nP)
                {
                    for(k in 1:p2) if(x$chain1$covNames2[k] == beta.names[i]) output[i,2] <- psrfBeta2[k]
                }
            }
            
            if(length(x$chain1$beta3.p) != 0){
                
                #beta3
                
                p3	= dim(x$chain1$beta3.p)[2]
                
                psrfBeta3 <- rep(NA, p3)
                for(j in 1:p3){
                    
                    #namPara = paste("beta_", j, sep = "")
                    
                    beta3 <- x$chain1$beta3[,j]
                    for(i in 2:nChain){
                        nam <- paste("chain", i, sep = "")
                        beta3 <- cbind(beta3, x[[nam]]$beta3[,j])
                    }
                    psrfBeta3[j] <- calcPSR(beta3)
                }
                for(i in 1:nP)
                {
                    for(k in 1:p3) if(x$chain1$covNames3[k] == beta.names[i]) output[i,3] <- psrfBeta3[k]
                }
            }
            
            if(nP > 0)
            {

                cat("\nRegression coefficients:\n")
                output.coef <- output
                print(round(output.coef, digits=digits))
            }


            ##
            cat("\nBaseline hazard function components:\n")
            
            if(class(x)[4] == "WB")
            {
                ##
                # alpha
                
                alpha <- x$chain1$alpha1.p
                for(i in 2:nChain){
                    nam <- paste("chain", i, sep = "")
                    alpha <- cbind(alpha, x[[nam]]$alpha1.p)
                }
                psrfAlpha1 <- calcPSR(alpha)
                
                alpha <- x$chain1$alpha2.p
                for(i in 2:nChain){
                    nam <- paste("chain", i, sep = "")
                    alpha <- cbind(alpha, x[[nam]]$alpha2.p)
                }
                psrfAlpha2 <- calcPSR(alpha)
                
                alpha <- x$chain1$alpha3.p
                for(i in 2:nChain){
                    nam <- paste("chain", i, sep = "")
                    alpha <- cbind(alpha, x[[nam]]$alpha3.p)
                }
                psrfAlpha3 <- calcPSR(alpha)
                
                # kappa
                
                kappa <- x$chain1$kappa1.p
                for(i in 2:nChain){
                    nam <- paste("chain", i, sep = "")
                    kappa <- cbind(kappa, x[[nam]]$kappa1.p)
                }
                psrfKappa1 <- calcPSR(kappa)
                
                kappa <- x$chain1$kappa2.p
                for(i in 2:nChain){
                    nam <- paste("chain", i, sep = "")
                    kappa <- cbind(kappa, x[[nam]]$kappa2.p)
                }
                psrfKappa2 <- calcPSR(kappa)
                
                kappa <- x$chain1$kappa3.p
                for(i in 2:nChain){
                    nam <- paste("chain", i, sep = "")
                    kappa <- cbind(kappa, x[[nam]]$kappa3.p)
                }
                psrfKappa3 <- calcPSR(kappa)
                
                bh_WB <- matrix(c(psrfKappa1, psrfKappa2, psrfKappa3, psrfAlpha1, psrfAlpha2, psrfAlpha3), 2, 3, byrow = T)
                dimnames(bh_WB) <- list(c("kappa", "alpha"), c("h1", "h2", "h3"))
                print(round(bh_WB, digits=digits))
                
            }
            
            if(class(x)[4] == "PEM")
            {
                ##
                
                ntime1  = length(x$chain1$time_lambda1)
                ntime2  = length(x$chain1$time_lambda2)
                ntime3  = length(x$chain1$time_lambda3)
                
                # lambda's
                
                psrfLam <- rep(NA, ntime1)
                
                for(j in 1:ntime1){
                    
                    lambda1 <- x$chain1$lambda1.fin[,j]
                    for(i in 2:nChain){
                        nam <- paste("chain", i, sep = "")
                        lambda1 <- cbind(lambda1, x[[nam]]$lambda1.fin[,j])
                    }
                    psrfLam[j] <- calcPSR(lambda1)
                }
                
                cat("\nlambda1: summary statistics", "\n")
                print(round(summary(psrfLam), digits=digits))
                
                psrfLam <- rep(NA, ntime2)
                
                for(j in 1:ntime2){
                    
                    lambda2 <- x$chain1$lambda2.fin[,j]
                    for(i in 2:nChain){
                        nam <- paste("chain", i, sep = "")
                        lambda2 <- cbind(lambda2, x[[nam]]$lambda2.fin[,j])
                    }
                    psrfLam[j] <- calcPSR(lambda2)
                }
                
                cat("\nlambda2: summary statistics", "\n")
                print(round(summary(psrfLam), digits=digits))
                
                psrfLam <- rep(NA, ntime3)
                
                for(j in 1:ntime3){
                    
                    lambda3 <- x$chain1$lambda3.fin[,j]
                    for(i in 2:nChain){
                        nam <- paste("chain", i, sep = "")
                        lambda3 <- cbind(lambda3, x[[nam]]$lambda3.fin[,j])
                    }
                    psrfLam[j] <- calcPSR(lambda3)
                }
                
                cat("\nlambda3: summary statistics", "\n")
                print(round(summary(psrfLam), digits=digits))
                
                
                # mu_lam
                
                mu <- x$chain1$mu_lam1.p
                for(i in 2:nChain){
                    nam <- paste("chain", i, sep = "")
                    mu <- cbind(mu, x[[nam]]$mu_lam1.p)
                }
                psrfMu1 <- calcPSR(mu)
                
                mu <- x$chain1$mu_lam2.p
                for(i in 2:nChain){
                    nam <- paste("chain", i, sep = "")
                    mu <- cbind(mu, x[[nam]]$mu_lam2.p)
                }
                psrfMu2 <- calcPSR(mu)
                
                mu <- x$chain1$mu_lam3.p
                for(i in 2:nChain){
                    nam <- paste("chain", i, sep = "")
                    mu <- cbind(mu, x[[nam]]$mu_lam3.p)
                }
                psrfMu3 <- calcPSR(mu)
                
                # sigSq_lam
                
                sig <- x$chain1$sigSq_lam1.p
                for(i in 2:nChain){
                    nam <- paste("chain", i, sep = "")
                    sig <- cbind(sig, x[[nam]]$sigSq_lam1.p)
                }
                psrfSig1 <- calcPSR(sig)

                
                sig <- x$chain1$sigSq_lam2.p
                for(i in 2:nChain){
                    nam <- paste("chain", i, sep = "")
                    sig <- cbind(sig, x[[nam]]$sigSq_lam2.p)
                }
                psrfSig2 <- calcPSR(sig)

                
                sig <- x$chain1$sigSq_lam3.p
                for(i in 2:nChain){
                    nam <- paste("chain", i, sep = "")
                    sig <- cbind(sig, x[[nam]]$sigSq_lam3.p)
                }
                psrfSig3 <- calcPSR(sig)
                
                # J
                
                J <- x$chain1$K1.p
                for(i in 2:nChain){
                    nam <- paste("chain", i, sep = "")
                    J <- cbind(J, x[[nam]]$K1.p)
                }
                psrfJ1 <- calcPSR(J)
                
                J <- x$chain1$K2.p
                for(i in 2:nChain){
                    nam <- paste("chain", i, sep = "")
                    J <- cbind(J, x[[nam]]$K2.p)
                }
                psrfJ2 <- calcPSR(J)
                
                J <- x$chain1$K3.p
                for(i in 2:nChain){
                    nam <- paste("chain", i, sep = "")
                    J <- cbind(J, x[[nam]]$K3.p)
                }
                psrfJ3 <- calcPSR(J)
                
                
                bh_PEM <- matrix(c(psrfMu1, psrfMu2, psrfMu3, psrfSig1, psrfSig2, psrfSig3, psrfJ1, psrfJ2, psrfJ3), 3, 3, byrow = T)
                dimnames(bh_PEM) <- list(c("mu", "sigmaSq", "K"), c("h1", "h2", "h3"))
                cat("\n")
                print(round(bh_PEM, digits=digits))
            }

        }
        
        if(class(x)[2] == "Surv")
        {
            beta.names <- c(x$chain1$covNames)
            nP         <- length(beta.names)
            output <- matrix(NA, nrow=nP, ncol=1)
            dimnames(output) <- list(beta.names, c("beta"))
            
            if(length(x$chain1$beta.p) != 0){
                
                #beta
                
                p	= dim(x$chain1$beta.p)[2]
                
                psrfBeta <- rep(NA, p)
                for(j in 1:p){
                    
                    beta <- x$chain1$beta[,j]
                    for(i in 2:nChain){
                        nam <- paste("chain", i, sep = "")
                        beta <- cbind(beta, x[[nam]]$beta[,j])
                    }
                    psrfBeta[j] <- calcPSR(beta)
                }
                
                for(i in 1:nP)
                {
                    for(k in 1:p) if(x$chain1$covNames[k] == beta.names[i]) output[i,1] <- psrfBeta[k]
                }
                
            }
            
            if(nP > 0)
            {
                cat("\nRegression coefficients:\n")
                output.coef <- output
                print(round(output.coef, digits=digits))
            }
            
            
            if(class(x)[4] == "WB")
            {
                ##
                # alpha
                
                alpha <- x$chain1$alpha.p
                for(i in 2:nChain){
                    nam <- paste("chain", i, sep = "")
                    alpha <- cbind(alpha, x[[nam]]$alpha.p)
                }
                psrfAlpha <- calcPSR(alpha)
                
                
                # kappa
                
                kappa <- x$chain1$kappa.p
                for(i in 2:nChain){
                    nam <- paste("chain", i, sep = "")
                    kappa <- cbind(kappa, x[[nam]]$kappa.p)
                }
                psrfKappa <- calcPSR(kappa)
                
                
                bh_WB <- matrix(c(psrfKappa, psrfAlpha), 2, 1, byrow = T)
                dimnames(bh_WB) <- list(c("kappa", "alpha"), c("h"))
                print(round(bh_WB, digits=digits))
            }
            
            if(class(x)[4] == "PEM")
            {
                ##
                ntime  = length(x$chain1$time_lambda)
                
                # lambda
                
                psrfLam <- rep(NA, ntime)
                
                for(j in 1:ntime){
                    
                    namPara = paste("beta_", j, sep = "")
                    
                    lambda <- x$chain1$lambda.fin[,j]
                    for(i in 2:nChain){
                        nam <- paste("chain", i, sep = "")
                        lambda <- cbind(lambda, x[[nam]]$lambda.fin[,j])
                    }
                    psrfLam[j] <- calcPSR(lambda)
                }
                
                cat("\n lambda: summary statistics", "\n")
                print(round(summary(psrfLam), 2))
                
                
                # mu_lam
                
                mu <- x$chain1$mu_lam.p
                for(i in 2:nChain){
                    nam <- paste("chain", i, sep = "")
                    mu <- cbind(mu, x[[nam]]$mu_lam.p)
                }
                psrfMu <- calcPSR(mu)
                
                # sigSq_lam
                
                sig <- x$chain1$sigSq_lam.p
                for(i in 2:nChain){
                    nam <- paste("chain", i, sep = "")
                    sig <- cbind(sig, x[[nam]]$sigSq_lam.p)
                }
                psrfSig <- calcPSR(sig)
                
                
                # J
                
                J <- x$chain1$K.p
                for(i in 2:nChain){
                    nam <- paste("chain", i, sep = "")
                    J <- cbind(J, x[[nam]]$K.p)
                }
                psrfJ <- calcPSR(J)
                
                
                bh_PEM <- matrix(c(psrfMu, psrfSig, psrfJ), 3, 1, byrow = T)
                dimnames(bh_PEM) <- list(c("mu", "sigmaSq", "K"), c("h"))
                cat("\n")
                print(round(bh_PEM, digits=digits))
            }

        }
        
        
    }
    else if(nChain == 1)
    {
        cat("Potential scale reduction factor cannot be calculated. \n")
        cat("The number of chains must be larger than 1. \n")	
    }
    
    cat("\n######\n")
    cat("Estimates\n")
    
    if(class(x)[2] == "ID")
    {
        ##
        cat("\nVariance of frailties, theta:\n")
        
        theta.p <- x$chain1$theta.p
        
        if(nChain > 1){
            for(i in 2:nChain){
                nam <- paste("chain", i, sep="")
                theta.p <- rbind(theta.p, x[[nam]]$theta.p)
            }
        }
        
        theta.pMed <- apply(theta.p, 2, median)
        theta.pSd <- apply(theta.p, 2, sd)
        theta.pUb <- apply(theta.p, 2, quantile, prob = 0.975)
        theta.pLb <- apply(theta.p, 2, quantile, prob = 0.025)
        
        tbl <- matrix(NA, 1, 4)
        dimnames(tbl) <- list("", c( "Estimate", "SD", "LL", "UL"))
        
        tbl[,1]	<- theta.pMed
        tbl[,2]	<- theta.pSd
        tbl[,3]	<- theta.pLb
        tbl[,4]	<- theta.pUb
        
        print(round(tbl, digits=digits))
        

        ##
        
        tbl_beta <- NULL
        
        if(length(x$chain1$beta1.p) != 0){
            
            
            p1	= dim(x$chain1$beta1.p)[2]
            beta.p <- x$chain1$beta1.p
            
            if(nChain > 1){
                for(i in 2:nChain){
                    nam <- paste("chain", i, sep="")
                    beta.p <- rbind(beta.p, x[[nam]]$beta1.p)
                }
            }
            
            
            beta.pMed <- apply(beta.p, 2, median)
            beta.pSd <- apply(beta.p, 2, sd)
            beta.pUb <- apply(beta.p, 2, quantile, prob = 0.975)
            beta.pLb <- apply(beta.p, 2, quantile, prob = 0.025)
            
            tbl1 <- matrix(NA, p1, 4)
            rownames(tbl1) <- x$chain1$covNames1
            
            tbl1[,1]	<- beta.pMed
            tbl1[,2]	<- beta.pSd
            tbl1[,3]	<- exp(beta.pLb)
            tbl1[,4]	<- exp(beta.pUb)
            tbl_beta    <- tbl1
        }
        
        
        
        if(length(x$chain1$beta2.p) != 0){
            
            p2	= dim(x$chain1$beta2.p)[2]
            beta.p <- x$chain1$beta2.p
            
            if(nChain > 1){
                for(i in 2:nChain){
                    nam <- paste("chain", i, sep="")
                    beta.p <- rbind(beta.p, x[[nam]]$beta2.p)
                }
            }
            
            
            beta.pMed <- apply(beta.p, 2, median)
            beta.pSd <- apply(beta.p, 2, sd)
            beta.pUb <- apply(beta.p, 2, quantile, prob = 0.975)
            beta.pLb <- apply(beta.p, 2, quantile, prob = 0.025)
            
            tbl2 <- matrix(NA, p2, 4)
            rownames(tbl2) <- x$chain1$covNames2
            
            tbl2[,1]	<- beta.pMed
            tbl2[,2]	<- beta.pSd
            tbl2[,3]	<- exp(beta.pLb)
            tbl2[,4]	<- exp(beta.pUb)
            tbl_beta     <- rbind(tbl_beta, tbl2)
        }
        
        if(length(x$chain1$beta3.p) != 0){
            
            p3	= dim(x$chain1$beta3.p)[2]
            beta.p <- x$chain1$beta3.p
            
            if(nChain > 1){
                for(i in 2:nChain){
                    nam <- paste("chain", i, sep="")
                    beta.p <- rbind(beta.p, x[[nam]]$beta3.p)
                }
            }
            
            
            beta.pMed <- apply(beta.p, 2, median)
            beta.pSd <- apply(beta.p, 2, sd)
            beta.pUb <- apply(beta.p, 2, quantile, prob = 0.975)
            beta.pLb <- apply(beta.p, 2, quantile, prob = 0.025)
            
            tbl3 <- matrix(NA, p3, 4)
            rownames(tbl3) <- x$chain1$covNames3
            
            tbl3[,1]	<- beta.pMed
            tbl3[,2]	<- beta.pSd
            tbl3[,3]	<- exp(beta.pLb)
            tbl3[,4]	<- exp(beta.pUb)
            tbl_beta     <- rbind(tbl_beta, tbl3)
        }
        

    }
    
    if(class(x)[2] == "Surv")
    {
        if(length(x$chain1$beta.p) != 0){
            
            p	= dim(x$chain1$beta.p)[2]
            beta.p <- x$chain1$beta.p
            
            if(nChain > 1){
                for(i in 2:nChain){
                    nam <- paste("chain", i, sep="")
                    beta.p <- rbind(beta.p, x[[nam]]$beta.p)
                }
            }
            
            
            beta.pMed <- apply(beta.p, 2, median)
            beta.pSd <- apply(beta.p, 2, sd)
            beta.pUb <- apply(beta.p, 2, quantile, prob = 0.975)
            beta.pLb <- apply(beta.p, 2, quantile, prob = 0.025)
            
            tbl_beta <- matrix(NA, p, 4)
            rownames(tbl_beta) <- x$chain1$covNames
            
            tbl_beta[,1]	<- beta.pMed
            tbl_beta[,2]	<- beta.pSd
            tbl_beta[,3]	<- exp(beta.pLb)
            tbl_beta[,4]	<- exp(beta.pUb)
        }
    }
    
    if(!is.null(tbl_beta))
    {
        cat("\nRegression coefficients:\n")
        colnames(tbl_beta) <- c( "Estimate", "SD", "LL", "UL")
        print(round(tbl_beta, digits=digits))
    }


    
    
    invisible()
}



####
## SUMMARY METHOD
####
##
summary.Freq <- function(object, digits=3, ...)
{
	obj <- object
  ##
  logEst  <- obj$estimate
  logSE   <- sqrt(diag(obj$Finv))
  results <- cbind(logEst, logEst - 1.96*logSE, logEst + 1.96*logSE)
  
  ##
  if(class(obj)[2] == "Surv")
  {
    ##
    #cat("\nRegression coefficients:\n")
    output.coef           <- results[-c(1:2),]
    dimnames(output.coef) <- list(unique(obj$myLabels[-c(1:2)]), c("beta", "LL", "UL"))
    ##
    #cat("\nBaseline hazard function components:\n")
    output.h0           <- results[c(1:2),]
    dimnames(output.h0) <- list(c("Weibull: log-kappa", "Weibull: log-alpha"), c("beta", "LL", "UL"))
    ##  
    value <- list(coef=output.coef, h0=output.h0, code=obj$code, logLike=obj$logLike, nP=nrow(results))
    class(value) <- c("summ.Freq", "Surv")
  }
  

  ##
  if(class(obj)[2] == "ID")
  {
    ##
    nP.0 <- ifelse(obj$frailty, 7, 6)
    nP.1 <- obj$nP[1]
    nP.2 <- obj$nP[2]
    nP.3 <- obj$nP[3]
    ##
    beta.names <- unique(obj$myLabels[-c(1:nP.0)])
    nP         <- length(beta.names)
    ##
    #cat("\nRegression coefficients:\n")
    output <- matrix(NA, nrow=nP, ncol=9)
    dimnames(output) <- list(beta.names, c("beta1", "LL", "UL", "beta2", "LL", "UL", "beta3", "LL", "UL"))
    for(i in 1:nP)
    {
      for(j in 1:nP.1) if(obj$myLabels[nP.0+j] == beta.names[i]) output[i,1:3] <- results[nP.0+j,]
      for(j in 1:nP.2) if(obj$myLabels[nP.0+nP.1+j] == beta.names[i]) output[i,4:6] <- results[nP.0+nP.1+j,]
      for(j in 1:nP.3) if(obj$myLabels[nP.0+nP.1+nP.2+j] == beta.names[i]) output[i,7:9] <- results[nP.0+nP.1+nP.2+j,]
    }
    output.coef <- output
    ##
    #cat("\nVariance of frailties:\n")
    output <- matrix(NA, nrow=1, ncol=3)
    dimnames(output) <- list(c("theta"), c("Estimate", "LL", "UL"))
    if(obj$frailty == TRUE)  output[1,] <- exp(results[7,])
    if(obj$frailty == FALSE) output[1,] <- rep(NA, 3)
    output.theta <- output
    ##
    #cat("\nBaseline hazard function components:\n")
    output <- matrix(NA, nrow=2, ncol=9)
    dimnames(output) <- list(c("Weibull: log-kappa", "Weibull: log-alpha"), c("h1-PM", "LL", "UL", "h2-PM", "LL", "UL", "h3-PM", "LL", "UL"))
    output[1,1:3] <- results[1,]
    output[1,4:6] <- results[3,]
    output[1,7:9] <- results[5,]
    output[2,1:3] <- results[2,]
    output[2,4:6] <- results[4,]
    output[2,7:9] <- results[6,] 
    output.h0 <- output
    ##
    value <- list(coef=output.coef, theta=output.theta, h0=output.h0, code=obj$code, logLike=obj$logLike, nP=nrow(results))
    if(class(obj)[5] == "semi-Markov")
    {
        class(value) <- c("summ.Freq", "ID", "semi-Markov")
    }
    if(class(obj)[5] == "Markov")
    {
        class(value) <- c("summ.Freq", "ID", "Markov")
    }
    
  }
  
  ##
  return(value)
}

summary.Bayes <- function(object, digits=3, ...)
{
    x <- object
	nChain = x$setup$nChain
    
    # convergence diagnostics
    
    psrf <- NULL
    
    if(nChain > 1){
        if(class(x)[2] == "ID")
        {
            theta <- x$chain1$theta.p
            for(i in 2:nChain){
                nam <- paste("chain", i, sep = "")
                theta <- cbind(theta, x[[nam]]$theta.p)
            }
            psrftheta <- matrix(calcPSR(theta), 1, 1)
            dimnames(psrftheta) <- list("", "")

            beta.names <- unique(c(x$chain1$covNames1, x$chain1$covNames2, x$chain1$covNames3))
            nP         <- length(beta.names)
            output <- matrix(NA, nrow=nP, ncol=3)
            dimnames(output) <- list(beta.names, c("beta1", "beta2", "beta3"))
            
            if(length(x$chain1$beta1.p) != 0){
                
                #beta1
                
                p1	= dim(x$chain1$beta1.p)[2]
                
                psrfBeta1 <- rep(NA, p1)
                for(j in 1:p1){
                    
                    #namPara = paste("beta_", j, sep = "")
                    
                    beta1 <- x$chain1$beta1[,j]
                    for(i in 2:nChain){
                        nam <- paste("chain", i, sep = "")
                        beta1 <- cbind(beta1, x[[nam]]$beta1[,j])
                    }
                    psrfBeta1[j] <- calcPSR(beta1)
                }
                
                for(i in 1:nP)
                {
                    for(k in 1:p1) if(x$chain1$covNames1[k] == beta.names[i]) output[i,1] <- psrfBeta1[k]
                }
                
            }
            
            if(length(x$chain1$beta2.p) != 0){
                
                #beta2
                
                p2	= dim(x$chain1$beta2.p)[2]
                
                psrfBeta2 <- rep(NA, p2)
                for(j in 1:p2){
                    
                    #namPara = paste("beta_", j, sep = "")
                    
                    beta2 <- x$chain1$beta2[,j]
                    for(i in 2:nChain){
                        nam <- paste("chain", i, sep = "")
                        beta2 <- cbind(beta2, x[[nam]]$beta2[,j])
                    }
                    psrfBeta2[j] <- calcPSR(beta2)
                }
                for(i in 1:nP)
                {
                    for(k in 1:p2) if(x$chain1$covNames2[k] == beta.names[i]) output[i,2] <- psrfBeta2[k]
                }
            }
            
            if(length(x$chain1$beta3.p) != 0){
                
                #beta3
                
                p3	= dim(x$chain1$beta3.p)[2]
                
                psrfBeta3 <- rep(NA, p3)
                for(j in 1:p3){
                    
                    #namPara = paste("beta_", j, sep = "")
                    
                    beta3 <- x$chain1$beta3[,j]
                    for(i in 2:nChain){
                        nam <- paste("chain", i, sep = "")
                        beta3 <- cbind(beta3, x[[nam]]$beta3[,j])
                    }
                    psrfBeta3[j] <- calcPSR(beta3)
                }
                for(i in 1:nP)
                {
                    for(k in 1:p3) if(x$chain1$covNames3[k] == beta.names[i]) output[i,3] <- psrfBeta3[k]
                }
            }
            
            psrfcoef <- NULL
            if(nP > 0)
            {
                psrfcoef <- output
            }
            
            ##
            if(class(x)[4] == "WB")
            {
                ##
                # alpha
                
                alpha <- x$chain1$alpha1.p
                for(i in 2:nChain){
                    nam <- paste("chain", i, sep = "")
                    alpha <- cbind(alpha, x[[nam]]$alpha1.p)
                }
                psrfAlpha1 <- calcPSR(alpha)
                
                alpha <- x$chain1$alpha2.p
                for(i in 2:nChain){
                    nam <- paste("chain", i, sep = "")
                    alpha <- cbind(alpha, x[[nam]]$alpha2.p)
                }
                psrfAlpha2 <- calcPSR(alpha)
                
                alpha <- x$chain1$alpha3.p
                for(i in 2:nChain){
                    nam <- paste("chain", i, sep = "")
                    alpha <- cbind(alpha, x[[nam]]$alpha3.p)
                }
                psrfAlpha3 <- calcPSR(alpha)
                
                # kappa
                
                kappa <- x$chain1$kappa1.p
                for(i in 2:nChain){
                    nam <- paste("chain", i, sep = "")
                    kappa <- cbind(kappa, x[[nam]]$kappa1.p)
                }
                psrfKappa1 <- calcPSR(kappa)
                
                kappa <- x$chain1$kappa2.p
                for(i in 2:nChain){
                    nam <- paste("chain", i, sep = "")
                    kappa <- cbind(kappa, x[[nam]]$kappa2.p)
                }
                psrfKappa2 <- calcPSR(kappa)
                
                kappa <- x$chain1$kappa3.p
                for(i in 2:nChain){
                    nam <- paste("chain", i, sep = "")
                    kappa <- cbind(kappa, x[[nam]]$kappa3.p)
                }
                psrfKappa3 <- calcPSR(kappa)
                
                bh <- matrix(c(psrfKappa1, psrfKappa2, psrfKappa3, psrfAlpha1, psrfAlpha2, psrfAlpha3), 2, 3, byrow = T)
                dimnames(bh) <- list(c("kappa", "alpha"), c("h1", "h2", "h3"))
                
                psrf <- list(theta=psrftheta, coef=psrfcoef, h0=bh)
            }
            
            if(class(x)[4] == "PEM")
            {
                ##
                
                ntime1  = length(x$chain1$time_lambda1)
                ntime2  = length(x$chain1$time_lambda2)
                ntime3  = length(x$chain1$time_lambda3)
                
                # lambda's
                
                psrfLam1 <- rep(NA, ntime1)
                
                for(j in 1:ntime1){
                    
                    lambda1 <- x$chain1$lambda1.fin[,j]
                    for(i in 2:nChain){
                        nam <- paste("chain", i, sep = "")
                        lambda1 <- cbind(lambda1, x[[nam]]$lambda1.fin[,j])
                    }
                    psrfLam1[j] <- calcPSR(lambda1)
                }
                
                psrfLam2 <- rep(NA, ntime2)
                
                for(j in 1:ntime2){
                    
                    lambda2 <- x$chain1$lambda2.fin[,j]
                    for(i in 2:nChain){
                        nam <- paste("chain", i, sep = "")
                        lambda2 <- cbind(lambda2, x[[nam]]$lambda2.fin[,j])
                    }
                    psrfLam2[j] <- calcPSR(lambda2)
                }
                
                psrfLam3 <- rep(NA, ntime3)
                
                for(j in 1:ntime3){
                    
                    lambda3 <- x$chain1$lambda3.fin[,j]
                    for(i in 2:nChain){
                        nam <- paste("chain", i, sep = "")
                        lambda3 <- cbind(lambda3, x[[nam]]$lambda3.fin[,j])
                    }
                    psrfLam3[j] <- calcPSR(lambda3)
                }

                # mu_lam
                
                mu <- x$chain1$mu_lam1.p
                for(i in 2:nChain){
                    nam <- paste("chain", i, sep = "")
                    mu <- cbind(mu, x[[nam]]$mu_lam1.p)
                }
                psrfMu1 <- calcPSR(mu)
                
                mu <- x$chain1$mu_lam2.p
                for(i in 2:nChain){
                    nam <- paste("chain", i, sep = "")
                    mu <- cbind(mu, x[[nam]]$mu_lam2.p)
                }
                psrfMu2 <- calcPSR(mu)
                
                mu <- x$chain1$mu_lam3.p
                for(i in 2:nChain){
                    nam <- paste("chain", i, sep = "")
                    mu <- cbind(mu, x[[nam]]$mu_lam3.p)
                }
                psrfMu3 <- calcPSR(mu)
                
                # sigSq_lam
                
                sig <- x$chain1$sigSq_lam1.p
                for(i in 2:nChain){
                    nam <- paste("chain", i, sep = "")
                    sig <- cbind(sig, x[[nam]]$sigSq_lam1.p)
                }
                psrfSig1 <- calcPSR(sig)
                
                
                sig <- x$chain1$sigSq_lam2.p
                for(i in 2:nChain){
                    nam <- paste("chain", i, sep = "")
                    sig <- cbind(sig, x[[nam]]$sigSq_lam2.p)
                }
                psrfSig2 <- calcPSR(sig)
                
                
                sig <- x$chain1$sigSq_lam3.p
                for(i in 2:nChain){
                    nam <- paste("chain", i, sep = "")
                    sig <- cbind(sig, x[[nam]]$sigSq_lam3.p)
                }
                psrfSig3 <- calcPSR(sig)
                
                # J
                
                J <- x$chain1$K1.p
                for(i in 2:nChain){
                    nam <- paste("chain", i, sep = "")
                    J <- cbind(J, x[[nam]]$K1.p)
                }
                psrfJ1 <- calcPSR(J)
                
                J <- x$chain1$K2.p
                for(i in 2:nChain){
                    nam <- paste("chain", i, sep = "")
                    J <- cbind(J, x[[nam]]$K2.p)
                }
                psrfJ2 <- calcPSR(J)
                
                J <- x$chain1$K3.p
                for(i in 2:nChain){
                    nam <- paste("chain", i, sep = "")
                    J <- cbind(J, x[[nam]]$K3.p)
                }
                psrfJ3 <- calcPSR(J)
                
                
                bh <- matrix(c(psrfMu1, psrfMu2, psrfMu3, psrfSig1, psrfSig2, psrfSig3, psrfJ1, psrfJ2, psrfJ3), 3, 3, byrow = T)
                dimnames(bh) <- list(c("mu", "sigmaSq", "K"), c("h1", "h2", "h3"))
                
                psrf <- list(theta=psrftheta, coef=psrfcoef, h0=bh, lambda1=psrfLam1, lambda2=psrfLam2, lambda3=psrfLam3)
                
            }
            

            
        }
        
        if(class(x)[2] == "Surv")
        {
            beta.names <- c(x$chain1$covNames)
            nP         <- length(beta.names)
            output <- matrix(NA, nrow=nP, ncol=1)
            dimnames(output) <- list(beta.names, c("beta"))
            
            if(length(x$chain1$beta.p) != 0){
                
                #beta
                
                p	= dim(x$chain1$beta.p)[2]
                
                psrfBeta <- rep(NA, p)
                for(j in 1:p){
                    
                    beta <- x$chain1$beta[,j]
                    for(i in 2:nChain){
                        nam <- paste("chain", i, sep = "")
                        beta <- cbind(beta, x[[nam]]$beta[,j])
                    }
                    psrfBeta[j] <- calcPSR(beta)
                }
                
                for(i in 1:nP)
                {
                    for(k in 1:p) if(x$chain1$covNames[k] == beta.names[i]) output[i,1] <- psrfBeta[k]
                }
                
            }
            
            psrfcoef <- NULL
            if(nP > 0)
            {
                psrfcoef <- output
            }
            
            if(class(x)[4] == "WB")
            {
                ##
                # alpha
                
                alpha <- x$chain1$alpha.p
                for(i in 2:nChain){
                    nam <- paste("chain", i, sep = "")
                    alpha <- cbind(alpha, x[[nam]]$alpha.p)
                }
                psrfAlpha <- calcPSR(alpha)
                
                
                # kappa
                
                kappa <- x$chain1$kappa.p
                for(i in 2:nChain){
                    nam <- paste("chain", i, sep = "")
                    kappa <- cbind(kappa, x[[nam]]$kappa.p)
                }
                psrfKappa <- calcPSR(kappa)
                
                
                bh <- matrix(c(psrfKappa, psrfAlpha), 2, 1, byrow = T)
                dimnames(bh) <- list(c("kappa", "alpha"), c("h"))
                
                psrf <- list(coef=psrfcoef, h0=bh)
                
            }
            
            if(class(x)[4] == "PEM")
            {
                ##
                ntime  = length(x$chain1$time_lambda)
                
                # lambda
                
                psrfLam <- rep(NA, ntime)
                
                for(j in 1:ntime){
                    
                    namPara = paste("beta_", j, sep = "")
                    
                    lambda <- x$chain1$lambda.fin[,j]
                    for(i in 2:nChain){
                        nam <- paste("chain", i, sep = "")
                        lambda <- cbind(lambda, x[[nam]]$lambda.fin[,j])
                    }
                    psrfLam[j] <- calcPSR(lambda)
                }
                
                # mu_lam
                
                mu <- x$chain1$mu_lam.p
                for(i in 2:nChain){
                    nam <- paste("chain", i, sep = "")
                    mu <- cbind(mu, x[[nam]]$mu_lam.p)
                }
                psrfMu <- calcPSR(mu)
                
                # sigSq_lam
                
                sig <- x$chain1$sigSq_lam.p
                for(i in 2:nChain){
                    nam <- paste("chain", i, sep = "")
                    sig <- cbind(sig, x[[nam]]$sigSq_lam.p)
                }
                psrfSig <- calcPSR(sig)
                
                
                # J
                
                J <- x$chain1$K.p
                for(i in 2:nChain){
                    nam <- paste("chain", i, sep = "")
                    J <- cbind(J, x[[nam]]$K.p)
                }
                psrfJ <- calcPSR(J)
                
                
                bh <- matrix(c(psrfMu, psrfSig, psrfJ), 3, 1, byrow = T)
                dimnames(bh) <- list(c("mu", "sigmaSq", "K"), c("h"))
                
                psrf <- list(coef=psrfcoef, h0=bh, lambda=psrfLam)
            }

            
        }
    }
    
    # estimates
    
    if(class(x)[2] == "ID")
    {
        ##
        theta.p <- x$chain1$theta.p
        
        if(nChain > 1){
            for(i in 2:nChain){
                nam <- paste("chain", i, sep="")
                theta.p <- rbind(theta.p, x[[nam]]$theta.p)
            }
        }
        
        theta.pMed <- apply(theta.p, 2, median)
        theta.pUb <- apply(theta.p, 2, quantile, prob = 0.975)
        theta.pLb <- apply(theta.p, 2, quantile, prob = 0.025)
        
        tbl_theta <- matrix(NA, 1, 3)
        dimnames(tbl_theta) <- list("", c( "theta", "LL", "UL"))
        
        tbl_theta[,1]	<- theta.pMed
        tbl_theta[,2]	<- theta.pLb
        tbl_theta[,3]	<- theta.pUb

        ##
        beta.names <- unique(c(x$chain1$covNames1, x$chain1$covNames2, x$chain1$covNames3))
        nP         <- length(beta.names)
        output <- matrix(NA, nrow=nP, ncol=9)
        dimnames(output) <- list(beta.names, c("exp(beta1)", "LL", "UL", "exp(beta2)", "LL", "UL", "exp(beta3)", "LL", "UL"))
        
        if(length(x$chain1$beta1.p) != 0){
            
            #beta1
            p1	= dim(x$chain1$beta1.p)[2]
            beta.p <- x$chain1$beta1.p
            
            if(nChain > 1){
                for(i in 2:nChain){
                    nam <- paste("chain", i, sep="")
                    beta.p <- rbind(beta.p, x[[nam]]$beta1.p)
                }
            }
            
            beta.pMed <- apply(exp(beta.p), 2, median)
            beta.pSd <- apply(exp(beta.p), 2, sd)
            beta.pUb <- apply(exp(beta.p), 2, quantile, prob = 0.975)
            beta.pLb <- apply(exp(beta.p), 2, quantile, prob = 0.025)
            
            tbl1 <- matrix(NA, p1, 3)
            rownames(tbl1) <- x$chain1$covNames1
            
            tbl1[,1]	<- beta.pMed
            tbl1[,2]	<- beta.pLb
            tbl1[,3]	<- beta.pUb
            
            for(i in 1:nP)
            {
                for(k in 1:p1) if(x$chain1$covNames1[k] == beta.names[i]) output[i,1:3] <- tbl1[k,]
            }
            
        }
        if(length(x$chain1$beta2.p) != 0){
            
            #beta2
            p2	= dim(x$chain1$beta2.p)[2]
            beta.p <- x$chain1$beta2.p
            
            if(nChain > 1){
                for(i in 2:nChain){
                    nam <- paste("chain", i, sep="")
                    beta.p <- rbind(beta.p, x[[nam]]$beta2.p)
                }
            }
            
            beta.pMed <- apply(exp(beta.p), 2, median)
            beta.pSd <- apply(exp(beta.p), 2, sd)
            beta.pUb <- apply(exp(beta.p), 2, quantile, prob = 0.975)
            beta.pLb <- apply(exp(beta.p), 2, quantile, prob = 0.025)
            
            tbl2 <- matrix(NA, p2, 3)
            rownames(tbl2) <- x$chain1$covNames2
            
            tbl2[,1]	<- beta.pMed
            tbl2[,2]	<- beta.pLb
            tbl2[,3]	<- beta.pUb
            
            for(i in 1:nP)
            {
                for(k in 1:p2) if(x$chain1$covNames2[k] == beta.names[i]) output[i,4:6] <- tbl2[k,]
            }
            
        }
        if(length(x$chain1$beta3.p) != 0){
            
            #beta3
            p3	= dim(x$chain1$beta3.p)[2]
            beta.p <- x$chain1$beta3.p
            
            if(nChain > 1){
                for(i in 2:nChain){
                    nam <- paste("chain", i, sep="")
                    beta.p <- rbind(beta.p, x[[nam]]$beta3.p)
                }
            }
            
            beta.pMed <- apply(exp(beta.p), 2, median)
            beta.pSd <- apply(exp(beta.p), 2, sd)
            beta.pUb <- apply(exp(beta.p), 2, quantile, prob = 0.975)
            beta.pLb <- apply(exp(beta.p), 2, quantile, prob = 0.025)
            
            tbl3 <- matrix(NA, p3, 3)
            rownames(tbl3) <- x$chain1$covNames3
            
            tbl3[,1]	<- beta.pMed
            tbl3[,2]	<- beta.pLb
            tbl3[,3]	<- beta.pUb
            
            for(i in 1:nP)
            {
                for(k in 1:p3) if(x$chain1$covNames3[k] == beta.names[i]) output[i,7:9] <- tbl3[k,]
            }
        }
        
        output.coef <- NULL
        if(nP > 0)
        {
            output.coef <- output
        }
        
        
        
        if(class(x)[4] == "WB")
        {
            ##
            alpha.p <- x$chain1$alpha1.p
            
            if(nChain > 1){
                for(i in 2:nChain){
                    nam <- paste("chain", i, sep="")
                    alpha.p <- rbind(alpha.p, x[[nam]]$alpha1.p)
                }
            }

            alpha.pMed <- apply(log(alpha.p), 2, median)
            alpha.pUb <- apply(log(alpha.p), 2, quantile, prob = 0.975)
            alpha.pLb <- apply(log(alpha.p), 2, quantile, prob = 0.025)
            
            tbl_a1 <- c(alpha.pMed,alpha.pLb, alpha.pUb)
            
            ##
            alpha.p <- x$chain1$alpha2.p
            
            if(nChain > 1){
                for(i in 2:nChain){
                    nam <- paste("chain", i, sep="")
                    alpha.p <- rbind(alpha.p, x[[nam]]$alpha2.p)
                }
            }
            
            alpha.pMed <- apply(log(alpha.p), 2, median)
            alpha.pUb <- apply(log(alpha.p), 2, quantile, prob = 0.975)
            alpha.pLb <- apply(log(alpha.p), 2, quantile, prob = 0.025)
            
            tbl_a2 <- c(alpha.pMed,alpha.pLb, alpha.pUb)
            
            
            ##
            alpha.p <- x$chain1$alpha3.p
            
            if(nChain > 1){
                for(i in 2:nChain){
                    nam <- paste("chain", i, sep="")
                    alpha.p <- rbind(alpha.p, x[[nam]]$alpha3.p)
                }
            }
            
            alpha.pMed <- apply(log(alpha.p), 2, median)
            alpha.pUb <- apply(log(alpha.p), 2, quantile, prob = 0.975)
            alpha.pLb <- apply(log(alpha.p), 2, quantile, prob = 0.025)
            
            tbl_a3 <- c(alpha.pMed,alpha.pLb, alpha.pUb)
            
            ##
            kappa.p <- x$chain1$kappa1.p
            
            if(nChain > 1){
                for(i in 2:nChain){
                    nam <- paste("chain", i, sep="")
                    kappa.p <- rbind(kappa.p, x[[nam]]$kappa1.p)
                }
            }
            
            kappa.pMed <- apply(log(kappa.p), 2, median)
            kappa.pUb <- apply(log(kappa.p), 2, quantile, prob = 0.975)
            kappa.pLb <- apply(log(kappa.p), 2, quantile, prob = 0.025)
            
            tbl_k1 <- c(kappa.pMed, kappa.pLb, kappa.pUb)
            
            ##
            kappa.p <- x$chain1$kappa2.p
            
            if(nChain > 1){
                for(i in 2:nChain){
                    nam <- paste("chain", i, sep="")
                    kappa.p <- rbind(kappa.p, x[[nam]]$kappa2.p)
                }							
            }	
            
            kappa.pMed <- apply(log(kappa.p), 2, median)
            kappa.pUb <- apply(log(kappa.p), 2, quantile, prob = 0.975)
            kappa.pLb <- apply(log(kappa.p), 2, quantile, prob = 0.025)
            
            tbl_k2 <- c(kappa.pMed, kappa.pLb, kappa.pUb)

            ##
            kappa.p <- x$chain1$kappa3.p
            
            if(nChain > 1){
                for(i in 2:nChain){
                    nam <- paste("chain", i, sep="")
                    kappa.p <- rbind(kappa.p, x[[nam]]$kappa3.p)
                }							
            }	
            
            kappa.pMed <- apply(log(kappa.p), 2, median)
            kappa.pUb <- apply(log(kappa.p), 2, quantile, prob = 0.975)
            kappa.pLb <- apply(log(kappa.p), 2, quantile, prob = 0.025)
            
            tbl_k3 <- c(kappa.pMed, kappa.pLb, kappa.pUb)

            bh  <- matrix(c(tbl_a1, tbl_a2, tbl_a3, tbl_k1, tbl_k2, tbl_k3), 2, 9, byrow = T)
            dimnames(bh) <- list(c("Weibull: log-kappa", "Weibull: log-alpha"), c("h1-PM", "LL", "UL", "h2-PM", "LL", "UL", "h3-PM", "LL", "UL"))
            
            value <- list(classFit=class(x), psrf=psrf, theta=tbl_theta, coef=output.coef, h0=bh)
            
        }
        if(class(x)[4] == "PEM")
        {
            ##
            mu_lam.p <- x$chain1$mu_lam1.p
            
            if(nChain > 1){
                for(i in 2:nChain){
                    nam <- paste("chain", i, sep="")
                    mu_lam.p <- rbind(mu_lam.p, x[[nam]]$mu_lam1.p)
                }
            }

            mu_lam.pMed <- apply(mu_lam.p, 2, median)
            mu_lam.pUb <- apply(mu_lam.p, 2, quantile, prob = 0.975)
            mu_lam.pLb <- apply(mu_lam.p, 2, quantile, prob = 0.025)
            
            tbl_m1 <- c(mu_lam.pMed, mu_lam.pLb, mu_lam.pUb)
     
            ##
            mu_lam.p <- x$chain1$mu_lam2.p
            
            if(nChain > 1){
                for(i in 2:nChain){
                    nam <- paste("chain", i, sep="")
                    mu_lam.p <- rbind(mu_lam.p, x[[nam]]$mu_lam2.p)
                }
            }
            
            mu_lam.pMed <- apply(mu_lam.p, 2, median)
            mu_lam.pUb <- apply(mu_lam.p, 2, quantile, prob = 0.975)
            mu_lam.pLb <- apply(mu_lam.p, 2, quantile, prob = 0.025)
            
            tbl_m2 <- c(mu_lam.pMed, mu_lam.pLb, mu_lam.pUb)
            
            ##
            mu_lam.p <- x$chain1$mu_lam3.p
            
            if(nChain > 1){
                for(i in 2:nChain){
                    nam <- paste("chain", i, sep="")
                    mu_lam.p <- rbind(mu_lam.p, x[[nam]]$mu_lam3.p)
                }
            }
            
            mu_lam.pMed <- apply(mu_lam.p, 2, median)
            mu_lam.pUb <- apply(mu_lam.p, 2, quantile, prob = 0.975)
            mu_lam.pLb <- apply(mu_lam.p, 2, quantile, prob = 0.025)
            
            tbl_m3 <- c(mu_lam.pMed, mu_lam.pLb, mu_lam.pUb)
            
            
            ##
            sigSq_lam.p <- x$chain1$sigSq_lam1.p
            
            if(nChain > 1){
                for(i in 2:nChain){
                    nam <- paste("chain", i, sep="")
                    sigSq_lam.p <- rbind(sigSq_lam.p, x[[nam]]$sigSq_lam1.p)
                }
            }
            
            sigSq_lam.pMed <- apply(sigSq_lam.p, 2, median)
            sigSq_lam.pUb <- apply(sigSq_lam.p, 2, quantile, prob = 0.975)
            sigSq_lam.pLb <- apply(sigSq_lam.p, 2, quantile, prob = 0.025)
            
            tbl_s1 <- c(sigSq_lam.pMed, sigSq_lam.pLb, sigSq_lam.pUb)
            
            ##
            sigSq_lam.p <- x$chain1$sigSq_lam2.p
            
            if(nChain > 1){
                for(i in 2:nChain){
                    nam <- paste("chain", i, sep="")
                    sigSq_lam.p <- rbind(sigSq_lam.p, x[[nam]]$sigSq_lam2.p)
                }
            }
            
            sigSq_lam.pMed <- apply(sigSq_lam.p, 2, median)
            sigSq_lam.pUb <- apply(sigSq_lam.p, 2, quantile, prob = 0.975)
            sigSq_lam.pLb <- apply(sigSq_lam.p, 2, quantile, prob = 0.025)
            
            tbl_s2 <- c(sigSq_lam.pMed, sigSq_lam.pLb, sigSq_lam.pUb)
            
            ##
            sigSq_lam.p <- x$chain1$sigSq_lam3.p
            
            if(nChain > 1){
                for(i in 2:nChain){
                    nam <- paste("chain", i, sep="")
                    sigSq_lam.p <- rbind(sigSq_lam.p, x[[nam]]$sigSq_lam3.p)
                }
            }
            
            sigSq_lam.pMed <- apply(sigSq_lam.p, 2, median)
            sigSq_lam.pUb <- apply(sigSq_lam.p, 2, quantile, prob = 0.975)
            sigSq_lam.pLb <- apply(sigSq_lam.p, 2, quantile, prob = 0.025)
            
            tbl_s3 <- c(sigSq_lam.pMed, sigSq_lam.pLb, sigSq_lam.pUb)
            
            ##
            J.p <- x$chain1$K1.p
            
            if(nChain > 1){
                for(i in 2:nChain){
                    nam <- paste("chain", i, sep="")
                    J.p <- rbind(J.p, x[[nam]]$K1.p)											
                }							
            }	

            J.pMed <- apply(J.p, 2, median)
            J.pUb <- apply(J.p, 2, quantile, prob = 0.975)		
            J.pLb <- apply(J.p, 2, quantile, prob = 0.025)
            
            tbl_j1 <- c(J.pMed, J.pLb, J.pUb)
            
            ##
            J.p <- x$chain1$K2.p
            
            if(nChain > 1){
                for(i in 2:nChain){
                    nam <- paste("chain", i, sep="")
                    J.p <- rbind(J.p, x[[nam]]$K2.p)											
                }							
            }	
            
            J.pMed <- apply(J.p, 2, median)
            J.pUb <- apply(J.p, 2, quantile, prob = 0.975)		
            J.pLb <- apply(J.p, 2, quantile, prob = 0.025)
            
            tbl_j2 <- c(J.pMed, J.pLb, J.pUb)
            
            ##
            J.p <- x$chain1$K3.p
            
            if(nChain > 1){
                for(i in 2:nChain){
                    nam <- paste("chain", i, sep="")
                    J.p <- rbind(J.p, x[[nam]]$K3.p)											
                }							
            }
            
            J.pMed <- apply(J.p, 2, median)
            J.pUb <- apply(J.p, 2, quantile, prob = 0.975)		
            J.pLb <- apply(J.p, 2, quantile, prob = 0.025)
            
            tbl_j3 <- c(J.pMed, J.pLb, J.pUb)
            
            bh <- matrix(c(tbl_m1, tbl_m2, tbl_m3, tbl_s1, tbl_s2, tbl_s3, tbl_j1, tbl_j2, tbl_j3), 3, 9, byrow = T)
            dimnames(bh) <- list(c("mu", "sigmaSq", "K"), c("h1-PM", "LL", "UL", "h2-PM", "LL", "UL", "h3-PM", "LL", "UL"))
            
            
            ##
            lambda.p <- x$chain1$lambda1.fin
            
            if(nChain > 1){
                for(i in 2:nChain){
                    nam <- paste("chain", i, sep="")
                    lambda.p <- rbind(lambda.p, x[[nam]]$lambda1.fin)
                }
            }
            
            lambda.pMed <- apply(lambda.p, 2, median)
            lambda.pUb <- apply(lambda.p, 2, quantile, prob = 0.975)
            lambda.pLb <- apply(lambda.p, 2, quantile, prob = 0.025)
            
            lambda1 <- cbind(x$chain1$time_lambda1, lambda.pMed, lambda.pLb, lambda.pUb)
            
            dimnames(lambda1) <- list(rep("", length(x$chain1$time_lambda1)), c("time", "lambda1-PM", "LL", "UL"))
            
            ##
            lambda.p <- x$chain1$lambda2.fin
            
            if(nChain > 1){
                for(i in 2:nChain){
                    nam <- paste("chain", i, sep="")
                    lambda.p <- rbind(lambda.p, x[[nam]]$lambda2.fin)
                }
            }
            
            lambda.pMed <- apply(lambda.p, 2, median)
            lambda.pUb <- apply(lambda.p, 2, quantile, prob = 0.975)
            lambda.pLb <- apply(lambda.p, 2, quantile, prob = 0.025)
            
            lambda2 <- cbind(x$chain1$time_lambda2, lambda.pMed, lambda.pLb, lambda.pUb)
            
            dimnames(lambda2) <- list(rep("", length(x$chain1$time_lambda2)), c("time", "lambda2-PM", "LL", "UL"))
            
            ##
            lambda.p <- x$chain1$lambda3.fin
            
            if(nChain > 1){
                for(i in 2:nChain){
                    nam <- paste("chain", i, sep="")
                    lambda.p <- rbind(lambda.p, x[[nam]]$lambda3.fin)
                }
            }
            
            lambda.pMed <- apply(lambda.p, 2, median)
            lambda.pUb <- apply(lambda.p, 2, quantile, prob = 0.975)
            lambda.pLb <- apply(lambda.p, 2, quantile, prob = 0.025)
            
            lambda3 <- cbind(x$chain1$time_lambda3, lambda.pMed, lambda.pLb, lambda.pUb)
            
            dimnames(lambda3) <- list(rep("", length(x$chain1$time_lambda3)), c("time", "lambda3-PM", "LL", "UL"))
            
            
            value <- list(classFit=class(x), psrf=psrf, theta=tbl_theta, coef=output.coef, h0=bh, lambda1=lambda1, lambda2=lambda2, lambda3=lambda3)
            
        }
        
        if(class(x)[3] == "Cor")
        {
            if(class(x)[5] == "MVN")
            {
                nS <- dim(x$chain1$Sigma_V.p)[3]
                Sigma <- array(NA, c(3,3, nS*nChain))
                Sigma[,,1:nS] <- x$chain1$Sigma_V.p
                
                if(nChain > 1){
                    for(i in 2:nChain){
                        nam <- paste("chain", i, sep="")
                        Sigma[,,(nS*(i-1)+1):(nS*i)] <- x[[nam]]$Sigma_V.p
                    }
                }
            }
            
            if(class(x)[5] == "DPM")
            {
                ##
                tau.p <- x$chain1$tau.p
                
                if(nChain > 1){
                    for(i in 2:nChain){
                        nam <- paste("chain", i, sep="")
                        tau.p <- rbind(tau.p, x[[nam]]$tau.p)
                    }
                }
                
                tau.pMed <- apply(tau.p, 2, median)
                tau.pUb <- apply(tau.p, 2, quantile, prob = 0.975)
                tau.pLb <- apply(tau.p, 2, quantile, prob = 0.025)
                
                tbl_tau <- matrix(NA, 1, 3)
                dimnames(tbl_tau) <- list("", c( "tau", "LL", "UL"))
                
                tbl_tau[,1]	<- tau.pMed
                tbl_tau[,2]	<- tau.pLb
                tbl_tau[,3]	<- tau.pUb
                
                
                nS <- dim(x$chain1$Sigma.p)[3]
                Sigma <- array(NA, c(3,3, nS*nChain))
                Sigma[,,1:nS] <- calVar_DPM_MVN(x$chain1)
                
                if(nChain > 1){
                    for(i in 2:nChain){
                        nam <- paste("chain", i, sep="")
                        Sigma[,,(nS*(i-1)+1):(nS*i)] <- calVar_DPM_MVN(x[[nam]])
                    }
                }
                value$tau <- tbl_tau
            }
            Sigma.Med <- apply(Sigma, c(1,2), median)
            Sigma.Sd <- apply(Sigma, c(1,2), sd)
            Sigma.Ub <- apply(Sigma, c(1,2), quantile, prob = 0.975)
            Sigma.Lb <- apply(Sigma, c(1,2), quantile, prob = 0.025)
            
            dimnames(Sigma.Med) <- list(c("", "", ""), c("Sigma_V-PM", "", ""))
            dimnames(Sigma.Sd) <- list(c("", "", ""), c("Sigma_V-SD", "", ""))
            dimnames(Sigma.Lb) <- list(c("", "", ""), c("Sigma_V-LL", "", ""))
            dimnames(Sigma.Ub) <- list(c("", "", ""), c("Sigma_V-UL", "", ""))
            
            value$Sigma.PM <- Sigma.Med
            value$Sigma.SD <- Sigma.Sd
            value$Sigma.UL <- Sigma.Ub
            value$Sigma.LL <- Sigma.Lb
        }
    }
    
    if(class(x)[2] == "Surv")
    {
        ##
        beta.names <- c(x$chain1$covNames)
        nP         <- length(beta.names)
        output <- matrix(NA, nrow=nP, ncol=3)
        dimnames(output) <- list(beta.names, c("exp(beta)", "LL", "UL"))
        
        
        if(length(x$chain1$beta.p) != 0){
            
            #beta
            p	= dim(x$chain1$beta.p)[2]
            beta.p <- x$chain1$beta.p
            
            if(nChain > 1){
                for(i in 2:nChain){
                    nam <- paste("chain", i, sep="")
                    beta.p <- rbind(beta.p, x[[nam]]$beta.p)
                }
            }
            
            beta.pMed <- apply(exp(beta.p), 2, median)
            beta.pSd <- apply(exp(beta.p), 2, sd)
            beta.pUb <- apply(exp(beta.p), 2, quantile, prob = 0.975)
            beta.pLb <- apply(exp(beta.p), 2, quantile, prob = 0.025)
            
            tbl <- matrix(NA, p, 3)
            rownames(tbl) <- x$chain1$covNames
            
            tbl[,1]	<- beta.pMed
            tbl[,2]	<- beta.pLb
            tbl[,3]	<- beta.pUb
            
            
            for(i in 1:nP)
            {
                for(k in 1:p) if(x$chain1$covNames[k] == beta.names[i]) output[i,1:3] <- tbl[k,]
            }
        }
        
        output.coef <- NULL
        if(nP > 0)
        {
            output.coef <- output
        }
        
        if(class(x)[4] == "WB")
        {
            ##
            alpha.p <- x$chain1$alpha.p
            
            if(nChain > 1){
                for(i in 2:nChain){
                    nam <- paste("chain", i, sep="")
                    alpha.p <- rbind(alpha.p, x[[nam]]$alpha.p)
                }
            }
            
            alpha.pMed <- apply(log(alpha.p), 2, median)
            alpha.pUb <- apply(log(alpha.p), 2, quantile, prob = 0.975)
            alpha.pLb <- apply(log(alpha.p), 2, quantile, prob = 0.025)
            
            tbl_a <- c(alpha.pMed,alpha.pLb, alpha.pUb)
            
            
            ##
            kappa.p <- x$chain1$kappa.p
            
            if(nChain > 1){
                for(i in 2:nChain){
                    nam <- paste("chain", i, sep="")
                    kappa.p <- rbind(kappa.p, x[[nam]]$kappa.p)
                }
            }
            
            kappa.pMed <- apply(log(kappa.p), 2, median)
            kappa.pUb <- apply(log(kappa.p), 2, quantile, prob = 0.975)
            kappa.pLb <- apply(log(kappa.p), 2, quantile, prob = 0.025)
            
            tbl_k <- c(kappa.pMed, kappa.pLb, kappa.pUb)
            
            bh  <- matrix(c(tbl_a, tbl_k), 2, 3, byrow = T)
            dimnames(bh) <- list(c("Weibull: log-kappa", "Weibull: log-alpha"), c("h-PM", "LL", "UL"))
            
            value <- list(coef=output.coef, h0=bh, psrf=psrf, classFit=class(x))
            
        }
        if(class(x)[4] == "PEM")
        {
            ##
            mu_lam.p <- x$chain1$mu_lam.p
            
            if(nChain > 1){
                for(i in 2:nChain){
                    nam <- paste("chain", i, sep="")
                    mu_lam.p <- rbind(mu_lam.p, x[[nam]]$mu_lam.p)
                }
            }
            
            mu_lam.pMed <- apply(mu_lam.p, 2, median)
            mu_lam.pUb <- apply(mu_lam.p, 2, quantile, prob = 0.975)
            mu_lam.pLb <- apply(mu_lam.p, 2, quantile, prob = 0.025)
            
            tbl_m <- c(mu_lam.pMed, mu_lam.pLb, mu_lam.pUb)
            
            ##
            sigSq_lam.p <- x$chain1$sigSq_lam.p
            
            if(nChain > 1){
                for(i in 2:nChain){
                    nam <- paste("chain", i, sep="")
                    sigSq_lam.p <- rbind(sigSq_lam.p, x[[nam]]$sigSq_lam.p)
                }
            }
            
            sigSq_lam.pMed <- apply(sigSq_lam.p, 2, median)
            sigSq_lam.pUb <- apply(sigSq_lam.p, 2, quantile, prob = 0.975)
            sigSq_lam.pLb <- apply(sigSq_lam.p, 2, quantile, prob = 0.025)
            
            tbl_s <- c(sigSq_lam.pMed, sigSq_lam.pLb, sigSq_lam.pUb)
            
            ##
            J.p <- x$chain1$K.p
            
            if(nChain > 1){
                for(i in 2:nChain){
                    nam <- paste("chain", i, sep="")
                    J.p <- rbind(J.p, x[[nam]]$K.p)
                }
            }
            
            J.pMed <- apply(J.p, 2, median)
            J.pUb <- apply(J.p, 2, quantile, prob = 0.975)
            J.pLb <- apply(J.p, 2, quantile, prob = 0.025)
            
            tbl_j <- c(J.pMed, J.pLb, J.pUb)
            
            bh <- matrix(c(tbl_m, tbl_s, tbl_j), 3, 3, byrow = T)
            dimnames(bh) <- list(c("mu", "sigmaSq", "K"), c("h-PM", "LL", "UL"))
            
            
            ##
            lambda.p <- x$chain1$lambda.fin
            
            if(nChain > 1){
                for(i in 2:nChain){
                    nam <- paste("chain", i, sep="")
                    lambda.p <- rbind(lambda.p, x[[nam]]$lambda.fin)
                }
            }
            
            lambda.pMed <- apply(lambda.p, 2, median)
            lambda.pUb <- apply(lambda.p, 2, quantile, prob = 0.975)
            lambda.pLb <- apply(lambda.p, 2, quantile, prob = 0.025)
            
            lambda <- cbind(x$chain1$time_lambda, lambda.pMed, lambda.pLb, lambda.pUb)
            
            dimnames(lambda) <- list(rep("", length(x$chain1$time_lambda)), c("time", "lambda-PM", "LL", "UL"))
            
            value <- list(coef=output.coef, h0=bh, psrf=psrf, lambda=lambda, classFit=class(x))
            
        }
        
        if(class(x)[3] == "Cor")
        {
            if(class(x)[5] == "Normal")
            {
                #sigmaV
                sigV <- 1/x$chain1$zeta.p
                
                if(nChain > 1){
                    for(i in 2:nChain){
                        nam <- paste("chain", i, sep="")
                        sigV <- rbind(sigV, 1/x[[nam]]$zeta.p)
                    }
                }
            }
            
            if(class(x)[5] == "DPM")
            {
                ##
                tau.p <- x$chain1$tau.p
                
                if(nChain > 1){
                    for(i in 2:nChain){
                        nam <- paste("chain", i, sep="")
                        tau.p <- rbind(tau.p, x[[nam]]$tau.p)
                    }
                }
                
                tau.pMed <- apply(tau.p, 2, median)
                tau.pUb <- apply(tau.p, 2, quantile, prob = 0.975)
                tau.pLb <- apply(tau.p, 2, quantile, prob = 0.025)
                
                tbl_tau <- matrix(NA, 1, 3)
                dimnames(tbl_tau) <- list("", c( "tau", "LL", "UL"))
                
                tbl_tau[,1]	<- tau.pMed
                tbl_tau[,2]	<- tau.pLb
                tbl_tau[,3]	<- tau.pUb
                
                
                nS <- dim(x$chain1$zeta.p)[1]
                sigV <- rep(NA, nS*nChain)
                sigV[1:nS] <- calVar_DPM_Normal(x$chain1)
                
                if(nChain > 1){
                    for(i in 2:nChain){
                        nam <- paste("chain", i, sep="")
                        sigV[(nS*(i-1)+1):(nS*i)] <- calVar_DPM_Normal(x[[nam]])
                    }
                }
                
                value$tau <- tbl_tau
            }
            
            sigVMed <- median(sigV)
            sigVSd <- sd(sigV)
            sigVUb <- quantile(sigV, prob = 0.975)
            sigVLb <- quantile(sigV, prob = 0.025)
            
            tbl_sigV <- matrix(NA, nrow=1, ncol=3)
            tbl_sigV[,1]	<- sigVMed
            tbl_sigV[,2]	<- sigVLb
            tbl_sigV[,3]	<- sigVUb
            dimnames(tbl_sigV) <- list("", c("sigma_V-PM", "LL", "UL"))
            
            value$sigma_V   <- tbl_sigV
        }


    }
    
    value$setup <- x$setup

# if(class(x)[3] == "Cor")
# {
#     class(value) <- c("summ.Bayes", as.vector(class(x)[2]), "Cor", as.vector(class(x)[4]), as.vector(class(x)[5]))
# }
# if(class(x)[3] == "Ind")
# {
#     class(value) <- c("summ.Bayes", as.vector(class(x)[2]), "Ind", as.vector(class(x)[4]))
# }

    class(value) <- "summ.Bayes"
 
 return(value)
 
}



####
## PRINT.SUMMARY METHOD
####
##
print.summ.Freq <- function(x, digits=3, ...)
{
	obj <- x
  ##
  if(class(obj)[2] == "Surv")
  {
      ##
      cat("\nAnalysis of independent univariate time-to-event data \n")
  }
  if(class(obj)[2] == "ID")
  {
      ##
      cat("\nAnalysis of independent semi-competing risks data \n")
      cat(class(obj)[3], "assumption for h3\n")
  }
  ##
  #cat("\nRegression coefficients:\n")
  #print(round(obj$coef, digits=digits))
  ##
  cat("\nHazard ratios:\n")
  print(round(exp(obj$coef), digits=digits))
  ##
  if(class(obj)[2] == "ID"){
    cat("\nVariance of frailties:\n")
    print(round(obj$theta, digits=digits))
  }
  ##
  cat("\nBaseline hazard function components:\n")
  print(round(obj$h0, digits=digits))
  ##
  invisible()
}

print.summ.Bayes <- function(x, digits=3, ...)
{
    nChain = x$setup$nChain
    
    if(x$classFit[2] == "ID")
    {
        if(x$classFit[3] == "Cor")
        {
            ##
            cat("\nAnalysis of cluster-correlated semi-competing risks data \n")
        }
        if(x$classFit[3] == "Ind")
        {
            ##
            cat("\nAnalysis of independent semi-competing risks data \n")
        }
        ##
        cat(x$setup$model, "assumption for h3\n")
    }
    if(x$classFit[2] == "Surv")
    {
        if(x$classFit[3] == "Cor")
        {
            ##
            cat("\nAnalysis of cluster-correlated univariate time-to-event data \n")
        }
        if(x$classFit[3] == "Ind")
        {
            ##
            cat("\nAnalysis of independent univariate time-to-event data \n")
        }
    }
    
    cat("\n#####\n")
    
    ##
    cat("\nHazard ratios:\n")
    print(round(x$coef, digits=digits))

    if(x$classFit[2] == "ID")
    {
        ##
        cat("\nVariance of frailties:\n")
        print(round(x$theta, digits=digits))
    }
    
    ##
    cat("\nBaseline hazard function components:\n")
    print(round(x$h0, digits=digits))
    
    if(x$classFit[3] == "Cor")
    {
        if(x$classFit[5] == "DPM")
        {
            ##
            cat("\nPrecision parameter of DPM prior:\n")
            print(round(x$tau, digits=digits))
        }
        
        if(x$classFit[2] == "ID")
        {
            ##
            cat("\nVariance-covariance matrix of cluster-specific random effects:\n")
            print(round(x$Sigma.PM, digits=digits))
        }
        if(x$classFit[2] == "Surv")
        {
            ##
            cat("\nVariance of cluster-specific random effects:\n")
            print(round(x$sigma_V, digits=digits))
        }
        
    }
    
    invisible()
}



####
## PLOT METHOD
####
##
plot.Freq <- function(x, tseq=c(0, 5, 10), plot=TRUE, plot.est="BS", xlab=NULL, ylab=NULL, ...)
{
    obj <- x
    T2seq <- tseq
    yLim <- NULL
    ##
    ## SEs based on the Delta method using log(-log(S0))
    
    ##
    if(class(obj)[2] == "Surv")
    {
        T2 <- seq(from=min(T2seq), to=max(T2seq), length=100)
        
        ##
        kappa    <- exp(obj$estimate[1])
        alpha    <- exp(obj$estimate[2])
        log_kappa    <- obj$estimate[1]
        log_alpha    <- obj$estimate[2]
        S0     <- exp(-(kappa*(T2)^alpha))
        ## Delta method based on log(-log(S0))
        #J            <- cbind(1/kappa, log(T2))
        J            <- cbind(1, exp(log_alpha)*log(T2))
        Var.loglogS0 <- J %*% obj$Finv[1:2,1:2] %*% t(J)
        se.loglogS0  <- sqrt(diag(Var.loglogS0))
        se.loglogS0[is.na(se.loglogS0)] <- 0
        LL <- S0^exp(-qnorm(0.025)*se.loglogS0)
        UL <- S0^exp(qnorm(0.025)*se.loglogS0)
        ##
        BS_tbl <- cbind(T2, S0, LL, UL)
        dimnames(BS_tbl) <- list(rep("", length(T2)), c("time", "S0", "LL", "UL"))
        
        ##
        h0  <- alpha*kappa*(T2)^(alpha-1)
        J   <- cbind(h0, h0*(1+alpha*log(T2)))
        Var.h0 <- J %*% obj$Finv[1:2,1:2] %*% t(J)
        se.h0  <- sqrt(diag(Var.h0))
        se.h0[is.nan(se.h0)] <- 0
        LLh0 <- h0 - qnorm(0.025)*se.h0
        ULh0 <- h0 + qnorm(0.025)*se.h0
        LLh0[LLh0 < 0] <- 0
     
     
        T2h <- T2
        if(T2[1] == 0)
        {
            T2h <- T2h[-1]
            h0 <- h0[-1]
            LLh0 <- LLh0[-1]
            ULh0 <- ULh0[-1]
        }
        
        BH_tbl <- cbind(T2h, h0, LLh0, ULh0)
        dimnames(BH_tbl) <- list(rep("", length(T2h)), c("time", "h0", "LL", "UL"))
        
        value <- list(h0=BH_tbl, S0=BS_tbl)
        
        ##
        if(is.null(yLim))
        {
            if(plot.est=="BS")
            {
                yLim <- seq(from=0, to=1, by=0.2)
            }
            if(plot.est=="BH")
            {
                grid <- (max(ULh0) - min(LLh0))/5
                yLim <- seq(from=min(LLh0), to=max(ULh0), by=grid)
            }
        }
        
        ##
        if(is.null(ylab))
        {
            if(plot.est=="BS")
            {
                ylab <- "Baseline survival"
            }
            if(plot.est=="BH")
            {
                ylab <- "Baseline hazard"
            }
        }
        
        ##
        if(is.null(xlab)) xlab <- "Time"
        
        ##
        if(plot == TRUE){
            if(plot.est == "BS")
            {
                ##
                plot(range(T2seq), range(yLim), xlab=xlab, ylab=ylab, type="n", main = expression(paste("Estimated ", S[0](t), "")), axes=FALSE)
                axis(1, at=T2seq)
                axis(2, at=yLim)
                lines(T2, S0, col="red", lwd=3)
                lines(T2, LL, col="red", lwd=3, lty=3)
                lines(T2, UL, col="red", lwd=3, lty=3)
            }
            if(plot.est == "BH")
            {
                ##
                plot(range(T2seq), range(yLim), xlab=xlab, ylab=ylab, type="n", main = expression(paste("Estimated ", h[0](t), "")), axes=FALSE)
                axis(1, at=T2seq)
                axis(2, at=round(yLim, 4))
                lines(T2h, h0, col="red", lwd=3)
                lines(T2h, LLh0, col="red", lwd=3, lty=3)
                lines(T2h, ULh0, col="red", lwd=3, lty=3)
            }
        }
        if(plot == FALSE) return(value)
    }
    
    ##
    if(class(obj)[2] == "ID")
    {
        ##
        T2 <- seq(from=min(T2seq), to=max(T2seq), length=100)
        ##
        
        kappa    <- exp(obj$estimate[1])
        alpha 	 <- exp(obj$estimate[2])
        log_alpha <- obj$estimate[2]
        S0.1   		<- exp(-kappa*(T2)^alpha)
        J            <- cbind(1, exp(log_alpha)*log(T2))
        Var.loglogS0 <- J %*% obj$Finv[1:2,1:2] %*% t(J)
        se.loglogS0  <- sqrt(diag(Var.loglogS0))
        LL.1         <- S0.1^exp(-qnorm(0.025)*se.loglogS0)
        UL.1         <- S0.1^exp(qnorm(0.025)*se.loglogS0)
        ##
        h0.1  <- alpha*kappa*(T2)^(alpha-1)
        J   <- cbind(h0.1, h0.1*(1+alpha*log(T2)))
        Var.h0.1 <- J %*% obj$Finv[1:2,1:2] %*% t(J)
        se.h0.1  <- sqrt(diag(Var.h0.1))
        se.h0.1[is.nan(se.h0.1)] <- 0
        LLh0.1 <- h0.1 - qnorm(0.025)*se.h0.1
        ULh0.1 <- h0.1 + qnorm(0.025)*se.h0.1
        LLh0.1[LLh0.1 < 0] <- 0
        
        ##
        kappa    <- exp(obj$estimate[3])
        alpha 	 <- exp(obj$estimate[4])
        log_alpha <- obj$estimate[4]
        S0.2   <- exp(-kappa*(T2)^alpha)
        J            <- cbind(1, exp(log_alpha)*log(T2))
        Var.loglogS0 <- J %*% obj$Finv[3:4,3:4] %*% t(J)
        se.loglogS0  <- sqrt(diag(Var.loglogS0))
        LL.2         <- S0.2^exp(-qnorm(0.025)*se.loglogS0)
        UL.2         <- S0.2^exp(qnorm(0.025)*se.loglogS0)
        ##
        h0.2  <- alpha*kappa*(T2)^(alpha-1)
        J   <- cbind(h0.2, h0.2*(1+alpha*log(T2)))
        Var.h0.2 <- J %*% obj$Finv[1:2,1:2] %*% t(J)
        se.h0.2  <- sqrt(diag(Var.h0.2))
        se.h0.2[is.nan(se.h0.2)] <- 0
        LLh0.2 <- h0.2 - qnorm(0.025)*se.h0.2
        ULh0.2 <- h0.2 + qnorm(0.025)*se.h0.2
        LLh0.2[LLh0.2 < 0] <- 0
        
        ##
        kappa    <- exp(obj$estimate[5])
        alpha 	 <- exp(obj$estimate[6])
        log_alpha <- obj$estimate[6]
        S0.3   <- exp(-kappa*(T2)^alpha)
        J            <- cbind(1, exp(log_alpha)*log(T2))
        Var.loglogS0 <- J %*% obj$Finv[5:6,5:6] %*% t(J)
        se.loglogS0  <- sqrt(diag(Var.loglogS0))
        LL.3         <- S0.3^exp(-qnorm(0.025)*se.loglogS0)
        UL.3         <- S0.3^exp(qnorm(0.025)*se.loglogS0)
        ##
        h0.3  <- alpha*kappa*(T2)^(alpha-1)
        J   <- cbind(h0.3, h0.3*(1+alpha*log(T2)))
        Var.h0.3 <- J %*% obj$Finv[1:2,1:2] %*% t(J)
        se.h0.3  <- sqrt(diag(Var.h0.3))
        se.h0.3[is.nan(se.h0.3)] <- 0
        LLh0.3 <- h0.3 - qnorm(0.025)*se.h0.3
        ULh0.3 <- h0.3 + qnorm(0.025)*se.h0.3
        LLh0.3[LLh0.3 < 0] <- 0

        T2h <- T2
        if(T2[1] == 0)
        {
            T2h <- T2h[-1]
            h0.1 <- h0.1[-1]
            LLh0.1 <- LLh0.1[-1]
            ULh0.1 <- ULh0.1[-1]
            
            h0.2 <- h0.2[-1]
            LLh0.2 <- LLh0.2[-1]
            ULh0.2 <- ULh0.2[-1]
            
            h0.3 <- h0.3[-1]
            LLh0.3 <- LLh0.3[-1]
            ULh0.3 <- ULh0.3[-1]
        }
        
        BH1_tbl <- cbind(T2h, h0.1, LLh0.1, ULh0.1)
        dimnames(BH1_tbl) <- list(rep("", length(T2h)), c("time", "h0.1", "LL.1", "UL.1"))
        BH2_tbl <- cbind(T2h, h0.2, LLh0.2, ULh0.2)
        dimnames(BH2_tbl) <- list(rep("", length(T2h)), c("time", "h0.2", "LL.2", "UL.2"))
        BH3_tbl <- cbind(T2h, h0.3, LLh0.3, ULh0.3)
        dimnames(BH3_tbl) <- list(rep("", length(T2h)), c("time", "h0.3", "LL.3", "UL.3"))
        
        BS1_tbl <- cbind(T2, S0.1, LL.1, UL.1)
        dimnames(BS1_tbl) <- list(rep("", length(T2)), c("time", "S0.1", "LL.1", "UL.1"))
        BS2_tbl <- cbind(T2, S0.2, LL.2, UL.2)
        dimnames(BS2_tbl) <- list(rep("", length(T2)), c("time", "S0.2", "LL.2", "UL.2"))
        BS3_tbl <- cbind(T2, S0.3, LL.3, UL.3)
        dimnames(BS3_tbl) <- list(rep("", length(T2)), c("time", "S0.3", "LL.3", "UL.3"))
        
        value <- list(h0.1=BH1_tbl, h0.2=BH2_tbl, h0.3=BH3_tbl, S0.1=BS1_tbl, S0.2=BS2_tbl, S0.3=BS3_tbl)
        
        ##
        if(is.null(yLim))
        {
            if(plot.est=="BS")
            {
                yLim <- seq(from=0, to=1, by=0.2)
            }
            if(plot.est=="BH")
            {
                grid <- (max(ULh0.1, ULh0.2, ULh0.3) - min(LLh0.1, LLh0.2, LLh0.3))/5
                yLim <- seq(from=min(LLh0.1, LLh0.2, LLh0.3), to=max(ULh0.1, ULh0.2, ULh0.3), by=grid)
            }
        }
        
        ##
        if(is.null(ylab))
        {
            if(plot.est=="BS")
            {
                ylab <- "Baseline survival"
            }
            if(plot.est=="BH")
            {
                ylab <- "Baseline hazard"
            }
        }
        
        ##
        if(is.null(xlab))
        {
            xlab <- c("Time", "Time", "Time")
            if(class(obj)[5] == "semi-Markov")
            {
                xlab[3] <- "Time since non-terminal event"
            }
        }
        
        ##
        if(plot == TRUE){
            if(plot.est == "BS")
            {
                ##
                par(mfrow=c(1,3))
                ##
                plot(range(T2seq), range(yLim), xlab=xlab[1], ylab=ylab, type="n", main = expression(paste("Estimated ", S[0][1](t), "")), axes=FALSE)
                axis(1, at=T2seq)
                axis(2, at=yLim)
                lines(T2, S0.1, col="blue", lwd=3)
                lines(T2, LL.1, col="blue", lwd=3, lty=3)
                lines(T2, UL.1, col="blue", lwd=3, lty=3)
                ##
                plot(range(T2seq), range(yLim), xlab=xlab[2], ylab=ylab, type="n", main = expression(paste("Estimated ", S[0][2](t), "")), axes=FALSE)
                axis(1, at=T2seq)
                axis(2, at=yLim)
                lines(T2, S0.2, col="red", lwd=3)
                lines(T2, LL.2, col="red", lwd=3, lty=3)
                lines(T2, UL.2, col="red", lwd=3, lty=3)
                ##
                plot(range(T2seq), range(yLim), xlab=xlab[3], ylab=ylab, type="n", main = expression(paste("Estimated ", S[0][3](t), "")), axes=FALSE)
                axis(1, at=T2seq)
                axis(2, at=yLim)
                lines(T2, S0.3, col="red", lwd=3)
                lines(T2, LL.3, col="red", lwd=3, lty=3)
                lines(T2, UL.3, col="red", lwd=3, lty=3)
            }
            if(plot.est == "BH")
            {
                ##
                par(mfrow=c(1,3))
                ##
                plot(range(T2seq), range(yLim), xlab=xlab[1], ylab=ylab, type="n", main = expression(paste("Estimated ", h[0][1](t), "")), axes=FALSE)
                axis(1, at=T2seq)
                axis(2, at=round(yLim, 4))
                lines(T2h, h0.1, col="blue", lwd=3)
                lines(T2h, LLh0.1, col="blue", lwd=3, lty=3)
                lines(T2h, ULh0.1, col="blue", lwd=3, lty=3)
                ##
                plot(range(T2seq), range(yLim), xlab=xlab[2], ylab=ylab, type="n", main = expression(paste("Estimated ", h[0][2](t), "")), axes=FALSE)
                axis(1, at=T2seq)
                axis(2, at=round(yLim, 4))
                lines(T2h, h0.2, col="red", lwd=3)
                lines(T2h, LLh0.2, col="red", lwd=3, lty=3)
                lines(T2h, ULh0.2, col="red", lwd=3, lty=3)
                ##
                plot(range(T2seq), range(yLim), xlab=xlab[3], ylab=ylab, type="n", main = expression(paste("Estimated ", h[0][3](t), "")), axes=FALSE)
                axis(1, at=T2seq)
                axis(2, at=round(yLim, 4))
                lines(T2h, h0.3, col="red", lwd=3)
                lines(T2h, LLh0.3, col="red", lwd=3, lty=3)
                lines(T2h, ULh0.3, col="red", lwd=3, lty=3)
            }
        }
        if(plot == FALSE) return(value)
    }
    ##
    invisible()
}


plot.Bayes <- function(x, tseq=c(0, 5, 10), plot=TRUE, plot.est="BS", xlab=NULL, ylab=NULL, ...)
{
    nChain = x$setup$nChain
    
    if(class(x)[2] == "ID")
    {
        if(class(x)[4] == "PEM")
        {
            time1 <- x$chain1$time_lambda1
            time2 <- x$chain1$time_lambda2
            time3 <- x$chain1$time_lambda3
            
            time1hz <- time1
            time2hz <- time2
            time3hz <- time3
            
            lambda1.fin	<- x$chain1$lambda1.fin
            lambda2.fin	<- x$chain1$lambda2.fin
            lambda3.fin	<- x$chain1$lambda3.fin
            
            if(nChain > 1){
                for(i in 2:nChain){
                    nam <- paste("chain", i, sep="")
                    lambda1.fin <- rbind(lambda1.fin, x[[nam]]$lambda1.fin)
                    lambda2.fin <- rbind(lambda2.fin, x[[nam]]$lambda2.fin)
                    lambda3.fin <- rbind(lambda3.fin, x[[nam]]$lambda3.fin)
                }
            }
            
            BH1Med  <-  apply(exp(lambda1.fin), 2, median)
            BH1Ub   <-  apply(exp(lambda1.fin), 2, quantile, prob = 0.975)
            BH1Lb   <-  apply(exp(lambda1.fin), 2, quantile, prob = 0.025)
            
            BH2Med  <-  apply(exp(lambda2.fin), 2, median)
            BH2Ub   <-  apply(exp(lambda2.fin), 2, quantile, prob = 0.975)
            BH2Lb   <-  apply(exp(lambda2.fin), 2, quantile, prob = 0.025)
            
            BH3Med  <-  apply(exp(lambda3.fin), 2, median)
            BH3Ub   <-  apply(exp(lambda3.fin), 2, quantile, prob = 0.975)
            BH3Lb   <-  apply(exp(lambda3.fin), 2, quantile, prob = 0.025)
            
            dif1    <- diff(c(0, time1hz))
            dif2    <- diff(c(0, time2hz))
            dif3    <- diff(c(0, time3hz))
            
            BS1     <- matrix(NA, dim(lambda1.fin)[1], dim(lambda1.fin)[2])
            for(i in 1:dim(lambda1.fin)[1])
            {
                BS1[i,] <- exp(-cumsum(exp(lambda1.fin[i,])* dif1) )
            }
            BS2     <- matrix(NA, dim(lambda2.fin)[1], dim(lambda2.fin)[2])
            for(i in 1:dim(lambda2.fin)[1])
            {
                BS2[i,] <- exp(-cumsum(exp(lambda2.fin[i,])* dif2) )
            }
            BS3     <- matrix(NA, dim(lambda3.fin)[1], dim(lambda3.fin)[2])
            for(i in 1:dim(lambda3.fin)[1])
            {
                BS3[i,] <- exp(-cumsum(exp(lambda3.fin[i,])* dif3) )
            }
            
            BS1Med  <-  apply(BS1, 2, median)
            BS1Ub   <-  apply(BS1, 2, quantile, prob = 0.975)
            BS1Lb   <-  apply(BS1, 2, quantile, prob = 0.025)
            
            BS2Med  <-  apply(BS2, 2, median)
            BS2Ub   <-  apply(BS2, 2, quantile, prob = 0.975)
            BS2Lb   <-  apply(BS2, 2, quantile, prob = 0.025)
            
            BS3Med  <-  apply(BS3, 2, median)
            BS3Ub   <-  apply(BS3, 2, quantile, prob = 0.975)
            BS3Lb   <-  apply(BS3, 2, quantile, prob = 0.025)
        }
        
        if(class(x)[4] == "WB")
        {
            time1 <- time2 <- time3 <- seq(from=min(tseq), to=max(tseq), length=100)
            nStore <- length(x$chain1$alpha1.p)
            numSpl <- nStore * nChain
            
            basehaz1 <- matrix(NA, numSpl, length(time1))
            basehaz2 <- matrix(NA, numSpl, length(time2))
            basehaz3 <- matrix(NA, numSpl, length(time3))
            basesurv1 <- matrix(NA, numSpl, length(time1))
            basesurv2 <- matrix(NA, numSpl, length(time2))
            basesurv3 <- matrix(NA, numSpl, length(time3))
            
            alpha1.p	<- x$chain1$alpha1.p
            alpha2.p	<- x$chain1$alpha2.p
            alpha3.p	<- x$chain1$alpha3.p
            kappa1.p	<- x$chain1$kappa1.p
            kappa2.p	<- x$chain1$kappa2.p
            kappa3.p	<- x$chain1$kappa3.p
            
            if(nChain > 1){
                for(i in 2:nChain){
                    nam <- paste("chain", i, sep="")
                    alpha1.p <- c(alpha1.p, x[[nam]]$alpha1.p)
                    alpha2.p <- c(alpha2.p, x[[nam]]$alpha2.p)
                    alpha3.p <- c(alpha3.p, x[[nam]]$alpha3.p)
                    
                    kappa1.p <- c(kappa1.p, x[[nam]]$kappa1.p)
                    kappa2.p <- c(kappa2.p, x[[nam]]$kappa2.p)
                    kappa3.p <- c(kappa3.p, x[[nam]]$kappa3.p)
                }
            }
            
            for(i in 1:numSpl){
                basehaz1[i, ] <- alpha1.p[i] * kappa1.p[i] * time1^(alpha1.p[i] - 1)
                basehaz2[i, ] <- alpha2.p[i] * kappa2.p[i] * time2^(alpha2.p[i] - 1)
                basehaz3[i, ] <- alpha3.p[i] * kappa3.p[i] * time3^(alpha3.p[i] - 1)
                
                basesurv1[i, ] <- exp(-kappa1.p[i] * time1^(alpha1.p[i]))
                basesurv2[i, ] <- exp(-kappa2.p[i] * time2^(alpha2.p[i]))
                basesurv3[i, ] <- exp(-kappa3.p[i] * time3^(alpha3.p[i]))
            }
            
            time1hz <- time1
            time2hz <- time2
            time3hz <- time3
            
            if(tseq[1] == 0){
                time1hz <- time1[-1]
                time2hz <- time2[-1]
                time3hz <- time3[-1]
                
                basehaz1 <- basehaz1[,-1]
                basehaz2 <- basehaz2[,-1]
                basehaz3 <- basehaz3[,-1]
            }
            
            BH1Med <- apply(basehaz1, 2, median)
            BH1Ub <- apply(basehaz1, 2, quantile, prob = 0.975)
            BH1Lb <- apply(basehaz1, 2, quantile, prob = 0.025)
            BH2Med <- apply(basehaz2, 2, median)
            BH2Ub <- apply(basehaz2, 2, quantile, prob = 0.975)
            BH2Lb <- apply(basehaz2, 2, quantile, prob = 0.025)
            BH3Med <- apply(basehaz3, 2, median)
            BH3Ub <- apply(basehaz3, 2, quantile, prob = 0.975)
            BH3Lb <- apply(basehaz3, 2, quantile, prob = 0.025)
            
            BS1Med <- apply(basesurv1, 2, median)
            BS1Ub <- apply(basesurv1, 2, quantile, prob = 0.975)
            BS1Lb <- apply(basesurv1, 2, quantile, prob = 0.025)
            BS2Med <- apply(basesurv2, 2, median)
            BS2Ub <- apply(basesurv2, 2, quantile, prob = 0.975)
            BS2Lb <- apply(basesurv2, 2, quantile, prob = 0.025)
            BS3Med <- apply(basesurv3, 2, median)
            BS3Ub <- apply(basesurv3, 2, quantile, prob = 0.975)
            BS3Lb <- apply(basesurv3, 2, quantile, prob = 0.025)
        }
        
        BH1_tbl <- cbind(time1hz, BH1Med, BH1Lb, BH1Ub)
        dimnames(BH1_tbl) <- list(rep("", length(time1hz)), c("time", "h0.1", "LL.1", "UL.1"))
        BH2_tbl <- cbind(time2hz, BH2Med, BH2Lb, BH2Ub)
        dimnames(BH2_tbl) <- list(rep("", length(time2hz)), c("time", "h0.2", "LL.2", "UL.2"))
        BH3_tbl <- cbind(time3hz, BH3Med, BH3Lb, BH3Ub)
        dimnames(BH3_tbl) <- list(rep("", length(time3hz)), c("time", "h0.3", "LL.3", "UL.3"))
        
        BS1_tbl <- cbind(time1, BS1Med, BS1Lb, BS1Ub)
        dimnames(BS1_tbl) <- list(rep("", length(time1)), c("time", "S0.1", "LL.1", "UL.1"))
        BS2_tbl <- cbind(time2, BS2Med, BS2Lb, BS2Ub)
        dimnames(BS2_tbl) <- list(rep("", length(time2)), c("time", "S0.2", "LL.2", "UL.2"))
        BS3_tbl <- cbind(time3, BS3Med, BS3Lb, BS3Ub)
        dimnames(BS3_tbl) <- list(rep("", length(time3)), c("time", "S0.3", "LL.3", "UL.3"))
        
        value <- list(h0.1=BH1_tbl, h0.2=BH2_tbl, h0.3=BH3_tbl, S0.1=BS1_tbl, S0.2=BS2_tbl, S0.3=BS3_tbl)
        
        if(plot == TRUE)
        {
            if(is.null(xlab))
            {
                xlab <- c("Time", "Time", "Time")
                if(x$setup$model == "semi-Markov")
                {
                   xlab[3] <- "Time since non-terminal event"
                }
            }
            
            if(plot.est == "BH")
            {
                if(is.null(ylab))
                {
                    ylab <- "Baseline hazard"
                }
                
                ygrid <- (max(BH1Ub, BH2Ub, BH3Ub) - 0)/5
                ylim <- seq(from=0, to=max(BH1Ub, BH2Ub, BH3Ub), by=ygrid)
                
                ##
                par(mfrow=c(1,3))
                ##
                plot(c(0, max(time1)), range(ylim), xlab=xlab[1], ylab=ylab, type="n", main = expression(paste("Estimated ", h[0][1](t), "")), axes=FALSE)
                if(class(x)[4] == "PEM")
                {
                    axis(1, at=c(0, max(time1)))
                }
                if(class(x)[4] == "WB")
                {
                    axis(1, at=tseq)
                }
                axis(2, at=round(ylim, 4))
                
                
                #if(time1hz[1] == 0)
                #{
                #    lines(time1hz, BH1Med, col="blue", lwd=3)
                #    lines(time1hz, BH1Ub, col="blue", lwd=3, lty=3)
                #    lines(time1hz, BH1Lb, col="blue", lwd=3, lty=3)
                #}else
                #{
                #    lines(unique(c(0, time1hz)), c(0, BH1Med), col="red", lwd=3)
                #    lines(unique(c(0, time1hz)), c(0, BH1Ub), col="red", lwd=3, lty=3)
                #    lines(unique(c(0, time1hz)), c(0, BH1Lb), col="red", lwd=3, lty=3)
                #}
                lines(time1hz, BH1Med, col="blue", lwd=3)
                lines(time1hz, BH1Ub, col="blue", lwd=3, lty=3)
                lines(time1hz, BH1Lb, col="blue", lwd=3, lty=3)

                ##
                plot(c(0, max(time2)), range(ylim), xlab=xlab[2], ylab=ylab, type="n", main = expression(paste("Estimated ", h[0][2](t), "")), axes=FALSE)
                if(class(x)[4] == "PEM")
                {
                    axis(1, at=c(0, max(time2)))
                }
                if(class(x)[4] == "WB")
                {
                    axis(1, at=tseq)
                }
                axis(2, at=round(ylim, 4))
                
                #if(time2hz[1] == 0)
                #{
                #    lines(time2hz, BH2Med, col="blue", lwd=3)
                #    lines(time2hz, BH2Ub, col="blue", lwd=3, lty=3)
                #    lines(time2hz, BH2Lb, col="blue", lwd=3, lty=3)
                #}else
                #{
                #    lines(unique(c(0, time2hz)), c(0, BH2Med), col="red", lwd=3)
                #    lines(unique(c(0, time2hz)), c(0, BH2Ub), col="red", lwd=3, lty=3)
                #    lines(unique(c(0, time2hz)), c(0, BH2Lb), col="red", lwd=3, lty=3)
                #}
                lines(time2hz, BH2Med, col="red", lwd=3)
                lines(time2hz, BH2Ub, col="red", lwd=3, lty=3)
                lines(time2hz, BH2Lb, col="red", lwd=3, lty=3)
                
                
                ##
                plot(c(0, max(time3)), range(ylim), xlab=xlab[3], ylab=ylab, type="n", main = expression(paste("Estimated ", h[0][3](t), "")), axes=FALSE)
                if(class(x)[4] == "PEM")
                {
                    axis(1, at=c(0, max(time3)))
                }
                if(class(x)[4] == "WB")
                {
                    axis(1, at=tseq)
                }
                axis(2, at=round(ylim, 4))
                
                #if(time3hz[1] == 0)
                #{
                #    lines(time3hz, BH3Med, col="blue", lwd=3)
                #    lines(time3hz, BH3Ub, col="blue", lwd=3, lty=3)
                #    lines(time3hz, BH3Lb, col="blue", lwd=3, lty=3)
                #}else
                #{
                #    lines(unique(c(0, time3hz)), c(0, BH3Med), col="red", lwd=3)
                #    lines(unique(c(0, time3hz)), c(0, BH3Ub), col="red", lwd=3, lty=3)
                #    lines(unique(c(0, time3hz)), c(0, BH3Lb), col="red", lwd=3, lty=3)
                #}
                lines(time3hz, BH3Med, col="red", lwd=3)
                lines(time3hz, BH3Ub, col="red", lwd=3, lty=3)
                lines(time3hz, BH3Lb, col="red", lwd=3, lty=3)

            }
            
            if(plot.est == "BS")
            {
                if(is.null(ylab))
                {
                    ylab <- "Baseline survival"
                }
                ylim <- seq(from=0, to=1, by=0.2)
                
                ##
                par(mfrow=c(1,3))
                ##
                plot(c(0, max(time1)), range(ylim), xlab=xlab[1], ylab=ylab, type="n", main = expression(paste("Estimated ", S[0][1](t), "")), axes=FALSE)
                if(class(x)[4] == "PEM")
                {
                    axis(1, at=c(0, max(time1)))
                }
                if(class(x)[4] == "WB")
                {
                    axis(1, at=tseq)
                }
                axis(2, at=ylim)
                
                if(time1[1] == 0)
                {
                    lines(time1, BS1Med, col="blue", lwd=3)
                    lines(time1, BS1Ub, col="blue", lwd=3, lty=3)
                    lines(time1, BS1Lb, col="blue", lwd=3, lty=3)
                }else
                {
                    lines(unique(c(0, time1)), c(1, BS1Med), col="red", lwd=3)
                    lines(unique(c(0, time1)), c(1, BS1Ub), col="red", lwd=3, lty=3)
                    lines(unique(c(0, time1)), c(1, BS1Lb), col="red", lwd=3, lty=3)
                }

                ##
                plot(c(0, max(time2)), range(ylim), xlab=xlab[2], ylab=ylab, type="n", main = expression(paste("Estimated ", S[0][2](t), "")), axes=FALSE)
                if(class(x)[4] == "PEM")
                {
                    axis(1, at=c(0, max(time2)))
                }
                if(class(x)[4] == "WB")
                {
                    axis(1, at=tseq)
                }
                axis(2, at=ylim)
                
                if(time2[1] == 0)
                {
                    lines(time2, BS2Med, col="blue", lwd=3)
                    lines(time2, BS2Ub, col="blue", lwd=3, lty=3)
                    lines(time2, BS2Lb, col="blue", lwd=3, lty=3)
                }else
                {
                    lines(unique(c(0, time2)), c(1, BS2Med), col="red", lwd=3)
                    lines(unique(c(0, time2)), c(1, BS2Ub), col="red", lwd=3, lty=3)
                    lines(unique(c(0, time2)), c(1, BS2Lb), col="red", lwd=3, lty=3)
                }

                ##
                plot(c(0, max(time3)), range(ylim), xlab=xlab[3], ylab=ylab, type="n", main = expression(paste("Estimated ", S[0][3](t), "")), axes=FALSE)
                if(class(x)[4] == "PEM")
                {
                    axis(1, at=c(0, max(time3)))
                }
                if(class(x)[4] == "WB")
                {
                    axis(1, at=tseq)
                }
                axis(2, at=ylim)
                
                if(time3[1] == 0)
                {
                    lines(time3, BS3Med, col="blue", lwd=3)
                    lines(time3, BS3Ub, col="blue", lwd=3, lty=3)
                    lines(time3, BS3Lb, col="blue", lwd=3, lty=3)
                }else
                {
                    lines(unique(c(0, time3)), c(1, BS3Med), col="red", lwd=3)
                    lines(unique(c(0, time3)), c(1, BS3Ub), col="red", lwd=3, lty=3)
                    lines(unique(c(0, time3)), c(1, BS3Lb), col="red", lwd=3, lty=3)
                }
                
            }
        }
        if(plot == FALSE)
        {
            return(value)
        }
        
    }
    
    if(class(x)[2] == "Surv")
    {
        if(class(x)[4] == "PEM")
        {
            time <- x$chain1$time_lambda
            timehz <- time
            
            lambda.fin	<- x$chain1$lambda.fin
            
            if(nChain > 1){
                for(i in 2:nChain){
                    nam <- paste("chain", i, sep="")
                    lambda.fin <- rbind(lambda.fin, x[[nam]]$lambda.fin)
                }
            }
            
            BHMed  <-  apply(exp(lambda.fin), 2, median)
            BHUb   <-  apply(exp(lambda.fin), 2, quantile, prob = 0.975)
            BHLb   <-  apply(exp(lambda.fin), 2, quantile, prob = 0.025)
            
            dif <- diff(c(0, timehz))
            
            BS     <- matrix(NA, dim(lambda.fin)[1], dim(lambda.fin)[2])
            for(i in 1:dim(lambda.fin)[1])
            {
                BS[i,] <- exp(-cumsum(exp(lambda.fin[i,])* dif) )
            }
            
            BSMed  <-  apply(BS, 2, median)
            BSUb   <-  apply(BS, 2, quantile, prob = 0.975)
            BSLb   <-  apply(BS, 2, quantile, prob = 0.025)
        }
        
        if(class(x)[4] == "WB")
        {
            time <- seq(from=min(tseq), to=max(tseq), length=100)
            nStore <- length(x$chain1$alpha.p)
            numSpl <- nStore * nChain
            
            basehaz <- matrix(NA, numSpl, length(time))
            basesurv <- matrix(NA, numSpl, length(time))
            
            alpha.p	<- x$chain1$alpha.p
            kappa.p	<- x$chain1$kappa.p
            
            if(nChain > 1){
                for(i in 2:nChain){
                    nam <- paste("chain", i, sep="")
                    alpha.p <- c(alpha.p, x[[nam]]$alpha.p)
                    kappa.p <- c(kappa.p, x[[nam]]$kappa.p)
                }
            }
            
            for(i in 1:numSpl){
                basehaz[i, ] <- alpha.p[i] * kappa.p[i] * time^(alpha.p[i] - 1)
                basesurv[i, ] <- exp(-kappa.p[i] * time^(alpha.p[i]))
            }
            
            timehz <- time
            
            if(tseq[1] == 0){
                timehz <- time[-1]
                basehaz <- basehaz[,-1]
            }
            
            BHMed <- apply(basehaz, 2, median)
            BHUb <- apply(basehaz, 2, quantile, prob = 0.975)
            BHLb <- apply(basehaz, 2, quantile, prob = 0.025)
            
            BSMed <- apply(basesurv, 2, median)
            BSUb <- apply(basesurv, 2, quantile, prob = 0.975)
            BSLb <- apply(basesurv, 2, quantile, prob = 0.025)
        }
        
        BH_tbl <- cbind(timehz, BHMed, BHLb, BHUb)
        dimnames(BH_tbl) <- list(rep("", length(timehz)), c("time", "h0", "LL", "UL"))
        
        BS_tbl <- cbind(time, BSMed, BSLb, BSUb)
        dimnames(BS_tbl) <- list(rep("", length(time)), c("time", "S0", "LL", "UL"))
        
        value <- list(h0=BH_tbl, S0=BS_tbl)
        
        
        if(plot == TRUE)
        {
            if(is.null(xlab))
            {
                xlab <- "Time"
            }
            
            if(plot.est == "BH")
            {
                if(is.null(ylab))
                {
                    ylab <- "Baseline hazard"
                }
                
                ygrid <- (max(BHUb) - 0)/5
                ylim <- seq(from=0, to=max(BHUb), by=ygrid)
                
                ##
                plot(c(0, max(time)), range(ylim), xlab=xlab, ylab=ylab, type="n", main = expression(paste("Estimated ", h[0](t), "")), axes=FALSE)
                if(class(x)[4] == "PEM")
                {
                    axis(1, at=c(0, max(time)))
                }
                if(class(x)[4] == "WB")
                {
                    axis(1, at=tseq)
                }
                axis(2, at=round(ylim, 4))
                
                #if(timehz[1] == 0)
                #{
                #    lines(timehz, BHMed, col="red", lwd=3)
                #    lines(timehz, BHUb, col="red", lwd=3, lty=3)
                #    lines(timehz, BHLb, col="red", lwd=3, lty=3)
                #}else
                #{
                #    lines(unique(c(0, timehz)), c(0, BHMed), col="red", lwd=3)
                #    lines(unique(c(0, timehz)), c(0, BHUb), col="red", lwd=3, lty=3)
                #    lines(unique(c(0, timehz)), c(0, BHLb), col="red", lwd=3, lty=3)
                #}
                lines(timehz, BHMed, col="red", lwd=3)
                lines(timehz, BHUb, col="red", lwd=3, lty=3)
                lines(timehz, BHLb, col="red", lwd=3, lty=3)

            }
            
            if(plot.est == "BS")
            {
                if(is.null(ylab))
                {
                    ylab <- "Baseline survival"
                }
                ylim <- seq(from=0, to=1, by=0.2)
                
                ##
                plot(c(0, max(time)), range(ylim), xlab=xlab, ylab=ylab, type="n", main = expression(paste("Estimated ", S[0](t), "")), axes=FALSE)
                if(class(x)[4] == "PEM")
                {
                    axis(1, at=c(0, max(time)))
                }
                if(class(x)[4] == "WB")
                {
                    axis(1, at=tseq)
                }
                axis(2, at=ylim)
                
                if(time[1] == 0)
                {
                    lines(time, BSMed, col="red", lwd=3)
                    lines(time, BSUb, col="red", lwd=3, lty=3)
                    lines(time, BSLb, col="red", lwd=3, lty=3)
                }else
                {
                    lines(unique(c(0, time)), c(1, BSMed), col="red", lwd=3)
                    lines(unique(c(0, time)), c(1, BSUb), col="red", lwd=3, lty=3)
                    lines(unique(c(0, time)), c(1, BSLb), col="red", lwd=3, lty=3)
                }
            }
        }
        if(plot == FALSE)
        {
            return(value)
        }
    }
    invisible()
}







