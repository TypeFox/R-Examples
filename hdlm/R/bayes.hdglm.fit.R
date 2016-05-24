bayes.hdglm.fit <- function(x,y,siglevel=0.05,bayesIters = NULL, bayesTune = 0.01) {

    #
    n <- nrow(x)
    p <- ncol(x)
    prob <- bayesTune
    if(is.null(bayesIters)) {
        iters <- 1000
    } else {
        iters <- bayesIters
    }

    # Number of workers:
    num <- getDoParWorkers()   
    if(num == 1) {
        out <- bayesianLatentFit(x=x,y=y,prob=prob,c=100,iter=iters,burnin=bayesIters, siglevel=siglevel)
    } else {
        iters <- ceiling(iters / num)
        outList <- foreach(icount(num), .inorder = FALSE) %dopar% {
            bayesianLatentFit(x=x,y=y,prob=prob,c=100,iter=iters,burnin=bayesIters, siglevel=siglevel)
        }
        DRAWS <- NULL
        FITTED <- NULL
        MODEL <- NULL
        for(j in 1:length(num)) {
            DRAWS <- rbind(DRAWS, outList[[j]]$beta)
            FITTED <- c(FITTED, outList[[j]]$fitted)
            MODEL <- rbind(MODEL, outList[[j]]$marginc)
        }
        MARGINC <- apply(MARGINC, 2, mean)
        out <- list(beta = DRAWS, fitted = FITTED, model = MARGINC)
    }

    #out <- blasso.vs(y, x, iters=10000, burn=100, thin=2, beta=rep(0,ncol(x)),
    #    sig2=1, tau=1, phi=0.5, tauprior=c(1,1), sig2prior=c(1,1))

    betaMCMC <- out$beta
    colnames(betaMCMC) <- colnames(x)
    
    point_estimator <- apply(betaMCMC, 2, mean)
    bound <- apply(betaMCMC, 2, quantile, probs = c(siglevel,1-siglevel), names=FALSE)
    pvalue <- 1 - apply(out$model,2,mean)

    # Set small coefficents to zero:
    index <- which(bound[1,] == 0 & bound[2,] == 0)
    point_estimator[index] <- 0

    fitted <- as.numeric(x %*% point_estimator > 0)
    resid <- fitted - y

    z <- list(coefficients=point_estimator, lower.bound=bound[1,], upper.bound=bound[2,], p.value=pvalue,
              effects=NULL, rank=dim(x), fitted.values=fitted, assign=NULL, residuals=resid,
              sigma.hat=NULL)
    return(z)
 }


