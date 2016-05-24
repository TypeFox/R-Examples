bayes.hdlm.fit <- function(x, y, siglevel=0.05, bayesIters=NULL, bayesTune=c(1,1)) {

    #
    n <- nrow(x)
    p <- ncol(x)
    if(is.null(bayesIters)) {
        iters <- 1000
    } else {
        iters <- bayesIters
    }

    # Number of workers:
    num <- getDoParWorkers()   
    if(num == 1) {
        out <- .C("blassoGibbs",
             X = as.double(x),
			 nn = as.integer(n),
			 pp = as.integer(p),
			 TT = as.integer(iters),
			 BB = as.integer(floor(iters*0.1)),
			 tthin = as.integer(2),
			 ttau = as.double(1),
			 ssig2 = as.double(var(y)),
			 pphi = as.double(0.01),
			 sig2prior = as.double(c(1,1)),
			 fits2 = as.integer(1),
			 tauprior = as.double(c(1,1)),
			 fittau = as.integer(1),
			 modunc = as.integer(TRUE),
			 phiprior = as.double(bayesTune),
			 fitphi = as.integer(1),
			 start = as.double(rep(0,p)),
			 bdraws = as.double(rep(0,p*iters)),
			 sig2draws = as.double(rep(0,iters)),
			 taudraws = as.double(rep(0,iters)),
			 phidraws = as.double(rep(0,iters)),
			 marginc = as.double(rep(0,p)),
			 rrb = as.integer(FALSE),
			 RB = as.double(rep(0,p)),
			 YtY = as.double(t(y)%*%y),
			 YtX = as.double(t(y)%*%x),
			 NOISY = as.integer(FALSE),
			 bprior = as.integer(1),
			 count = as.integer(0),
			 NAOK = TRUE,
			 PACKAGE="hdlm")
    } else {
        iters <- ceiling(iters / num)
        outList <- foreach(icount(num), .inorder = FALSE) %dopar% {
            .C("blassoGibbs",
             X = as.double(x),
			 nn = as.integer(n),
			 pp = as.integer(p),
			 TT = as.integer(iters),
			 BB = as.integer(floor(iters*0.1)),
			 tthin = as.integer(2),
			 ttau = as.double(1),
			 ssig2 = as.double(var(y)),
			 pphi = as.double(0.01),
			 sig2prior = as.double(c(1,1)),
			 fits2 = as.integer(1),
			 tauprior = as.double(c(1,1)),
			 fittau = as.integer(1),
			 modunc = as.integer(TRUE),
			 phiprior = as.double(c(1,1)),
			 fitphi = as.integer(1),
			 start = as.double(rep(0,p)),
			 bdraws = as.double(rep(0,p*iters)),
			 sig2draws = as.double(rep(0,iters)),
			 taudraws = as.double(rep(0,iters)),
			 phidraws = as.double(rep(0,iters)),
			 marginc = as.double(rep(0,p)),
			 rrb = as.integer(FALSE),
			 RB = as.double(rep(0,p)),
			 YtY = as.double(t(y)%*%y),
			 YtX = as.double(t(y)%*%x),
			 NOISY = as.integer(FALSE),
			 bprior = as.integer(1),
			 count = as.integer(0),
			 NAOK = TRUE,
			 PACKAGE="hdlm")
        }
        DRAWS <- NULL
        SIGDRAWS <- NULL
        MARGINC <- NULL
        for(j in 1:length(num)) {
            DRAWS <- rbind(DRAWS, outList[[j]]$bdraws)
            SIGDRAWS <- c(SIGDRAWS, outList[[j]]$sig2draws)
            MARGINC <- rbind(MARGINC, outList[[j]]$marginc)
        }
        MARGINC <- apply(MARGINC, 2, mean)
        out <- list(bdraws = DRAWS, sig2draws = SIGDRAWS, marginc = MARGINC)
    }

    #out <- blasso.vs(y, x, iters=10000, burn=100, thin=2, beta=rep(0,ncol(x)),
    #    sig2=1, tau=1, phi=0.5, tauprior=c(1,1), sig2prior=c(1,1))

    betaMCMC <- matrix(out$bdraws, ncol=p, byrow=T) 
    
    sigma_hat <- mean(sqrt(out$sig2draws))
    
    point_estimator <- apply(betaMCMC, 2, mean)
    bound <- apply(betaMCMC, 2, quantile, probs = c(siglevel,1-siglevel), names=FALSE)
    pvalue <- 1 - out$marginc

    # Set small coefficents to zero:
    index <- which(bound[1,] == 0 & bound[2,] == 0)
    point_estimator[index] <- 0

    fitted <- x %*% point_estimator
    resid <- fitted - y

    z <- list(coefficients=point_estimator, lower.bound=bound[1,], upper.bound=bound[2,], p.value=pvalue,
              effects=NULL, rank=dim(x), fitted.values=fitted, assign=NULL, residuals=resid,
              sigma.hat=sigma_hat)
    return(z)
 }


