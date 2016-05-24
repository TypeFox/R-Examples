getInitialKnots <- function(n, lambda, l){
    delta <- l+1
    linterval <- floor((n-2*(l+1)) / lambda)
    if(linterval < delta){
        stop("The number of knots is too larger. The knots have to be
             at least 'l+1' points away from eath other.")
    }
    else{
        knots <- rep(0, lambda)
        knots[1] <- sample(1:linterval+(l+1), 1)
        if(lambda > 1){
            for(i in 2:lambda)
                knots[i] <- sample((knots[i-1]+delta):(linterval*(i)+l+1), 1)
        }
        knots
  }
}


getDesignMatrix <- function(x, xknots, l, l0){
    .Call("curvePolynomialGetDM", l0, l, x, xknots)
    
}

fitLM <- function(y, dm){
    s <- lm(y ~ dm)
    rss <- sum((s$residuals)^2)
    fitted <- s$fitted.values
    betas <- as.vector(s$coefficients)
    list(rss=rss, fitted=fitted, betas=betas)
}


birthKnot <- function(y, x, n, knots, l, l0, sigma2, rss, fitted, betas,
                  gamma.shape, gamma.rate){
    pool <- (l+2):(n-l-1)
    candidates <- pool[-match(knots, pool)]
    lcandidate <- length(candidates)
    if(lcandidate > 0){
        choose <- rep(0, lcandidate)
        for(i in 1:lcandidate){
            choose[i] <- all(abs(candidates[i] - knots) >= (l+1))
        }
        candidates <- candidates[choose > 0]
        lcandidate <- length(candidates)
    }
    if(lcandidate > 0){
        choose <- sample(1:lcandidate, 1)
        newknot <- candidates[choose]
        knots.new <- c(knots, newknot)
        knots.new <- sort(knots.new)

        dm.new <- getDesignMatrix(x, x[knots.new], l, l0)
        rlm.new <- fitLM(y, dm.new)
        
        rss.new <- rlm.new$rss
        k <- length(knots)
        Zk <- 2*(l+1) + k*(2*l+1)
        alpha <- min(1, 1/sqrt(n) * ((2*gamma.rate + rss)/(2*gamma.rate + rss.new))^(n/2 + gamma.shape) * (n-Zk)/n)
                
        if(runif(1) < alpha){
            knots <- knots.new
            rss <- rss.new
            fitted <- rlm.new$fitted
            betas <- rlm.new$betas
         }
    }
    list(knots=knots, rss=rss, fitted=fitted, betas=betas)
}

deathKnot <- function(y, x, n, knots, l, l0, sigma2, rss, fitted, betas,
                  gamma.shape, gamma.rate){
    
    if(length(knots) > 1){
        choose <- sample(1:length(knots), 1)
        knots.new <- knots[-choose]

        dm.new <- getDesignMatrix(x, x[knots.new], l, l0)
        rlm.new <- fitLM(y, dm.new)
  
        rss.new <- rlm.new$rss
        k <- length(knots.new)
        Zk <- 2*(l+1) + k*(2*l+1)
        alpha <- min(1, sqrt(n) * ((2*gamma.rate + rss)/(2*gamma.rate + rss.new))^(n/2 + gamma.shape) * n/(n-Zk))
        
        if(runif(1) < alpha){
            knots <- knots.new
            rss <- rss.new
            fitted <- rlm.new$fitted
            betas <- rlm.new$betas
        }
    }
    list(knots=knots, rss=rss, fitted=fitted, betas=betas)
}

moveKnot <- function(y, x, n, knots, l, l0, sigma2, rss, fitted, betas,
                 gamma.shape, gamma.rate){
    choose <- sample(1:length(knots), 1)
    old <- knots[choose]
    candidates <- max(old - (l+1), l+2) : min(old + (l+1), n-l-1)
    candidates <- candidates[-which(candidates==old)]
    lcandidate <- length(candidates)
    if(lcandidate > 0){
        choose <- rep(0, lcandidate)
        for(i in 1:lcandidate){
            choose[i] <- all(abs(candidates[i] - c(1, knots, n)) >= (l+1))
        }
        candidates <- candidates[choose > 0]
        lcandidate <- length(candidates)
        if(lcandidate > 0){
            new <- candidates[sample(1:lcandidate, 1)]
            knots.new <- c(knots[-which(knots==old)], new)
            knots.new <- sort(knots.new)
            dm.new <- getDesignMatrix(x, x[knots.new], l, l0)
            rlm.new <- fitLM(y, dm.new)
            rss.new <- rlm.new$rss
            alpha <- min(1, ((2*gamma.rate + rss)/(2*gamma.rate + rss.new))^(n/2 + gamma.shape))
                         
            if(runif(1) < alpha){
                knots <- knots.new
                rss <- rss.new
                fitted <- rlm.new$fitted
                betas <- rlm.new$betas
            }
        }
    }
    list(knots=knots, rss=rss, fitted=fitted, betas=betas)
}

curve.polynomial.rjmcmc <- function(y, x, lambda, l, l0, c=0.4,
                                    gamma.shape=1e-3, gamma.rate=1e-3,
                                    maxit=10000, err=1e-8, verbose=TRUE){

    if(length(y) != length(x))
        stop("The length of 'y' has to be the same as that of 'x'.")
    if(l0 > l)
        stop("'l0' has to be less than or equal to 'l'.")
    if(c > 0.5)
        stop("'c' has be less than or equal to 0.5.")
    
    ox <- order(x)
    x <- x[ox]
    y <- y[ox]
    
    #step 1
    #choose lambda knots locations uniformly along the range (between two
    #boundaries) at least l+1 points away from each other
    n <- length(x)
    knots <- getInitialKnots(n, lambda, l)
    k <- length(knots)
    dm <- getDesignMatrix(x, x[knots], l, l0)
    rlm <- fitLM(y, dm)
    rss.old <- rlm$rss
    fitted.save <- matrix(0, nrow=n, ncol=maxit+1)
    fitted.save[,1] <- rlm$fitted
    betas.save <- list(rlm$betas)
    knots.save <- list(knots)
    sigma2.save <- rss.old / (n-(l+1)-k*(l-l0+1))
    iter <- 1
    flag <- 0
    while(flag == 0){
        iter <- iter + 1
        #step 2
        #set k equal to the number of interior knots in the present model
        k <- length(knots.save[iter-1])

        #step 3
        #generate u uniformly on [0,1]
        u <- runif(1)
    
        #step 4
        p.kplus1 <- dpois(k+1, lambda)
        p.k <- dpois(k, lambda)
        bk <- c * min(1, p.kplus1 / p.k)
        dk <- c * min(1, p.k / p.kplus1)
        if(u <= bk){
            newmove <- birthKnot(y, x, n, knots.save[[iter-1]], l, l0,
                             sigma2.save[iter-1], rss.old,
                             fitted.save[,iter-1], betas.save[[iter-1]],
                             gamma.shape, gamma.rate)
        }
        else{
            if(u <= bk + dk){
                newmove <- deathKnot(y, x, n, knots.save[[iter-1]], l, l0,
                                 sigma2.save[iter-1], rss.old,
                                 fitted.save[,iter-1], betas.save[[iter-1]],
                                 gamma.shape, gamma.rate)
           
            }
            else{
                newmove <- moveKnot(y, x, n, knots.save[[iter-1]], l, l0,
                                sigma2.save[iter-1], rss.old,
                                fitted.save[,iter-1], betas.save[[iter-1]],
                                gamma.shape, gamma.rate)
            }
        }
        knots.save <- c(knots.save, list(newmove$knots))
        rss.new <- newmove$rss
        fitted.save[,iter] <- newmove$fitted
        betas.save <- c(betas.save, list(newmove$betas))
        
       
        #step 5
        #draw sigma2
        sigma2.new <- 1 / rgamma(1, shape=(n/2 - 2) + gamma.shape,
                                rate=gamma.rate + rss.new/2)
        sigma2.save <- c(sigma2.save, sigma2.new)
        #step 6
        #check stopping rule
        if(rss.new!=rss.old && abs(rss.new - rss.old) / rss.old < err){
            cat(paste("The mean-squared error converged at iteration ",
                      iter, ".\n", sep=""))
            flag <- 1
        }
        rss.old <- rss.new
        
        if (verbose && iter %% 1000 == 0){
            cat(paste(iter, " iterations", " have finished.\n", sep=""))
        }
        if(iter >= maxit){
            print("The iteration number exceeded the limit")
            flag <- 1
        }
    }
    list(knots.save=knots.save, betas.save=betas.save,
         fitted.save=fitted.save[,1:iter], sigma2.save=sigma2.save)
    
}
