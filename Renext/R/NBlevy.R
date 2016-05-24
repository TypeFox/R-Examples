##=======================================================================
## Estimation of the two parameters gamma and prob ob a Negative Binomial
## count process.
##=======================================================================

"NBlevy" <- function(N,
                     gamma = NA,
                     prob = NA,
                     w = rep(1, length(N)),
                     sum.w = sum(w),
                     interval = c(0.01, 1000),
                     optim = TRUE,
                     plot = FALSE,
                     ...) { 
  
  
    if (sum.w < sum(w))
        stop("sum(w) must be less or equatl to 'sum.w'")
    else if (sum.w > sum(w)) 
        stop("incomplete durations in 'w' are NOT IMPLEMENTED YET")
    
    n.N <- length(N)
    
    S.N <- sd(N)
    Nbar <- mean(N)
    Rate <- sum(N) / sum.w
    
    CV.w <- sd(w) / mean(w)
    
    ## I the duration are near identical
    if (CV.w < 0.05) {
        p.ini <- Nbar / S.N / S.N
        gam.ini <- Nbar * p.ini / (1 - p.ini)
    }
    
    ##==========================================================================
    ## Compute log-likelihood as well as the score vector the observed
    ## information matrix for the vector [ gamma, p]
    ## ==========================================================================
    
    "NB.loglik" <- function(gamprob) {
    
        gamw <- gamprob[1] * w
        p <- gamprob[2]
        q <- 1 - p
        
        ll <- sum(lgamma(gamw + N) - lgamma(gamw) - lgamma(N + 1) +
                      gamw * log(p) + N * log(q))
        
        score <- c(sum(w * (digamma(gamw + N) - digamma(gamw) + log(p) )),
                   sum(gamw / p - N / q))
        
        info.vec <- -c(sum(w * w * (trigamma(gamw + N) - trigamma(gamw))),
                       rep(sum(w) / p, 2),
                       - sum(gamw / p / p + N / q / q))
        
        info <- matrix(info.vec, ncol = 2, nrow = 2)
        
        list(par = gamprob,
             ll = ll,
             score = score,
             info = info)
    }
    
    ## concentrated log-likelihood
    "NB.loglik.conc" <- function(gam) {
        gamw <- gam*w
        R <- sum(N)/sum(w)
        p <- gam/(gam + R)
        ##cat(gam, p, "\n")
        sum(lgamma(gamw + N) - lgamma(gamw) - lgamma(N + 1) + gamw * log(p) +
                N*log(1 - p)) 
    }
    ## derivative of the concentrated log-likelihood
    "NB.dloglik.conc" <- function(gam) {
        gamw <- gam * w
        R <- sum(N) / sum(w)
        p <- gam / (gam + R) 
        sum(w * (digamma(gamw + N) - digamma(gamw) + log(p)) ) 
    }
    
    if (optim) {
        ## find gamma with an 1D optim...
        res <- optimize(f = NB.loglik.conc, interval = interval, maximum = TRUE)
        gamma.hat <- res$maximum
        
    } else {
        ## ... or  with a zero finding.
        res <- uniroot(f = NB.dloglik.conc, interval = interval)
        gamma.hat <- res$root
        ## print(res)
    }
    
    p.hat <- gamma.hat / (gamma.hat + Rate)
    
    Lres <- NB.loglik(gamprob = c(gamma.hat, p.hat))
    
    pnames <- c("gamma", "prob")
    
    estimate = c(gamma.hat, p.hat)
    names(estimate) <- pnames
    
    colnames(Lres$info) <- pnames
    rownames(Lres$info) <- pnames
    
    C.hat <- solve(Lres$info)
    colnames(C.hat) <- pnames
    rownames(C.hat) <- pnames
    
    rate.hat <- gamma.hat * (1 - p.hat) / p.hat
    rate.vec <- c((1 - p.hat) / p.hat, -gamma.hat / p.hat / p.hat)
    sd.rate <- drop(t(rate.vec) %*% C.hat %*% rate.vec)
    
    Res <- list(gamma = gamma.hat,
                prob = p.hat,
                estimate = estimate,
                sd   = sqrt(diag(C.hat)),
                cov = C.hat,
                rate = rate.hat,
                sd.rate = sd.rate,
                loglik = Lres$ll,
                score = Lres$score,
                info = Lres$info)
    
    if (plot) {
        
        col <- col2rgb("SteelBlue3") / 256
        col1 <- rgb(col[1], col[2], col[3], 0.5)
        col <- col2rgb("orangered") / 256
        col2 <- rgb(col[1], col[2], col[3], 0.5)
        
        Esp.N <- gamma.hat * w * (1 - p.hat) / p.hat
        mat <- cbind(N, Esp.N)
        colnames(mat) <- c("Obs", "Esp")
        rmat <- range(mat)
        
        plot(x = 1L:n.N, y = N,
             ylim = rmat,
             xlab = "",
             type = "o", lwd = 2, col = col1, bg = col1, pch = 21)
        
        points(x = 1:n.N, y = Esp.N,
               type = "o", lwd = 2, col = col2, bg = col2, pch = 22)
        
        legend("topleft",
               legend = c("Obs.", "Esp"),
               col = c(col1, col2),
               lwd = 3, lty = "solid", pch = c(21, 22))
    }
    
    Res
    
}

##================================================================
## Tests
##================================================================

if (FALSE) {
    nint <- 100
    gam <- 6
    prob <- 0.20
    
    w <- rgamma(nint, shape = 3, scale = 1/5)
    ## w <- rep(1, nint)
    N <- rnbinom(nint, size = w*gam, prob = prob)
    mu <- w*gam*(1-prob)/prob
    
    Res <- NBlevy(N = N, w = w)
}








