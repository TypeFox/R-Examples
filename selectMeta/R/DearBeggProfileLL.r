DearBeggProfileLL <- function(z, res0, lam, conf.level = 0.95, maxiter = 500){

    y <- res0$y
    u <- res0$u
    teststat <- abs(y) / u
    n <- length(y)
    k <- 1 + floor(n / 2)
    
    c <- - 0.5 * qchisq(conf.level, df = 1)
    
    ## log-likelihood l(what(theta), theta, sigmahat(theta))
    size <- 10 * (k + 1)
    inipop <- matrix(runif(size * (k + 1)), ncol = k + 1, nrow = size, byrow = TRUE)
    for (i in 1:nrow(inipop)){inipop[i, ] <- c(sort(inipop[i, 1:k]), runif(1, 0, 20))}

    d0 <- DEoptim::DEoptim(fn = DearBeggToMinimizeProfile, lower = c(rep(0, k), 0), upper = c(rep(1, k), 50), 
        control = DEoptim.control(strategy = 2, bs = FALSE, NP = size, trace = FALSE, itermax = maxiter, CR = 0.9, F = 0.8, 
        initialpop = inipop), z, y, u, lam) 
        
    w <- as.numeric(d0$optim$bestmem)[1:k] 
    sigma <- as.numeric(d0$optim$bestmem)[k + 1] 
    hij <- Hij(z, sigma, y, u, teststat)$Hij  
    l1 <- DearBeggLoglik(w, z, sigma, y, u, hij, lam)$LL

    ## loglik at maximum
    hij2 <- Hij(res0$theta, sigma = res0$sigma, y = res0$y, u = res0$u, teststat = abs(res0$y) / res0$u)$Hij
    l2 <- DearBeggLoglik(w = res0$w, theta = res0$theta, sigma = res0$sigma, y = res0$y, u = res0$u, hij = hij2, lam = lam)$LL
    
    res <- l1 - l2 - c
    return(res)
}
