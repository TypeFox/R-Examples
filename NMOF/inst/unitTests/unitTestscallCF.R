## -*- truncate-lines: t; -*-
## pricing with the cf
test.callCF <- function() {

                                        # HESTON

    S <- 100; X <- 100; tau <- 1; r <- 0.02; q <- 0.08;
    v0 <- 0.2^2; vT <- 0.2^2            ## variance, not volatility
    rho <- -0.3; k <- 0.2; sigma <- 0.3 ## stoch. vol: Heston/Bates
    temp <- callCF(cf = cfHeston, S = S, X = X, tau = tau,
                   r = r, q = q, v0 = v0, vT = vT, rho = rho, k = k,
                   sigma = sigma, implVol = FALSE)
    checkEquals(round(temp,3), 4.268)

    S <- 100; X <- 100; tau <- 1; r <- 0.00; q <- 0.08;
    v0 <- 0.2^2; vT <- 0.2^2            # variance, not volatility
    rho <- -0.3; k <- 0.2; sigma <- 0.3 # stoch. vol: Heston/Bates
    temp <- callCF(cf = cfHeston, S=S, X=X, tau=tau, r=r, q = q,
                   v0 = v0, vT = vT, rho = rho, k = k,
                   sigma = sigma, implVol = FALSE)
    checkEquals(round(temp,3), 3.621)

    S <- 100; X <- 100; tau <- 1; r <- 0.05; q <- 0.00;
    v0 <- 0.2^2; vT <- 0.2^2            ## variance, not volatility
    rho <- -0.3; k <- 0.2; sigma <- 0.3 ## stoch. vol: Heston/Bates
    temp <- callCF(cf = cfHeston, S=S, X=X, tau=tau, r=r, q = q,
                   v0 = v0, vT = vT, rho = rho, k = k,
                   sigma = sigma, implVol = FALSE)
    checkEquals(round(temp,3), 10.055)

    S <- 100; X <- 90; tau <- 0.1; r <- 0.05; q <- 0.00;
    v0 <- 0.2^2; vT <- 0.2^2            ## variance, not volatility
    rho <- -0.3; k <- 0.2; sigma <- 0.3 ## stoch. vol: Heston/Bates
    temp <- callCF(cf = cfHeston, S=S, X=X, tau=tau, r=r, q = q,
                   v0 = v0, vT = vT, rho = rho, k = k,
                   sigma = sigma, implVol = FALSE)
    checkEquals(round(temp,3), 10.586)

    
                                        # BSM

    S <- 100; X <- 100; tau <- 1; r <- 0.02; q <- 0.08; v <- 0.2^2
    temp1 <- callCF(cf = cfBSM, S = S, X = X, tau = tau, r = r, q = q,
                    v = v, implVol = FALSE)
    temp2 <- vanillaOptionEuropean(S=S, X=X, tau=tau, r=r, q=q,
                                   v=v, greeks = FALSE)
    checkEquals(round(temp1,4), round(temp2, 4))

    S <- 100; X <- 100; tau <- 1; r <- 0.00; q <- 0.08; v <- 0.2^2
    temp1 <- callCF(cf = cfBSM, S = S, X = X, tau = tau, r = r, q = q,
                    v = v, implVol = FALSE)
    temp2 <- vanillaOptionEuropean(S=S, X=X, tau=tau, r=r, q=q, v=v, greeks = FALSE)
    checkEquals(round(temp1,4), round(temp2,4))

    S <- 100; X <- 100; tau <- 1; r <- 0.02; q <- 0.00; v <- 0.2^2
    temp1 <- callCF(cf = cfBSM, S = S, X = X, tau = tau, r = r, q = q,
                    v = v, implVol = FALSE)
    temp2 <- vanillaOptionEuropean(S=S, X=X, tau=tau, r=r, q=q, v=v, greeks = FALSE)
    checkEquals(round(temp1,4), round(temp2,4))

    S <- 100; X <- 100; tau <- 1; r <- 0.00; q <- 0.00; v <- 0.2^2
    temp1 <- callCF(cf = cfBSM, S = S, X = X, tau = tau, r = r, q = q,
                    v = v, implVol = FALSE)
    temp2 <- vanillaOptionEuropean(S=S, X=X, tau=tau, r=r, q=q, v=v, greeks = FALSE)
    checkEquals(round(temp1,4), round(temp2,4))

    S <- 100; X <- 90; tau <- 1/12; r <- 0.00; q <- 0.00; v <- 0.2^2
    temp1 <- callCF(cf = cfBSM, S = S, X = X, tau = tau, r = r, q = q,
                    v = v, implVol = FALSE)
    temp2 <- vanillaOptionEuropean(S=S, X=X, tau=tau, r=r, q=q, v=v, greeks = FALSE)
    checkEquals(round(temp1,4), round(temp2,4))

    S <- 100; X <- 90; tau <- 1/12; r <- 0.00; q <- 0.00; v <- 0.2^2
    temp1 <- callCF(cf = cfBSM, S = S, X = X, tau = tau, r = r, q = q,
                    v = v, implVol = TRUE)
    checkEquals(round(temp1$impliedVol^2,4), v)


                                        # BATES

    ## http://www.wilmott.com/messageview.cfm?catid=4&threadid=63014
    S <- 100; X <- 100; tau <- 1; r <- 0.04; q <- 0.02;
    v0 <- 0.2^2; vT <- 0.2^2            ## variance, not volatility
    rho <- -0.5; k <- 5; sigma <- 0.25  ## stoch. vol: Heston/Bates
    lambda <- 0.2; muJ <- 0.15; vJ <- 0.1^2 ## jumps: Merton/Bates/Heston
    temp <- callCF(cf = cfBates, S = S, X = X, tau = tau,
                   r = r, q = q,
                   v0 = v0, vT = vT, rho = rho, k = k, sigma = sigma,
                   lambda = lambda, muJ = muJ, vJ = vJ, implVol = FALSE)
    checkEquals(round(temp,5), 9.19049)

    S <- 100; X <- 110; tau <- 0.5; r <- 0.04; q <- 0.02;
    v0 <- 0.2^2; vT <- 0.2^2            ## variance, not volatility
    rho <- -0.5; k <- 5; sigma <- 0.25  ## stoch. vol: Heston/Bates
    lambda <- 0.2; muJ <- 0.15; vJ <- 0.1^2 ## jumps: Merton/Bates/Heston
    temp <- callCF(cf = cfBates, S = S, X = X, tau = tau,
                   r = r, q = q,
                   v0 = v0, vT = vT, rho = rho, k = k, sigma = sigma,
                   lambda = lambda, muJ = muJ, vJ = vJ, implVol = FALSE)
    checkEquals(round(temp,5), 2.66063)

    S <- 100; X <- 100; tau <- 0.5; r <- 0.00; q <- 0.00;
    v0 <- 0.2^2; vT <- 0.2^2            ## variance, not volatility
    rho <- -0.5; k <- 5; sigma <- 0.25  ## stoch. vol: Heston/Bates
    lambda <- 0.2; muJ <- 0.15; vJ <- 0.1^2 ## jumps: Merton/Bates/Heston
    temp <- callCF(cf = cfBates, S = S, X = X, tau = tau,
                   r = r, q = q,
                   v0 = v0, vT = vT, rho = rho, k = k, sigma = sigma,
                   lambda = lambda, muJ = muJ, vJ = vJ, implVol = TRUE)
    checkTrue(length(temp) == 2L)
    checkTrue(all(sapply(temp, is.finite)))

    
                                        # MERTON

    temp <- callCF(cf = cfMerton, S = S, X = X, tau = tau,
                   r = r, q = q, v = v, lambda = lambda,
                   muJ = muJ, vJ = vJ, implVol = FALSE)

    ## boundary cases
    S <- 100; X <- 100; tau <- 1; r <- 0.02; q <- 0.08; v <- 0.2^2
    temp1 <- callCF(cf = cfBSM, S = S, X = X, tau = tau, r = r, q = q,
                    v = v, implVol = TRUE)
    temp2 <- vanillaOptionEuropean(S=S, X=X, tau=tau, r=r, q=q,
                                   v=v, greeks = TRUE)
    v0 <- 0.2^2; vT <- 0.2^2
    rho <- 0; k <- 0; sigma <- 0
    lambda <- 0.; muJ <- 0.15; vJ <- 0.1^2
    temp3 <- callCF(cf = cfBates, S = S, X = X, tau = tau,
                   r = r, q = q,
                   v0 = v0, vT = vT, rho = rho, k = k, sigma = sigma,
                   lambda = lambda, muJ = muJ, vJ = vJ, implVol = TRUE)
    temp4 <- callCF(cf = cfHeston, S = S, X = X,
                    tau = tau, r = r, q = q,
                    v0 = v0, vT = vT, rho = rho, k = k,
                    sigma = sigma, implVol = TRUE)
    temp5 <- callCF(cf = cfMerton, S = S, X = X, tau = tau,
                    r = r, q = q, v = v, lambda = lambda,
                    muJ = muJ, vJ = vJ, implVol = TRUE)
    ## --> prices
    checkEquals(round(temp1[[1L]],5), round(temp2[[1L]],5))
    checkEquals(round(temp1[[1L]],5), round(temp3[[1L]],5))
    checkEquals(round(temp1[[1L]],5), round(temp4[[1L]],5))
    checkEquals(round(temp1[[1L]],5), round(temp5[[1L]],5))

    ## --> implied vols
    checkEquals(round(temp1[[2L]],5), round(temp3[[2L]],5))
    checkEquals(round(temp1[[2L]],5), round(temp4[[2L]],5))
    checkEquals(round(temp1[[2L]],5), round(temp5[[2L]],5))


                                        # VARIANCE GAMMA

    S <- 100; X <- 100; tau <- 1; r <- 0.1; q <- 0.0
    nu <- 0.2; theta <- -0.14; sigma <- 0.12

    temp <- callCF(cf = cfVG, S = S, X = X, tau = tau,
                   r = r, q = q,
                   nu = nu, theta=theta, sigma = sigma, implVol = FALSE)
    checkEquals(round(temp,2), 11.37)

}
