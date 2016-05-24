## -*- truncate-lines: t; -*-

## - option pricing
## - testFunctions
## - restartOpt

test.callHestoncf <- function() {
    S <- 100; X <- 100; tau <- 1; r <- 0.02; q <- 0.01;
    v0 <- 0.2^2; vT <- 0.2^2
    rho <- -0.5; k <- 0.5; sigma <- 0.5
    result <- callHestoncf(S = S, X = X, tau = tau, r = r, q = q,
                           v0 = v0, vT = vT, rho = rho, k = k,
                           sigma = sigma)
    checkEquals(round(result[[1L]], 3), 7.119)

    S <- 100; X <- 100; tau <- 1; r <- 0.02; q <- 0.01; v0 <- 0.2^2; vT <- 0.2^2
    rho <- -0.5; k <- 0.5; sigma <- 0.01
    result <- callHestoncf(S = S, X = X, tau = tau, r = r, q = q,
                           v0 = v0, vT = vT, rho = rho, k = k, sigma = sigma)
    checkEquals(round(result[[1L]], 3), 8.347)

    S <- 100; X <- 90; tau <- 1; r <- 0.02; q <- 0.01; v0 <- 0.2^2; vT <- 0.2^2
    rho <- -0.5; k <- 0.5; sigma <- 1
    result <- callHestoncf(S = S, X = X, tau = tau, r = r, q = q,
                           v0 = v0, vT = vT, rho = rho, k = k, sigma = sigma)
    checkEquals(round(result[[1L]], 3), 13.362)
}



## EUROPEAN
test.EuropeanCall <- function() {
    S0 <- 10; X <- 10; r <- 0.02; tau <- 1; sigma <- 0.20; M = 101
    res <- EuropeanCall(S0 = S0, X = X, r = r, tau = tau, sigma = sigma, M = M)
    checkEquals(round(res,2), 0.89)

    S0 <- 10; X <- 6; r <- 0.02; tau <- 1; sigma <- 0.20; M = 101
    res <- EuropeanCall(S0 = S0, X = X, r = r, tau = tau, sigma = sigma, M = M)
    checkEquals(round(res,2), 4.12)

    S0 <- 10; X <- 10; r <- 0.00; tau <- 1; sigma <- 0.20; M = 101
    res <- EuropeanCall(S0 = S0, X = X, r = r, tau = tau, sigma = sigma, M = M)
    checkEquals(round(res,2), 0.80)

    S0 <- 10; X <- 10; r <- 0.02; tau <- 1/12; sigma <- 0.20; M = 101
    res <- EuropeanCall(S0 = S0, X = X, r = r, tau = tau, sigma = sigma, M = M)
    checkEquals(round(res,2), 0.24)

    S0 <- 10; X <- 10; r <- 0.02; tau <- 1/12; sigma <- 0.80; M = 101
    res <- EuropeanCall(S0 = S0, X = X, r = r, tau = tau, sigma = sigma, M = M)
    checkEquals(round(res,2), 0.93)

    S0 <- 10; X <- 10; r <- 0.02; tau <- 1/12; sigma <- 0.02; M = 101
    res <- EuropeanCall(S0 = S0, X = X, r = r, tau = tau, sigma = sigma, M = M)
    checkEquals(round(res,2), 0.03)

}
test.EuropeanCallBE <- function() {
    ## EuropeanCall and EuropeanCallBE should give the same results
    S0 <- 10; X <- 10; r <- 0.02; tau <- 1; sigma <- 0.20; M = 101
    res <- EuropeanCall(S0 = S0, X = X, r = r, tau = tau, sigma = sigma, M = M)
    res2 <- EuropeanCallBE(S0 = S0, X = X, r = r, tau = tau, sigma = sigma, M = M)
    checkEquals(res,res2)

    S0 <- 10; X <- 6; r <- 0.02; tau <- 1; sigma <- 0.20; M = 101
    res <- EuropeanCall(   S0 = S0, X = X, r = r, tau = tau, sigma = sigma, M = M)
    res2 <- EuropeanCallBE(S0 = S0, X = X, r = r, tau = tau, sigma = sigma, M = M)
    checkEquals(res,res2)

    S0 <- 10; X <- 10; r <- 0.00; tau <- 1; sigma <- 0.20; M = 101
    res <- EuropeanCall(S0 = S0, X = X, r = r, tau = tau, sigma = sigma, M = M)
    res2 <- EuropeanCallBE(S0 = S0, X = X, r = r, tau = tau, sigma = sigma, M = M)
    checkEquals(res,res2)

    S0 <- 10; X <- 10; r <- 0.02; tau <- 1/12; sigma <- 0.20; M = 101
    res <- EuropeanCall(S0 = S0, X = X,
                        r = r, tau = tau, sigma = sigma, M = M)
    res2 <- EuropeanCallBE(S0 = S0, X = X,
                           r = r, tau = tau, sigma = sigma, M = M)
    checkEquals(res,res2)

    S0 <- 10; X <- 10; r <- 0.02; tau <- 1/12; sigma <- 0.80; M = 101
    res <- EuropeanCall(S0 = S0, X = X,
                        r = r, tau = tau, sigma = sigma, M = M)
    res2 <- EuropeanCallBE(S0 = S0, X = X,
                           r = r, tau = tau, sigma = sigma, M = M)
    checkEquals(res,res2)

    S0 <- 10; X <- 10; r <- 0.02; tau <- 1/12; sigma <- 0.02; M = 101
    res <- EuropeanCall(S0 = S0, X = X,
                        r = r, tau = tau, sigma = sigma, M = M)
    res2 <- EuropeanCallBE(S0 = S0, X = X,
                           r = r, tau = tau, sigma = sigma, M = M)
    checkEquals(res,res2)
}

## TESTfunctions
test.testFunctions <- function() {
    x <- rep(0,10L)
    checkEqualsNumeric(tfAckley(x), 0)
    checkEqualsNumeric(tfGriewank(x), 0)
    checkEqualsNumeric(tfRastrigin(x), 0)

    x <- rep(1,10L)
    checkEqualsNumeric(tfRosenbrock(x), 0)

    x <- rep(420.9687, 10)
    checkEqualsNumeric(tfSchwefel(x), -418.9829*10,
                       tolerance = 1e-4)

    x <- c(-0.0244, 0.2106)
    checkEqualsNumeric(tfTrefethen(x), -3.306868,
                       tolerance = 1e-4)
}


## EUROPEAN BSM
test.vanillaOptionEuropean <- function() {
    ## PRICES
    S <- 100; X <- 100; tau <- 1; r <- 0.02; q <- 0.00; vol <- 0.3
    x <- vanillaOptionEuropean(S, X, tau, r, q, vol^2, type = "call")$value
    checkEquals(round(x,3), 12.822)

    S <- 100; X <- 100; tau <- 1; r <- 0.00; q <- 0.00; vol <- 0.3
    x <- vanillaOptionEuropean(S, X, tau, r, q, vol^2, type = "call")$value
    y <- vanillaOptionEuropean(S, X, tau, r, q, vol^2, type = "put")$value
    checkEquals(x, y)

    S <- 100; X <- 95; tau <- 0.5; r <- 0.03; q <- 0.06; vol <- 0.1
    x <- vanillaOptionEuropean(S, X, tau, r, q, vol^2, type = "put")$value
    checkEquals(round(x,3), 1.305)

    S <- 30; X <- 32; tau <- 0.5; r <- 0.03; q <- 0.00; vol <- 0.2;
    tauD <- c(0.1,0.2,0.3); D <- c(1,2,3)
    x <- vanillaOptionEuropean(S, X, tau, r, q, vol^2,
        tauD = tauD, D = D, type = "put")$value
    checkEquals(round(x,3), 7.523)

    S <- 30; X <- 32; tau <- 0.5; r <- 0.03; q <- 0.00; vol <- 0.2;
    tauD <- c(0.1,0.2,0.3,1); D <- c(1,2,3,5)
    x <- vanillaOptionEuropean(S, X, tau, r, q, vol^2,
        tauD = tauD, D = D, type = "put")$value
    checkEquals(round(x,3), 7.523)

    ## dividends
    S <- 30; X <- 30; tau <- 0.5; r <- 0.03; vol <- 0.1;
    tauD <- c(0.1,0.2,0.3); D <- c(1,2,3)
    x1 <- vanillaOptionEuropean(S, X, tau, r, q=0, vol^2,
                               tauD = tauD, D = D, type = "put")$value
    x2 <- vanillaOptionEuropean(S - sum(exp(-r*tauD)*D),
                               X, tau, r, q=0, vol^2,type = "put")$value
    tauD <- c(0.1,0.2,0.3,1); D <- c(1,2,3,20) ## div beyond expiry
    x3 <- vanillaOptionEuropean(S, X, tau, r, q=0, vol^2,
                               tauD = tauD, D = D, type = "put")$value
    checkEquals(x1,x2)
    checkEquals(x2,x3)
    checkEquals(x1,x3)
    
    ## ERRORS IN INPUTS
    ## ... q and D specified
    S <- 30; X <- 30; tau <- 0.5; r <- 0.03; q <- 0.06; vol <- 0.1;
    tauD <- c(0.1,0.2,0.3); D <- c(1,2,3)
    checkException(vanillaOptionEuropean(S, X, tau, r, q, vol^2,
            tauD = tauD, D = D, type = "put")$value, silent = TRUE)


    
    ## passing VOL instead of VARIANCE
    S <- 100; X <- 100; tau <- 1; r <- 0.05; q <- 0.072
    v <- 0.22^2  ## variance, not volatility
    vol <- 0.22
    p1 <- vanillaOptionEuropean(S=S, X = X, tau=tau, r=r, q=q, v=v,     ## with variance
                          type = "call", greeks = FALSE) 
    p2 <- vanillaOptionEuropean(S=S, X = X, tau=tau, r=r, q=q, vol=vol, ## with vol
                          type = "call", greeks = FALSE)
    p3 <- vanillaOptionEuropean(S=S, X = X, tau=tau, r=r, q=q, vol=vol, ## vol ignored
                                type = "call", greeks = FALSE, v = 0.2^2)
    p4 <- vanillaOptionEuropean(S=S, X = X, tau=tau, r=r, q=q,
                                type = "call", greeks = FALSE, v = 0.2^2)
    checkEquals(p1,p2)
    checkEquals(p3,p4)

    
    ## GREEXS
    S <- 30; X <- 30; tau <- 0.5; r <- 0.03; q <- 0.06; vol <- 0.1;
    x <- vanillaOptionEuropean(S, X, tau, r, q, vol^2, type = "put")
    checkEqualsNumeric(round(x$delta,6),  round(-0.55330738389122,6))
    checkEqualsNumeric(round(x$gamma,6),  round(0.1796756,6))
    checkEqualsNumeric(round(x$vega,6),   round(8.085402,6))
    checkEqualsNumeric(round(x$theta,6),  round(-1.274542,6))
    checkEqualsNumeric(round(x$rho,6),    round(-8.83253,6))
    checkEqualsNumeric(round(x$rhoDiv,6), round(8.299611,6))

    ## --- non-BSM ---
    S <- 100; X <- 100; tau <- 1; r <- 0.02; q <- 0.08;
    v0 <- 0.2^2; vT <- 0.2^2            ## variance, not volatility
    rho <- -0.3; k <- 0.2; sigma <- 0.3 ## stoch. vol: Heston/Bates
    temp1 <- callCF(cf = cfHeston, S = S, X = X, tau = tau,
                    r = r, q = q, v0 = v0, vT = vT, rho = rho,
                    k = k, sigma = sigma, implVol = FALSE)
    temp2 <- vanillaOptionEuropean(S = 100, X = 100, tau, r, q = q,
                                   tauD = 0, D = 0, type = "call",
                                   greeks = FALSE, model = "heston",
                                   v0 = v0, vT = vT, rho = rho, k = k,
                                   sigma = sigma)
    checkEqualsNumeric(temp1, temp2)

    ## forward difference
    fd <- function(cf, S, X, tau, r, q, ..., implVol = FALSE,
                   uniroot.control = list(), 
                   uniroot.info = FALSE, what = "S") {
        h <- 1e-6
        r1 <- r
        S1 <- S
        if (what == "r")
            r1 <- r + h
        if (what == "S")
            S1 <- S + h
        
        (callCF(cf = cf, S=S1, X=X, tau=tau, r=r1, q = q, v0 =
                v0, vT = vT, rho = rho, k = k, sigma = sigma,
                implVol = FALSE) -
         callCF(cf = cf, S=S, X=X, tau=tau,
                r=r, q = q, v0 = v0, vT = vT, rho = rho, k = k,
                sigma = sigma, implVol = FALSE))/h
    }

    for (X in seq(80, 120, by = 5))
        for (S in seq(80, 120, by = 5))
            for (tau in seq(0.5, 3, by = 0.1)) {
                D1 <- vanillaOptionEuropean(S = S, X = X, tau = tau,
                                            r = r, q = q,
                                            tauD = 0, D = 0, type = "call",
                                            greeks = TRUE, model = "heston",
                                            v0 = v0, vT = vT, rho = rho, k = k,
                                            sigma = sigma)$delta
                D2 <- fd(cf = cfHeston, S = S, X = X, tau = tau,
                         r = r, q = q,
                         v0 = v0, vT = vT, rho = rho, k = k,
                         sigma = sigma, implVol = FALSE, what = "S")
                checkTrue(abs(D1-D2) < 0.015)
            }    
}


## AMERICAN BSM
test.vanillaOptionAmerican <- function() {
    # PRICES
    S <- 100; X <- 100; tau <- 1; r <- 0.02; q <- 0.00; vol <- 0.3
    x <- vanillaOptionAmerican(S, X, tau, r, q, vol^2,
        type = "call", M = 101)$value
    checkEquals(round(x,3), 12.850)

    S <- 100; X <- 100; tau <- 1; r <- 0.00; q <- 0.00; vol <- 0.3
    x <- vanillaOptionAmerican(S, X, tau, r, q, vol^2,
        type = "call", M = 101)$value
    y <- vanillaOptionAmerican(S, X, tau, r, q, vol^2,
        type = "put", M = 101)$value
    checkEquals(x, y)

    S <- 100; X <- 95; tau <- 0.5; r <- 0.03; q <- 0.06; vol <- 0.1
    x <- vanillaOptionAmerican(S, X, tau, r, q, vol^2,
        type = "put", M = 101)$value
    checkEquals(round(x,3), 1.303)

    S <- 30; X <- 32; tau <- 0.5; r <- 0.03; q <- 0.00; vol <- 0.2;
    tauD <- c(0.1,0.2,0.3); D <- c(1,2,3)
    x <- vanillaOptionAmerican(S, X, tau, r, q, vol^2,
        tauD = tauD, D = D, type = "call")$value
    checkEquals(round(x,3), 0.153)

    S <- 30; X <- 32; tau <- 0.5; r <- 0.03; q <- 0.00; vol <- 0.2;
    tauD <- c(0.1,0.2,0.3,1); D <- c(1,2,3,5)
    x <- vanillaOptionAmerican(S, X, tau, r, q, vol^2,
        tauD = tauD, D = D, type = "call")$value
    checkEquals(round(x,3), 0.153)

    ## ERRORS IN INPUTS
    ## ... q and D specified
    S <- 30; X <- 30; tau <- 0.5; r <- 0.03; q <- 0.06; vol <- 0.1;
    tauD <- c(0.1,0.2,0.3); D <- c(1,2,3)
    checkException(vanillaOptionAmerican(S, X, tau, r, q, vol^2,
            tauD = tauD, D = D, type = "put")$value, silent = TRUE)

    # GREEKS
    # TODO...
}


test.vanillaOptionImpliedVol <- function() {
    # EUROPEAN
    S <- 100; X <- 100; tau <- 1; r <- 0.02; q <- 0.00; vol <- 0.3
    p <- vanillaOptionEuropean(S, X, tau, r, q, vol^2, type = "call")$value
    ivol <- vanillaOptionImpliedVol(exercise = "european", price = p,
        S = S, X = X, tau = tau, r = r, q = q, tauD = 0, D = 0,
        type = "call")
    checkEquals(round(ivol,4), vol)

    S <- 100; X <- 100; tau <- 1; r <- 0.00; q <- 0.00; vol <- 0.3
    p <- vanillaOptionEuropean(S, X, tau, r, q, vol^2, type = "call")$value
    ivol <- vanillaOptionImpliedVol(exercise = "european", price = p,
        S = S, X = X, tau = tau, r = r, q = q, tauD = 0, D = 0,
        type = "call")
    checkEquals(round(ivol,4), vol)

    S <- 100; X <- 100; tau <- 1; r <- 0.00; q <- 0.00; vol <- 0.3
    p <- vanillaOptionEuropean(S, X, tau, r, q, vol^2, type = "put")$value
    ivol <- vanillaOptionImpliedVol(exercise = "european", price = p,
        S = S, X = X, tau = tau, r = r, q = q, tauD = 0, D = 0,
        type = "put")
    checkEquals(round(ivol,4), vol)

    S <- 100; X <- 95; tau <- 0.5; r <- 0.03; q <- 0.06; vol <- 0.1
    p <- vanillaOptionEuropean(S, X, tau, r, q, vol^2, type = "put")$value
    ivol <- vanillaOptionImpliedVol(exercise = "european", price = p,
        S = S, X = X, tau = tau, r = r, q = q, tauD = 0, D = 0,
        type = "put")
    checkEquals(round(ivol,4), vol)

    S <- 30; X <- 32; tau <- 0.5; r <- 0.03; q <- 0.00; vol <- 0.2;
    tauD <- c(0.1,0.2,0.3); D <- c(1,2,3)
    p <- vanillaOptionEuropean(S, X, tau, r, q, vol^2,
        tauD = tauD, D = D, type = "put")$value
    ivol <- vanillaOptionImpliedVol(exercise = "european", price = p,
        S = S, X = X, tau = tau, r = r, q = q, tauD = tauD, D = D,
        type = "put")
    checkEquals(round(ivol,4), vol)

    S <- 30; X <- 32; tau <- 0.5; r <- 0.03; q <- 0.00; vol <- 0.2;
    tauD <- c(0.1,0.2,0.3,1); D <- c(1,2,3,5)
    p <- vanillaOptionEuropean(S, X, tau, r, q, vol^2,
        tauD = tauD, D = D, type = "put")$value
    ivol <- vanillaOptionImpliedVol(exercise = "european", price = p,
        S = S, X = X, tau = tau, r = r, q = q, tauD = tauD, D = D,
        type = "put")
    checkEquals(round(ivol,4), vol)

    # AMERICAN
    S <- 100; X <- 100; tau <- 1; r <- 0.02; q <- 0.00; vol <- 0.3
    p <- vanillaOptionAmerican(S, X, tau, r, q, vol^2,
        type = "call", M = 101)$value
    ivol <- vanillaOptionImpliedVol(exercise = "american", price = p,
        S = S, X = X, tau = tau, r = r, q = q, tauD = 0, D = 0,
        type = "call", uniroot.control = list(interval = c(0.01, 0.8)))
    checkEquals(round(ivol,4), vol)

    S <- 100; X <- 100; tau <- 1; r <- 0.00; q <- 0.00; vol <- 0.3
    p <- vanillaOptionAmerican(S, X, tau, r, q, vol^2,
        type = "call", M = 101)$value
    ivol <- vanillaOptionImpliedVol(exercise = "american", price = p,
        S = S, X = X, tau = tau, r = r, q = q, tauD = 0, D = 0,
        type = "call")
    checkEquals(round(ivol,4), vol)

    S <- 100; X <- 100; tau <- 1; r <- 0.00; q <- 0.00; vol <- 0.3
    p <- vanillaOptionAmerican(S, X, tau, r, q, vol^2,
        type = "put", M = 101)$value
    ivol <- vanillaOptionImpliedVol(exercise = "american", price = p,
        S = S, X = X, tau = tau, r = r, q = q, tauD = 0, D = 0,
        type = "put",uniroot.control = list(interval = c(0.01, 0.5)))
    checkEquals(round(ivol,4), vol)

    S <- 100; X <- 95; tau <- 0.5; r <- 0.03; q <- 0.06; vol <- 0.1
    p <- vanillaOptionAmerican(S, X, tau, r, q, vol^2,
        type = "put", M = 101)$value
    ivol <- vanillaOptionImpliedVol(exercise = "american", price = p,
        S = S, X = X, tau = tau, r = r, q = q, tauD = 0, D = 0,
        type = "put")
    checkEquals(round(ivol,4), vol)

    S <- 30; X <- 32; tau <- 0.5; r <- 0.03; q <- 0.00; vol <- 0.2;
    tauD <- c(0.1,0.2,0.3); D <- c(1,2,3)
    p <- vanillaOptionAmerican(S, X, tau, r, q, vol^2,
        tauD = tauD, D = D, type = "call")$value
    ivol <- vanillaOptionImpliedVol(exercise = "american", price = p,
        S = S, X = X, tau = tau, r = r, q = q, tauD = tauD, D = D,
        type = "call")
    checkEquals(round(ivol,4), vol)

    S <- 30; X <- 32; tau <- 0.5; r <- 0.03; q <- 0.00; vol <- 0.2;
    tauD <- c(0.1,0.2,0.3,1); D <- c(1,2,3,5)
    p <- vanillaOptionAmerican(S, X, tau, r, q, vol^2,
        tauD = tauD, D = D, type = "call")$value
    ivol <- vanillaOptionImpliedVol(exercise = "american", price = p,
        S = S, X = X, tau = tau, r = r, q = q, tauD = tauD, D = D,
        type = "call")
    checkEquals(round(ivol,4), vol)
}
