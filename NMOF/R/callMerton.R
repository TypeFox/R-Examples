## function call = callMerton(S,X,tau,r,q,v,lambda,muJ,vJ,N)
## % callMerton.m -- version 2010-10-24
## % S      = spot
## % X      = strike
## % tau    = time to mat
## % r      = riskfree rate
## % q      = dividend yield
## % v      = variance (volatility squared)
## % lambda = intensity of poisson process
## % muJ    = mean jump size
## % vJ     = variance of jump process
## % N      = number of jumps to be included in sum
## lambda2  = lambda*(1+muJ); call = 0;
## for n=0:N
##     v_n = v + n*vJ/tau;
##     r_n = r - lambda*muJ+ n*log(1+muJ)/tau;
##     call = call + ( exp(-lambda2*tau) * (lambda2*tau)^n ) * ...
##         callBSM(S,X,tau,r_n,q,v_n)/ exp( sum(log(1:n)) );
## end

callMerton <- function(S, X, tau, r, q, v, lambda, muJ, vJ, N,
                       implVol = FALSE) {

    callBSM <- function(S,X,tau,r,q,v) {
        d1 <- (log(S/X) + (r - q + v / 2)*tau) / (sqrt(v)*sqrt(tau))
        d2 <- d1 - sqrt(v)*sqrt(tau)
        S * exp(-q * tau) * pnorm(d1) -  X * exp(-r * tau) * pnorm(d2)
    }

    lambda2  <- lambda * (1 + muJ)
    result <- 0
    
    for (n in 0:N) {
        v_n <- v + n*vJ/tau
        r_n <- r - lambda * muJ + n * log(1 + muJ)/tau
        result <- result + ( exp(-lambda2 * tau) * (lambda2 * tau)^n ) * 
            callBSM(S,X,tau,r_n,q,v_n)/factorial(n)
    }

    if (implVol) {
        diffPrice <- function(vol,call,S,X,tau,r,q){
            d1 <- (log(S/X)+(r - q + vol^2/2)*tau)/(vol*sqrt(tau))
            d2 <- d1 - vol*sqrt(tau)
            callBSM <- S * exp(-q * tau) * pnorm(d1) -
                X * exp(-r * tau) * pnorm(d2)
            call - callBSM
        }
        impliedVol <- uniroot(diffPrice, interval = c(0.0001, 2),
                              call = result, S = S, X = X,
                              tau = tau, r = r, q = q)[[1L]]
        result <- list(value = result, impliedVol = impliedVol)
    }
    result
}

if (FALSE) {

    S <- 100
    X <- 100
    tau <- 1
    r <- 0.0075
    q <- 0.00
    v <- 0.2^2
    lambda <- 0.0
    muJ    <- -0.2
    vJ     <- 0.6^2
    N      <- 100

    ## require("compiler")
    ## system.time(for (i in 1:200)
    ##             callMerton(S, X, tau, r, q, v, lambda, muJ, vJ, N, implVol = FALSE))
    ## callMerton2 <- cmpfun(callMerton, options = list(optimize = 3))
    ## system.time(for (i in 1:200)
    ##             callMerton2(S, X, tau, r, q, v, lambda, muJ, vJ, N, implVol = FALSE))
    callMerton(S, X, tau, r, q, v, lambda, muJ, vJ, N, implVol = TRUE)
    callCF(cf = cfMerton, S = S, X = X, tau = tau, r = r, q = q,
           v = v, lambda = lambda, muJ = muJ, vJ = vJ, implVol = FALSE)
    vanillaOptionEuropean(S,X,tau,r,q,v, greeks = FALSE)
}
