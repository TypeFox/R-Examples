callCF <- function(cf, S, X, tau, r, q = 0, ..., implVol = FALSE,
                   uniroot.control = list(), uniroot.info = FALSE) {
    ucon <- list(interval = c(1e-5, 2),
                 tol = .Machine$double.eps^0.25,
                 maxiter = 1000)
    ucon[names(uniroot.control)] <- uniroot.control

    P1 <- function(om, S, X, tau, r, q, ...)
        Re(exp(-1i * log(X) * om) * cf(om - 1i, S, tau, r, q, ...) /
                (1i * om * S * exp((r-q) * tau)))
    P2 <- function(om,S,X,tau,r,q, ...)
        Re(exp(-1i * log(X) * om) * cf(om, S, tau, r, q, ...) /
                (1i * om))
    vP1 <- 0.5 + 1/pi * integrate(P1, lower = 1e-8, upper = Inf,
                                  S, X, tau, r, q, ...)$value
    vP2 <- 0.5 + 1/pi * integrate(P2, lower = 1e-8, upper = Inf,
                                  S, X, tau, r, q, ...)$value
    result <- exp(-q * tau) * S * vP1 - exp(-r * tau) * X * vP2;

    if (implVol) {
        diffPrice <- function(vol, call, S, X, tau, r, q) {
            d1 <- (log(S/X) + (r - q + vol^2/2)*tau)/(vol*sqrt(tau))
            d2 <- d1 - vol*sqrt(tau)
            callBSM <- S * exp(-q * tau) * pnorm(d1) -
                X * exp(-r * tau) * pnorm(d2)
            call - callBSM
        }
        impliedVol <- uniroot(diffPrice,
                              interval = ucon$interval,
                              call = result, S = S, X = X,
                              tau = tau, r = r, q = q,
                              tol = ucon$tol,
                              maxiter = ucon$maxiter)
        result <- list(value = result, impliedVol = impliedVol[[1L]])
        if (uniroot.info)
            result <- c(result, uniroot = impliedVol)
    }
    result
}

cfHeston <- function(om, S, tau, r, q, v0, vT, rho, k, sigma) {
    if (sigma < 1e-8)
        sigma <- 1e-8
    d <- sqrt((rho * sigma * 1i * om - k)^2 + sigma^2 *
              (1i * om + om ^ 2))
    g <- (k - rho * sigma * 1i * om - d) /
        (k - rho * sigma * 1i * om + d)
    cf1 <- 1i * om * (log(S) + (r - q) * tau)
    cf2 <- vT*k/(sigma^2)*((k - rho * sigma * 1i * om - d) *
            tau - 2 * log((1 - g * exp(-d * tau)) / (1 - g)))
    cf3 <- v0 / sigma^2 * (k - rho * sigma * 1i * om - d) *
        (1 - exp(-d * tau)) / (1 - g * exp(-d * tau))
    exp(cf1 + cf2 + cf3)
}

cfBSM <- function(om, S, tau, r, q, v)
    exp(1i * om * log(S) + 1i * tau * (r - q) * om -
            0.5 * tau * v * (1i * om + om ^ 2))

cfBates <- function(om, S, tau, r, q,
                    v0, vT, rho, k, sigma, lambda, muJ, vJ) {
    if (sigma < 1e-8)
        sigma <- 1e-8
    sigma <- max(sigma,1e-4)
    om1i <- om * 1i
    d <- sqrt( (rho*sigma*om1i - k)^2 + sigma^2 * (om1i + om^2) )
    g <- (k - rho*sigma*om1i - d) / (k - rho*sigma*om1i + d)
    cf1 <- om1i * (log(S) + (r - q) * tau)
    cf2 <- vT*k / (sigma^2) * ((k - rho*sigma*om1i - d) * tau -
                               2 * log((1 - g * exp(-d * tau)) / (1 - g)))
    cf3 <- v0/sigma^2*(k - rho*sigma*om1i - d)*(1 - exp(-d*tau)) /
        (1-g*exp(-d*tau))
    cf4 <- -lambda * muJ * om1i * tau + lambda * tau *
        ((1+muJ)^(om1i) * exp( vJ*(om1i/2) * (om1i-1) )-1)
    exp(cf1 + cf2 + cf3 + cf4)
}

cfMerton <- function(om, S, tau, r, q, v, lambda, muJ, vJ) {
    om1i <- om * 1i
    cf1 <- om1i*log(S) + om1i*tau*(r-q-0.5*v-lambda*muJ) -
        0.5*(om^2)*v*tau
    cf2 <- lambda * tau * (exp(om1i * log(1 + muJ) -
                               0.5 * om1i * vJ - 0.5 * vJ * om^2) - 1)
    exp(cf1 + cf2)
}

cfVG <- function(om, S, tau, r, q, nu, theta, sigma) {
    om1i <- om * 1i
    sigma2 <- sigma^2
    w <- log(1 - theta*nu - 0.5*nu*sigma2)/nu
    temp <- om1i*log(S) + om1i*(r-q+w)*tau
    temp <- exp(temp)
    temp / ((1 - om1i*theta*nu + 0.5*sigma2*nu*om^2)^(tau/nu))
}
