dComovVol <-
function(z, omega.m = 0.27, omega.lambda = 0.73, H.0 = 71) {
    D.H <- 3.e5/H.0
    omega.k <- 1 - omega.m - omega.lambda
    E <- function(x) sqrt(omega.m*(1+x)^3 + omega.k*(1+x)^2 + omega.lambda)
    dcomovVol <- D.H*(1+z)^2*
        D.A(z, omega.m=omega.m, omega.lambda=omega.lambda, H.0=H.0)^2/E(z)
    return(dcomovVol)
}

