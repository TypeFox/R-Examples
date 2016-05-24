Local_LL_all <- function(x, w, phi){
    n <- length(x)
    dx <- diff(x)
    ll <- sum(w * phi) - sum(dx * J00(phi[1:(n - 1)], phi[2:n]))
    grad <- matrix(w, ncol = 1)
    grad[1:(n - 1)] <- grad[1:(n - 1)] - (dx * J10(phi[1:(n - 1)], phi[2:n]))
    grad[2:n] <- grad[2:n] - (dx * J10(phi[2:n], phi[1:(n - 1)]))
    tmp <- c(dx * J20(phi[1:(n - 1)], phi[2:n]), 0) + c(0, dx * J20(phi[2:n], phi[1:(n - 1)]))
    tmp <- tmp + mean(tmp) * 1e-12
    mhess2 <- matrix(0, nrow = n, ncol = n)
    mhess3 <- mhess2
    mhess1 <- tmp
    tmp <- c(0, dx * J11(phi[1:(n - 1)], phi[2:n]))
    tmp.up <- diag(tmp[2:n], nrow = n - 1, ncol = n - 1)
    mhess2[1:(n - 1), 2:n] <- tmp.up
    mhess3[2:n, 1:(n - 1)] <- diag(tmp[2:n], nrow = n - 1, ncol = n - 1)
    mhess <- diag(mhess1) + mhess2 + mhess3
    phi_new <- phi + solve(mhess) %*% grad
    dirderiv <- t(grad) %*% (phi_new - phi)
    return(list(ll = ll, phi_new = phi_new, dirderiv = dirderiv))
}
