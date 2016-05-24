betaIV <- function(L, V, X, Y, kend, k, intercept) {
    Vm <- matrix(V[(kend + 1):(nrow(V) - 1), -((kend + 1):(kend + k))],
                 nrow = nrow(V) - kend - 1)
    Swx <- Vm[, -ncol(Vm)]
    Sxw <- t(Swx)
    Sww <- matrix(V[(kend + 1):(nrow(V) - 1), (kend + 1):(nrow(V) - 1)],
                  nrow = nrow(V) - kend - 1)
    Swy <- Vm[, ncol(Vm)]
    Lm <- L[-((kend + 1):(kend + k))]
    Mx <- Lm[1:(length(Lm) - 1)]
    My <- Lm[length(Lm)]
    
    part1 <- Sxw %*% solve(Sww)
    b1 <- solve(part1 %*% Swx) %*% part1 %*% Swy
    biv <- b1
    resid <- Y - X %*% biv
    
    if (intercept) {
      b0 <- My - sum(b1 * Mx)
      biv <- matrix(rbind(b0, b1), ncol = 1, dimnames = NULL)
      resid <- Y - b0 - X %*% b1
    }
    
    list(beta=biv, resid=resid)
}
