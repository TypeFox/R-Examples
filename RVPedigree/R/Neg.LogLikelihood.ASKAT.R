Neg.LogLikelihood.ASKAT <- function(delta, S, Ut.y, Ut.x, n) {
    W    <- diag(1/(delta+S))
    beta <- solve(t(Ut.x) %*% W %*% Ut.x) %*% t(Ut.x) %*% W %*% Ut.y
    s.g  <- mean((Ut.y-Ut.x %*% beta)^2 / (delta+S))
    LL   <- n * log(s.g) + sum(log(S+delta))
    return(LL)
}
