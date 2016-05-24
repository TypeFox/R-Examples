logLike.poisson <-
function (param, y, X, w, variant)
{
    J <- NCOL(w)
    beta <- vector2matrix(betav = param, variant = variant, J = J)
    eta <- sapply(1:J, function(j) {
        Xj <- cbind(X[, , j])
        betaj <- beta[, j, drop = FALSE]
        Xj %*% betaj
    })
    lambda <- exp(eta)
    ff <- sapply(1:J, function(j) dpois(y, lambda[, j]))
    gg <- rowSums(w * ff)
    return(sum(log(gg)))
}
