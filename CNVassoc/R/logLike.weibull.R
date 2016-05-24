logLike.weibull <-
function (param, y, cens, X, w, variant)
{
    J <- NCOL(w)
    K <- sum(ifelse(variant, J, 1))
    beta <- vector2matrix(betav = param[1:K], variant = variant, J = J)
    alpha <- param[K + 1]
    eta <- sapply(1:J, function(j) {
        Xj <- cbind(X[, , j])
        betaj <- beta[, j, drop = FALSE]
        Xj %*% betaj
    })
    lambda <- exp(eta) ##
    scale <- lambda^(-1/alpha) ##
    ff <- lambda
    ff[cens == 1, ] <- sapply(1:J, function(j) dweibull(y[cens == 1], alpha, scale[cens == 1 ,j])) ## not censored
    ff[cens == 0, ] <- sapply(1:J, function(j) pweibull(y[cens == 0], alpha, scale[cens == 0 ,j], lower.tail = FALSE)) ## censored
    gg <- rowSums(w * ff)
    return(sum(log(gg)))
}
