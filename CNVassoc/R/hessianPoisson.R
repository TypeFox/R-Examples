hessianPoisson<-
function (beta, y, w, X, variant)
{
    n <- NROW(w)
    J <- NCOL(w)
    K <- NCOL(X)
    y <- sapply(1:J, function(j) y)
    eta <- sapply(1:J, function(j) {
        Xj <- cbind(X[, , j])
        betaj <- beta[, j, drop = FALSE]
        Xj %*% betaj
    })
    lambda <- exp(eta)  ##
    h <- w * exp(-exp(eta)) * exp(y * eta)

    g <- rowSums(h)
    dlogP <- function(k) {
        if (variant[k])
            colSums(dh(k)/g)
        else sum(rowSums(dh(k))/g)
    }
    dh <- function(k) h * (y - lambda) * X[, k, ]    #######
    S <- NULL
    for (k in 1:K) {
        Sk <- dlogP(k)
        names(Sk) <- rep(paste("var", k, sep = ""), length(Sk))
        S <- c(S, Sk)
    }

    d2h <- function(k, kprime, j, jprime) X[, k, j] * dh(kprime)[,
        jprime] * (y - lambda)[, j] - X[, k, j] * X[, kprime, jprime] *
        h[, j] * lambda[, jprime]                          #######
    d2logP <- function(k, kprime, variantk, variantkprime) {
        if (variantk & variantkprime) {
            res <- matrix(0, J, J)
            for (j in 1:J) {
                for (jprime in 1:J) {
                  dhk <- dh(k)[, j]
                  dhkprime <- dh(kprime)[, jprime]
                  d2hkkprime <- d2h(k, kprime, j, jprime)
                  if (j != jprime)
                    res[j, jprime] <- sum(-dhk * dhkprime/g^2)
                  if (j == jprime)
                    res[j, jprime] <- sum((d2hkkprime * g - dhk *
                      dhkprime)/g^2)
                }
            }
        }
        if (variantk & !variantkprime) {
            res <- matrix(NA, J, 1)
            for (j in 1:J) {
                dhk <- dh(k)[, j]
                dhkprime <- rowSums(dh(kprime))
                d2hkkprime <- d2h(k, kprime, j, j)
                res[j, 1] <- sum((d2hkkprime * g - dhk * dhkprime)/g^2)
            }
            res <- cbind(res)
        }
        if (!variantk & variantkprime) {
            res <- matrix(NA, 1, J)
            for (j in 1:J) {
                dhk <- rowSums(dh(k))
                dhkprime <- dh(kprime)[, j]
                d2hkkprime <- d2h(k, kprime, j, j)
                res[1, j] <- sum((d2hkkprime * g - dhk * dhkprime)/g^2)
            }
            res <- rbind(res)
        }
        if (!variantk & !variantkprime) {
            dhk <- rowSums(dh(k))
            dhkprime <- rowSums(dh(kprime))
            d2hkkprime <- rowSums(d2h(k, kprime, 1:J, 1:J))
            res <- cbind(sum((d2hkkprime * g - dhk * dhkprime)/g^2))
        }
        res
    }
    H <- NULL
    for (k in 1:K) {
        Hfilak <- NULL
        for (kprime in 1:K) {
            res <- d2logP(k, kprime, variant[k], variant[kprime])
            colnames(res) <- rep(paste("var", kprime, sep = ""),
                NCOL(res))
            rownames(res) <- rep(paste("var", k, sep = ""), NROW(res))
            Hfilak <- cbind(Hfilak, res)
        }
        H <- rbind(H, Hfilak)
    }
    return(list(S = S, H = H))
}
