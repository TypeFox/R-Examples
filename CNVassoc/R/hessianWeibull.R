hessianWeibull <-
function (beta, alpha, y, cens, w, X, variant)
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
    lambda <- exp(eta)
    h <- w * ifelsem(cens == 1, lambda * alpha * y^(alpha-1) * exp(-lambda * y^alpha), exp(-lambda * y^alpha))
    g <- rowSums(h)
    dlogP <- function(k, which.param) {
        vari <- if (which.param == "beta") variant[k] else FALSE
        dhk <- dh(k, which.param = which.param)
        if (vari)
            colSums(dhk/g)
        else
            sum(rowSums(dhk)/g)
    }
    dh <- function(k, which.param) {
        if (which.param == "beta")
            res <- h * X[, k, ] * ifelsem(cens == 1, 1 - lambda * y^alpha, - y^alpha * lambda)
        if (which.param == "alpha")
            res <- - h * ifelsem(cens == 1, lambda * y^alpha * log(y) - log(y) - 1/alpha, lambda * y^alpha * log(y))
        res
    }
    S <- NULL
    for (k in 1:K) {
        Sk <- dlogP(k, "beta")
        names(Sk) <- rep(paste("var", k, sep = ""), length(Sk))
        S <- c(S, Sk)
    }
    S <- c(S, alpha = dlogP(1, "alpha"))
    d2h <- function(k, kprime, which.param1, which.param2) {
        if (which.param1 == "beta" & which.param2 == "beta")
            res <- X[, k, ] * ifelsem(cens == 1,
                dh(kprime, "beta") * (1 - lambda * y^alpha) - h * lambda * y^alpha * X[, kprime, ],
                -y^alpha * lambda * (X[, kprime, ] * h + dh(kprime,"beta")))

        if (which.param1 == "beta" & which.param2 == "alpha")
            res <- X[, k, ] * ifelsem(cens == 1,
                dh(1,"alpha") * (1 - lambda * y^alpha) - h * lambda * y^alpha * log(y),
                -lambda * y^alpha * (dh(1, "alpha") + h * log(y)))

        if (which.param1 == "alpha" & which.param2 == "alpha")
            res <- ifelsem(cens == 1,
                -dh(1, "alpha") * (y^alpha * log(y) * lambda - log(y) - 1/alpha) - h * (y^alpha * log(y)^2 * lambda + 1/alpha^2),                
                -lambda * log(y) * (dh(1, "alpha") * y^alpha + h * y^alpha * log(y)))
        res
    }
    d2logP <- function(k, kprime, which.param1, which.param2) {
        variantk <- variant[k]
        variantkprime <- variant[kprime]
        if (which.param1 == "beta" & which.param2 == "beta") {
            if (variantk & variantkprime) {
                res <- matrix(0, J, J)
                for (j in 1:J) {
                  for (jprime in 1:J) {
                    dhk <- dh(k, "beta")[, j]
                    dhkprime <- dh(kprime, "beta")[, jprime]
                    d2hkkprime <- d2h(k, kprime, "beta", "beta")[,
                      j]
                    if (j != jprime)
                      res[j, jprime] <- sum(-dhk * dhkprime/g^2)
                    if (j == jprime)
                      res[j, jprime] <- sum((d2hkkprime * g -
                        dhk * dhkprime)/g^2)
                  }
                }
            }
            if (variantk & !variantkprime) {
                res <- matrix(NA, J, 1)
                for (j in 1:J) {
                  dhk <- dh(k, "beta")[, j]
                  dhkprime <- rowSums(dh(kprime, "beta"))
                  d2hkkprime <- d2h(k, kprime, "beta", "beta")[,
                    j]
                  res[j, 1] <- sum((d2hkkprime * g - dhk * dhkprime)/g^2)
                }
                res <- cbind(res)
            }
            if (!variantk & variantkprime) {
                res <- matrix(NA, 1, J)
                for (j in 1:J) {
                  dhk <- rowSums(dh(k, "beta"))
                  dhkprime <- dh(kprime, "beta")[, j]
                  d2hkkprime <- d2h(k, kprime, "beta", "beta")[,
                    j]
                  res[1, j] <- sum((d2hkkprime * g - dhk * dhkprime)/g^2)
                }
                res <- rbind(res)
            }
            if (!variantk & !variantkprime) {
                dhk <- rowSums(dh(k, "beta"))
                dhkprime <- rowSums(dh(kprime, "beta"))
                d2hkkprime <- rowSums(d2h(k, kprime, "beta",
                  "beta"))
                res <- cbind(sum((d2hkkprime * g - dhk * dhkprime)/g^2))
            }
        }
        if (which.param1 == "beta" & which.param2 == "alpha") {
            if (variantk) {
                res <- matrix(0, J, 1)
                for (j in 1:J) {
                  dhbeta <- dh(k, "beta")[, j]
                  dhalpha <- rowSums(dh(1, "alpha"))
                  d2hbetaalpha <- d2h(k, 1, "beta", "alpha")[,
                    j]
                  res[j, 1] <- sum((d2hbetaalpha * g - dhbeta *
                    dhalpha)/g^2)
                }
            }
            if (!variantk) {
                dhbeta <- rowSums(dh(k, "beta"))
                dhalpha <- rowSums(dh(1, "alpha"))
                d2hbetaalpha <- rowSums(d2h(k, kprime, "beta",
                  "alpha"))
                res <- sum((d2hbetaalpha * g - dhbeta * dhalpha)/g^2)
            }
        }
        if (which.param1 == "alpha" & which.param2 == "alpha") {
            dhalpha <- rowSums(dh(1, "alpha"))
            d2halpha2 <- rowSums(d2h(k, kprime, "alpha", "alpha"))
            res <- sum((d2halpha2 * g - dhalpha^2)/g^2)
        }
        res
    }
    Hbeta <- NULL
    for (k in 1:K) {
        Hfilak <- NULL
        for (kprime in 1:K) {
            res <- d2logP(k, kprime, "beta", "beta")
            colnames(res) <- rep(paste("var", kprime, sep = ""),
                NCOL(res))
            rownames(res) <- rep(paste("var", k, sep = ""), NROW(res))
            Hfilak <- cbind(Hfilak, res)
        }
        Hbeta <- rbind(Hbeta, Hfilak)
    }
    Hbetaalpha <- NULL
    for (k in 1:K) Hbetaalpha <- rbind(Hbetaalpha, d2logP(k, 1, "beta", "alpha"))
    Halpha2 <- d2logP(1, 1, "alpha", "alpha")
    H <- rbind(cbind(Hbeta, Hbetaalpha), cbind(t(Hbetaalpha), Halpha2))
    colnames(H)[ncol(H)] <- "alpha"
    rownames(H)[nrow(H)] <- "alpha"
    return(list(S = S, H = H))
}





# fit amb EMWeibull
#beta=fit$beta
#alpha=fit$alpha
#y=dsim$resp
#cens=dsim$cens
#w=pp
#X=X
#variant=fit$variant
#param <- matrix2vector(fit$beta, fit$variant)
#param <- c(param, fit$alpha)
#logLike.weibull(param, y, cens, X, w, variant)
#oo<-optim(par=param,logLike.weibull, y=y, cens=cens, X=X, w=w, variant=variant,hessian=TRUE,control=list(fnscale=-1))









