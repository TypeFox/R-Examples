intcox.derivs <-
function (data, covar, lambda0u, lambda0v, beta)        # derivatives and likelihood
{
    eps <- 10^(-10)                                     # for numerical stability
    cn <- data$cens
    cens.nr <- (1:length(cn))[cn == 0]
    ncens.nr <- (1:length(cn))[cn == 3]
    e <- exp(beta %*% t(covar))                         # often used terms
    qu <- lambda0u * e                                  #
    qu.n <- qu[cn == 3]                                 #
    qu.c <- qu[cn == 0]                                 #
    qv <- lambda0v * e[cn == 3]                         #
    equ <- exp(-qu)                                     #
    equ.n <- equ[cn == 3]                               #
    equ.c <- equ[cn == 0]                               #
    eqv <- exp(-qv)                                     #
    likeli <- sum(log(equ.n - eqv)) - sum(qu.c)         # Likelihood
    l1u.n <- -qu.n * equ.n/(equ.n - eqv)                # d L/d F0
    l1u.c <- -qu.c                                      #
    l1u <- NULL                                         #
    l1u[ncens.nr] <- l1u.n                              #
    l1u[cens.nr] <- l1u.c                               #
    l1v <- qv * eqv/(equ.n - eqv)                       #
    l1 <- c(l1u, l1v)                                   #
    temp.l2 <- (-qu.n * equ.n + qv * eqv)/(equ.n - eqv) # d L/d beta
    if (sum(3 - cn) == 0) {
        l2 <- as.vector(apply(as.matrix(temp.l2 * covar[cn ==
            3, ]), 2, sum))
    }
    else {
        l2 <- as.vector(apply(as.matrix(temp.l2 * covar[cn ==
            3, ]), 2, sum) - apply(as.matrix(qu.c * covar[cn ==
            0, ]), 2, sum))
    }
    g1u.n <- qu.n * equ.n/(equ.n - eqv) * (1 - qu.n + qu.n *
        equ.n/(equ.n - eqv))                            # diagonal elements of -d^2 L/d F0^2
    g1u.c <- qu.c
    g1u <- NULL
    g1u[ncens.nr] <- g1u.n
    g1u[cens.nr] <- g1u.c
    g1v <- -qv * eqv/(equ.n - eqv) * (1 - qv - qv * eqv/(equ.n -
        eqv))
    g1 <- c(g1u, g1v)
    g1 <- g1 + eps                                      # for numerical stability
    temp.g2 <- (-qu.n * equ.n * (1 - qu.n) + qv * eqv * (1 -
        qv) - (qu.n * equ.n - qv * eqv)^2/(equ.n - eqv))/(equ.n -
        eqv)                                            # diagonal elements of -d^2 L/d beta^2
    if (sum(3 - cn) == 0) {
        g2 <- -as.vector(apply(as.matrix(temp.g2 * covar[cn ==
            3, ] * covar[cn == 3, ]), 2, sum))
    }
    else {
        g2 <- -as.vector(apply(as.matrix(temp.g2 * covar[cn ==
            3, ] * covar[cn == 3, ]), 2, sum) - apply(as.matrix(qu.c *
            covar[cn == 0, ] * covar[cn == 0, ]), 2, sum))
    }
    g2 <- g2 + eps                                      # for numerical stability
    derivs.ret <- list(l1 = l1, l2 = l2, g1 = g1, g2 = g2, likeli = likeli)
    return(derivs.ret)
}
