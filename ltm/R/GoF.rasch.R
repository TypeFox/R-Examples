GoF.rasch <-
function (object, B = 49, ...) {
    if (!inherits(object, "rasch"))
        stop("Use only with 'rasch' objects.\n")
    nas <- any(na.ind <- is.na(object$X))
    pearson.chi <- function (object) {
        fits <- fitted(object)
        X <- fits[, -ncol(fits), drop = FALSE]
        Obs <- observedFreqs(object, X)
        Exp <- fits[, ncol(fits)]
        if (any(ind <- Exp == 0))
            Exp[ind] <- 0.001
        if (nas)
            sum((Obs - Exp)^2/Exp, na.rm = TRUE)
        else
            sum((Obs - Exp)^2/Exp, na.rm = TRUE) + sum(Obs) - sum(Exp)
    }
    rmvlogis <- function (betas) {
        z <- rnorm(n)
        pr <- plogis(cbind(1, z) %*% t(betas))
        X <- matrix(0, n, p)
        for (i in 1:p)
            X[, i] <- ifelse(runif(n) < pr[, i], 1, 0)
        X
    }
    constraint <- object$constraint
    betas <- object$coef
    betas <- c(betas[, -2], betas[1, 2])
    if (!is.null(constraint))
        betas <- betas[-constraint[, 1]]
    Var.betas <- vcov(object)
    n <- nrow(object$X)
    p <- ncol(object$X)
    Tobs <- pearson.chi(object)
    Ts <- numeric(B)
    for (i in 1:B) {
        betas. <- mvrnorm(1, betas, Var.betas)
        betas. <- betas.rasch(betas., constraint, p)
        X <- rmvlogis(betas.)
        if (nas)
            X[na.ind] <- NA
        Ts[i] <- pearson.chi(rasch(X, constraint = constraint, start.val = c(betas.[, 1], betas.[1, 2]), 
                                    control = object$control))
    }
    p.val <- (1 + sum(Ts >= Tobs)) / (B + 1)
    out <- list(Tobs = Tobs, p.value = p.val, B = B, call = object$call)
    class(out) <- "GoF.rasch"
    out
}
