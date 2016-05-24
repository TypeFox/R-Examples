shrink.mda <-
function (object, sub.df = NULL, tot.df = NULL, ...) 
{
    if (is.null(sub.df) & is.null(tot.df)) {
        warning("No shrinkage parameters (sub.df or tot.df) given")
        return(object)
    }
    dimension <- object$dimension
    lambda <- object$values[seq(dimension)]
    theta.mod <- object$theta.mod
    theta <- object$means
    alpha <- sqrt(lambda)
    sqima <- sqrt(1 - lambda)
    theta <- scale(theta, FALSE, alpha/sqima)
    sub.prior <- unlist(object$sub.prior)
    prior <- object$prior
    Rj <- sapply(object$assign.theta, length)
    dp <- sub.prior * rep(prior, Rj)
    cl <- rep(seq(Rj), Rj)
    P <- diag(dp) + transformPenalty(prior = dp, cl = cl, df = sub.df, 
        tot.df = tot.df)
    K <- t(theta) %*% P %*% theta
    TT <- chol((K + t(K))/2)
    Tinv <- backsolve(TT, diag(length(lambda)))
    M <- t(Tinv) %*% (lambda * Tinv)
    ed <- svd(M)
    thetan <- ed$v
    lambda <- ed$d
    discr.eigen <- lambda/(1 - lambda)
    pe <- (100 * cumsum(discr.eigen))/sum(discr.eigen)
    dimension <- min(dimension, sum(lambda > .Machine$double.eps))
    if (dimension == 0) {
        warning("degenerate problem; no discrimination")
        return(structure(list(dimension = 0, fit = object$fit,
                              call = object$call), 
                         class = "fda"))
    }
    thetan <- thetan[, seq(dimension), drop = FALSE]
    pe <- pe[seq(dimension)]
    alpha <- sqrt(lambda[seq(dimension)])
    sqima <- sqrt(1 - lambda[seq(dimension)])
    dm <- dimnames(theta)
    vnames <- dm[[2]][seq(dimension)]
    means <- scale(theta %*% Tinv %*% thetan, FALSE, sqima/alpha)
    dimnames(means) <- list(dm[[1]], vnames)
    names(lambda) <- c(vnames, rep("", length(lambda) - dimension))
    names(pe) <- vnames
    theta.mod <- theta.mod %*% Tinv %*% thetan
    object$confusion <- object$deviance <- NULL
    incl.names <- c("percent.explained", "values", "means", "theta.mod", 
        "dimension")
    rl <- list(pe, lambda, means, theta.mod, dimension)
    names(rl) <- incl.names
    object$sub.df <- sub.df
    object$tot.df <- tot.df
    object[incl.names] <- rl
    object$weights <- NULL
    update(object, sub.df = sub.df, tot.df = tot.df, weights = object, 
        ...)
}

