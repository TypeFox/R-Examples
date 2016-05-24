mda.fit <-
function (x, g, weights, theta, assign.theta, Rj, sub.df=NULL, tot.df=NULL, 
    dimension = R - 1, eps = .Machine$double.eps, method = polyreg, 
    ...) 
{
    this.call <- match.call()
    n <- nrow(x)
    cnames <- names(weights)
    J <- length(cnames)
    R <- sum(Rj)
    wtots <- lapply(weights, function(x) apply(x, 2, sum))
    sub.prior <- lapply(wtots, function(x) x/sum(x))
    dp <- unlist(wtots)
    subclass.names <- names(dp)
    dp <- dp/sum(dp)
    if (missing(theta)) 
        theta <- contr.helmert(R)
    theta <- contr.fda(dp, theta)
    if (!(is.null(sub.df) & is.null(tot.df))) {
        Q <- diag(dp) + transformPenalty(prior = dp, cl = rep(seq(J), 
            Rj), df = sub.df, tot.df = tot.df)
        theta <- fix.theta(theta, Q)
    }
    Theta <- matrix(0, n, R - 1)
    obs.weights <- double(n)
    for (j in seq(J)) {
        Theta[g == j, ] <- weights[[j]] %*% theta[assign.theta[[j]], 
            , drop = FALSE]
        obs.weights[g == j] <- weights[[j]] %*% rep(1, Rj[j])
    }
    fit <- method(x, Theta, obs.weights, ...)
    Theta <- Theta * obs.weights
    ssm <- t(Theta) %*% fitted(fit)/n
    ed <- svd(ssm, nu = 0)
    thetan <- ed$v
    lambda <- ed$d
    lambda[lambda > 1 - eps] <- 1 - eps
    discr.eigen <- lambda/(1 - lambda)
    pe <- (100 * cumsum(discr.eigen))/sum(discr.eigen)
    dimension <- min(dimension, sum(lambda > eps))
    if (dimension == 0) {
        warning("degenerate problem; no discrimination")
        return(structure(list(dimension = 0, fit = fit, call = this.call), 
            class = "fda"))
    }
    thetan <- thetan[, seq(dimension), drop = FALSE]
    pe <- pe[seq(dimension)]
    alpha <- sqrt(lambda[seq(dimension)])
    sqima <- sqrt(1 - lambda[seq(dimension)])
    vnames <- paste("v", seq(dimension), sep = "")
    means <- scale(theta %*% thetan, FALSE, sqima/alpha)
    dimnames(means) <- list(subclass.names, vnames)
    names(lambda) <- c(vnames, rep("", length(lambda) - dimension))
    names(pe) <- vnames
    list(percent.explained = pe, values = lambda, means = means, 
        theta.mod = thetan, dimension = dimension, sub.prior = sub.prior, 
        fit = fit, call = this.call)
}

