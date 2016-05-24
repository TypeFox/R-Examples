"residuals.grouped" <-
function(object, standardized = FALSE, B = 100, ...){
    if(!inherits(object, "grouped"))
        stop("Use only with 'grouped' objects.\n")
    dets <- object$details
    rtrunc <- switch(dets$distribution,
        "normal" = function(n, mu, sigma, L, U) mu + sigma * qnorm(runif(n, pnorm(L, mu, sigma), pnorm(U, mu, sigma))), 
        "t" = function(n, mu, sigma, L, U) mu + sigma * qt(runif(n, pt((L - mu) / sigma, df), pt((U - mu) / sigma, df)), df),
        "logistic" = function(n, mu, sigma, L, U) mu + sigma * qlogis(runif(n, plogis(L, mu, sigma), plogis(U, mu, sigma))))
    n <- dets$n
    link <- dets$link
    y <- dets$y
    qa <- switch(link, "identity" = y[, 1], "log" = log(y[, 1]), "logit" = qlogis(y[, 1]))
    qb <- switch(link, "identity" = y[, 2], "log" = log(y[, 2]), "logit" = qlogis(y[, 2]))
    X <- dets$X
    px <- ncol(X)
    V <- vcov(object)
    df <- object$df
    sigma <- numeric(B)
    mat.res <- matrix(0.0, n, B)
    for(b in 1:B){
        theta <- mvrnorm(1, object$coefficients, V)
        sigma[b] <- theta[px + 1]
        mu <- c(X %*% theta[1:px])
        z <- rtrunc(n, mu, sigma[b], qa, qb)
        mat.res[, b] <- z - mu
    }
    res <- if(standardized){
        mean.res <- rowMeans(mat.res)
        var.btn <- rowSums((mat.res - mean.res)^2) / (B - 1)
	  mat.res <- t(t(mat.res) / sigma)
        sigma <- mean(sigma) + (1 + 1/B) * var.btn
        mean.res / sigma
    } else rowMeans(mat.res)
    lis <- list(residuals = res, mat.res = mat.res, nam.res = rownames(X), B = B, standardized = standardized, 
                    y = object$details$y, fitted = fitted(object))
    class(lis) <- "resid.grouped"
    lis
}

