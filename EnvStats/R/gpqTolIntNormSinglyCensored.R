gpqTolIntNormSinglyCensored <-
function (n, n.cen, p, probs, nmc, method = "mle", censoring.side = "left", 
    seed = NULL, names = TRUE) 
{
    if (length(n) != 1 || !is.numeric(n) || n != trunc(n) || 
        n < 3) 
        stop("'n' must be an integer greater than or equal to 3")
    if (length(n.cen) != 1 || !is.numeric(n.cen) || n.cen != 
        trunc(n.cen) || n.cen < 1 || n.cen > (n - 2)) 
        stop("'n.cen' must be an integer between 1 and 'n-2'")
    if (length(p) != 1 || !is.numeric(p) || p < .Machine$double.eps || 
        p > (1 - .Machine$double.eps)) 
        stop("'p' must be a numeric scalar strictly greater than 0 and strictly less than 1")
    if (!is.numeric(probs) || any(probs < .Machine$double.eps) || 
        any(probs > (1 - .Machine$double.eps))) 
        stop(paste("'probs' must be a numeric vector and all values of", 
            "'probs' must be strictly greater than 0 and strictly less than 1"))
    if (length(nmc) != 1 || !is.numeric(nmc) || nmc != trunc(nmc) || 
        nmc < 10) 
        stop("'nmc' must be an integer greater than or equal to 10")
    method <- match.arg(method, c("mle", "bcmle", "qq.reg", "qq.reg.w.cen.level", 
        "impute.w.qq.reg", "impute.w.qq.reg.w.cen.level", "impute.w.mle", 
        "iterative.impute.w.qq.reg", "half.cen.level", "m.est"))
    censoring.side <- match.arg(censoring.side, c("left", "right"))
    if (!is.null(seed)) 
        set.seed(seed)
    censored <- rep(c(TRUE, FALSE), c(n.cen, n - n.cen))
    if (censoring.side == "right") 
        censored <- rev(censored)
    cen.index <- (1:n)[censored]
    cen.level.index <- ifelse(censoring.side == "left", n.cen + 
        1, n - n.cen)
    z.p <- qnorm(p)
    t.vec <- numeric(nmc)
    for (i in 1:nmc) {
        x <- sort(rnorm(n))
        x[cen.index] <- x[cen.level.index]
        params <- enormSinglyCensored(x, censored, censoring.side = censoring.side, 
            method = method, ci = FALSE)$parameters
        t.vec[i] <- (z.p - params["mean"])/params["sd"]
    }
    quantile(x = t.vec, probs = probs, names = names)
}
