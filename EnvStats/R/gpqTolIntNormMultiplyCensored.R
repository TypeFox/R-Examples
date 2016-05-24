gpqTolIntNormMultiplyCensored <-
function (n, cen.index, p, probs, nmc, method = "mle", censoring.side = "left", 
    seed = NULL, names = TRUE) 
{
    if (length(n) != 1 || !is.numeric(n) || n != trunc(n) || 
        n < 3) 
        stop("'n' must be an integer greater than or equal to 3")
    n.cen <- length(cen.index)
    if (n.cen < 1 || n.cen > (n - 2)) 
        stop("The length of 'cen.index' must be between 1 and 'n-2'")
    if (!is.numeric(cen.index) || !all(cen.index == trunc(cen.index)) || 
        any(cen.index < 1) || any(cen.index > n) || !all(cen.index == 
        sort(cen.index))) 
        stop(paste("'cen.index' must be a sorted vector of integers", 
            "with possible values between 1 and 'n'"))
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
    method <- match.arg(method, c("mle", "qq.reg", "impute.w.qq.reg", 
        "half.cen.level"))
    censoring.side <- match.arg(censoring.side, c("left", "right"))
    if (!is.null(seed)) 
        set.seed(seed)
    censored <- rep(FALSE, n)
    censored[cen.index] <- TRUE
    obs.index <- (1:n)[!censored]
    if (censoring.side == "left") {
        outer.mat <- outer(cen.index, obs.index, "<")
        apply.fcn <- function(x, n) {
            if (sum(x)) 
                ret.val <- min((1:n)[x])
            else ret.val <- n
        }
    }
    else {
        outer.mat <- outer(cen.index, obs.index, ">")
        apply.fcn <- function(x, n) {
            if (sum(x)) 
                ret.val <- max((1:n)[x])
            else ret.val <- 1
        }
    }
    set.index <- obs.index[apply(outer.mat, 1, apply.fcn, n = n - 
        n.cen)]
    z.p <- qnorm(p)
    t.vec <- numeric(nmc)
    for (i in 1:nmc) {
        x <- sort(rnorm(n))
        x[cen.index] <- x[set.index]
        params <- enormMultiplyCensored(x, censored, censoring.side = censoring.side, 
            method = method, ci = FALSE)$parameters
        t.vec[i] <- (z.p - params["mean"])/params["sd"]
    }
    quantile(x = t.vec, probs = probs, names = names)
}
