sfSinglyCensoredGofTest <-
function (x, censored, censoring.side = c("left", "right"), distribution = c("norm", 
    "lnorm", "lnormAlt"), est.arg.list = NULL) 
{
    if (!is.vector(x, mode = "numeric") || is.factor(x)) 
        stop("'x' must be a numeric vector")
    if (!((is.vector(censored, mode = "numeric") && !is.factor(censored)) || 
        is.vector(censored, mode = "logical"))) 
        stop("'censored' must be a logical or numeric vector")
    if (length(censored) != length(x)) 
        stop("'censored' must be the same length as 'x'")
    data.name <- deparse(substitute(x))
    censoring.name <- deparse(substitute(censored))
    if ((bad.obs <- sum(!(ok <- is.finite(x) & is.finite(as.numeric(censored))))) > 
        0) {
        is.not.finite.warning(x)
        is.not.finite.warning(as.numeric(censored))
        x <- x[ok]
        censored <- censored[ok]
        warning(paste(bad.obs, "observations with NA/NaN/Inf in 'x' and 'censored' removed."))
    }
    if (is.numeric(censored)) {
        if (!all(censored == 0 | censored == 1)) 
            stop(paste("When 'censored' is a numeric vector, all values of", 
                "'censored' must be 0 (not censored) or 1 (censored)."))
        censored <- as.logical(censored)
    }
    n.cen <- sum(censored)
    if (n.cen == 0) {
        warning(paste("No censored values indicated by 'censored',", 
            "so the function 'sfGofTest' was called."))
        ret.list <- sfGofTest(x = x, distribution = distribution, 
            est.arg.list = est.arg.list)
        ret.list$data.name <- data.name
        ret.list$bad.obs <- bad.obs
        return(ret.list)
    }
    x.no.cen <- x[!censored]
    if (length(unique(x.no.cen)) < 2) 
        stop("'x' must contain at least 2 non-missing, uncensored, distinct values.")
    N <- length(x)
    n <- N - n.cen
    if (N < 5) 
        stop("'x' must contain at least 5 non-missing values.")
    if (N < 20 || N > 5000) 
        warning(paste("Too few or too many observations.  This approximation only works", 
            "if the number of observations is between 20 and 5000."))
    censoring.side <- match.arg(censoring.side)
    distribution <- match.arg(distribution)
    if (any(distribution == c("lnorm", "lnormAlt")) && any(x <= 
        0)) 
        stop("All values of 'x' must be positive for a lognormal distribution")
    T1 <- unique(x[censored])
    if (length(T1) > 1) 
        stop("More than one censoring level.")
    if (censoring.side == "left") {
        if (T1 > min(x.no.cen)) 
            stop(paste("For singly left-censored data,", "all uncensored observations must be bigger than", 
                "or equal to the censoring level."))
    }
    else {
        if (T1 < max(x.no.cen)) 
            stop(paste("For singly right-censored data,", "all uncensored observations must be less than", 
                "or equal to the censoring level."))
    }
    est.fcn <- paste("e", distribution, "SinglyCensored", sep = "")
    ret.list <- do.call(est.fcn, c(list(x = x, censored = censored, 
        censoring.side = censoring.side), est.arg.list))
    nrl <- names(ret.list)
    names(ret.list)[match("parameters", nrl)] <- "distribution.parameters"
    names(ret.list)[match("method", nrl)] <- "estimation.method"
    ret.list$data.name <- data.name
    ret.list$censoring.name <- censoring.name
    ret.list$bad.obs <- bad.obs
    ret.list$dist.abb <- distribution
    new.x <- switch(distribution, norm = sort(x.no.cen), lnorm = , 
        lnormAlt = sort(log(x.no.cen)), )
    m.tilda <- qnorm(ppoints(N, a = 3/8))
    ss.m <- sum(m.tilda^2)
    c.vec <- m.tilda/sqrt(ss.m)
    if (censoring.side == "right") 
        index <- 1:n
    else index <- (n.cen + 1):N
    W <- cor(c.vec[index], new.x)^2
    w <- log(1 - W)
    nu <- log(N)
    u.mu <- log(nu) - nu
    mu <- -1.2725 + 1.0521 * u.mu
    u.sigma <- log(nu) + 2/nu
    sigma <- 1.0308 - 0.26758 * u.sigma
    z <- (w - mu)/sigma
    u <- log(N)
    D <- 0.76676 * u + 0.015814 * u^2
    R <- numeric(3)
    alpha <- c(0.9, 0.95, 0.99)
    names(R) <- format(alpha)
    R[c("0.90", "0.95")] <- .royston.93.array[c("0.90", "0.95"), 
        "W", "A"] + .royston.93.array[c("0.90", "0.95"), "W", 
        "B"] * .royston.93.array[c("0.90", "0.95"), "W", "C"]^u
    R["0.99"] <- .royston.93.array["0.99", "W", "A"] + .royston.93.array["0.99", 
        "W", "B"] * u
    delta <- n.cen/N
    qna <- qnorm(alpha)
    Z <- qna + D * R^(-log(delta))
    fit.coef <- lm(Z ~ qna)$coef
    mu.z <- fit.coef[1]
    sigma.z <- fit.coef[2]
    z.prime <- (z - mu.z)/sigma.z
    names(z.prime) <- NULL
    p <- 1 - pnorm(z.prime)
    sep.string <- paste("\n", space(33), sep = "")
    ret.list <- c(ret.list, list(statistic = W, parameters = c(N, 
        delta), z.value = z.prime, p.value = p, alternative = paste("True cdf does not equal the", 
        sep.string, ret.list$distribution, " Distribution.", 
        sep = ""), method = paste("Shapiro-Francia GOF", "(Singly Censored Data)", 
        sep = sep.string), data = x, censored = censored))
    names(ret.list$statistic) <- "W"
    names(ret.list$parameters) <- c("N", "DELTA")
    ret.list <- ret.list[c("distribution", "dist.abb", "distribution.parameters", 
        "n.param.est", "estimation.method", "statistic", "sample.size", 
        "censoring.side", "censoring.levels", "percent.censored", 
        "parameters", "z.value", "p.value", "alternative", "method", 
        "data", "data.name", "censored", "censoring.name", "bad.obs")]
    oldClass(ret.list) <- "gofCensored"
    ret.list
}
