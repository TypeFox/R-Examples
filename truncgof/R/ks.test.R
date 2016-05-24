"ks.test" <-
function (x, distn, fit, H = NA, 
    alternative = c("two.sided", "less", "greater"), sim = 100, tol = 1e-04, estfun = NA) 
{   
    if (!is.function(try(get(distn), silent = TRUE)))
       stop("'distn' must be a character string of a distribution function")
    if (is.na(H)) H <- -Inf 
    est <- try(eval(parse(text = estfun)), silent = TRUE)
    if (!is.na(est[1])) {
        if (class(est) == "try-error")
            stop(paste("can not evaluat 'estfun':\n", est[1]))
        if (!is.vector(est) || !is.list(est))
            stop("value of 'estfun' is not a vector of 'mode(list)'")
        if (!all(names(est) %in% names(formals(distn))))
            stop("values of 'estfun' dosn't match arguments in 'distn'")
    }
    ks.p <- function(n, zH, z, j) sqrt(n)/(1-zH) * max(zH + j/n * (1-zH) - z)
    ks.m <- function(n, zH, z, j) sqrt(n)/(1-zH) * max(z - (zH + (j-1)/n * (1-zH)))
    ks   <- function(n, zH, z, j) max(sqrt(n)/(1-zH) * max(zH + j/n * (1-zH) - z), 
        sqrt(n)/(1-zH) * max(z - (zH + (j-1)/n * (1-zH))))
    AL <- match.arg(alternative)
    TS <- switch(AL, "two.sided" = ks, "less" = ks.m, "greater" = ks.p)
    val <- mctest(x, distn, fit, H, sim, tol, TS, estfun)
    STATISTIC <- val$TS
    names(STATISTIC) <- switch(AL, two.sided = "KS", less = "KS-", greater = "KS+")
    res <- list(method = "Kolmogorov-Smirnov Test", statistic = STATISTIC, 
        p.value = val$p.value, alternative = AL, method = "", 
        data.name = deparse(substitute(x)), sim.no = val$sim, treshold = H)
    class(res) <- "mchtest"
    return(res)  
}
