"w2.test" <-
function (x, distn, fit, H = NA, sim = 100, tol = 1e-04, estfun = NA)
{
    if (!is.function(try(get(distn), silent = TRUE)))
       stop("'distn' must be a character string of a distribution function")
    if (is.na(H)) H <- -Inf
    est <- try(eval(parse(text = estfun)), silent = TRUE)    
    if (!is.na(est[1])) {
        if (class(est) == "try-error")
            stop(paste("can not evaluate 'estfun':\n", est[1]))
        if (!is.vector(est) || !is.list(est))
            stop("value of 'estfun' is not a vector of 'mode(list)'")
        if (!all(names(est) %in% names(formals(distn))))
            stop("values of 'estfun' dosn't match arguments of 'distn'")
    }
    w2 <- function(n, zH, z, j) 
          n/3  + n*zH/(1-zH) + 1/(n*(1-zH))*sum((1-2*j)*z) + 1/(1-zH)^2 * sum((z-zH)^2)
    val <- mctest(x, distn, fit, H, sim, tol, STATISTIC = w2, estfun)
    STATISTIC <- val$TS
    names(STATISTIC) <- "W"
    res <- list(method = "Cramer-von Mises test", statistic = STATISTIC, 
        p.value = val$p.value, data.name = deparse(substitute(x)), 
        sim.no = val$sim, treshold = H)
    class(res) <- "mchtest"
    return(res)
}
