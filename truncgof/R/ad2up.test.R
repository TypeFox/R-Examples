"ad2up.test" <-
function (x, distn, fit, H = NA, sim = 100, tol = 1e-04, estfun = NA)
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
    ad2up <- function(n, zH, z, j)
             -2*n*log(1-zH) + 2*sum(log(1-z)) + (1-zH)/n * sum((1+2*(n-j))/(1-z))
    val <- mctest(x, distn, fit, H, sim, tol, STATISTIC = ad2up, estfun)
    STATISTIC <- val$TS
    names(STATISTIC) <- "AD2up"
    res <- list(method = "Quadratic Class Anderson-Darling Upper Tail Test", statistic = STATISTIC, 
        p.value = val$p.value, data.name = deparse(substitute(x)), 
        sim.no = val$sim, treshold = H)
    class(res) <- "mchtest"
    return(res)
}
