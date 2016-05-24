egammaAltCensored.jackknife <-
function (x, censored, censoring.side, est.fcn) 
{
    N <- length(x)
    jack.vec <- numeric(N)
    new.N <- N - 1
    for (i in 1:N) {
        new.x <- x[-i]
        new.censored <- censored[-i]
        new.n.cen <- sum(new.censored)
        if (new.n.cen == 0) 
            jack.vec[i] <- egammaAlt(new.x)$parameters["mean"]
        else {
            jack.vec[i] <- do.call(est.fcn, list(x = new.x, censored = new.censored, 
                censoring.side = censoring.side, ci = FALSE))$parameters["mean"]
        }
    }
    jack.vec
}
