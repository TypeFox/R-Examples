elnormAltMultiplyCensored.jackknife <-
function (x, censored, N, cen.levels, K, c.vec, n.cen, censoring.side, 
    est.fcn, ci.type, conf.level, ...) 
{
    jack.vec <- numeric(N)
    new.N <- N - 1
    for (i in 1:N) {
        new.x <- x[-i]
        new.censored <- censored[-i]
        new.n.cen <- sum(new.censored)
        if (new.n.cen == 0) 
            jack.vec[i] <- elnormAlt(new.x)$parameters[1]
        else {
            new.x.cen <- new.x[new.censored]
            new.c.vec <- table(new.x.cen)
            new.cen.levels <- sort(unique(new.x.cen))
            new.K <- length(new.cen.levels)
            jack.vec[i] <- do.call(est.fcn, list(x = new.x, censored = new.censored, 
                N = new.N, cen.levels = new.cen.levels, K = new.K, 
                c.vec = new.c.vec, n.cen = new.n.cen, censoring.side = censoring.side, 
                ci = FALSE, ...))$parameters[1]
        }
    }
    jack.vec
}
