elnormAltSinglyCensored.jackknife <-
function (x, censored, N, T1, censoring.side, est.fcn, ci.type, 
    conf.level, ...) 
{
    jack.vec <- numeric(N)
    new.N <- N - 1
    for (i in 1:N) {
        new.x <- x[-i]
        new.censored <- censored[-i]
        new.n.cen <- sum(new.censored)
        if (new.n.cen == 0) 
            jack.vec[i] <- elnormAlt(new.x)$parameters[1]
        else jack.vec[i] <- do.call(est.fcn, list(x = new.x, 
            censored = new.censored, N = new.N, T1 = T1, n.cen = new.n.cen, 
            censoring.side = censoring.side, ci = FALSE, ...))$parameters[1]
    }
    jack.vec
}
