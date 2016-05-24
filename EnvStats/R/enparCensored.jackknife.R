enparCensored.jackknife <-
function (x, censored, censoring.side, correct.se, left.censored.min, 
    right.censored.max, est.fcn) 
{
    N <- length(x)
    jack.vec <- numeric(N)
    new.N <- N - 1
    for (i in 1:N) {
        new.x <- x[-i]
        new.censored <- censored[-i]
        new.n.cen <- sum(new.censored)
        if (new.n.cen == 0) 
            jack.vec[i] <- mean(new.x)
        else {
            new.x.cen <- new.x[new.censored]
            jack.vec[i] <- do.call(est.fcn, list(x = new.x, censored = new.censored, 
                censoring.side = censoring.side, correct.se = correct.se, 
                left.censored.min = left.censored.min, right.censored.max = right.censored.max, 
                ci = FALSE))$parameters["mean"]
        }
    }
    jack.vec
}
