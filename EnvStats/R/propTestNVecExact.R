propTestNVecExact <-
function (p, p0, alpha, power, alternative, round.up, n.min, 
    n.max, tol.alpha, tol, maxiter) 
{
    dum.list <- propTestPower(n.or.n1 = n.min, p.or.p1 = p, p0.or.p2 = p0, 
        alpha = alpha, sample.type = "one.sample", alternative = alternative, 
        approx = FALSE, warn = FALSE, return.exact.list = TRUE)
    if (!is.na(dum.list$power) && (dum.list$power >= power & 
        dum.list$alpha <= alpha)) {
        ret.val <- c(n = n.min, unlist(dum.list))
    }
    else {
        fcn.for.root <- function(n, power, p0, p, alpha, alternative) {
            dum.list <- propTestPower(n.or.n1 = ceiling(n), p.or.p1 = p, 
                p0.or.p2 = p0, alpha = alpha, sample.type = "one.sample", 
                alternative = alternative, approx = FALSE, warn = FALSE, 
                return.exact.list = TRUE)
            test.power <- dum.list$power
            if (is.na(test.power)) 
                test.power <- alpha
            test.power - power
        }
        f.lower <- fcn.for.root(n = n.min, power, p0, p, alpha, 
            alternative)
        f.upper <- fcn.for.root(n = n.max, power, p0, p, alpha, 
            alternative)
        old.n <- uniroot(fcn.for.root, power = power, p0 = p0, 
            p = p, alpha = alpha, alternative = alternative, 
            lower = n.min, upper = n.max, f.lower = f.lower, 
            f.upper = f.upper, tol = tol, maxiter = maxiter)$root
        old.n <- ceiling(old.n)
        new.n.vec <- n.min:old.n
        dum.list <- propTestPower(n.or.n1 = new.n.vec, p.or.p1 = p, 
            p0.or.p2 = p0, alpha = alpha, sample.type = "one.sample", 
            alternative = alternative, approx = FALSE, warn = FALSE, 
            return.exact.list = TRUE)
        dum.power <- dum.list$power
        dum.power[is.na(dum.power)] <- alpha
        index <- (dum.power >= power) & (abs(dum.list$alpha - 
            alpha) <= tol.alpha)
        if (!any(index)) {
            while (!any(index) & max(new.n.vec) < n.max) {
                old.n <- max(new.n.vec) + 1
                new.n.vec <- old.n:min(old.n + 100, n.max)
                dum.list <- propTestPower(n.or.n1 = new.n.vec, 
                  p.or.p1 = p, p0.or.p2 = p0, alpha = alpha, 
                  sample.type = "one.sample", alternative = alternative, 
                  approx = FALSE, warn = FALSE, return.exact.list = TRUE)
                dum.power <- dum.list$power
                dum.power[is.na(dum.power)] <- alpha
                index <- (dum.list$power >= power) & (abs(dum.list$alpha - 
                  alpha) <= tol.alpha)
            }
        }
        index <- ((1:length(new.n.vec))[index])[1]
        ret.val <- c(n = new.n.vec[index], sapply(dum.list, function(x) x[index]))
    }
    ret.val
}
