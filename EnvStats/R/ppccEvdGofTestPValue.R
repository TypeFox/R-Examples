ppccEvdGofTestPValue <-
function (r, n, tol = 1e-07) 
{
    stat <- 1000 * (1 - r)
    n.vec <- as.numeric(dimnames(.ppcc.evd.gof.crit.val.mat)[[1]])
    sig.vals.vec <- as.numeric(dimnames(.ppcc.evd.gof.crit.val.mat)[[2]])
    if (!is.na(index <- match(n, n.vec))) {
        crit.vals <- .ppcc.evd.gof.crit.val.mat[index, ]
    }
    else {
        crit.vals <- numeric(6)
        for (i in 1:6) crit.vals[i] <- exp(approx(log(n.vec), 
            log(.ppcc.evd.gof.crit.val.mat[, i]), xout = log(n))$y)
    }
    if (stat - crit.vals[1] > tol) 
        p.val <- "< 0.005"
    else if (crit.vals[6] - stat > tol) 
        p.val <- "> 0.5"
    else p.val <- exp(approx(log(crit.vals), log(sig.vals.vec), 
        xout = log(stat))$y)
    p.val
}
