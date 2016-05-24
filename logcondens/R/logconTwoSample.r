logconTwoSample <- function(x, y, which = c("MLE", "smooth"), M = 999, n.grid = 500, display = TRUE, seed0 = 1977){

    n1 <- length(x)
    n2 <- length(y)
    fac <- sqrt(n1 * n2 / (n1 + n2))
    test.stat <- matrix(NA, nrow = M, ncol = 2)
    
    set.seed(seed0)
    
    for (m in 1:M){
        per <- sample(c(x, y))
        x.per <- per[1:n1]
        y.per <- per[(n1 + 1):(n1 + n2)]
        res1 <- logConDens(x.per, smoothed = FALSE)
        res2 <- logConDens(y.per, smoothed = FALSE)

        test.stat[m, ] <- maxDiffCDF(res1, res2, which = which, n.grid = n.grid)$test.stat
        if (display){print(paste("Iteration ", m, " of ", M, " done", sep = ''))}
    }

    # calculate p-value
    test.stat[is.na(test.stat)] <- -Inf
    test.stat <- apply(test.stat, 2, sort)
    test.stat.orig <- maxDiffCDF(logConDens(x, smoothed = FALSE), logConDens(y, smoothed = FALSE), which = which, n.grid = n.grid)$test.stat

    ps <- rep(NA, 2)
    for (i in 1:2){ps[i] <- (1 + sum(test.stat[, i] >= test.stat.orig[i])) / (1 + M)}    
    test.stat[test.stat == -Inf] <- NA
    res <- list("p.value" = ps, "test.stat.orig" = test.stat.orig * fac, "test.stats" = test.stat * fac)
    return(res)
}







#
