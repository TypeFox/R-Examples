poisson.mtest <- 
function(x, R = 999) {
    # parametric bootstrap mean distance test of Poisson distribution
    n <- length(x)
    lambda <- mean(x)
    bootobj <- boot(x, statistic = poisson.m, R = R, sim = "parametric", 
            ran.gen = function(x, y) {rpois(n, lambda)})
    p <- 1 - mean(bootobj$t < bootobj$t0)
    names(bootobj$t0) <- "test statistic"
    names(lambda) <- "mean"
    e <- list(
        method = paste("Mean distance test of Poisson distribution", sep = ""),
        statistic = bootobj$t0, 
        p.value = p, 
        data.name = paste("sample size ", n, ", replicates ", R, sep=""),
        estimate = lambda)
    class(e) <- "htest"        
    e           
}

poisson.m<- 
function(x) {
    # mean distance statistic for Poissonity
    n <- length(x)
    stat <- 0
    e <- .C("poisMstat", 
            x = as.integer(x),
            nx = as.integer(n), 
            stat = as.double(stat), 
            PACKAGE = "energy")$stat
    e
}
