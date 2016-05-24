qlnormMix <-
function (p, meanlog1 = 0, sdlog1 = 1, meanlog2 = 0, sdlog2 = 1, 
    p.mix = 0.5) 
{
    names.p <- names(p)
    arg.mat <- cbind.no.warn(p = as.vector(p), meanlog1 = as.vector(meanlog1), 
        sdlog1 = as.vector(sdlog1), meanlog2 = as.vector(meanlog2), 
        sdlog2 = as.vector(sdlog2), p.mix = as.vector(p.mix))
    na.index <- is.na.matrix(arg.mat)
    if (all(na.index)) 
        q <- rep(NA, nrow(arg.mat))
    else {
        q <- numeric(nrow(arg.mat))
        q[na.index] <- NA
        q.no.na <- q[!na.index]
        for (i in c("p", "meanlog1", "sdlog1", "meanlog2", "sdlog2", 
            "p.mix")) assign(i, arg.mat[!na.index, i])
        if (any(p < 0 | p > 1)) 
            stop("All non-missing values of 'p' must be between 0 and 1.")
        if (any(c(sdlog1, sdlog2) < .Machine$double.eps)) 
            stop("All non-missing values of 'sdlog1' and 'sdlog2' must be positive.")
        if (any(p.mix < 0 | p.mix > 1)) 
            stop("All non-missing values of 'p.mix' must be between 0 and 1.")
        q.no.na[p == 0] <- 0
        q.no.na[p == 1] <- Inf
        index <- (1:length(q.no.na))[0 < p & p < 1]
        if (any(index)) {
            o.fcn <- function(q, meanlog1, sdlog1, meanlog2, 
                sdlog2, p.mix, p) {
                (plnormMix(q, meanlog1, sdlog1, meanlog2, sdlog2, 
                  p.mix) - p)^2
            }
            for (i in index) {
                q.no.na[i] <- nlminb(start = (1 - p.mix[i]) * 
                  qlnorm(p[i], meanlog1[i], sdlog1[i]) + p.mix[i] * 
                  qlnorm(p[i], meanlog2[i], sdlog2[i]), o.fcn, 
                  lower = 0, meanlog1 = meanlog1[i], sdlog1 = sdlog1[i], 
                  meanlog2 = meanlog2[i], sdlog2 = sdlog2[i], 
                  p.mix = p.mix[i], p = p[i])$par
            }
        }
        q[!na.index] <- q.no.na
    }
    if (!is.null(names.p)) 
        names(q) <- rep(names.p, length = length(q))
    q
}
