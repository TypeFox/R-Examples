qnormMix <-
function (p, mean1 = 0, sd1 = 1, mean2 = 0, sd2 = 1, p.mix = 0.5) 
{
    names.p <- names(p)
    arg.mat <- cbind.no.warn(p = as.vector(p), mean1 = as.vector(mean1), 
        sd1 = as.vector(sd1), mean2 = as.vector(mean2), sd2 = as.vector(sd2), 
        p.mix = as.vector(p.mix))
    na.index <- is.na.matrix(arg.mat)
    if (all(na.index)) 
        q <- rep(NA, nrow(arg.mat))
    else {
        q <- numeric(nrow(arg.mat))
        q[na.index] <- NA
        q.no.na <- q[!na.index]
        for (i in c("p", "mean1", "sd1", "mean2", "sd2", "p.mix")) assign(i, 
            arg.mat[!na.index, i])
        if (any(p < 0 | p > 1)) 
            stop("All non-missing values of 'p' must be between 0 and 1.")
        if (any(c(sd1, sd2) < .Machine$double.eps)) 
            stop("All non-missing values of 'sd1' and 'sd2' must be positive.")
        if (any(p.mix < 0 | p.mix > 1)) 
            stop("All non-missing values of 'p.mix' must be between 0 and 1.")
        q.no.na[p == 0] <- -Inf
        q.no.na[p == 1] <- Inf
        index <- (1:length(q.no.na))[0 < p & p < 1]
        if (any(index)) {
            o.fcn <- function(q, mean1, sd1, mean2, sd2, p.mix, 
                p) {
                (pnormMix(q, mean1, sd1, mean2, sd2, p.mix) - 
                  p)^2
            }
            for (i in index) {
                q.no.na[i] <- nlminb(start = (1 - p.mix[i]) * 
                  qnorm(p[i], mean1[i], sd1[i]) + p.mix[i] * 
                  qnorm(p[i], mean2[i], sd2[i]), o.fcn, mean1 = mean1[i], 
                  sd1 = sd1[i], mean2 = mean2[i], sd2 = sd2[i], 
                  p.mix = p.mix[i], p = p[i])$par
            }
        }
        q[!na.index] <- q.no.na
    }
    if (!is.null(names.p)) 
        names(q) <- rep(names.p, length = length(q))
    q
}
