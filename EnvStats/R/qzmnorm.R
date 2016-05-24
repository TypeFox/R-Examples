qzmnorm <-
function (p, mean = 0, sd = 1, p.zero = 0.5) 
{
    names.p <- names(p)
    arg.mat <- cbind.no.warn(p = as.vector(p), mean = as.vector(mean), 
        sd = as.vector(sd), p.zero = as.vector(p.zero))
    na.index <- is.na.matrix(arg.mat)
    if (all(na.index)) 
        q <- rep(NA, nrow(arg.mat))
    else {
        q <- numeric(nrow(arg.mat))
        q[na.index] <- NA
        q.no.na <- q[!na.index]
        for (i in c("p", "mean", "sd", "p.zero")) assign(i, arg.mat[!na.index, 
            i])
        if (any(p < 0) || any(p > 1)) 
            stop("All non-missing values of 'p' must be between 0 and 1.")
        if (any(sd < .Machine$double.eps)) 
            stop("All non-missing values of 'sd' must be positive.")
        if (any(p.zero <= 0 | p.zero >= 1)) 
            stop(paste("All non-missing values of 'p.zero' must be", 
                "greater than 0 and less than 1."))
        p0 <- pzmnorm(0, mean, sd, p.zero)
        q.no.na[p == 0] <- -Inf
        q.no.na[p == p0] <- 0
        q.no.na[p == 1] <- Inf
        o.fcn <- function(q, mean, sd, p.zero, p) {
            ((1 - p.zero) * pnorm(q, mean, sd) - p)^2
        }
        index <- (1:length(q.no.na))[0 < p & p < p0]
        if (any(index)) {
            for (i in index) {
                q.no.na[i] <- nlminb(start = min((1 - p.zero[i]) * 
                  qnorm(p[i], mean[i], sd[i]), -.Machine$double.eps), 
                  o.fcn, upper = 0, mean = mean[i], sd = sd[i], 
                  p.zero = p.zero[i], p = p[i])$par
            }
        }
        index <- (1:length(q.no.na))[p0 < p & p < 1]
        if (any(index)) {
            for (i in index) {
                q.no.na[i] <- nlminb(start = max((1 - p.zero[i]) * 
                  qnorm(p[i], mean[i], sd[i]), 1e+07 * .Machine$double.eps), 
                  o.fcn, lower = 0, mean = mean[i], sd = sd[i], 
                  p.zero = p.zero[i], p = p[i] - p.zero[i])$par
            }
        }
        q[!na.index] <- q.no.na
    }
    if (!is.null(names.p)) 
        names(q) <- rep(names.p, length = length(q))
    q
}
