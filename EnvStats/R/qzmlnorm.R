qzmlnorm <-
function (p, meanlog = 0, sdlog = 1, p.zero = 0.5) 
{
    names.p <- names(p)
    arg.mat <- cbind.no.warn(p = as.vector(p), meanlog = as.vector(meanlog), 
        sdlog = as.vector(sdlog), p.zero = as.vector(p.zero))
    na.index <- is.na.matrix(arg.mat)
    if (all(na.index)) 
        q <- rep(NA, nrow(arg.mat))
    else {
        q <- numeric(nrow(arg.mat))
        q[na.index] <- NA
        q.no.na <- q[!na.index]
        for (i in c("p", "meanlog", "sdlog", "p.zero")) assign(i, 
            arg.mat[!na.index, i])
        if (any(p < 0) || any(p > 1)) 
            stop("All non-missing values of 'p' must be between 0 and 1.")
        if (any(sdlog < .Machine$double.eps)) 
            stop("All non-missing values of 'sdlog' must be positive.")
        if (any(p.zero <= 0 | p.zero >= 1)) 
            stop(paste("All non-missing values of 'p.zero' must be", 
                "greater than 0 and less than 1."))
        q.no.na[p <= p.zero] <- 0
        q.no.na[p == 1] <- Inf
        index <- (1:length(q.no.na))[p.zero < p & p < 1]
        if (any(index)) {
            o.fcn <- function(q, meanlog, sdlog, p.zero, p) {
                ((1 - p.zero) * plnorm(q, meanlog, sdlog) - p)^2
            }
            for (i in index) {
                q.no.na[i] <- nlminb(start = (1 - p.zero[i]) * 
                  qlnorm(p[i], meanlog[i], sdlog[i]), o.fcn, 
                  lower = 0, meanlog = meanlog[i], sdlog = sdlog[i], 
                  p.zero = p.zero[i], p = p[i] - p.zero[i])$par
            }
        }
        q[!na.index] <- q.no.na
    }
    if (!is.null(names.p)) 
        names(q) <- rep(names.p, length = length(q))
    q
}
