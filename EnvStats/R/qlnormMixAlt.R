qlnormMixAlt <-
function (p, mean1 = exp(1/2), cv1 = sqrt(exp(1) - 1), mean2 = exp(1/2), 
    cv2 = sqrt(exp(1) - 1), p.mix = 0.5) 
{
    names.p <- names(p)
    arg.mat <- cbind.no.warn(p = as.vector(p), mean1 = as.vector(mean1), 
        cv1 = as.vector(cv1), mean2 = as.vector(mean2), cv2 = as.vector(cv2), 
        p.mix = as.vector(p.mix))
    na.index <- is.na.matrix(arg.mat)
    if (all(na.index)) 
        q <- rep(NA, nrow(arg.mat))
    else {
        q <- numeric(nrow(arg.mat))
        q[na.index] <- NA
        q.no.na <- q[!na.index]
        for (i in c("p", "mean1", "cv1", "mean2", "cv2", "p.mix")) assign(i, 
            arg.mat[!na.index, i])
        if (any(p < 0 | p > 1)) 
            stop("All non-missing values of 'p' must be between 0 and 1.")
        if (any(c(mean1, mean2, cv1, cv2) < .Machine$double.eps)) 
            stop("All non-missing values of 'mean1', 'mean2', 'cv1', and 'cv2' must be positive.")
        if (any(p.mix < 0 | p.mix > 1)) 
            stop("All non-missing values of 'p.mix' must be between 0 and 1.")
        q.no.na[p == 0] <- 0
        q.no.na[p == 1] <- Inf
        index <- (1:length(q.no.na))[0 < p & p < 1]
        if (any(index)) {
            o.fcn <- function(q, mean1, cv1, mean2, cv2, p.mix, 
                p) {
                (plnormMixAlt(q, mean1, cv1, mean2, cv2, p.mix) - 
                  p)^2
            }
            for (i in index) {
                q.no.na[i] <- nlminb(start = (1 - p.mix[i]) * 
                  qlnormAlt(p[i], mean1[i], cv1[i]) + p.mix[i] * 
                  qlnormAlt(p[i], mean2[i], cv2[i]), o.fcn, lower = 0, 
                  mean1 = mean1[i], cv1 = cv1[i], mean2 = mean2[i], 
                  cv2 = cv2[i], p.mix = p.mix[i], p = p[i])$par
            }
        }
        q[!na.index] <- q.no.na
    }
    if (!is.null(names.p)) 
        names(q) <- rep(names.p, length = length(q))
    q
}
