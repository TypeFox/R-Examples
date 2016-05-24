plnormMixAlt <-
function (q, mean1 = exp(1/2), cv1 = sqrt(exp(1) - 1), mean2 = exp(1/2), 
    cv2 = sqrt(exp(1) - 1), p.mix = 0.5) 
{
    names.q <- names(q)
    arg.mat <- cbind.no.warn(q = as.vector(q), mean1 = as.vector(mean1), 
        cv1 = as.vector(cv1), mean2 = as.vector(mean2), cv2 = as.vector(cv2), 
        p.mix = as.vector(p.mix))
    for (i in c("q", "mean1", "cv1", "mean2", "cv2", "p.mix")) assign(i, 
        arg.mat[, i])
    na.index <- is.na.matrix(arg.mat)
    if (all(na.index)) 
        p <- rep(NA, length(q))
    else {
        if (any(c(mean1[!na.index], mean2[!na.index], cv1[!na.index], 
            cv2[!na.index]) < .Machine$double.eps)) 
            stop("All values of 'mean1', 'mean2', 'cv1', and 'cv2' must be positive.")
        if (any(p.mix[!na.index] < 0 | p.mix[!na.index] > 1)) 
            stop("All values of 'p.mix' must be between 0 and 1.")
        p <- (1 - p.mix) * plnormAlt(q, mean1, cv1) + p.mix * 
            plnormAlt(q, mean2, cv2)
    }
    if (!is.null(names.q)) 
        names(p) <- rep(names.q, length = length(p))
    else names(p) <- NULL
    p
}
