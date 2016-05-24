confint.expectreg <-
function (object, parm = NULL, level = 0.95, ...) 
{
    if (is.null(object$covmat)) 
        stop("No covariance matrix calculated.")
    res = list()
    if (any(object$design[, 1] != 1)) 
        center = FALSE
    else center = TRUE
    if (is.null(parm)) 
        for (i in 1:length(object$asymmetries)) {
            res[[i]] = matrix(NA, nrow = nrow(object$design), 
                ncol = 2)
            colnames(res[[i]]) = c(paste(eval((1 - level)/2), 
                "%"), paste(eval((1 + level)/2), "%"))
            for (j in 1:nrow(object$design)) {
                deviation = qnorm((1 + level)/2) * sqrt(t(object$design[j, 
                  ]) %*% object$covmat[[i]] %*% object$design[j, 
                  ])
                res[[i]][j, ] = c(fitted(object)[j, i] - deviation, 
                  fitted(object)[j, i] + deviation)
            }
        }
    else {
        nb = NULL
        for (i in 1:length(object$coefficients)) nb = c(nb, nrow(object$coefficients[[i]]))
        PB = NULL
        Koff = NULL
        for (k in 1:length(parm)) {
            co = which(names(object$covariates) == parm[k])
            partbasis = (sum(nb[0:(co - 1)]) + 1):(sum(nb[0:co]))
            if (center) 
                partbasis = partbasis + 1
            PB = c(PB, partbasis)
            Koff = rbind(Koff, object$coefficients[[co]])
        }
        if (center) {
            PB = c(1, PB)
            Koff = rbind(object$intercept, Koff)
        }
        B = object$design[, PB, drop = FALSE]
        fitti = B %*% Koff
        for (i in 1:length(object$asymmetries)) {
            res[[i]] = matrix(NA, nrow = nrow(B), ncol = 2)
            colnames(res[[i]]) = c(paste(eval((1 - level)/2), 
                "%"), paste(eval((1 + level)/2), "%"))
            for (j in 1:nrow(B)) {
                deviation = qnorm((1 + level)/2) * sqrt(t(B[j, 
                  ]) %*% object$covmat[[i]][PB, PB] %*% B[j, 
                  ])
                res[[i]][j, ] = c(fitti[j, i] - deviation, fitti[j, 
                  i] + deviation)
            }
        }
    }
    res
}
