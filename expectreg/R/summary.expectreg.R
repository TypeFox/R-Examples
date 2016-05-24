summary.expectreg <-
function (object, ...) 
{
    if (is.null(object$covmat)) 
        stop("No covariance matrix calculated.")
    if (any(object$design[, 1] != 1)) 
        center = FALSE
    else center = TRUE
    res = matrix(NA, nrow = 1, ncol = 2 * length(object$asymmetries))
    cnam = NULL
    for (i in 1:length(object$asymmetries)) {
        cnam = c(cnam, object$asymmetries[i], "sd")
        res[1, 2 * i - 1] = object$intercept[i]
        if (center) 
            res[1, 2 * i] = sqrt(object$covmat[[i]][1, 1])
        else res[1, 2 * i] = 0
    }
    rnam = "intercept"
    nb = NULL
    for (i in 1:length(object$coefficients)) nb = c(nb, nrow(object$coefficients[[i]]))
    for (k in 1:length(object$effects)) if (object$effects[k] == 
        "parametric" || object$effects[k] == "ridge") {
        if (nb[k] == 1) 
            rnam = c(rnam, names(object$covariates)[k])
        else for (n in 1:nb[k]) rnam = c(rnam, paste(names(object$covariates)[k], 
            n))
        partbasis = (sum(nb[0:(k - 1)]) + 1):(sum(nb[0:k]))
        if (center) 
            partbasis = partbasis + 1
        tmp = matrix(NA, nrow = nb[k], ncol = 2 * length(object$asymmetries))
        for (i in 1:length(object$asymmetries)) {
            tmp[, 2 * i - 1] = object$coef[[k]][, i]
            if (nb[k] == 1) 
                tmp[, 2 * i] = sqrt(object$covmat[[i]][partbasis, 
                  partbasis])
            else tmp[, 2 * i] = sqrt(diag(object$covmat[[i]][partbasis, 
                partbasis]))
        }
        res = rbind(res, tmp)
    }
    colnames(res) = cnam
    rownames(res) = rnam
    res
}
