"summary.resid.grouped" <-
function(object, K = 2, observed = NULL, ...){
    if(!inherits(object, "resid.grouped"))
        stop("Use only with 'resid.grouped' objects.\n")
    if(!object$standardized)
        warning("for outlier detection you should use standardized residuals.\n")
    y <- object$y
    nam.P <- paste("P(|e|>", K, ")", sep = "")
    dat <- cbind("fitted" = round(object$fitted, 3), 
                 "resid" = round(object$residuals, 3),
                 round(rowMeans(abs(object$mat.res) > K), 3))
    if(!is.null(observed))
        dat <- cbind("y" = round(observed, 3), dat)
    colnames(dat)[ncol(dat)] <- nam.P
    rownames(dat) <- object$nam.res
    dat
}

