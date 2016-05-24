vcov.zlm <-
function (object, include.const = FALSE, ...) 
{
    Xmat = as.matrix(model.frame(object)[, -1, drop = FALSE])
    if (ncol(Xmat) < 1) 
        stop("Needs at least one non-constant regressor")
    regnames = colnames(Xmat)
    Xmat = Xmat - matrix(colMeans(Xmat), nrow(Xmat), ncol(Xmat), 
        byrow = TRUE)
    xxinv = chol2inv(chol(crossprod(Xmat)))
    outmat = ((object$coef2moments[[2]] - object$coefficients[[2]]^2)/xxinv[[1]]) * 
        xxinv
    if (include.const) {
        outmat = rbind(rep(NA, nrow(outmat) + 1), cbind(rep(NA, 
            ncol(outmat)), outmat))
        regnames = c("(Intercept)", regnames)
    }
    colnames(outmat) <- rownames(outmat) <- regnames
    return(outmat)
}
