summary.lqa <-
function (object, dispersion = NULL, ...) 
{
### computation of dispersion parameter:

    est.disp <- FALSE
    df.r <- object$df.residual
    if (is.null (dispersion)) 
        dispersion <- if (object$family$family %in% c("poisson", "binomial")) 
            1
        else if (df.r > 0) {
            est.disp <- TRUE
            if (any (object$weights == 0)) 
                warning ("observations with zero weight not used for calculating dispersion")
            sum ((object$weights * object$residuals^2)[object$weights >  0]) / df.r
            }
            else {
              est.disp <- TRUE
              NaN
            }

    coef.p <- object$coefficients  
    coef.table <- as.matrix (coef.p)
    dimnames(coef.table) <- list (names (coef.p), "Estimate")

    ans <- c (object, list (deviance.resid = residuals (object, type = "deviance"), dispersion = dispersion))
    ans$coefficients <- coef.table
    class(ans) <- "summary.lqa"
    return(ans)
}

